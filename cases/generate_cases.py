#!/usr/bin/env python3
"""
generate_cases.py -- build OpenFOAM case folders from caseMatrix.csv

One row of the CSV = one simulation. The CSV is the single source of truth for
the whole parametric study; nothing about a case should be hand-edited inside
the generated folder, because regenerating would silently discard it.

Layout expected:

    cases/
      caseMatrix.csv
      generate_cases.py
      templates/
        modelD/            <- template case, tokens written as @column_name@
        modelP/
        modelT/
      10-modelD-droplet/   <- generated output
      20-modelP-pool/
      30-modelT-thermal/

Two substitution mechanisms, because OpenFOAM needs both:

  1. constant/caseParameters -- an OpenFOAM dictionary written into every case.
     Template dictionaries do  #include "$FOAM_CASE/constant/caseParameters"
     and then reference values as  $I_A,  $dt_s  etc. Preferred: values stay
     live, and you can tweak one case by hand without regenerating.

  2. @token@ text substitution across every text file in the template. Needed
     where OpenFOAM's own $-expansion does not reach -- most importantly inside
     the verbatim  #{ ... #}  blocks of codedFvOption sources.

Usage:
    python3 generate_cases.py --list
    python3 generate_cases.py --dry-run
    python3 generate_cases.py --only "M0*"
    python3 generate_cases.py --force

Author: generated for the MHD-WAAM thesis. Python 3.8+, standard library only.
"""

import argparse
import csv
import fnmatch
import json
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

# --------------------------------------------------------------------------
# configuration
# --------------------------------------------------------------------------

MODEL_DIRS = {
    "D": "10-modelD-droplet",
    "P": "20-modelP-pool",
    "T": "30-modelT-thermal",
}

REQUIRED_COLUMNS = {
    "case_id", "model", "enabled", "I_A", "U_V", "travel_speed_m_s",
    "wire_feed_m_s", "wire_dia_m", "B_mT", "B_dir", "B_freq_Hz",
    "cell_size_m", "deposit_len_m", "dt_s", "sym",
}

# Files that must never be token-substituted or treated as text.
BINARY_SUFFIXES = {".png", ".jpg", ".jpeg", ".pdf", ".gz", ".zip", ".stl",
                   ".vtu", ".vtk", ".catemp", ".tar"}

# Geometry defaults per model, in metres. Any of the keys x_min_m, x_max_m,
# y_min_m, y_max_m, z_min_m, z_max_m, nx, ny, nz may also appear as CSV columns,
# in which case the CSV value wins. Add a column only for the cases that need it.
#
#   D  small box around the wire tip -- no travel, no substrate to speak of
#   P  substrate plus air, spanning the deposit with a margin at each end
#   T  as P but taller, to clear five stacked layers
MODEL_GEOM = {
    "D": {"pad_x_m": 0.004, "span_x": False, "half_width_m": 0.004,
          "substrate_z_m": 0.002, "air_z_m": 0.010},
    "P": {"pad_x_m": 0.005, "span_x": True, "half_width_m": 0.015,
          "substrate_z_m": 0.004, "air_z_m": 0.012},
    "T": {"pad_x_m": 0.005, "span_x": True, "half_width_m": 0.015,
          "substrate_z_m": 0.004, "air_z_m": 0.020},
}

GEOM_KEYS = ("x_min_m", "x_max_m", "y_min_m", "y_max_m", "z_min_m", "z_max_m")

# Rough planning figure for this solver class: VOF + energy + enthalpy-porosity
# + electric potential. Used only to warn before you queue something that will
# not fit in a 5 GB WSL allocation.
GB_PER_MILLION_CELLS = 3.5
RAM_WARN_GB = 4.0

TOKEN_RE = re.compile(r"@([A-Za-z_][A-Za-z0-9_]*)@")


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def truthy(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def as_float(row, key, default=None):
    raw = (row.get(key) or "").strip()
    if raw == "":
        if default is None:
            die(f"case {row.get('case_id')}: column '{key}' is empty and has no default")
        return default
    try:
        return float(raw)
    except ValueError:
        die(f"case {row.get('case_id')}: column '{key}' = '{raw}' is not a number")


def git_describe(path):
    """Record which commit of the thesis repo produced these cases."""
    try:
        out = subprocess.run(
            ["git", "-C", str(path), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=5,
        )
        if out.returncode == 0:
            return out.stdout.strip()
    except Exception:
        pass
    return "unknown"


def is_probably_text(path):
    if path.suffix.lower() in BINARY_SUFFIXES:
        return False
    try:
        with open(path, "rb") as fh:
            chunk = fh.read(4096)
        return b"\0" not in chunk
    except OSError:
        return False


def has_results(case_dir):
    """True if the folder holds run output we must not silently destroy."""
    if not case_dir.exists():
        return False
    for child in case_dir.iterdir():
        if child.is_dir():
            if child.name.startswith("processor"):
                return True
            if child.name in {"postProcessing", "VTK"}:
                return True
            # a numeric directory other than 0 is a written time step
            try:
                if float(child.name) > 0:
                    return True
            except ValueError:
                pass
    return False


# --------------------------------------------------------------------------
# derived quantities
# --------------------------------------------------------------------------

def derive(row):
    """Everything computable from the CSV, so the CSV stays human-editable.

    Keep intent in the spreadsheet; keep arithmetic here. That way a reviewer
    reading caseMatrix.csv sees decisions, not derived clutter.
    """
    d = dict(row)

    # --- wire ---
    d["Rw_m"] = as_float(row, "wire_dia_m") / 2.0

    # --- applied magnetic field, mT -> T, resolved onto axes ---
    b_t = as_float(row, "B_mT") / 1000.0
    direction = (row.get("B_dir") or "none").strip().lower()
    components = {
        "x": (b_t, 0.0, 0.0),
        "y": (0.0, b_t, 0.0),
        "z": (0.0, 0.0, b_t),
        "none": (0.0, 0.0, 0.0),
        "": (0.0, 0.0, 0.0),
    }
    if direction not in components:
        die(f"case {row['case_id']}: B_dir='{direction}' must be one of x, y, z, none")
    d["Bx_T"], d["By_T"], d["Bz_T"] = components[direction]
    d["B_T"] = b_t
    d["B_omega_rad_s"] = 2.0 * 3.141592653589793 * as_float(row, "B_freq_Hz", 0.0)
    d["B_is_alternating"] = 1 if as_float(row, "B_freq_Hz", 0.0) > 0 else 0

    # --- symmetry: a half domain is only legitimate for B=0 or a y-directed field ---
    sym = (row.get("sym") or "half").strip().lower()
    if sym not in {"half", "full"}:
        die(f"case {row['case_id']}: sym='{sym}' must be 'half' or 'full'")
    d["sym"] = sym

    # --- domain box: per-model defaults, overridable by CSV columns ---
    model = row["model"].strip().upper()
    geom = MODEL_GEOM.get(model, MODEL_GEOM["P"])
    pad = geom["pad_x_m"]
    halfw = geom["half_width_m"]
    length = as_float(row, "deposit_len_m")
    cell = as_float(row, "cell_size_m")

    # Model D does not travel, so its box is centred on the wire tip rather
    # than spanning the deposit.
    if geom["span_x"]:
        d["x_min_m"], d["x_max_m"] = -pad, length + pad
    else:
        d["x_min_m"], d["x_max_m"] = -pad, pad
    d["y_min_m"] = 0.0 if sym == "half" else -halfw
    d["y_max_m"] = halfw
    d["z_min_m"] = -geom["substrate_z_m"]
    d["z_max_m"] = geom["air_z_m"]

    for key in GEOM_KEYS:                       # CSV wins if the column exists
        raw = (row.get(key) or "").strip()
        if raw != "":
            d[key] = float(raw)

    for axis, lo, hi in (("nx", "x_min_m", "x_max_m"),
                         ("ny", "y_min_m", "y_max_m"),
                         ("nz", "z_min_m", "z_max_m")):
        raw = (row.get(axis) or "").strip()
        d[axis] = int(raw) if raw != "" else max(1, round((d[hi] - d[lo]) / cell))

    d["n_cells"] = d["nx"] * d["ny"] * d["nz"]
    d["est_ram_gb"] = round(d["n_cells"] / 1e6 * GB_PER_MILLION_CELLS, 2)

    # --- timing: fall back to a single traverse per layer ---
    end_time = as_float(row, "end_time_s", 0.0)
    if end_time <= 0:
        layers = max(1, int(as_float(row, "n_layers", 1)))
        end_time = length / as_float(row, "travel_speed_m_s") * layers
    d["end_time_s"] = round(end_time, 6)

    d["arc_power_W"] = (as_float(row, "eta_arc", 0.7)
                        * as_float(row, "U_V") * as_float(row, "I_A"))
    return d


# --------------------------------------------------------------------------
# writing
# --------------------------------------------------------------------------

def write_case_parameters(case_dir, values):
    """OpenFOAM dictionary that template dicts #include and reference as $name."""
    target = case_dir / "constant"
    target.mkdir(parents=True, exist_ok=True)
    lines = [
        "/*--------------------------------*- C++ -*----------------------------------*\\",
        "  GENERATED BY generate_cases.py -- DO NOT EDIT BY HAND.",
        "  Edit caseMatrix.csv and regenerate, or your change will be lost.",
        "\\*---------------------------------------------------------------------------*/",
        "FoamFile",
        "{",
        "    version     2.0;",
        "    format      ascii;",
        "    class       dictionary;",
        "    object      caseParameters;",
        "}",
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //",
        "",
    ]
    for key in sorted(values):
        val = values[key]
        if isinstance(val, str) and (val.strip() == "" or " " in val):
            continue  # free text (notes) -- not a valid dictionary entry
        lines.append(f"{key:<24} {val};")
    lines += [
        "",
        f"Bext                     ({values['Bx_T']} {values['By_T']} {values['Bz_T']});",
        "",
        "// ************************************************************************* //",
        "",
    ]
    (target / "caseParameters").write_text("\n".join(lines))


def substitute_tree(case_dir, values, strict=True):
    """Replace @token@ everywhere in the copied template. Returns unknown tokens."""
    unknown = set()
    lookup = {k: str(v) for k, v in values.items()}

    def repl(match):
        key = match.group(1)
        if key in lookup:
            return lookup[key]
        unknown.add(key)
        return match.group(0)

    for path in case_dir.rglob("*"):
        if not path.is_file() or not is_probably_text(path):
            continue
        try:
            text = path.read_text()
        except UnicodeDecodeError:
            continue
        if "@" not in text:
            continue
        new = TOKEN_RE.sub(repl, text)
        if new != text:
            path.write_text(new)
    return unknown


def build_case(row, templates_dir, out_root, args, repo_root):
    case_id = row["case_id"].strip()
    model = row["model"].strip().upper()
    if model not in MODEL_DIRS:
        die(f"case {case_id}: model='{model}' must be one of {sorted(MODEL_DIRS)}")

    template = templates_dir / f"model{model}"
    if not template.is_dir():
        die(f"case {case_id}: no template at {template}")

    dest = out_root / MODEL_DIRS[model] / case_id
    values = derive(row)

    warn = ""
    if values["est_ram_gb"] > RAM_WARN_GB:
        warn = f"  <-- ~{values['est_ram_gb']} GB, will not fit locally"

    print(f"  {case_id:<20} {model}  {values['n_cells']:>9,} cells  "
          f"{values['sym']:<4}  t_end={values['end_time_s']:<6}{warn}")

    if args.dry_run:
        return "dry-run"

    if dest.exists():
        if has_results(dest) and not args.force:
            print(f"      SKIPPED -- folder holds run output. Use --force to overwrite.")
            return "skipped"
        shutil.rmtree(dest)

    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(template, dest)

    write_case_parameters(dest, values)
    unknown = substitute_tree(dest, values)
    if unknown:
        print(f"      WARNING unresolved tokens: {', '.join(sorted(unknown))}")

    meta = {
        "case_id": case_id,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "thesis_git_commit": git_describe(repo_root),
        "generator": Path(__file__).name,
        "csv_row": dict(row),
        "derived": {k: v for k, v in values.items() if k not in row},
    }
    (dest / "case_meta.json").write_text(json.dumps(meta, indent=2))
    return "built"


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def main():
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--csv", default=here / "caseMatrix.csv", type=Path)
    p.add_argument("--templates", default=here / "templates", type=Path)
    p.add_argument("--out", default=here, type=Path)
    p.add_argument("--only", default="*", help="glob on case_id, e.g. 'M0*'")
    p.add_argument("--group", default=None, help="only rows with this group value")
    p.add_argument("--dry-run", action="store_true", help="report, change nothing")
    p.add_argument("--force", action="store_true", help="overwrite cases holding results")
    p.add_argument("--all", action="store_true", help="include rows with enabled=0")
    p.add_argument("--list", action="store_true", help="print the matrix and exit")
    args = p.parse_args()

    if not args.csv.is_file():
        die(f"no CSV at {args.csv}")

    with open(args.csv, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        die("CSV has no rows")

    missing = REQUIRED_COLUMNS - set(rows[0])
    if missing:
        die(f"CSV is missing required columns: {', '.join(sorted(missing))}")

    seen = set()
    for r in rows:
        cid = r["case_id"].strip()
        if cid in seen:
            die(f"duplicate case_id '{cid}'")
        seen.add(cid)

    if args.list:
        print(f"{'case_id':<20} {'M':<2} {'grp':<12} {'B':<10} {'sym':<5} notes")
        print("-" * 100)
        for r in rows:
            b = "-" if r["B_dir"] in ("none", "") else f"{r['B_mT']}mT {r['B_dir']}"
            if float(r.get("B_freq_Hz") or 0) > 0:
                b += f" @{r['B_freq_Hz']}Hz"
            flag = " " if truthy(r["enabled"]) else "x"
            print(f"{flag}{r['case_id']:<19} {r['model']:<2} {r.get('group',''):<12} "
                  f"{b:<10} {r['sym']:<5} {r.get('notes','')[:40]}")
        print(f"\n{sum(1 for r in rows if truthy(r['enabled']))} of {len(rows)} enabled "
              f"('x' marks disabled)")
        return

    selected = [
        r for r in rows
        if (args.all or truthy(r["enabled"]))
        and fnmatch.fnmatch(r["case_id"].strip(), args.only)
        and (args.group is None or r.get("group", "").strip() == args.group)
    ]
    if not selected:
        die("no rows selected -- check --only / --group / the enabled column")

    print(f"\n{'DRY RUN -- ' if args.dry_run else ''}generating {len(selected)} case(s) "
          f"into {args.out}")
    print("cell counts assume a UNIFORM mesh -- they are an upper bound and will "
          "drop once\nblockMeshDict grading and AMR are added.\n")
    tally = {}
    for row in selected:
        result = build_case(row, args.templates, args.out, args, here.parent)
        tally[result] = tally.get(result, 0) + 1

    print("\n" + "  ".join(f"{k}: {v}" for k, v in sorted(tally.items())))
    total_cells = sum(derive(r)["n_cells"] for r in selected)
    print(f"total cells across selection: {total_cells:,}")
    if not args.dry_run:
        print("\nNext: cd into a case, run blockMesh, then the solver.")


if __name__ == "__main__":
    main()
