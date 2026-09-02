#!/usr/bin/env python3
"""Make sigma_e_gas_ratio live, re-centre its sweep, and state the arc origin.

THE BUG THIS FIXES
    caseMatrix.csv carries a sigma_e_gas_ratio column, and S00/S01 vary it to
    study the conductivity floor. But transportProperties is copied VERBATIM,
    so nothing reads that column: S00 and S01 would reproduce the baseline
    exactly and you would conclude "no sensitivity to the floor" from a study
    that never ran. A silent false negative.

THREE CHANGES
    1. transportProperties gas Elec_conductivity becomes @sigma_e_gas@, filled
       by substitute_tree (which runs at line 552, right after the file is
       installed at 551).
    2. derive() computes sigma_e_gas = metal_sigma * ratio, reading the metal
       value FROM the property set so the two can never disagree. The metal
       conductivity stays hardcoded in transportProperties where its full
       provenance is documented; it is not duplicated into the CSV.
    3. The sweep is re-centred. It bracketed 1e-3; Gate B measured the floor
       at 1e-6, so 1e-2/1e-4 both sit ABOVE the new baseline and probe one
       direction only. Becomes 1e-4 / 1e-8, two orders each way.

    Also states the arc origin explicitly, which Phase 3's heat source needs.

USAGE
    python3 patch_sigma_gas.py            # apply
    python3 patch_sigma_gas.py --revert
"""
import csv, shutil, sys
from pathlib import Path

GEN = Path.home() / "thesis/cases/generate_cases.py"
CSV = Path.home() / "thesis/cases/caseMatrix.csv"
MT  = Path.home() / "thesis/materials/make_transport.py"

BASELINE = "1e-6"
SWEEP = {"S00-sigmafloor-1e2": "1e-4", "S01-sigmafloor-1e4": "1e-8"}
RENAME = {"S00-sigmafloor-1e2": "S00-sigmafloor-1e4",
          "S01-sigmafloor-1e4": "S01-sigmafloor-1e8"}

# ---------------------------------------------------------------- generator
GEN_HELPER_ANCHOR = "def derive(row):"
GEN_HELPER = '''def metal_sigma_e():
    """The metal's electrical conductivity, read from the property set itself.

    Deliberately NOT duplicated into caseMatrix.csv. The value carries a long
    provenance comment in transportProperties (Desai + Matthiessen offset,
    cross-checked to 1% against a measured alloy); copying the bare number into
    a spreadsheet would separate it from its justification and let the two
    drift. The gas floor is then always a ratio of the number the metal uses.
    """
    for line in PROPERTIES_SRC.read_text().splitlines():
        s = line.split("//")[0].strip()
        if s.startswith("Elec_conductivity"):      # metal block precedes gas
            return float(s.split()[1].rstrip(";"))
    die(f"no Elec_conductivity found in {PROPERTIES_SRC}")


'''

GEN_DERIVE_ANCHOR = '    d["Rw_m"] = as_float(row, "wire_dia_m") / 2.0'
GEN_DERIVE_ADD = '''
    # --- electrical conductivity of the gas phase ---
    # The CSV holds a RATIO; the absolute value is derived from the metal's, so
    # S00/S01 actually vary what they claim to. Fills @sigma_e_gas@ in
    # transportProperties via substitute_tree.
    d["sigma_e_gas"] = metal_sigma_e() * as_float(row, "sigma_e_gas_ratio")
    # --- arc start point, for the Goldak heat source (Phase 3) ---
    # The domain is built with x_min = -pad_x and z_min = -substrate_z, so the
    # deposit corner sits at the ORIGIN by construction. Stated explicitly here
    # so the heat source reads a value instead of re-deriving mesh convention.
    d["arc_x0_m"] = 0.0
    d["arc_y0_m"] = 0.0          # symmetry plane
    d["arc_z0_m"] = 0.0          # substrate top surface'''

# ---------------------------------------------------------------- properties
MT_OLD = '''    // [MEAS-ish] the tutorial's own floor, ratio 1e-6 to the metal. Held
    // div(muH) at 1e-12 across the jump in Gate B. Do not lower without
    // re-checking divergence.
    Elec_conductivity   3.65;'''
MT_NEW = '''    // Filled per case by generate_cases.py from sigma_e_gas_ratio in
    // caseMatrix.csv, as ratio * the metal value above. Baseline ratio 1e-6 --
    // the MTHD tutorial's own floor, which held div(muH) at 1e-12 across the
    // jump in Gate B. S00/S01 sweep it two orders either side.
    // Cold argon is effectively an insulator; this floor exists only to stop
    // the induction equation going singular as sigma -> 0, so smaller is more
    // physical until numerics complain. That is what the sweep tests.
    Elec_conductivity   @sigma_e_gas@;'''


def revert():
    n = 0
    for p in (GEN, CSV, MT):
        bak = p.with_suffix(p.suffix + ".sigbak")
        if bak.exists():
            shutil.copy2(bak, p); bak.unlink(); n += 1
            print(f"  reverted {p.name}")
    print(f"  {n} file(s) restored" if n else "  nothing to revert")


def apply():
    gen, mt = GEN.read_text(), MT.read_text()
    for name, text, anchors in (("generate_cases.py", gen,
                                 [GEN_HELPER_ANCHOR, GEN_DERIVE_ANCHOR]),
                                ("make_transport.py", mt, [MT_OLD])):
        for a in anchors:
            if text.count(a) != 1:
                sys.exit(f"ERROR: {name}: expected one match for:\n{a}\n"
                         f"found {text.count(a)}. Nothing modified.")
    if "sigma_e_gas" in gen:
        sys.exit("ERROR: generate_cases.py already patched. --revert first.")

    rows = list(csv.DictReader(CSV.open()))
    if "sigma_e_gas_ratio" not in rows[0]:
        sys.exit("ERROR: no sigma_e_gas_ratio column")

    for p in (GEN, CSV, MT):
        shutil.copy2(p, p.with_suffix(p.suffix + ".sigbak"))

    gen = gen.replace(GEN_HELPER_ANCHOR, GEN_HELPER + GEN_HELPER_ANCHOR, 1)
    gen = gen.replace(GEN_DERIVE_ANCHOR, GEN_DERIVE_ANCHOR + GEN_DERIVE_ADD, 1)
    GEN.write_text(gen)
    MT.write_text(mt.replace(MT_OLD, MT_NEW, 1))

    changed = []
    for r in rows:
        cid = r["case_id"].strip()
        old = r["sigma_e_gas_ratio"]
        new = SWEEP.get(cid, BASELINE)
        if old != new:
            r["sigma_e_gas_ratio"] = new
            changed.append((cid, old, new))
        if cid in RENAME:
            r["case_id"] = RENAME[cid]
            r["notes"] = r["notes"].replace(
                f"Conductivity floor sensitivity {old}",
                f"Conductivity floor sensitivity {new}")
    # twin references must follow a rename
    for r in rows:
        t = (r.get("noise_twin_of") or "").strip()
        if t in RENAME:
            r["noise_twin_of"] = RENAME[t]

    with CSV.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys())
        w.writeheader(); w.writerows(rows)

    print(f"  patched generate_cases.py, make_transport.py, caseMatrix.csv")
    print(f"  backups: *.sigbak\n")
    print(f"  {len(changed)} row(s) changed:")
    for cid, o, n in changed[:6]:
        print(f"    {cid:<24} {o} -> {n}")
    if len(changed) > 6:
        print(f"    ... and {len(changed)-6} more at the {BASELINE} baseline")
    print(f"\n  renamed: {', '.join(f'{k} -> {v}' for k, v in RENAME.items())}")
    print("\n  Next:")
    print("    python3 ~/thesis/materials/make_transport.py   # regenerate with the token")
    print("    cd ~/thesis/cases && python3 generate_cases.py --only R00-SPSL-paper --force")


if __name__ == "__main__":
    revert() if "--revert" in sys.argv else apply()
