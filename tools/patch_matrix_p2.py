#!/usr/bin/env python3
"""Revise caseMatrix.csv against what Phase 2 actually measured.

The matrix was built in Phase 0 around an estimated "20-30 mT crossover".
Phase 2 replaced that estimate with six measured thresholds and one finding
that changes what the magnetic group is even for. This brings the matrix into
line. All changes are CSV-only -- no new columns, nothing downstream to rewire.

  Marangoni ......... 0.48 mT     gravity ............ 2.4 mT
  Ha = 1 ............ 6.2 mT      arc pressure ....... 7.1 mT
  arc SELF-field ... 13.5 mT      N = 1 (bulk brake) . 269 mT

SIX CHANGES

1. M01-M12 DISABLED pending Phase 6.
   Current-driven stirring beats motional braking by 200-500x at every usable
   field strength, and the shipped solver has ONLY the braking term. Running
   these before current injection spends twelve cluster runs confirming a null
   already predicted on paper.

2. NEW M-null-003 at 3 mT -- a negative control.
   Below Ha = 1, so the model MUST show nothing. Every other magnetic case
   expects an effect, which leaves nothing to distinguish a real result from a
   bug. A case that must come out empty is the cheapest such check there is.

3. TWO MORE noise twins on M00.
   Phase 1: "one twin pair is not a noise floor, it is one sample of a
   distribution." Headline numbers need three to five perturbed runs and a
   standard deviation. Every comparison here had exactly one.

4. R00-SPSL-paper-noise ENABLED.
   The audit has warned on every run that paperfidelity has no enabled twin.
   Phase 5 needs it, and enabling it later means re-running R00 as well.

5. Magnetic notes REWRITTEN against the measured thresholds, and against the
   mesh. At 350 um cells you need B < 11 mT for five cells across a Hartmann
   layer, and below 6.2 mT there is no layer. EVERY magnetic case sits in the
   unresolved window. Survivable -- post-Phase-6 the effect is bulk stirring,
   linear in B, needing no layer resolved -- but it has to be recorded, because
   "did you resolve the Hartmann layer?" is the first question anyone will ask.

6. T-model solver column: additiveFoam -> TBD-thermal.
   Those rows name a solver you do not have and have not chosen. The route is
   a Phase 8 decision; the column should not pretend otherwise meanwhile.

USAGE
    python3 patch_matrix_p2.py           # apply
    python3 patch_matrix_p2.py --revert
"""
import csv, shutil, sys
from pathlib import Path

CSV = Path.home() / "thesis/cases/caseMatrix.csv"
BAK = CSV.with_suffix(".csv.p2bak")

# --- 5. notes carrying the Phase 2 numbers ---------------------------------
NOTES = {
 "M00-B000":
   "Reduced-mesh baseline at B=0. Reference for the whole magnetic group.",
 "M-null-003":
   "NEGATIVE CONTROL. 3mT sits below Ha=1 (6.2mT) so NO effect is the correct "
   "result. If this one shows an effect something is wrong.",
 "M01-By-010":
   "Transverse-Y 10mT = 0.74x the arc self-field (13.5mT). Ha=1.6 and the "
   "Hartmann layer is 1.85mm ~ 5 cells: the only marginally resolved case.",
 "M02-By-030":
   "Transverse-Y 30mT = 2.2x self-field. Ha=4.9 and the layer is 617um < 2 "
   "cells: UNRESOLVED. Bulk stirring only. Previously mis-noted as near the "
   "Marangoni crossover which is actually 0.48mT.",
 "M03-By-050":
   "Transverse-Y 50mT = 3.7x self-field. Ha=8.1 and the layer is 370um ~ 1 "
   "cell: unresolved. Damps u_x and u_z ie the main Marangoni loop.",
 "M04-By-100":
   "Transverse-Y 100mT = 7.4x self-field. Ha=16.2 and the layer is 185um < 1 "
   "cell: badly unresolved. Keep for the stirring trend not for layers.",
 "M05-Bx-030":
   "Longitudinal-X 30mT. J x B acts along y which DESTROYS the XZ symmetry so "
   "a full domain is mandatory. Damps the transverse circulation.",
 "M06-Bx-050":
   "Longitudinal-X 50mT = 3.7x self-field. Symmetry broken so full domain.",
 "M07-Bz-030":
   "Axial-Z 30mT. Axial current gives zero force but the arc current also "
   "spreads radially and J_r x B_z is AZIMUTHAL giving swirl. Symmetry broken.",
 "M08-Bz-050":
   "Axial-Z 50mT. Azimuthal swirl about the arc axis. Full domain.",
 "M09-Bz-100":
   "Axial-Z 100mT = 7.4x self-field. Strongest swirl case. Full domain.",
 "M10-Bz050-f010":
   "Alternating Bz 50mT at 10Hz. Skin depth 83mm >> 3mm pool so the field "
   "fully penetrates and quasi-static holds. Needs a time-varying H boundary "
   "condition via codedFixedValue - no solver work.",
 "M11-Bz050-f050":
   "Alternating Bz 50mT at 50Hz. Skin depth 37mm >> pool: full penetration.",
 "M12-Bz050-f070":
   "Alternating Bz 50mT at 70Hz. Skin depth 31mm >> pool: full penetration. "
   "Highest frequency in the sweep and still fully penetrating.",
}

P6 = "[PHASE 6] "        # prefix marking cases blocked on current injection


def revert():
    if not BAK.exists():
        print("  nothing to revert"); return
    shutil.copy2(BAK, CSV); BAK.unlink(); print(f"  reverted {CSV.name}")


def apply():
    rows = list(csv.DictReader(CSV.open()))
    cols = list(rows[0].keys())
    by_id = {r["case_id"].strip(): r for r in rows}

    if "M-null-003" in by_id:
        sys.exit("ERROR: already patched (M-null-003 exists). --revert first.")
    for need in ("M00-B000", "M00-B000-noise", "M01-By-010",
                 "R00-SPSL-paper-noise"):
        if need not in by_id:
            sys.exit(f"ERROR: expected row '{need}' not found. Nothing changed.")

    shutil.copy2(CSV, BAK)
    log = []

    # 1. disable the magnetic sweep until current injection exists
    for r in rows:
        cid = r["case_id"].strip()
        if cid.startswith("M") and cid != "M00-B000" and "noise" not in cid:
            if r["enabled"].strip() == "1":
                r["enabled"] = "0"
                log.append(f"    disabled  {cid}")

    # 5. rewrite notes, tagging what is blocked on Phase 6
    for cid, note in NOTES.items():
        if cid in by_id:
            blocked = cid.startswith("M") and cid not in ("M00-B000",)
            by_id[cid]["notes"] = (P6 if blocked else "") + note

    # 2. negative control, cloned from M01 so every other column matches
    null = dict(by_id["M01-By-010"])
    null.update({"case_id": "M-null-003", "enabled": "0", "B_mT": "3",
                 "B_dir": "y", "notes": P6 + NOTES["M-null-003"]})
    rows.insert(rows.index(by_id["M01-By-010"]), null)
    log.append("    ADDED     M-null-003  (3 mT, below Ha=1, must show nothing)")

    # 3. two more twins so M00 has a distribution, not a single sample
    base = by_id["M00-B000-noise"]
    for tag, tol in (("2", "1e-09"), ("3", "5e-08")):
        t = dict(base)
        t.update({"case_id": f"M00-B000-noise{tag}", "p_rgh_tol": tol,
                  "notes": (f"Noise-floor twin {tag} of M00-B000. Identical "
                            f"physics with p_rgh tolerance {tol}. Phase 1 "
                            f"lesson: one pair is one SAMPLE not a floor - "
                            f"3 to 5 are needed for a standard deviation.")})
        rows.insert(rows.index(base) + 1, t)
        log.append(f"    ADDED     M00-B000-noise{tag}  (p_rgh_tol {tol})")

    # 4. paperfidelity needs its floor before Phase 5, not during
    by_id["R00-SPSL-paper-noise"]["enabled"] = "1"
    log.append("    ENABLED   R00-SPSL-paper-noise")

    # 6. stop naming a solver that has not been chosen
    for r in rows:
        if r["solver"].strip() == "additiveFoam":
            r["solver"] = "TBD-thermal"
            log.append(f"    solver    {r['case_id'].strip()} -> TBD-thermal")

    with CSV.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader(); w.writerows(rows)

    print(f"  patched {CSV.name}   (backup: {BAK.name})\n")
    for line in log:
        print(line)
    live = sum(1 for r in rows if r["enabled"].strip() == "1")
    print(f"\n  {len(rows)} rows, {live} enabled (was 17)")
    print("""
  Still outstanding, needing NEW COLUMNS rather than CSV edits:
    * surface tension sensitivity -- Phase 2 called the pure-Al 0.85 N/m the
      single largest uncertainty in the replication, and there is no case
    * wall thermal BC -- heat reaches the substrate underside within the run,
      so fixed-300 vs adiabatic is a real and currently untested assumption

  Next:
    python3 generate_cases.py --list
    python3 generate_cases.py --dry-run --all""")


if __name__ == "__main__":
    revert() if "--revert" in sys.argv else apply()
