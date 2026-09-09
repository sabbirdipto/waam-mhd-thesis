#!/usr/bin/env python3
"""
hartmann_check.py -- Gate A. Compare a computed Hartmann profile with theory.

Hartmann flow: fully developed flow between parallel plates at y = +/-L with a
transverse magnetic field B. The Lorentz force flattens the parabolic
Poiseuille profile into a plug with thin boundary layers, and the analytical
solution is exact:

    u(y) / u(0) = [cosh(Ha) - cosh(Ha*y/L)] / [cosh(Ha) - 1]

    Ha = B * L * sqrt(sigma / mu)          (Hartmann number)

Ha -> 0 recovers the parabola; large Ha gives a flat core. Because we normalise
by the centreline value, the pressure gradient cancels and the comparison tests
the MHD coupling alone -- which is the point.

Usage
-----
    # in the case directory, after the run
    postProcess -latestTime -func sample

    ./hartmann_check.py postProcessing/sample/<time>/line_centreProfile_U_H.csv

Options let you override what is read from constant/transportProperties:
    --sigma  electrical conductivity of the fluid   (default 1)
    --mu     dynamic viscosity, rho*nu              (default 1)
    --B      applied flux density, mu_mag * H       (default: from the data)

Python 3.8+, standard library only. Optional matplotlib for a plot.
"""

import argparse
import csv
import math
import sys
from pathlib import Path


def read_profile(path):
    """Return (y, ux) sorted by y, from an OpenFOAM 'sets' csv."""
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        sys.exit(f"{path}: no rows")

    cols = {k.strip().lower(): k for k in rows[0]}

    def pick(*cands):
        for c in cands:
            if c in cols:
                return cols[c]
        return None

    ycol = pick("y")
    uxcol = pick("u_x", "ux", "u_0")
    if ycol is None or uxcol is None:
        sys.exit(f"{path}: need a 'y' and a 'U_x' column; found {list(rows[0])}")

    hy = pick("h_y", "hy", "h_1")
    data = []
    for r in rows:
        try:
            data.append((float(r[ycol]), float(r[uxcol]),
                         float(r[hy]) if hy else None))
        except (TypeError, ValueError):
            continue
    data.sort()
    if len(data) < 5:
        sys.exit(f"{path}: only {len(data)} usable rows")
    return data


def analytic(y, L, Ha):
    """Normalised Hartmann profile, u(y)/u(0)."""
    if Ha < 1e-9:                       # Poiseuille limit
        return 1.0 - (y / L) ** 2
    return (math.cosh(Ha) - math.cosh(Ha * y / L)) / (math.cosh(Ha) - 1.0)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("csv", type=Path)
    p.add_argument("--sigma", type=float, default=None,
                   help="electrical conductivity. Default: read from the case")
    p.add_argument("--mu", type=float, default=None,
                   help="dynamic viscosity rho*nu. Default: read from the case")
    p.add_argument("--B", type=float, default=None,
                   help="applied flux density. Default: computed as mu_mag*H from "
                        "the sampled H column -- far stricter than fitting")
    p.add_argument("--mu-mag", type=float, default=None, dest="mumag",
                   help="magnetic permeability, to turn sampled H into B. "
                        "Default: read from the case")
    p.add_argument("--L", type=float, default=None,
                   help="channel half-width (wall position). Default: inferred as "
                        "half-span of cell centres PLUS half a cell, which is the "
                        "true wall on a uniform mesh")
    p.add_argument("--plot", action="store_true")
    args = p.parse_args()

    # ------------------------------------------------------------------
    # Resolve material properties FROM THE CASE. Assuming unity turned a
    # passing benchmark (1.3869e-02) into a 63% FAIL, because hartmann1d
    # hides elec_conductivity_air_ = 1000 behind an #include.
    # ------------------------------------------------------------------
    def _scalars(path):
        """Every 'key value;' pair in a foam dict. Comments stripped."""
        out = {}
        if not path.exists():
            return out
        for raw in path.read_text().splitlines():
            line = raw.split("//")[0].strip()
            if not line.endswith(";"):
                continue
            bits = line[:-1].split()
            if len(bits) == 2:
                try:
                    out[bits[0]] = float(bits[1])
                except ValueError:
                    pass
        return out

    def _all_values(path, key):
        """Every value assigned to `key`, so disagreement can be detected."""
        vals = []
        if not path.exists():
            return vals
        for raw in path.read_text().splitlines():
            line = raw.split("//")[0].strip()
            if not line.endswith(";"):
                continue
            bits = line[:-1].split()
            if len(bits) == 2 and bits[0] == key:
                try:
                    vals.append(float(bits[1]))
                except ValueError:
                    pass
        return vals

    case = args.csv.resolve().parent
    for _ in range(6):
        if (case / "constant").is_dir():
            break
        case = case.parent
    else:
        case = None

    src = {}
    if case is not None:
        ref = _scalars(case / "constant" / "propertiesReference")
        tp = case / "constant" / "transportProperties"

        if args.sigma is None and "elec_conductivity_air_" in ref:
            args.sigma = ref["elec_conductivity_air_"]
            src["sigma"] = "constant/propertiesReference"
        if args.mumag is None and "magnetic_permeability_air_" in ref:
            args.mumag = ref["magnetic_permeability_air_"]
            src["mu_mag"] = "constant/propertiesReference"

        if args.mu is None:
            nus, rhos = _all_values(tp, "nu"), _all_values(tp, "rho")
            if nus and rhos and len(set(nus)) == 1 and len(set(rhos)) == 1:
                args.mu = rhos[0] * nus[0]
                src["mu"] = "constant/transportProperties (rho*nu)"
            elif nus and rhos:
                sys.exit(
                    "ERROR: the phases in constant/transportProperties disagree "
                    f"(nu={sorted(set(nus))}, rho={sorted(set(rhos))}).\n"
                    "  Which phase is flowing is not something this script will "
                    "guess -- in hartmann1d 'metal' is marked \"Unused Phase\".\n"
                    "  Pass --mu explicitly.")

    missing = [n for n, v in (("--sigma", args.sigma), ("--mu", args.mu)) if v is None]
    if missing:
        sys.exit(
            f"ERROR: could not determine {', '.join(missing)} from the case"
            + (f" at {case}" if case else " (no constant/ directory found)")
            + ".\n  Refusing to assume 1: doing so scores a passing Hartmann "
              "benchmark\n  as a 63% FAIL. Pass the value(s) explicitly.")
    if args.mumag is None:
        args.mumag = 1.0
        src["mu_mag"] = "assumed 1 (not found in the case)"

    for name, val in (("sigma", args.sigma), ("mu", args.mu),
                      ("mu_mag", args.mumag)):
        print(f"  {name:<16} : {val:<12g} [{src.get(name, 'given on command line')}]")


    data = read_profile(args.csv)
    ys = [d[0] for d in data]
    us = [d[1] for d in data]
    # The sampled points are CELL CENTRES, so the outermost sample sits half a
    # cell inside the wall. Using half the sampled span as L puts the analytical
    # zero at the first cell centre and manufactures a huge spurious error there
    # while the core still matches -- a very recognisable signature.
    y0 = (max(ys) + min(ys)) / 2.0      # recentre in case the line is offset
    span = (max(ys) - min(ys)) / 2.0
    dy = (max(ys) - min(ys)) / (len(ys) - 1) if len(ys) > 1 else 0.0
    L = args.L if args.L is not None else span + dy / 2.0
    Lhow = "given" if args.L is not None else f"half-span {span:.4f} + half-cell {dy/2:.4f}"

    # centreline value: interpolate at y = y0
    mid = min(range(len(ys)), key=lambda i: abs(ys[i] - y0))
    u0 = us[mid]
    if abs(u0) < 1e-30:
        sys.exit("centreline velocity is zero -- did the case converge?")

    # Ha must be IMPOSED, not fitted. A fitted Ha absorbs a wrong profile by
    # sliding along the Hartmann family and turns the test into a tautology.
    fitted = False
    if args.B is not None:
        B = args.B
        how = f"B = {B:g} from --B"
    else:
        hy = [d[2] for d in data if d[2] is not None]
        if hy:
            B = args.mumag * (sum(hy) / len(hy))
            how = f"B = mu_mag*<H_y> = {B:.6g} from the sampled H column"
        else:
            B = None
            how = "FITTED -- weak test, see warning"
    if B is not None:
        Ha = B * L * math.sqrt(args.sigma / args.mu)
    else:
        fitted = True
        best, Ha = None, 0.0
        h = 0.0
        while h <= 200.0:
            err = sum((us[i] / u0 - analytic(ys[i] - y0, L, h)) ** 2
                      for i in range(len(ys)))
            if best is None or err < best:
                best, Ha = err, h
            h += 0.05

    dev, worst = 0.0, 0.0
    for y, u, _ in data:
        a = analytic(y - y0, L, Ha)
        d = abs(u / u0 - a)
        if d > dev:
            dev, worst = d, y

    print(f"\n  samples          : {len(data)}")
    print(f"  half-width L     : {L:.4f}   ({Lhow})")
    print(f"  centreline u(0)  : {u0:.6g}")
    print(f"  Hartmann number  : {Ha:.3f}   ({how})")
    print(f"  max |u/u0 - theory| : {dev:.4e}   at y = {worst:.4f}")

    if dev < 0.02:
        v = "PASS  - profile matches theory to better than 2%"
    elif dev < 0.05:
        v = "MARGINAL - within 5%; check mesh resolution in the Hartmann layer"
    else:
        v = "*** FAIL - investigate before building on this module ***"
    print(f"  VERDICT          : {v}")
    if fitted:
        print("  WARNING          : Ha was FITTED, not imposed. A fitted Ha slides")
        print("                     along the Hartmann family and can absorb a wrong")
        print("                     profile. Pass --B, or sample H, for a real test.")
    print()

    print(f"  {'y':>10} {'computed':>12} {'theory':>12} {'diff':>12}")
    print("  " + "-" * 48)
    step = max(1, len(data) // 15)
    for i in range(0, len(data), step):
        y, u, _ = data[i]
        a = analytic(y - y0, L, Ha)
        print(f"  {y:10.4f} {u/u0:12.6f} {a:12.6f} {u/u0-a:12.3e}")
    print()

    if args.plot:
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
        except ImportError:
            print("  (matplotlib not installed -- skipping plot)")
            return
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.plot([u / u0 for u in us], ys, "o", ms=3, label="computed")
        fine = [min(ys) + i * (max(ys) - min(ys)) / 400 for i in range(401)]
        ax.plot([analytic(y - y0, L, Ha) for y in fine], fine, "-",
                label=f"theory, Ha={Ha:.2f}")
        ax.set_xlabel("$u_x / u_x(0)$")
        ax.set_ylabel("y")
        ax.legend()
        ax.grid(alpha=.3)
        fig.tight_layout()
        out = args.csv.parent / "hartmann_check.png"
        fig.savefig(out, dpi=150)
        print(f"  plot written to {out}\n")


if __name__ == "__main__":
    main()
