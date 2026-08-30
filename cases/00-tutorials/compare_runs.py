#!/usr/bin/env python3
"""
compare_runs.py -- compare an OpenFOAM field between two or more case runs.

Purpose: distinguish "this solution is chaotic" from "the parallel path is
broken". Chaotic cases amplify round-off, so a tiny SERIAL perturbation
produces differences of the same order as serial-vs-parallel. If a serial
perturbation stays at round-off while parallel explodes, the parallel path
has a real bug.

Usage:
    ./compare_runs.py T 0.001 0.01 -- P2D-ser P2D-ser2 P2D-par

    ./compare_runs.py <field> <time> [<time> ...] -- <caseA> <caseB> [...]

The FIRST case listed is the reference; every other case is compared to it.

Handles ascii uniform and nonuniform scalar fields, gzipped or not.
Python 3.8+, standard library only.
"""

import gzip
import re
import sys
from pathlib import Path

NONUNIFORM = re.compile(
    r"internalField\s+nonuniform\s+List<scalar>\s*\n?\s*(\d+)\s*\n?\s*\((.*?)\)\s*;",
    re.S,
)
UNIFORM = re.compile(r"internalField\s+uniform\s+([-\d.eE+]+)\s*;")


def read_field(case, time, field):
    """Return a list of cell values, or None with a reason."""
    base = Path(case) / time
    for cand in (base / field, base / (field + ".gz")):
        if not cand.is_file():
            continue
        opener = gzip.open if cand.suffix == ".gz" else open
        with opener(cand, "rt", errors="replace") as fh:
            txt = fh.read()
        m = NONUNIFORM.search(txt)
        if m:
            vals = [float(v) for v in m.group(2).split()]
            n = int(m.group(1))
            if len(vals) != n:
                return None, f"{cand}: header says {n} values, parsed {len(vals)}"
            return vals, None
        u = UNIFORM.search(txt)
        if u:
            return [float(u.group(1))], None
        return None, f"{cand}: could not parse (binary writeFormat?)"
    return None, f"{base}/{field} not found"


def compare(ref, other):
    """Return (max_abs, rel, n_over_1, worst) for two equal-length lists."""
    pairs = list(zip(ref, other))
    diffs = [abs(a - b) for a, b in pairs]
    max_abs = max(diffs)
    span = max(ref) - min(ref)
    rel = max_abs / span if span > 0 else 0.0
    n_over_1 = sum(1 for d in diffs if d > 1.0)
    worst = diffs.index(max_abs)
    return max_abs, rel, n_over_1, worst


def verdict(rel):
    if rel < 1e-9:
        return "identical to round-off"
    if rel < 1e-6:
        return "round-off level"
    if rel < 1e-3:
        return "small - sensitive but consistent"
    return "LARGE - chaotic amplification or a bug"


def main():
    argv = sys.argv[1:]
    if "--" not in argv or len(argv) < 4:
        print(__doc__)
        sys.exit(1)
    split = argv.index("--")
    field = argv[0]
    times = argv[1:split]
    cases = argv[split + 1:]
    if len(cases) < 2:
        print("need at least two cases")
        sys.exit(1)

    ref_case = cases[0]
    print(f"\nfield: {field}    reference case: {ref_case}\n")
    header = f"{'time':<10}{'case':<14}{'n':<8}{'max|d|':>12}{'relative':>12}  {'>1K':>6}  verdict"
    print(header)
    print("-" * len(header))

    for t in times:
        ref, err = read_field(ref_case, t, field)
        if err:
            print(f"{t:<10}{ref_case:<14}  SKIP: {err}")
            continue
        print(f"{t:<10}{ref_case:<14}{len(ref):<8}{'--':>12}{'(reference)':>12}")
        for c in cases[1:]:
            other, err = read_field(c, t, field)
            if err:
                print(f"{'':<10}{c:<14}  SKIP: {err}")
                continue
            if len(other) != len(ref):
                print(f"{'':<10}{c:<14}  SKIP: {len(other)} cells vs {len(ref)}")
                continue
            ma, rel, n1, worst = compare(ref, other)
            print(f"{'':<10}{c:<14}{len(other):<8}{ma:12.4e}{rel:12.4e}  {n1:>6}  {verdict(rel)}")
        print()

    print("How to read this:")
    print("  If the PERTURBED SERIAL case differs as much as the PARALLEL case,")
    print("  the solution is chaotic and parallel is fine.")
    print("  If perturbed serial stays at round-off while parallel is large,")
    print("  the parallel path has a genuine bug.\n")


if __name__ == "__main__":
    main()
