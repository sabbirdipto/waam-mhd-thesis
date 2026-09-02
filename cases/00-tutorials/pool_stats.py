#!/usr/bin/env python3
"""Integral flow statistics from OpenFOAM time directories.

In a chaotic solution a single-cell extremum is not a basis for comparing two
cases. This reports several measures so you can pick the one with the tightest
noise floor for YOUR case.

Measured on a Plate2D MHD/no-MHD/noise-twin triple, the floors were:

    melt volume (cells with alpha>0.5)   0.5 %     <- tightest
    mean |U| in metal                    2.6 %     <- best velocity measure
    max  |U| in metal                    7.5 %
    alpha-weighted mean |U|              9.9 %
    max  |U| whole domain               25   %     <- useless

alpha-weighting is WORSE than restricting to alpha>0.5, because it includes the
partial interface cells -- and the VOF interface is where round-off amplifies
fastest. Excluding the interface is what makes a measure robust.

Time-averaging over the quasi-steady window (--last N) reduces the floor
further at no compute cost, since the time directories already exist.

Usage:
    ./pool_stats.py caseA caseB ...            # latest time only
    ./pool_stats.py --last 5 caseA caseB ...   # average over the last 5 writes
"""
import gzip, re, sys, math
from pathlib import Path

VEC = re.compile(r"internalField\s+nonuniform\s+List<vector>\s*\n?\s*(\d+)\s*\n?\s*\((.*?)\)\s*;", re.S)
SCA = re.compile(r"internalField\s+nonuniform\s+List<scalar>\s*\n?\s*(\d+)\s*\n?\s*\((.*?)\)\s*;", re.S)
UNI = re.compile(r"internalField\s+uniform\s+([-\d.eE+]+)\s*;")

def read(p):
    for c in (p, Path(str(p) + ".gz")):
        if c.is_file():
            op = gzip.open if c.suffix == ".gz" else open
            with op(c, "rt", errors="replace") as fh:
                return fh.read()
    return None

def vectors(txt):
    m = VEC.search(txt)
    if not m: return None
    v = [float(x) for x in m.group(2).replace("(", " ").replace(")", " ").split()]
    return [(v[i], v[i+1], v[i+2]) for i in range(0, len(v) - 2, 3)]

def scalars(txt, n):
    m = SCA.search(txt)
    if m: return [float(x) for x in m.group(2).split()]
    u = UNI.search(txt)
    return [float(u.group(1))] * n if u else None

args = sys.argv[1:]
nlast = 1
series = False
while args and args[0].startswith("--"):
    if args[0] == "--last":
        nlast = int(args[1]); args = args[2:]
    elif args[0] == "--series":
        series = True; args = args[1:]
    else:
        sys.exit(f"unknown option {args[0]}")

if series:
    # one row per (case, time): paste or send this straight to a plotting tool
    print("case,time,melt_cells,mean_U_metal,max_U_metal,alpha_wtd_mean,max_U_all")
    for case in args:
        c = Path(case)
        ts = sorted((d for d in c.iterdir() if d.is_dir()
                     and re.fullmatch(r"[0-9.eE+-]+", d.name) and float(d.name) > 0),
                    key=lambda d: float(d.name))
        for t in ts:
            ut, at = read(t / "U"), read(t / "alpha.metal")
            if ut is None: continue
            U = vectors(ut)
            if U is None: continue
            a = scalars(at, len(U)) if at else [1.0] * len(U)
            mag = [math.sqrt(x*x+y*y+z*z) for x, y, z in U]
            metal = [m for m, av in zip(mag, a) if av > 0.5]
            ws = sum(a)
            print(f"{c.name},{t.name},{len(metal)},"
                  f"{(sum(metal)/len(metal) if metal else 0):.6g},"
                  f"{(max(metal) if metal else 0):.6g},"
                  f"{(sum(m*av for m, av in zip(mag, a))/ws if ws else 0):.6g},"
                  f"{max(mag):.6g}")
    sys.exit(0)

for case in args:
    c = Path(case)
    times = sorted((d for d in c.iterdir() if d.is_dir() and re.fullmatch(r"[0-9.eE+-]+", d.name)
                    and float(d.name) > 0), key=lambda d: float(d.name))
    if not times:
        print(f"{c.name}: no time directories"); continue
    use = times[-nlast:]
    acc = {"maxall": [], "maxmetal": [], "meanmetal": [], "wmean": [], "ncell": []}
    for t in use:
        ut, at = read(t / "U"), read(t / "alpha.metal")
        if ut is None: continue
        U = vectors(ut)
        if U is None: continue
        a = scalars(at, len(U)) if at else [1.0] * len(U)
        mag = [math.sqrt(x*x + y*y + z*z) for x, y, z in U]
        metal = [m for m, av in zip(mag, a) if av > 0.5]
        wsum = sum(a)
        acc["maxall"].append(max(mag))
        acc["ncell"].append(len(metal))
        if metal:
            acc["maxmetal"].append(max(metal))
            acc["meanmetal"].append(sum(metal) / len(metal))
        acc["wmean"].append(sum(m * av for m, av in zip(mag, a)) / wsum if wsum else 0.0)

    if not acc["maxall"]:
        print(f"{c.name}: could not parse any time directory (binary?)"); continue

    def av(k): return sum(acc[k]) / len(acc[k]) if acc[k] else float("nan")
    span = use[0].name if nlast == 1 else f"{use[0].name}..{use[-1].name}"
    print(f"\n{c.name}   (t = {span}, {len(acc['maxall'])} write(s))")
    print(f"  melt volume (cells) : {av('ncell'):9.1f}        <-- tightest floor, ~0.5%")
    print(f"  mean|U| in metal    : {av('meanmetal'):9.4f} m/s   <-- USE THIS for velocity, ~2.6%")
    print(f"  max|U| in metal     : {av('maxmetal'):9.4f} m/s")
    print(f"  alpha-weighted mean : {av('wmean'):9.4f} m/s   (includes noisy interface cells)")
    print(f"  max|U| whole domain : {av('maxall'):9.4f} m/s   (often in GAS -- unreliable)")
