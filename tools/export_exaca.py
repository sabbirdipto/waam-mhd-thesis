#!/usr/bin/env python3
"""export_exaca.py -- turn exacaRecorder output into an ExaCA input package.

Run from inside a case directory (it works on partial data while the run is
still going, so the pipeline can be checked early):

    python3 ~/thesis/tools/export_exaca.py
    python3 ~/thesis/tools/export_exaca.py --xmin 0 --xmax 10     # cut a slab, mm
    python3 ~/thesis/tools/export_exaca.py --out /some/dir

Reads   postProcessing/exacaRecorder/events_p*_from*.csv   (every rank, every
        run segment -- a restart writes a new segment file)

Writes  <out>/exaca_<case>.csv   x,y,z,tm,tl,cr        ExaCA temperature file
        <out>/gr_<case>.csv      x,y,z,tl,cr,G,R       solidification map
        <out>/README.txt         definitions, grid, units, caveats
        default <out> = ~/thesis/results/exaca_<case>

What it does, in order
  1. merge all ranks and segments
  2. drop events whose melt time is unknown (tm = -1: the cell was already
     molten when a restarted segment began)
  3. keep only cells in the UNIFORM 0.5 mm core of the mesh (x -12..32 mm,
     y -6..6 mm, all z) and verify every kept centre sits on the regular
     lattice (k + 1/2) * 0.5 mm -- ExaCA reads a regular grid, and cells in the
     graded outer band would break that
  4. optional x-range cut, for a smaller ExaCA domain
  5. sort by position then time, so each cell's melt events are consecutive
     and chronological
  6. write the files and print a summary, including a per-x table that shows
     where the bead reaches steady state
"""
import argparse, csv, glob, math, os, statistics as st, sys

H = 0.5e-3                          # thermal grid spacing [m]
CORE = dict(x=(-12e-3, 32e-3), y=(-6e-3, 6e-3), z=(-4e-3, 5e-3))
TL, TS, DTCR = 906.0, 815.0, 10.0   # liquidus, solidus, cr window [K]


def on_lattice(v, h=H, tol=1e-6):
    k = v / h - 0.5
    return abs(k - round(k)) < tol


def idx(v, h=H):
    """integer lattice index k of a cell centre at (k + 1/2) * h.
    Not round(v/h): Python rounds halves to even, which merges neighbouring cells."""
    return round(v / h - 0.5)


def med(v):
    return st.median(v) if v else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xmin", type=float, help="mm")
    ap.add_argument("--xmax", type=float, help="mm")
    ap.add_argument("--out")
    a = ap.parse_args()

    case = os.path.basename(os.getcwd())
    files = sorted(glob.glob("postProcessing/exacaRecorder/events_p*_from*.csv"))
    if not files:
        sys.exit("no postProcessing/exacaRecorder/events_p*.csv here -- run from the case directory")

    rows, bad = [], 0
    for f in files:
        with open(f) as fh:
            rd = csv.reader(fh)
            next(rd, None)
            for r in rd:
                if len(r) != 7:
                    bad += 1
                    continue
                try:
                    rows.append([float(v) for v in r])
                except ValueError:
                    bad += 1
    n_raw = len(rows)

    unknown = [r for r in rows if r[3] < 0]
    rows = [r for r in rows if r[3] >= 0]

    def in_core(r):
        return all(CORE[k][0] < r[i] < CORE[k][1] for i, k in enumerate("xyz"))
    outside = [r for r in rows if not in_core(r)]
    rows = [r for r in rows if in_core(r)]
    off = [r for r in rows if not all(on_lattice(r[i]) for i in range(3))]
    if off:
        sys.exit(f"{len(off)} events inside the core are NOT on the 0.5 mm lattice "
                 f"(first: {off[0][:3]}). The mesh is not what this script assumes -- stop.")

    if a.xmin is not None:
        rows = [r for r in rows if r[0] >= a.xmin * 1e-3]
    if a.xmax is not None:
        rows = [r for r in rows if r[0] <= a.xmax * 1e-3]
    bad_order = sum(1 for r in rows if r[4] <= r[3])
    rows.sort(key=lambda r: (idx(r[0]), idx(r[1]), idx(r[2]), r[3]))
    if not rows:
        sys.exit("no events left after filtering")

    out = a.out or os.path.expanduser(f"~/thesis/results/exaca_{case}")
    os.makedirs(out, exist_ok=True)
    f_exa = os.path.join(out, f"exaca_{case}.csv")
    f_gr = os.path.join(out, f"gr_{case}.csv")

    with open(f_exa, "w") as fh:
        fh.write("x,y,z,tm,tl,cr\n")
        for x, y, z, tm, tl, cr, G in rows:
            fh.write(f"{x:.6e},{y:.6e},{z:.6e},{tm:.9e},{tl:.9e},{cr:.6e}\n")
    with open(f_gr, "w") as fh:
        fh.write("x,y,z,tl,cr,G,R\n")
        for x, y, z, tm, tl, cr, G in rows:
            R = cr / G if G > 1.0 else float("nan")
            fh.write(f"{x:.6e},{y:.6e},{z:.6e},{tl:.9e},{cr:.6e},{G:.6e},{R:.6e}\n")

    cells = {}
    for r in rows:
        cells.setdefault((idx(r[0]), idx(r[1]), idx(r[2])), []).append(r)
    remelted = sum(1 for v in cells.values() if len(v) > 1)
    col = list(zip(*rows))
    x, y, z, tm, tl, cr, G = col
    R = [c / g for c, g in zip(cr, G) if g > 1.0]
    lo = [min(v) for v in (x, y, z)]
    hi = [max(v) for v in (x, y, z)]
    ext = [(h - l) + H for l, h in zip(lo, hi)]      # cell-centre span + one cell

    print(f"case {case}")
    print(f"  {len(files)} file(s), {n_raw} raw events, {bad} malformed lines")
    print(f"  dropped: {len(unknown)} with unknown tm (restart), {len(outside)} outside the uniform core")
    print(f"  kept   : {len(rows)} events in {len(cells)} cells ({remelted} cells melted more than once)")
    print(f"  tl <= tm: {bad_order}   (must be 0)")
    print(f"  domain (cell centres) mm: x {lo[0]*1e3:.2f}..{hi[0]*1e3:.2f}  "
          f"y {lo[1]*1e3:.2f}..{hi[1]*1e3:.2f}  z {lo[2]*1e3:.2f}..{hi[2]*1e3:.2f}")
    for name, v, u in (("tm", tm, "s"), ("tl", tl, "s"), ("cr", cr, "K/s"),
                       ("G", G, "K/m"), ("R", R, "m/s")):
        if v:
            print(f"  {name:2s} {u:4s} min {min(v):10.4g}  median {med(v):10.4g}  max {max(v):10.4g}")
    for um in (5, 10):
        n = 1
        for e in ext:
            n *= math.ceil(e / (um * 1e-6))
        print(f"  ExaCA cells for this box at {um:2d} um: {n/1e6:8.1f} M")

    print("  per-x steady-state check (2 mm bins):")
    print("      x mm   events   median cr K/s   median G K/m   median tl s")
    b0 = math.floor(lo[0] * 1e3 / 2) * 2
    b = b0
    while b <= hi[0] * 1e3:
        sel = [r for r in rows if b <= r[0] * 1e3 < b + 2]
        if sel:
            print(f"   {b:5.0f}..{b+2:<4.0f}{len(sel):7d}   {med([r[5] for r in sel]):12.4g}"
                  f"   {med([r[6] for r in sel]):12.4g}   {med([r[4] for r in sel]):10.4f}")
        b += 2

    with open(os.path.join(out, "README.txt"), "w") as fh:
        fh.write(f"""ExaCA input from OpenFOAM case {case}
{'=' * (26 + len(case))}

FILES
  exaca_{case}.csv   ExaCA temperature file: x,y,z,tm,tl,cr
  gr_{case}.csv      solidification map:    x,y,z,tl,cr,G,R

UNITS            x,y,z m | tm,tl s | cr K/s | G K/m | R = cr/G m/s

GRID             regular, spacing 500 um (HeatTransferCellSize = 500 um).
                 Points are cell centres at (k + 1/2) * 0.5 mm; every point in
                 the file was checked to lie on that lattice.
                 Rows sorted by position, then time: a cell that melted more
                 than once has consecutive rows in chronological order.

COORDINATES      z = 0 is the top of the 4 mm plate, +z is the build direction.
                 y = 0 is the bead centreline. The torch starts at x = -10 mm at
                 t = 0 and travels +x at 20 mm/s.

DEFINITIONS      liquid metal : metal fraction > 0.5 and T >= {TL:.0f} K (liquidus)
                 tm : time the cell became liquid metal. Interpolated liquidus
                      crossing if it melted in place; arrival time if metal
                      arrived already molten (droplet)
                 tl : time it cooled back through the liquidus, interpolated
                 cr : {DTCR:.0f} K / (time taken to cool from {TL:.0f} K to {TL-DTCR:.0f} K),
                      i.e. the mean cooling rate over the first {DTCR:.0f} K of
                      undercooling -- the quantity ExaCA's linear undercooling
                      uses. A dip below the liquidus that re-melts before
                      {TL-DTCR:.0f} K does not count as solidification.
                 G  : |grad T| at tl (unreliable at the free surface, where
                      the neighbouring cell is gas)
                 Recorded every solver time step (~1e-5 s), not from snapshots.

MATERIAL         Al-5Mg (5A05 plate, ER5356 wire). Solidus {TS:.0f} K, liquidus
                 {TL:.0f} K, freezing range {TL-TS:.0f} K.
                 ExaCA needs an interfacial response function for Al-5Mg; it is
                 not part of this package.

PROCESS          GMAW WAAM, 135 A, 16 V, travel 20 mm/s, wire 1.2 mm at 150 mm/s,
                 single bead on a 50 x 30 x 4 mm plate (Zhao et al. 2021).

CAVEATS          Model state "M2": energy-consistent arc (1113.7 W into the pool,
                 722 W carried by droplets) with Zhao's heat-source dimensions.
                 Penetration matches Zhao (~1 mm liquid) but the fused zone is
                 ~3.6 mm wide against his ~7.8 mm bead, so solidification
                 conditions are those of a narrower pool than the experiment.
                 Treat this as pipeline data. The calibrated 0.25 mm runs will
                 replace it with the same format.
                 The first ~5 mm of track is a start transient; use the per-x
                 table below to choose a steady region.
                 Cells inside the bounding box that never appear here never
                 melted: plate below z = 0, gas above the bead. Only cells
                 listed here have a solidification history.

SUMMARY
  events {len(rows)} in {len(cells)} cells, {remelted} remelted
  domain (cell centres) mm: x {lo[0]*1e3:.2f}..{hi[0]*1e3:.2f}  y {lo[1]*1e3:.2f}..{hi[1]*1e3:.2f}  z {lo[2]*1e3:.2f}..{hi[2]*1e3:.2f}
  cr K/s median {med(list(cr)):.4g}   G K/m median {med(list(G)):.4g}
""")
    print(f"  wrote {f_exa}")
    print(f"        {f_gr}")
    print(f"        {os.path.join(out, 'README.txt')}")


if __name__ == "__main__":
    main()
