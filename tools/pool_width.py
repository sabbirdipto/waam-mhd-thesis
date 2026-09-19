#!/usr/bin/env python3
"""How wide is the melted substrate, against how wide the bead actually is?

THE QUESTION THIS SETTLES

    The bead is 6.5 mm wide against Zhao's 8.2, and the height error follows
    from the width error (a parabolic cap at the right area and the right
    width has the right height automatically). Three hypotheses compete for
    the width deficit, and they call for completely different work:

      FORCE    the substrate melts 8 mm wide but the liquid does not spread
               out to fill it, because 663 Pa of flattening load is missing.
               Fix: implement Eqs. 22-24.

      THERMAL  the substrate only melts 6.5 mm wide, so the bead is already
               exactly as wide as the molten region it can wet. Fix: nothing
               mechanical helps; the heat input distribution is wrong.

      MESH     the interface and the fusion depth are under-resolved at
               0.5 mm against Zhao's 0.25 mm. Fix: refine.

    One measurement separates FORCE from THERMAL, and it needs no new run:

        w_fuse   the transverse extent of MELTED SUBSTRATE at the plate
                 surface -- how wide a puddle the arc has actually made
        w_toe    the transverse extent of the deposited bead, from
                 bead_profile.py

    w_fuse near 8 mm with w_toe at 6.5 says the metal is not spreading into
    available liquid: FORCE. w_fuse near 6.5 says the bead already fills its
    puddle and no force will widen it: THERMAL. There is no third reading.

WHY THE SUBSTRATE, SPECIFICALLY, AND NOT 'THE POOL'

    pool_stats_waam.py counts every metal cell above the solidus, which
    includes the deposited bead -- metal that arrived at 1506 K and is still
    cooling. That number is dominated by the deposit and says almost nothing
    about how far the arc melted the plate sideways.

    What decides spreading is the FUSION ZONE: substrate that was solid and
    is now liquid. Liquid metal running outward over cold solid plate freezes
    at the contact line and stops; running over its own melt it keeps going.
    So the measurement is restricted to cells BELOW z = 0 in the topmost
    substrate layer, which is the plate surface the bead is spreading across.

    'w melt' is reported beside it -- the widest molten extent at any depth --
    because a pool that is wide at the surface and narrow below is a different
    object from one that is wide throughout, and the Lorentz force is supposed
    to change exactly that.

THE EDGE IS INTERPOLATED, NOT COUNTED

    Taking the outermost cell CENTRE above the solidus under-reports the
    fusion width by half a cell each side: 0.5 mm on an 8 mm width is 6%,
    against a 15% gate. The solidus crossing is therefore linearly
    interpolated in T between the last cell above it and the first below.
    Walking outward from the hottest column rather than inward from the
    domain edge, so a detached hot spot elsewhere in the plane cannot bridge
    the gap and report one pool spanning both.

BEFORE FIRST USE, in the case directory:
    postProcess -func writeCellCentres -time 0
    postProcess -func writeCellVolumes -time 0

USAGE
    cd <case>
    python3 pool_width.py            # every written time, summary per time
    python3 pool_width.py 0.44       # per-x profile at one time
"""
import gzip, re, sys
from pathlib import Path

CASE = Path.cwd()
TSOL, TLIQ = 815.0, 906.0
SNAP_TOL = 1e-7
ZHAO_W = 8.2


def slurp(path):
    if path.exists():
        return path.read_text()
    gz = Path(str(path) + ".gz")
    if gz.exists():
        with gzip.open(gz, "rt") as f:
            return f.read()
    return None


def read_field(path):
    txt = slurp(path)
    if txt is None:
        return None
    m = re.search(r"internalField\s+nonuniform\s+List<(scalar|vector)>\s*\n?\s*(\d+)\s*\(", txt)
    if m:
        kind, n = m.group(1), int(m.group(2))
        tail = txt[m.end():]
        if kind == "scalar":
            return [float(v) for v in tail.split(None, n)[:n]]
        f = [float(t.strip("()")) for t in tail.split(None, 3*n)[:3*n]]
        return list(zip(f[0::3], f[1::3], f[2::3]))
    m = re.search(r"internalField\s+uniform\s+([^;]+);", txt)
    if m:
        tok = m.group(1).strip()
        return [tuple(float(x) for x in tok.strip("()").split())] if tok.startswith("(") \
               else [float(tok)]
    return None


def fit(f, n):
    if f is None:
        return None
    if len(f) == 1:
        return f*n
    return f if len(f) == n else None


def times():
    out = []
    for d in CASE.iterdir():
        if d.is_dir():
            try:
                out.append((float(d.name), d))
            except ValueError:
                pass
    return sorted(out)


def planes(coords, tol=SNAP_TOL):
    """Distinct cell-centre planes, merged on a tolerance.

    Exact equality is the wrong test: centres on one physical plane are
    written to six significant figures and come back differing in the ninth
    decimal, which splits a plane in two.
    """
    vals = sorted({round(c, 12) for c in coords})
    reps, look = [], {}
    for v in vals:
        if reps and v - reps[-1] <= tol:
            look[v] = len(reps) - 1
        else:
            reps.append(v)
            look[v] = len(reps) - 1
    return reps, look


def edge(ys, Ts, jpk, direction, level):
    """Interpolated y where T falls through `level`, walking out from the peak."""
    j, n = jpk, len(ys)
    while 0 <= j + direction < n and Ts[j + direction] is not None \
            and Ts[j + direction] >= level:
        j += direction
    k = j + direction
    if not (0 <= k < n) or Ts[k] is None:
        return ys[j]
    dT = Ts[j] - Ts[k]
    if dT <= 0.0:
        return ys[j]
    return ys[j] + (Ts[j] - level)/dT*(ys[k] - ys[j])


def width(ys, Ts, level):
    """Full transverse extent where T >= level, or 0 if nothing is that hot."""
    hot = [j for j, t in enumerate(Ts) if t is not None and t >= level]
    if not hot:
        return 0.0
    jpk = max(hot, key=lambda j: Ts[j])
    return (edge(ys, Ts, jpk, +1, level) - edge(ys, Ts, jpk, -1, level))*1e3


def main():
    want = float(sys.argv[1]) if len(sys.argv) > 1 else None
    ts = times()
    if not ts:
        sys.exit("ERROR: no time directories here. Run this inside the case.")

    C = None
    for _, d in ts:
        cx, cy, cz = (read_field(d / n) for n in ("Cx", "Cy", "Cz"))
        if cx and cy and cz and len(cx) == len(cz) > 1:
            C, where = list(zip(cx, cy, cz)), d.name
            break
    if C is None:
        sys.exit("ERROR: no cell centres. Run\n"
                 "  postProcess -func writeCellCentres -time 0")
    n = len(C)

    xs, xl = planes(c[0] for c in C)
    ys, yl = planes(c[1] for c in C)

    # The topmost SUBSTRATE layer: the highest cell-centre plane below z = 0.
    # This is the plate surface the bead spreads across. Cells above it are
    # deposit, and their temperature says nothing about how far the arc melted
    # the plate sideways.
    zsub = max(c[2] for c in C if c[2] < 0.0)
    surf = [i for i in range(n) if abs(C[i][2] - zsub) < SNAP_TOL]
    print(f"  mesh: {n:,} cells from {where}/  "
          f"({len(xs)} x-planes, {len(ys)} y-planes)")
    print(f"  substrate surface layer at z = {zsub*1e3:.3f} mm, "
          f"{len(surf):,} cells")
    print(f"  solidus {TSOL:g} K; Zhao bead width {ZHAO_W:g} mm\n")

    def measure(d):
        T = fit(read_field(d / "T"), n)
        a = fit(read_field(d / "alpha.metal"), n)
        if T is None or a is None:
            return None
        # Substrate surface: T by (x, y), metal only.
        grid = {}
        for i in surf:
            if a[i] > 0.5:
                grid[(xl[round(C[i][0], 12)], yl[round(C[i][1], 12)])] = T[i]
        # Any depth, metal only -- the widest molten extent anywhere.
        deep = {}
        for i in range(n):
            if a[i] > 0.5 and C[i][2] < 0.0:
                k = (xl[round(C[i][0], 12)], yl[round(C[i][1], 12)])
                deep[k] = max(deep.get(k, 0.0), T[i])
        out = []
        for ix in range(len(xs)):
            Ts = [grid.get((ix, jy)) for jy in range(len(ys))]
            Td = [deep.get((ix, jy)) for jy in range(len(ys))]
            wf = width(ys, Ts, TSOL)
            wl = width(ys, Ts, TLIQ)
            wm = width(ys, Td, TSOL)
            if wf > 0 or wm > 0:
                out.append((xs[ix]*1e3, wf, wl, wm))
        return out

    if want is not None:
        hit = [(t, d) for t, d in ts if abs(t - want) < 1e-9]
        if not hit:
            sys.exit("ERROR: no time " + str(want) + ". Have: "
                     + ", ".join(d.name for _, d in ts))
        rows = measure(hit[0][1])
        print(f"  --- fusion width along the track at t = {hit[0][1].name} s ---\n")
        print(f"  {'x':>8}  {'w fuse':>8}  {'w liq':>8}  {'w melt':>8}")
        print(f"  {'mm':>8}  {'mm':>8}  {'mm':>8}  {'mm':>8}")
        print("  " + "-"*40)
        for x, wf, wl, wm in rows:
            print(f"  {x:8.2f}  {wf:8.3f}  {wl:8.3f}  {wm:8.3f}")
        footer()
        return

    print(f"  {'time':>14}  {'w fuse max':>10}  {'w liq max':>9}  "
          f"{'w melt max':>10}  {'x at max':>8}  {'vs Zhao':>8}")
    print(f"  {'s':>14}  {'mm':>10}  {'mm':>9}  {'mm':>10}  {'mm':>8}  {'%':>8}")
    print("  " + "-"*70)
    for t, d in ts:
        rows = measure(d)
        if not rows:
            continue
        best = max(rows, key=lambda r: r[1])
        if best[1] <= 0 and max(r[3] for r in rows) <= 0:
            print(f"  {d.name:>14}  {'no fusion':>10}")
            continue
        print(f"  {d.name:>14}  {best[1]:10.3f}  {best[2]:9.3f}  "
              f"{max(r[3] for r in rows):10.3f}  {best[0]:8.2f}  "
              f"{100*(best[1]/ZHAO_W - 1):+8.1f}")
    footer()


def footer():
    print(f"""
  HOW TO READ THIS

  w fuse   Transverse extent of MELTED SUBSTRATE at the plate surface --
           substrate metal above the {TSOL:g} K solidus, in the topmost cell layer
           below z = 0. This is the puddle the bead can spread across, and it
           is the number to put against Zhao's {ZHAO_W:g} mm bead width.

  w liq    The same, but above the liquidus. The gap between w fuse and w liq
           is the mushy rim. If that rim is a large fraction of w fuse, the
           pool edge is slush rather than liquid, and slush does not flow --
           note that laserbeamFoam also zeroes the Marangoni force below the
           liquidus, so the rim gets no outward drive either.

  w melt   Widest molten extent at ANY depth in the substrate. Wider than
           w fuse means the pool is undercut -- broader below the surface than
           at it. That is the shape the Lorentz force is supposed to produce,
           so it is worth having a before-and-after number for it.

  THE READING

    w fuse near {ZHAO_W:g} mm, bead w toe near 6.5 mm
        The substrate melts as wide as Zhao's and the metal is NOT spreading
        into it. A force problem. Eqs. 22-24 are the right next job.

    w fuse near 6.5 mm, matching w toe
        The bead already fills its puddle. No body force widens a bead beyond
        the metal it can wet, so neither the Lorentz force nor mesh refinement
        addresses this, and the heat input distribution -- Goldak b, the
        eta_1 split, or both -- is what needs re-examining.

    w fuse near ZERO
        The arc is not melting the plate at all and the deposit is sitting on
        solid metal. Everything downstream of that is a wetting problem, not a
        welding one.""")


if __name__ == "__main__":
    main()
