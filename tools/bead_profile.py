#!/usr/bin/env python3
"""Bead geometry along the track: crown height, width, cross-sectional area.

WHY THIS SCRIPT EXISTS

  The one validation gate we are failing is the bead crown: 2.60 mm against
  Zhao's 1.55 mm (simulated) and 1.80 mm (experimental), i.e. +68% against a
  15% tolerance. That 2.60 mm was read off a ParaView Plot Over Line, and it
  has two problems, one of measurement and one of location.

  THE MEASUREMENT PROBLEM. A Plot Over Line through the bead reports where
  the alpha = 0.5 contour crosses, and alpha = 0.5 is a threshold. We already
  learned once (check_droplet.py rev 2) that thresholding a VOF field on a
  smeared interface loses real material: 1.5 mm^3 of aluminium appeared to
  evaporate when a droplet flattened and a ring of cells fell from alpha 0.55
  to 0.45. The same failure applies to a height. This script never
  thresholds. For each vertical column of cells above the plate it computes

      h(x,y) = sum over the column of  alpha_c * V_c / (dx * dy)

  which is the EQUIVALENT METAL HEIGHT of that column -- the height the metal
  in that column would have if it were collected into a sharp-topped block of
  the same footprint. For a column through the crown apex this is exactly the
  crown height, and it is exact whether the interface is one cell thick or
  four. It is the height analogue of 'aV', not of 'a>.5'.

  THE LOCATION PROBLEM, which is the more serious of the two. Travel speed is
  120 cm/min = 20 mm/s. At t = 0.5 s the arc had advanced 10 mm -- roughly two
  pool lengths. The first droplets in any weld land on a plate that has not
  melted yet, so they do not wet and spread, they pile; that pile stays at the
  start of the track for the rest of the run, and it is the tallest thing in
  the domain. Zhao's 1.55 / 1.80 mm is a steady-state cross-section from a
  macrograph cut well down a finished track, where the bead has stopped
  changing. Measuring our start-transient pile against their developed bead
  compares two different objects.

  This script therefore reports the bead as a FUNCTION OF x and identifies the
  developed window automatically, rather than quoting one number from wherever
  the line happened to be drawn. Two boundaries define that window:

    - the trailing edge of the pool, x_tail = the smallest x at which any
      metal cell is still above the solidus. Everything behind that has
      solidified and will not change again, so it is the only region whose
      geometry is final.
    - the end of the start transient, taken as 5 mm ahead of where deposition
      begins. Five millimetres is a little over one pool length and about
      seven droplet spacings (0.665 mm each).

  If fewer than three x-stations fall between them, the bead is not developed
  yet and the script says so instead of quoting a number.

WHAT IT REPORTS

  A         Cross-sectional area of deposited metal above z = 0, mm^2. The
            primary geometric quantity, because it is a pure sum(alpha*V) and
            therefore conserved -- it must equal (wire feed rate)/(travel
            speed) once the bead is developed, independent of what shape the
            free surface takes. If A is right and h is wrong, the problem is
            surface shape (arc pressure, surface tension). If A itself is
            wrong, the problem is mass delivery.
  h crown   h(x, y = 0), the equivalent metal height on the centreline. THIS
            is the number to compare against Zhao 1.55 / 1.80 mm.
  h max     max over y of h(x,y). Equal to 'h crown' for a symmetric bead;
            larger means the crown has moved off centre, which matters because
            we already know the Marangoni circulation in this case is
            one-directional rather than the symmetric double cell it should be.
  w toe     Full width where h exceeds 10% of the local maximum. This is the
            closest model analogue of the toe-to-toe width an experimentalist
            measures on a macrograph, and it is what compares against Zhao's
            8.2 mm.
  w half    Full width at half maximum. Reported beside 'w toe' for the same
            reason 'a>.5' is reported beside 'aV': the gap between two
            definitions of the same quantity is information, and hiding it
            invites quoting whichever one flatters.
  ripple    Peak-to-peak variation of A(x) across the developed window, as a
            percentage of its mean. Droplets land every 0.665 mm; if that
            periodicity survives into the solidified bead the deposit is
            humping rather than flowing, which is a real defect and a real
            difference from Zhao's smooth bead. A few percent is discretisation
            noise. Tens of percent is humping.

HALF vs FULL DOMAIN
  Read from constant/caseParameters, the same file the solver reads. On a half
  domain (sym half) the mesh spans y >= 0 only, so areas and widths are
  doubled here and the fact is printed. Getting this wrong halves the crown
  comparison silently, which is exactly the class of error that produced the
  756 W / 1512 W energy discrepancy.

BEFORE FIRST USE, in the case directory:
    postProcess -func writeCellCentres -time 0
    postProcess -func writeCellVolumes -time 0

  These write Cx, Cy, Cz and V into 0/. They describe the mesh, not the
  solution, and only need regenerating when the mesh changes -- but Allrun
  recreates 0/ from initial/, so they DO need regenerating after every clean.

USAGE
    cd <case>
    python3 bead_profile.py             # summary of the developed window, per time
    python3 bead_profile.py 1.5         # full x-profile and transverse section at t
    python3 bead_profile.py 1.5 --csv   # also write bead_profile_1.5.csv
"""
import gzip, math, re, sys
from pathlib import Path

CASE = Path.cwd()
TSOL, TLIQ = 815.0, 906.0      # Al-5Mg, materials/transportProperties
TRANSIENT_MM = 5.0             # start-transient margin, ~1 pool length
MIN_H_MM = 0.05                # a column with less than this is bare plate


# --------------------------------------------------------------------------
# Field reading. Identical machinery to pool_stats_waam.py -- see the notes
# there on why the parser reads a declared token count rather than hunting
# for the closing parenthesis, and why '.gz' cannot go through with_suffix.
# --------------------------------------------------------------------------
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
            toks = tail.split(None, n)[:n]
            if len(toks) != n:
                sys.exit(f"ERROR: {path.name} declared {n} values, found {len(toks)}")
            return [float(v) for v in toks]
        toks = tail.split(None, 3*n)[:3*n]
        f = [float(t.strip("()")) for t in toks]
        return list(zip(f[0::3], f[1::3], f[2::3]))
    m = re.search(r"internalField\s+uniform\s+([^;]+);", txt)
    if m:
        tok = m.group(1).strip()
        return [tuple(float(x) for x in tok.strip("()").split())] if tok.startswith("(") \
               else [float(tok)]
    return None


def fit(f, n):
    """Make a parsed field safe to index per cell, or return None.

    A field that nothing has touched yet is written as 'internalField uniform
    300', and read_field returns that as a ONE-element list -- correct, and
    the cheapest possible representation, but fatal the moment something
    indexes it by cell number. That is exactly how this script died on its
    first contact with a real case: T in 0/ is uniform 300, the loop reached
    cell 1, and IndexError.

    Broadcasting is the right repair rather than skipping the time, because a
    uniform field is not missing data -- every cell really does hold that
    value. A field of any OTHER wrong length is a genuine mismatch (a stale
    mesh, a partially written file) and is rejected as None instead, so a
    caller can skip it rather than silently read a shifted array.
    """
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


def case_params():
    """sym, arc_x0_m and arc travel speed, from the file the solver reads.

    None of these are guessed. A wrong 'sym' silently halves every width and
    area; a wrong arc_x0 puts the arc in the wrong place and mislabels which
    part of the bead has solidified.
    """
    f = CASE / "constant" / "caseParameters"
    p = {"sym": "full", "x0": 0.0, "v": None}
    if not f.exists():
        return p, "defaults (no caseParameters)"
    txt = f.read_text()
    s = re.search(r"^\s*sym\s+(\w+)\s*;", txt, re.M)
    if s:
        p["sym"] = s.group(1)
    for key, name in (("x0", "arc_x0_m"), ("v", "arcV_mps")):
        m = re.search(rf"^\s*{name}\s+([0-9eE+.-]+)\s*;", txt, re.M)
        if m:
            p[key] = float(m.group(1))
    if p["v"] is None:                       # older cases spell it differently
        for name in ("arc_speed_mps", "travel_speed_mps", "arcV_m_per_s"):
            m = re.search(rf"^\s*{name}\s+([0-9eE+.-]+)\s*;", txt, re.M)
            if m:
                p["v"] = float(m.group(1))
                break
    return p, "constant/caseParameters"


# --------------------------------------------------------------------------
# Mesh structure. Cell sizes are derived from the spacing between neighbouring
# cell-centre planes rather than assumed uniform, so a graded mesh gives the
# same answer -- the same reason pool_stats_waam.py refuses to guess volumes.
# --------------------------------------------------------------------------
def stations(coords):
    """Sorted unique centre coordinates, and the cell extent at each.

    For an interior plane the extent is the half-distance to each neighbour
    summed, which is exact for a hex mesh however it is graded. At the two
    ends there is only one neighbour, so the one-sided spacing is used; those
    planes are domain boundaries and never carry bead.
    """
    xs = sorted({round(c, 9) for c in coords})
    if len(xs) < 2:
        sys.exit("ERROR: mesh has a single cell plane in one direction.")
    d = []
    for i, x in enumerate(xs):
        lo = xs[i] - xs[i-1] if i > 0 else xs[1] - xs[0]
        hi = xs[i+1] - xs[i] if i < len(xs)-1 else xs[-1] - xs[-2]
        d.append(0.5*lo + 0.5*hi)
    return xs, {x: dd for x, dd in zip(xs, d)}


def build_columns(C, V):
    """Map every cell above z = 0 to its (x-station, y-station) column.

    Returns the station lists, their extents, and for each cell an (i, j)
    index pair -- or None for cells at or below the plate surface, which are
    substrate and never part of the deposit.
    """
    xs, dx = stations(c[0] for c in C)
    ys, dy = stations(c[1] for c in C)
    xi = {x: i for i, x in enumerate(xs)}
    yj = {y: j for j, y in enumerate(ys)}
    idx = []
    for c in C:
        if c[2] > 0.0:
            idx.append((xi[round(c[0], 9)], yj[round(c[1], 9)]))
        else:
            idx.append(None)
    return xs, ys, [dx[x] for x in xs], [dy[y] for y in ys], idx


# --------------------------------------------------------------------------
# Geometry of one written time
# --------------------------------------------------------------------------
def profile(a, T, C, V, xs, ys, dx, dy, idx, fac):
    """Per-x-station area, and the full h(x,y) column-height map.

    h is in mm, A in mm^2, both already multiplied by the half/full domain
    factor. x_tail is the smallest x carrying a cell above the solidus, i.e.
    the trailing edge of the pool; everything behind it has solidified.
    """
    nx, ny = len(xs), len(ys)
    h = [[0.0]*ny for _ in range(nx)]
    A = [0.0]*nx
    x_tail = None
    for c in range(len(V)):
        p = idx[c]
        av = a[c]*V[c]
        if p is not None and av > 0.0:
            i, j = p
            A[i] += av
            h[i][j] += av/(dx[i]*dy[j])
        if T is not None and a[c] > 0.5 and T[c] > TSOL:
            x = C[c][0]
            if x_tail is None or x < x_tail:
                x_tail = x
    A = [v*1e6/dxi*fac for v, dxi in zip(A, dx)]      # m^3/m -> mm^2, doubled if half
    h = [[v*1e3 for v in row] for row in h]           # m -> mm
    return A, h, x_tail


def _interp_at(ys, hrow, y0):
    """h interpolated to an arbitrary y, linearly between column centres.

    The centreline matters here. On a full domain the y stations straddle
    y = 0 rather than landing on it, so 'the station nearest zero' is half a
    cell off centre; on a graded mesh it could be more. Interpolating costs
    nothing and removes the question.
    """
    if y0 <= ys[0]:
        return hrow[0]
    if y0 >= ys[-1]:
        return hrow[-1]
    for j in range(len(ys)-1):
        if ys[j] <= y0 <= ys[j+1]:
            f = (y0 - ys[j])/(ys[j+1] - ys[j])
            return hrow[j] + f*(hrow[j+1] - hrow[j])
    return hrow[0]


def _crossing(ys, hrow, lvl, jpk, direction):
    """Where h falls through lvl, walking outward from the peak column.

    Walking OUTWARD FROM THE PEAK rather than inward from the domain edge is
    what makes this robust to a second, detached blob of metal somewhere else
    in the plane -- a satellite droplet in flight, say, or the far end of the
    track at a different height. Inward-from-the-edge would jump the gap and
    report a width spanning both.

    The crossing is linearly interpolated between the last column above the
    level and the first below it. Taking the cell centre instead would
    under-report by half a cell on each side: 0.5 mm on an 8.2 mm width is
    6%, against a 15% gate, so it is not a rounding detail.
    """
    j = jpk
    n = len(ys)
    while 0 <= j + direction < n and hrow[j + direction] >= lvl:
        j += direction
    k = j + direction
    if not (0 <= k < n):
        return ys[j]                      # bead runs off the mesh; report the edge
    dh = hrow[j] - hrow[k]
    if dh <= 0.0:
        return ys[j]
    return ys[j] + (hrow[j] - lvl)/dh*(ys[k] - ys[j])


def widths(hrow, ys, fac):
    """Toe-to-toe (10% of peak) and half-maximum widths of one cross-section.

    On a HALF domain the mesh spans y >= 0 and the symmetry plane is y = 0, so
    the full width is twice the distance from the plane to the outer crossing
    -- not twice the span between the innermost and outermost columns, which
    would be the same thing only by accident and different whenever the bead
    does not reach the last column.
    """
    hmax = max(hrow) if hrow else 0.0
    if hmax <= MIN_H_MM:
        return 0.0, 0.0, 0.0, 0.0
    jpk = max(range(len(ys)), key=lambda j: hrow[j])
    out = []
    for lvl in (0.10*hmax, 0.50*hmax):
        yp = _crossing(ys, hrow, lvl, jpk, +1)
        if fac == 2.0:
            out.append(2.0*yp*1e3)
        else:
            yn = _crossing(ys, hrow, lvl, jpk, -1)
            out.append((yp - yn)*1e3)
    return out[0], out[1], _interp_at(ys, hrow, 0.0), hmax


def developed(A, h, xs, x_tail, fac):
    """The window of x that has solidified and is past the start transient.

    Returns (i0, i1) inclusive, or None. Deposition start is measured, not
    assumed: it is the first station carrying more than MIN_H_MM of metal on
    any column.
    """
    start = None
    for i in range(len(xs)):
        if max(h[i]) > MIN_H_MM:
            start = xs[i]
            break
    if start is None or x_tail is None:
        return None, start, x_tail
    lo, hi = start + TRANSIENT_MM*1e-3, x_tail
    win = [i for i in range(len(xs)) if lo <= xs[i] <= hi and max(h[i]) > MIN_H_MM]
    if len(win) < 3:
        return None, start, x_tail
    return (win[0], win[-1]), start, x_tail


# --------------------------------------------------------------------------
def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    want = float(args[0]) if args else None
    to_csv = "--csv" in sys.argv

    ts = times()
    if not ts:
        sys.exit("ERROR: no time directories here. Run this inside the case.")

    C = V = None
    for _, d in ts:
        cx, cy, cz, v = (read_field(d / n) for n in ("Cx", "Cy", "Cz", "V"))
        if cx and cy and cz and v and len(cx) == len(v) > 1:
            C, V, where = list(zip(cx, cy, cz)), v, d.name
            break
    if C is None:
        sys.exit("ERROR: no cell centres / volumes in any time directory.\n"
                 "  postProcess -func writeCellCentres -time 0\n"
                 "  postProcess -func writeCellVolumes -time 0")

    par, src = case_params()
    fac = 2.0 if par["sym"] == "half" else 1.0
    xs, ys, dx, dy, idx = build_columns(C, V)

    print(f"  mesh: {len(V):,} cells from {where}/  "
          f"({len(xs)} x-stations, {len(ys)} y-stations)")
    print(f"  domain: {par['sym']}  ->  areas and widths x{fac:g}   [{src}]")
    if par["v"]:
        print(f"  arc: x0 = {par['x0']*1e3:+.1f} mm, v = {par['v']*1e3:.1f} mm/s")
    print(f"  developed window = [deposition start + {TRANSIENT_MM:g} mm, pool trailing edge]")
    print()

    # ---------------- single time: full profile ----------------
    if want is not None:
        hit = [(t, d) for t, d in ts if abs(t - want) < 1e-9]
        if not hit:
            sys.exit(f"ERROR: no time directory {want}. Have: "
                     + ", ".join(d.name for _, d in ts))
        t, d = hit[0]
        n = len(V)
        a = fit(read_field(d / "alpha.metal"), n)
        T = fit(read_field(d / "T"), n)
        if a is None:
            sys.exit(f"ERROR: alpha.metal at t = {d.name} does not match the mesh "
                     f"({n:,} cells).")
        A, h, x_tail = profile(a, T, C, V, xs, ys, dx, dy, idx, fac)
        win, start, x_tail = developed(A, h, xs, x_tail, fac)

        hdr = (f"  {'x':>8}  {'A':>7}  {'h crown':>8}  {'h max':>7}  "
               f"{'w toe':>7}  {'w half':>7}  {'':>4}")
        print(f"  --- profile along the track at t = {d.name} s ---\n")
        print(hdr)
        print(f"  {'mm':>8}  {'mm^2':>7}  {'mm':>8}  {'mm':>7}  {'mm':>7}  {'mm':>7}")
        print("  " + "-"*len(hdr.strip()))
        rows = []
        for i, x in enumerate(xs):
            if max(h[i]) <= MIN_H_MM:
                continue
            wt, wh, hc, hm = widths(h[i], ys, fac)
            tag = "  <<" if win and win[0] <= i <= win[1] else ""
            print(f"  {x*1e3:8.2f}  {A[i]:7.3f}  {hc:8.3f}  {hm:7.3f}  "
                  f"{wt:7.3f}  {wh:7.3f}{tag}")
            rows.append((x*1e3, A[i], hc, hm, wt, wh))
        print("\n  '<<' marks the developed window: solidified and past the "
              "start transient.")
        if x_tail is not None:
            print(f"  pool trailing edge at x = {x_tail*1e3:.2f} mm; "
                  f"deposition starts at x = {start*1e3:.2f} mm")

        if win:
            i0, i1 = win
            jc = min(range(len(ys)), key=lambda j: abs(ys[j]))
            print(f"\n  --- transverse section, averaged over x = "
                  f"{xs[i0]*1e3:.2f} to {xs[i1]*1e3:.2f} mm ---\n")
            print(f"  {'y':>8}  {'h':>8}")
            print(f"  {'mm':>8}  {'mm':>8}")
            print("  " + "-"*18)
            for j, y in enumerate(ys):
                hbar = sum(h[i][j] for i in range(i0, i1+1))/(i1-i0+1)
                if hbar > 1e-4:
                    print(f"  {y*1e3:8.2f}  {hbar:8.3f}")

        if to_csv:
            out = CASE / f"bead_profile_{d.name}.csv"
            with out.open("w") as f:
                f.write("x_mm,A_mm2,h_crown_mm,h_max_mm,w_toe_mm,w_half_mm\n")
                for r in rows:
                    f.write(",".join(f"{v:.5f}" for v in r) + "\n")
            print(f"\n  wrote {out.name}")
        footer()
        return

    # ---------------- all times: developed-window summary ----------------
    hdr = (f"  {'time':>7}  {'x arc':>7}  {'window':>13}  {'n':>3}  "
           f"{'A':>7}  {'h crown':>8}  {'h max':>7}  {'w toe':>7}  {'w half':>7}  {'ripple':>7}")
    print(hdr)
    print(f"  {'s':>7}  {'mm':>7}  {'mm':>13}  {'':>3}  "
          f"{'mm^2':>7}  {'mm':>8}  {'mm':>7}  {'mm':>7}  {'mm':>7}  {'%':>7}")
    print("  " + "-"*len(hdr.strip()))

    for t, d in ts:
        a = fit(read_field(d / "alpha.metal"), len(V))
        if a is None:
            continue
        T = fit(read_field(d / "T"), len(V))
        A, h, x_tail = profile(a, T, C, V, xs, ys, dx, dy, idx, fac)
        win, start, x_tail = developed(A, h, xs, x_tail, fac)
        xarc = (par["x0"] + par["v"]*t)*1e3 if par["v"] is not None else float("nan")
        if not win:
            why = "no deposit" if start is None else "not developed"
            print(f"  {d.name:>7}  {xarc:7.2f}  {why:>13}")
            continue
        i0, i1 = win
        nst = i1 - i0 + 1
        Aw = [A[i] for i in range(i0, i1+1)]
        Abar = sum(Aw)/nst
        ripple = 100.0*(max(Aw) - min(Aw))/Abar if Abar > 0 else float("nan")
        wt = wh = hc = hm = 0.0
        for i in range(i0, i1+1):
            a1, a2, a3, a4 = widths(h[i], ys, fac)
            wt += a1; wh += a2; hc += a3; hm += a4
        print(f"  {d.name:>7}  {xarc:7.2f}  "
              f"{xs[i0]*1e3:5.1f} to {xs[i1]*1e3:5.1f}  {nst:3d}  "
              f"{Abar:7.3f}  {hc/nst:8.3f}  {hm/nst:7.3f}  "
              f"{wt/nst:7.3f}  {wh/nst:7.3f}  {ripple:7.1f}")
    footer()


def footer():
    print(f"""
  HOW TO READ THIS

  A         Deposited cross-section above z = 0, mm^2. Pure sum(alpha*V), so
            it is conserved and can be checked against mass in independently:
            once the bead is developed, A must equal the wire volume rate
            divided by the travel speed. If A is right and 'h crown' is wrong,
            the metal is all there and the FREE SURFACE SHAPE is wrong --
            which is the arc-pressure hypothesis. If A itself is wrong, the
            problem is upstream, in delivery.

  h crown   Equivalent metal height on the centreline. This is the number to
            put against Zhao 1.55 mm (simulated) and 1.80 mm (experimental).
            It is an alpha-weighted column integral, not an alpha = 0.5
            contour crossing, so a smeared interface does not inflate or
            deflate it.

  h max     Same quantity maximised over y. If this exceeds 'h crown' by more
            than a few percent the crown is off-centre, which would be a
            direct geometric consequence of the one-directional Marangoni
            circulation already seen in the U_x profile.

  w toe     Full width at 10% of peak height -- the model analogue of the
            toe-to-toe width measured on a macrograph. Compare to Zhao 8.2 mm.
  w half    Full width at half maximum, reported beside it so the choice of
            width definition stays visible rather than being made silently.

  ripple    Peak-to-peak spread of A(x) across the developed window, percent
            of mean. Droplets land every 0.665 mm (v/f). A few percent is
            discretisation. Tens of percent means the droplet periodicity is
            surviving into the solidified bead -- humping, not flow -- and
            that is a physical difference from Zhao's smooth bead, not a
            post-processing artefact.

  window    The x range that has BOTH solidified (behind the pool trailing
            edge) AND cleared the start transient (more than {TRANSIENT_MM:g} mm past
            first deposition). Only this region is comparable to a
            steady-state macrograph. 'not developed' on every row means the
            run is too short for a geometry comparison to mean anything, and
            no number quoted from it should enter the thesis.""")


if __name__ == "__main__":
    main()
