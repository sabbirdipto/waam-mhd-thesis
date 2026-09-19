#!/usr/bin/env python3
"""Compare two bead_profile.py runs station-by-station. Ablation and mesh tests.

WHY A TOOL RATHER THAN READING TWO TABLES SIDE BY SIDE

    The window-mean rows in bead-summary.txt cannot be compared directly, and
    A1' vs A2 showed exactly how that misleads. Their developed windows began
    one station apart -- -7.8 vs -7.2 mm -- because deposition started at a
    marginally different x. The extra station in A1' is the thinnest in the
    window, so dropping it RAISES the mean. Read off the summary rows, arc
    pressure looked like it made the crown 4.8% TALLER. Compared station by
    station over the common x range it is +2.2%, and the sign of the aspect
    ratio flips from worse to better. The window boundary was carrying more of
    the difference than the physics was.

    So: intersect the two developed windows, compare only where both have
    solidified bead, and say how many stations that was.

THE MESH CASE, WHICH IS WHY THIS INTERPOLATES

    Two runs at the same cell size share an x grid and match exactly. A mesh
    test does not: 0.5 mm puts stations at -7.25, -6.75, -6.25 and 0.35 mm
    puts them at -7.19, -6.84, -6.49. Nothing lines up.

    Both profiles are therefore linearly interpolated onto a common grid, and
    the grid is taken from the COARSER run. Interpolating the fine result down
    is honest -- it is a local average of data that exists. Interpolating the
    coarse result UP would invent structure between its stations and then
    report the invention as a difference. The direction matters.

WHAT IT REPORTS

    A          cross-section, mm^2. Should be nearly unchanged by anything
               except a delivery bug: it is set by wire feed over travel speed
               and is conserved. A real shift here is the first thing to
               explain, before any shape claim is made.
    h crown    centreline height, the Table 5 number.
    w toe      full width at 10% of peak, the macrograph analogue.
    h/w        THE aspect ratio. Zhao sim 1.55/8.20 = 0.1890.
    h/(3A/2w)  the parabolic ratio. Zhao's own numbers give 0.999, so this is
               the sharpest shape target available and it needs no figure
               reading. 1.0 is a parabolic cap; above 1.0 is more peaked.

USAGE
    python3 compare_bead.py REF.txt TEST.txt
    python3 compare_bead.py A.txt B.txt --labels "A1' 0.5mm" "M1 0.35mm"

    Both arguments are per-x profile files -- the output of
        python3 bead_profile.py <time>
    not the summary. The developed window is taken from the '<<' markers.
"""
import re, sys

ZHAO_H, ZHAO_W, ZHAO_A = 1.550, 8.200, 8.480

ROW = re.compile(
    r"^\s*(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)(\s*<<)?\s*$"
)


def read(path):
    """(x, A, h_crown, w_toe) rows, and the subset flagged as developed.

    The '<<' marker is what bead_profile.py prints for stations that have BOTH
    solidified and cleared the start transient. Rows without it are still
    molten or still in the lead-in, and a geometry number taken from them is
    not comparable to a steady-state macrograph -- which is the entire reason
    the marker exists. They are parsed anyway so the profile can be
    interpolated across the window edges, but only marked rows set the range.
    """
    rows, dev = [], []
    for ln in open(path):
        m = ROW.match(ln.rstrip("\n"))
        if not m:
            continue
        x, A, hc, hm, wt, wh = (float(g) for g in m.groups()[:6])
        rows.append((x, A, hc, wt))
        if m.group(7):
            dev.append(x)
    if not rows:
        sys.exit(f"ERROR: no profile rows in {path}.\n"
                 f"  This wants the per-x profile -- 'python3 bead_profile.py "
                 f"<time>' -- not bead-summary.txt.")
    rows.sort()
    return rows, (min(dev), max(dev)) if dev else None


def interp(rows, x):
    """Linear interpolation of (A, h, w) to an arbitrary x. None outside."""
    if x < rows[0][0] or x > rows[-1][0]:
        return None
    for i in range(len(rows) - 1):
        x0, x1 = rows[i][0], rows[i + 1][0]
        if x0 <= x <= x1:
            f = 0.0 if x1 == x0 else (x - x0) / (x1 - x0)
            return tuple(rows[i][k] + f * (rows[i + 1][k] - rows[i][k])
                         for k in (1, 2, 3))
    return None


def spacing(rows):
    xs = [r[0] for r in rows]
    d = sorted(xs[i + 1] - xs[i] for i in range(len(xs) - 1))
    return d[len(d) // 2]                      # median: robust to graded ends


def main():
    # --labels takes two VALUES, and they are not paths. Filtering only on a
    # leading '--' leaves them in the positional list, which made the first
    # test of this script print its own usage instead of a comparison.
    argv, labs = sys.argv[1:], ["REF", "TEST"]
    if "--labels" in argv:
        i = argv.index("--labels")
        labs = (argv[i + 1:i + 3] + ["REF", "TEST"])[:2]
        argv = argv[:i] + argv[i + 3:]
    args = [a for a in argv if not a.startswith("--")]
    if len(args) != 2:
        sys.exit(__doc__.split("USAGE")[1])

    (ra, wa), (rb, wb) = read(args[0]), read(args[1])
    if wa is None or wb is None:
        sys.exit("ERROR: one file has no '<<' developed stations. Nothing to "
                 "compare -- the run is too short for a geometry claim.")

    # Common grid from the COARSER run. See the module docstring: averaging a
    # fine result down is honest, inventing structure in a coarse one is not.
    da, db = spacing(ra), spacing(rb)
    step = max(da, db)
    lo, hi = max(wa[0], wb[0]), min(wa[1], wb[1])
    if hi - lo < 2 * step:
        sys.exit(f"ERROR: developed windows barely overlap "
                 f"({lo:.2f} to {hi:.2f} mm, {step:.2f} mm stations).\n"
                 f"  {labs[0]}: {wa[0]:.2f} to {wa[1]:.2f}\n"
                 f"  {labs[1]}: {wb[0]:.2f} to {wb[1]:.2f}\n"
                 f"  Run the shorter one further before comparing.")

    xs, x = [], lo
    while x <= hi + 1e-9:
        xs.append(round(x, 4)); x += step

    print(f"  {labs[0]}  ({da:.3f} mm stations, window {wa[0]:.2f} to {wa[1]:.2f})")
    print(f"  {labs[1]}  ({db:.3f} mm stations, window {wb[0]:.2f} to {wb[1]:.2f})")
    print(f"  compared on {len(xs)} points, {lo:.2f} to {hi:.2f} mm, "
          f"step {step:.3f} mm (from the coarser run)\n")

    hdr = (f"  {'x':>7} | {'A ref':>6} {'A test':>6} {'d%':>6} | "
           f"{'h ref':>6} {'h test':>6} {'d%':>6} | "
           f"{'w ref':>6} {'w test':>6} {'d%':>6}")
    print(hdr)
    print(f"  {'mm':>7} | {'mm^2':>6} {'mm^2':>6} {'':>6} | "
          f"{'mm':>6} {'mm':>6} {'':>6} | {'mm':>6} {'mm':>6} {'':>6}")
    print("  " + "-" * (len(hdr) - 2))

    acc = {k: [0.0, 0.0] for k in "Ahw"}
    n = 0
    for x in xs:
        a, b = interp(ra, x), interp(rb, x)
        if a is None or b is None:
            continue
        n += 1
        for k, i in zip("Ahw", (0, 1, 2)):
            acc[k][0] += a[i]; acc[k][1] += b[i]
        print(f"  {x:7.2f} | {a[0]:6.3f} {b[0]:6.3f} {100*(b[0]/a[0]-1):+6.1f} | "
              f"{a[1]:6.3f} {b[1]:6.3f} {100*(b[1]/a[1]-1):+6.1f} | "
              f"{a[2]:6.3f} {b[2]:6.3f} {100*(b[2]/a[2]-1):+6.1f}")

    if not n:
        sys.exit("ERROR: no overlapping points after interpolation.")
    m = {k: (v[0]/n, v[1]/n) for k, v in acc.items()}
    print("  " + "-" * (len(hdr) - 2))
    print(f"  {'MEAN':>7} | {m['A'][0]:6.3f} {m['A'][1]:6.3f} "
          f"{100*(m['A'][1]/m['A'][0]-1):+6.1f} | "
          f"{m['h'][0]:6.3f} {m['h'][1]:6.3f} {100*(m['h'][1]/m['h'][0]-1):+6.1f} | "
          f"{m['w'][0]:6.3f} {m['w'][1]:6.3f} {100*(m['w'][1]/m['w'][0]-1):+6.1f}")

    def shape(A, h, w):
        return h/w, h/(1.5*A/w)

    ra_hw, ra_par = shape(m['A'][0], m['h'][0], m['w'][0])
    rb_hw, rb_par = shape(m['A'][1], m['h'][1], m['w'][1])
    zh_hw = ZHAO_H/ZHAO_W

    print(f"""
  SHAPE

                       {labs[0]:>14}  {labs[1]:>14}     Zhao sim
    h/w                {ra_hw:>14.4f}  {rb_hw:>14.4f}  {zh_hw:>11.4f}
    vs Zhao            {ra_hw/zh_hw:>13.2f}x  {rb_hw/zh_hw:>13.2f}x  {'1.00x':>11}
    h/(3A/2w)          {ra_par:>14.3f}  {rb_par:>14.3f}  {ZHAO_H/(1.5*ZHAO_A/ZHAO_W):>11.3f}

    change in h/w      {100*(rb_hw/ra_hw-1):+13.1f}%
    gap to Zhao closed {100*(ra_hw-rb_hw)/(ra_hw-zh_hw):+13.1f}%

  'gap closed' is the change as a fraction of the distance to Zhao. 100% means
  the test run reaches Zhao's aspect ratio; a negative number means it moved
  away. Read it beside the A column: if A itself shifted by more than a couple
  of percent the two runs are not depositing the same metal and no shape
  comparison between them is clean, whatever the shape numbers say.""")


if __name__ == "__main__":
    main()
