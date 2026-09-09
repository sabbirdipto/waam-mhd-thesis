#!/usr/bin/env python3
"""Did the droplet mass actually arrive, and did all of it stay? Mass only.

WHAT CHANGED IN THIS REVISION (rev 2, after the 0.1 s droplet run)

  1. THE OLD 'metal all' COLUMN MEASURED A THRESHOLD, NOT A MASS.
     It summed the WHOLE volume of every cell with alpha > 0.5. That is not
     the metal volume; it is the volume of the cells that happen to be more
     than half metal. The two differ by whatever is in the interface cells,
     and the interface is exactly where a droplet lives.

     The 0.1 s run showed the consequence. After droplet 1 the column read
     18005.5 mm^3; five milliseconds later, 18004.0. One and a half cubic
     millimetres of aluminium apparently evaporated. Nothing of the sort
     happened -- the droplet flattened on impact, its interface smeared over
     more cells, and a ring of cells fell from alpha = 0.55 to alpha = 0.45
     and stopped being counted at all, while the metal they contained was
     still there.

     The honest measure is the integral the VOF method actually conserves:

         V_metal = sum over cells of  alpha_c * V_c

     This is what isoAdvector/MULES bound and conserve to machine precision,
     so it is the number that can legitimately be held to a conservation
     test. It is reported as 'aV'. The thresholded figure is kept beside it
     as 'a>.5' precisely so the gap between them stays visible: that gap is
     interface smearing, and it growing over time is a real signal about
     advection quality even though it is not mass loss.

  2. THE SUBSTRATE IS NO LONGER HARDCODED AT 18,000 mm^3.
     That figure was a full-domain, 60 x 50 x 6 mm assumption, and it was
     already wrong once for the half domain. The initial time IS the
     substrate, so the script measures it and reports DEPOSITED metal --
     the difference -- directly, instead of asking the reader to subtract
     five significant figures in their head.

  3. TIMES ARE PRINTED AS THE DIRECTORY NAME.
     '%7.3f' rendered the 0.0035 s verification run as '0.003' and '0.004',
     collapsing the injection step into the step before it. The directory
     name is the exact time OpenFOAM wrote; there is no reason to reformat
     it and several ways to do so wrongly.

  4. Gzipped fields are read, and the parser no longer walks parentheses
     one character at a time. See the notes in pool_stats_waam.py.

USAGE
    cd <case>
    python3 check_droplet.py           # every written time
    python3 check_droplet.py 0.04      # one time
"""
import gzip, math, re, sys
from pathlib import Path

CASE = Path.cwd()
TSOL, TLIQ = 815.0, 906.0


def slurp(path):
    """Text of an OpenFOAM field, gzipped or not.

    OpenFOAM appends '.gz' to the whole filename -- 'alpha.metal.gz' -- so
    Path.with_suffix is the wrong tool: it would replace '.metal'.
    """
    if path.exists():
        return path.read_text()
    gz = Path(str(path) + ".gz")
    if gz.exists():
        with gzip.open(gz, "rt") as f:
            return f.read()
    return None


def read_field(path):
    """internalField as floats, or 3-tuples for a vector.

    The list length is declared in the file, so the closing paren need not be
    hunted for. str.split(None, k) does at most k splits and returns k+1
    fields; the first k are clean tokens and the trailing ')' is left in the
    discarded remainder.
    """
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
        if len(toks) != 3*n:
            sys.exit(f"ERROR: {path.name} declared {n} vectors, found {len(toks)//3}")
        f = [float(t.strip("()")) for t in toks]
        return list(zip(f[0::3], f[1::3], f[2::3]))
    m = re.search(r"internalField\s+uniform\s+([^;]+);", txt)
    if m:
        tok = m.group(1).strip()
        return [tuple(float(x) for x in tok.strip("()").split())] if tok.startswith("(") \
               else [float(tok)]
    return None


def times():
    out = []
    for d in CASE.iterdir():
        if d.is_dir():
            try:
                out.append((float(d.name), d))
            except ValueError:
                pass
    return sorted(out)


def nominal_droplet():
    """Volume of one droplet, mm^3, from the case the solver reads.

    Same source and same arithmetic as updateDropletSource.H: Zhao Eq. 30
    radius, halved with the domain if the case is symmetric. Returned as
    None rather than guessed if the entry is absent, because a wrong
    denominator here would silently rescale the droplet count.
    """
    f = CASE / "constant" / "caseParameters"
    if not f.exists():
        return None
    txt = f.read_text()
    m = re.search(r"^\s*droplet_radius_m\s+([0-9eE+.-]+)\s*;", txt, re.M)
    if not m:
        return None
    r = float(m.group(1))
    s = re.search(r"^\s*sym\s+(\w+)\s*;", txt, re.M)
    sym = 0.5 if (s and s.group(1) == "half") else 1.0
    return (4.0/3.0)*math.pi*r**3*sym*1e9


def main():
    want = float(sys.argv[1]) if len(sys.argv) > 1 else None
    ts = times()
    if not ts:
        sys.exit("ERROR: no time directories here. Run this inside the case.")

    mesh = None
    for _, d in ts:
        cz, v = read_field(d / "Cz"), read_field(d / "V")
        if cz and v and len(cz) == len(v) > 1:
            mesh = (cz, v)
            break
    if mesh is None:
        sys.exit("ERROR: no Cz / V. Run  postProcess -func writeCellCentres -time 0\n"
                 "                   and postProcess -func writeCellVolumes -time 0")
    Cz, V = mesh
    n = len(V)

    # The substrate, measured rather than assumed. The first written time is
    # the initial condition: nothing has been deposited yet, so whatever
    # metal exists then is the plate.
    a_init = read_field(ts[0][1] / "alpha.metal")
    if not a_init or len(a_init) != n:
        sys.exit(f"ERROR: cannot read alpha.metal at t = {ts[0][1].name}")
    SUB = sum(a_init[c]*V[c] for c in range(n))*1e9

    VD = nominal_droplet()
    print(f"  substrate at t = {ts[0][1].name}: {SUB:,.1f} mm^3  (measured, sum alpha*V)")
    if VD:
        print(f"  nominal droplet: {VD:.3f} mm^3  [caseParameters droplet_radius_m]")
    else:
        print("  nominal droplet: unknown (no droplet_radius_m) -- 'drops' omitted")
    print()

    hdr = (f"  {'time':>8}  {'bead a>.5':>9}  {'bead aV':>8}  {'T bead':>8}  "
           f"{'deposited':>9}  {'drops':>6}  {'T metal':>8}  {'>sol':>5}  {'>liq':>5}")
    print(hdr)
    print(f"  {'s':>8}  {'mm^3':>9}  {'mm^3':>8}  {'K':>8}  "
          f"{'mm^3':>9}  {'':>6}  {'K':>8}  {'cells':>5}  {'cells':>5}")
    print("  " + "-"*len(hdr.strip()))

    for t, d in ts:
        if want is not None and abs(t - want) > 1e-9:
            continue
        T = read_field(d / "T")
        a = read_field(d / "alpha.metal")
        if not T or not a or len(T) != n or len(a) != n:
            continue

        # ABOVE z = 0 is where a droplet lands and where the bead grows. The
        # substrate is entirely below it, so any metal here arrived.
        up_thr = sum(V[c] for c in range(n) if Cz[c] > 0.0 and a[c] > 0.5)*1e9
        up_aV  = sum(a[c]*V[c] for c in range(n) if Cz[c] > 0.0)*1e9
        upcells = [c for c in range(n) if Cz[c] > 0.0 and a[c] > 0.5]
        T_bead = max((T[c] for c in upcells), default=float("nan"))

        dep = sum(a[c]*V[c] for c in range(n))*1e9 - SUB
        drops = f"{dep/VD:6.2f}" if VD else "     -"

        met = [c for c in range(n) if a[c] > 0.5]
        T_met = max((T[c] for c in met), default=float("nan"))
        nsol = sum(1 for c in met if T[c] > TSOL)
        nliq = sum(1 for c in met if T[c] > TLIQ)

        print(f"  {d.name:>8}  {up_thr:9.2f}  {up_aV:8.2f}  {T_bead:8.1f}  "
              f"{dep:9.3f}  {drops}  {T_met:8.1f}  {nsol:5d}  {nliq:5d}")

    print(f"""
  HOW TO READ THIS

  bead a>.5   Volume of the CELLS that are more than half metal above z = 0.
              This is the old 'metal>0' number, kept only for comparison.
  bead aV     sum(alpha*V) over the same region -- the actual metal there.
              a>.5 BELOW aV means the bead is spread thinly across many
              partly-filled cells; the gap between the two columns is
              interface smearing, and it should not grow without bound.

  deposited   sum(alpha*V) over the whole domain, MINUS the substrate
              measured at t = {ts[0][1].name}. This is the conservation test, and it
              is the only mass column that can fail honestly: VOF conserves
              sum(alpha*V) to machine precision, so a real fall here means
              mass is leaving through a boundary or being clipped.

  drops       'deposited' in units of one nominal droplet. Should climb in
              near-integer steps: 1.00, 2.00, 3.00. A steady shortfall (0.95
              each) means cells are being clipped at placement -- the
              injector reports this itself as 'placed X of Y nominal'. A
              value that FALLS is mass loss.

  T bead      Hottest cell in the bead. Injected at droplet_T_K and cooling
              fast: a droplet meeting a plate of its own alloy drops to the
              mean of the two temperatures almost instantly, because with
              equal thermal effusivities the contact temperature is
              (T_drop + T_plate)/2 with no material contrast to slow it.

  >sol/>liq   Metal cells above solidus / liquidus. For a real weld pool
              these must be non-zero at EVERY row. Non-zero only on the rows
              next to an injection means each droplet freezes before the
              next arrives, which is deposition without a pool.""")


if __name__ == "__main__":
    main()
