#!/usr/bin/env python3
"""Measure the weld pool in waam-arc-01. Works on a GRADED mesh.

WHAT CHANGED IN THIS REVISION (rev 3, after the 0.1 s droplet run)

  1. 'max T' USED TO REPORT THE ARGON.
     The old column maximised T over EVERY cell in the domain. With droplets
     running, the gas above the impact site gets hot: at t = 0.040 s this
     column read 1277.8 K while the hottest METAL cell was 872.5 K. A 400 K
     error, and it fell on the one number the runaway gate watches.
     It is now split into 'T metal' and 'T gas', the same way 'U metal' and
     'U gas' were split for exactly the same reason one revision earlier.
     The gas is not noise -- it is where vaporisation and arc pressure will
     act later -- but it is not the pool, and it must not be read as the pool.

  2. Fields may be GZIPPED.
     writeCompression is now on, because a 0.5 s run at writeInterval 0.01
     writes 50 times x ~15 MB. OpenFOAM reads 'T.gz' transparently and so,
     now, does this script.

  3. The PARSER no longer walks the file character by character.
     The old version matched parentheses one char at a time over a ~1.9 MB
     string to find the end of the list. That is ~1-2 s per field in Python;
     across 3 fields x 50 times it was five minutes of pure bookkeeping.
     The list length is DECLARED in the file, so the closing paren does not
     need finding: read exactly that many tokens and stop.

WHAT IT REPORTS, per written time
    T metal   hottest metal cell. THIS is the runaway gate.
    T gas     hottest argon cell. Expected to exceed T metal near a droplet;
              it is the plume, not the pool.
              ('mean T' was dropped in rev 4 to make room. It measured the
              volume-weighted metal temperature, which over this domain is
              dominated by the cold bulk -- it moved 300.3 -> 323.1 K across
              the entire 0.5 s run -- and it says the same thing as 'E kept'
              in less interpretable units.)
    bottom    hottest cell on the substrate underside. While this sits at
              300 K the adiabatic bottom wall is free. Once it climbs, the
              missing convective/radiative condition (Zhao Eq. 29, h = 80)
              starts to matter, because heat that should be leaving is not.
    far met   hottest METAL cell in the outermost y layer, against farField.
              THIS is the domain-independence readout: while it stays at
              300 K the far boundary condition cannot matter, because there
              is no gradient there for it to act on.
    far gas   hottest ARGON cell in the same layer. Split out in rev 4 after
              the 0.5 s run, where the single combined column jumped from
              300.0 to 321.4 at t = 0.12 and then wandered 308 -> 311 -> 347
              -> 374 -> 426 -> 415. Conduction into a far boundary is smooth
              and monotonic; that trace is not conduction, it is the argon
              plume being blown 25 mm sideways at 4-10 m/s and arriving as
              hot gas. Reading it as a metal temperature would retire the
              domain-independence claim on the strength of the shielding
              gas. It is still a real signal -- hot argon piling against a
              wall it cannot leave through drives spurious buoyancy -- but
              it is a gas-domain question, not a substrate one.
    pool      cells above the SOLIDUS, and their real volume and depth.
    U metal   max |U| in the metal. Zero whenever nothing is molten -- the
              Darcy source drives velocity to machine zero in solid cells,
              so this column is a liquid-fraction detector as much as a
              velocity one.
    E kept    thermal energy in the metal as a percentage of ARC energy
              delivered. With droplets this EXCEEDS 100 and should: each
              droplet carries ~16.5 J of sensible enthalpy and up to ~5.9 J
              of latent heat, against 1512 W x 33.3 ms = 50 J of arc energy
              per droplet interval. Roughly a third again on top.

BEFORE FIRST USE, in the case directory:
    postProcess -func writeCellCentres -time 0
    postProcess -func writeCellVolumes -time 0

    These write Cx, Cy, Cz and V into 0/. They describe the mesh, not the
    solution, so they only need generating once per mesh.

USAGE
    cd ~/thesis/cases/00-tutorials/waam-arc-01
    python3 pool_stats_waam.py
"""
import gzip, math, re, sys
from pathlib import Path

CASE = Path.cwd()

TSOL, TLIQ = 815.0, 906.0      # materials/transportProperties
RHO, CP    = 2650.0, 940.0     # cp from closing the energy balance on the
                               # uniform run: 1000/1.064 = 940 J/kg/K
T0         = 300.0


def _arc_power():
    """Arc power actually delivered, read from the case.

    This was hardcoded at 756.0 -- the HALF-domain figure -- and a full-domain
    run delivers 1512 W, so 'E kept' came out at 204% instead of 102% and
    looked like an energy conservation failure when nothing was wrong. Read
    it from constant/caseParameters, the same file the solver reads.
    """
    f = CASE / "constant" / "caseParameters"
    if not f.exists():
        return 756.0, "assumed half domain (no caseParameters)"
    txt = f.read_text()
    m = re.search(r"^\s*arc_power_W\s+([0-9eE+.-]+)\s*;", txt, re.M)
    if not m:
        return 756.0, "assumed (no arc_power_W)"
    p = float(m.group(1))
    s = re.search(r"^\s*sym\s+(\w+)\s*;", txt, re.M)
    half = bool(s and s.group(1) == "half")
    return p*(0.5 if half else 1.0), f"caseParameters, {'half' if half else 'full'} domain"


ARC_W, ARC_SRC = _arc_power()


def slurp(path):
    """Text of an OpenFOAM field, whether or not writeCompression is on.

    OpenFOAM appends '.gz' to the WHOLE name -- 'alpha.metal.gz', not
    'alpha.gz' -- so this cannot use Path.with_suffix, which would eat the
    '.metal' part and look for a file that never existed.
    """
    if path.exists():
        return path.read_text()
    gz = Path(str(path) + ".gz")
    if gz.exists():
        with gzip.open(gz, "rt") as f:
            return f.read()
    return None


def read_field(path):
    """internalField as a list of floats, or of 3-tuples for a vector.

    The count is declared in the file, immediately before the opening paren:

        internalField   nonuniform List<scalar>
        143616
        (
        3.0000000e+02
        ...

    so the parser reads exactly that many whitespace tokens and stops. It
    never looks for the closing paren, which is what made the old version
    slow. str.split(None, k) performs at most k splits and returns k+1
    fields, the last being the unparsed remainder -- so taking the first k
    gives clean tokens with the trailing ')' and everything after it left
    untouched in the discarded tail.
    """
    txt = slurp(path)
    if txt is None:
        return None
    m = re.search(r"internalField\s+nonuniform\s+List<(scalar|vector)>\s*\n?\s*(\d+)\s*\(",
                  txt)
    if m:
        kind, n = m.group(1), int(m.group(2))
        tail = txt[m.end():]
        if kind == "scalar":
            toks = tail.split(None, n)[:n]
            if len(toks) != n:
                sys.exit(f"ERROR: {path.name} declared {n} values, found {len(toks)}")
            return [float(v) for v in toks]
        # vectors are written '(a b c)', so three whitespace tokens each,
        # the first carrying a leading '(' and the third a trailing ')'
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
    sys.exit(f"ERROR: could not parse internalField in {path}")


def times():
    out = []
    for d in CASE.iterdir():
        if d.is_dir():
            try:
                out.append((float(d.name), d))
            except ValueError:
                pass
    return sorted(out)


def find_mesh():
    """Cell centres and volumes, from OpenFOAM. No fallback by design."""
    for _, d in times():
        cx, cy, cz, v = (read_field(d / n) for n in ("Cx", "Cy", "Cz", "V"))
        if cx and cy and cz and v and len(cx) == len(v) > 1:
            return list(zip(cx, cy, cz)), v, d.name
    sys.exit(
        "ERROR: no cell centres / volumes found in any time directory.\n\n"
        "  This script will not guess the mesh. On a graded mesh the cells\n"
        "  differ in volume by a factor of ~16, so assuming uniformity would\n"
        "  give confidently wrong depths and pool volumes.\n\n"
        "  Generate them once for this mesh:\n"
        "      postProcess -func writeCellCentres -time 0\n"
        "      postProcess -func writeCellVolumes -time 0")


def main():
    ts = times()
    if not ts:
        sys.exit("ERROR: no time directories here. Run this inside the case.")

    C, V, where = find_mesh()
    n = len(V)
    print(f"  arc:  {ARC_W:g} W  [{ARC_SRC}]")
    print(f"  mesh: {n:,} cells, geometry read from {where}/  "
          f"(cell volume {min(V)*1e9:.3g} to {max(V)*1e9:.3g} mm^3, "
          f"ratio {max(V)/min(V):.0f})")

    # Sanity check the geometry against the physics we know: at the first
    # time, every cell below z = 0 must be metal and every cell above must
    # not be. Catches a mismatched mesh before any number is reported.
    a0 = read_field(ts[0][1] / "alpha.metal")
    if a0 and len(a0) == n:
        bad = sum(1 for c in range(n) if (C[c][2] < 0.0) != (a0[c] > 0.5))
        if bad:
            sys.exit(f"ERROR: {bad} cells disagree with the substrate at z < 0. "
                     f"The cell centres do not match this solution.")
        print(f"  geometry check passed on all {n:,} cells\n")
    else:
        print("  (could not cross-check geometry against alpha.metal)\n")

    zmin = min(c[2] for c in C)

    # The outermost cell layer against each patch, by cell-centre coordinate.
    # Taken from the mesh rather than from the nominal domain size, because on
    # a graded mesh the last centre sits half of a LARGE cell inside the face
    # and that offset is not something to guess.
    ymax = max(c[1] for c in C)
    bottom_cells = [c for c in range(n) if C[c][2] < zmin + 1e-9]
    far_cells    = [c for c in range(n) if C[c][1] > ymax - 1e-9]

    hdr = (f"  {'time':>7}  {'T metal':>8}  {'T gas':>8}  {'bottom':>8}  "
           f"{'far met':>8}  {'far gas':>8}  {'pool':>6}  {'volume':>8}  {'depth':>6}  "
           f"{'U metal':>7}  {'U gas':>7}  {'E kept':>7}")
    print(hdr)
    print(f"  {'s':>7}  {'K':>8}  {'K':>8}  {'K':>8}  "
          f"{'K':>8}  {'K':>8}  {'cells':>6}  {'mm^3':>8}  {'mm':>6}  "
          f"{'m/s':>7}  {'m/s':>7}  {'%':>7}")
    print("  " + "-"*len(hdr.strip()))

    for t, d in ts:
        T = read_field(d / "T")
        if T is None or len(T) != n:
            continue
        a = read_field(d / "alpha.metal")
        if a is None or len(a) != n:
            a = a0
        U = read_field(d / "U")

        metal  = [c for c in range(n) if a[c] > 0.5]
        gas    = [c for c in range(n) if a[c] <= 0.5]
        molten = [c for c in metal if T[c] > TSOL]

        vol = sum(V[c] for c in molten)*1e9
        depth = -min(C[c][2] for c in molten)*1e3 if molten else 0.0

        Ein = RHO*CP*sum((T[c] - T0)*V[c] for c in metal)
        frac = 100.0*Ein/(ARC_W*t) if t > 0 else float("nan")

        # 'bottom' needs no phase split: the plate underside at z = zmin is
        # entirely substrate, so every cell in that layer is metal already.
        # The far-y layer is NOT -- it spans the full z range, so most of it
        # is argon, and the argon there is downstream of the plume.
        botT = max((T[c] for c in bottom_cells), default=float("nan"))
        farM = max((T[c] for c in far_cells if a[c] > 0.5), default=float("nan"))
        farG = max((T[c] for c in far_cells if a[c] <= 0.5), default=float("nan"))

        # Split EVERY extensive maximum by phase, not just velocity. Buoyant
        # argon over the arc is hotter and an order of magnitude faster than
        # anything in the metal, and reporting either as a pool quantity is
        # wrong in the flattering direction -- it looks like a vigorous,
        # superheated pool when the metal has not melted at all.
        def maxT_of(cells):
            return max((T[c] for c in cells), default=float("nan"))

        def umax_of(cells):
            if not (U and len(U) == n):
                return float("nan")
            return max((math.sqrt(U[c][0]**2 + U[c][1]**2 + U[c][2]**2)
                        for c in cells), default=0.0)

        print(f"  {d.name:>7}  {maxT_of(metal):8.1f}  {maxT_of(gas):8.1f}  "
              f"{botT:8.1f}  {farM:8.1f}  {farG:8.1f}  "
              f"{len(molten):6d}  {vol:8.2f}  {depth:6.2f}  "
              f"{umax_of(metal):7.4f}  {umax_of(gas):7.4f}  {frac:7.1f}")

    print(f"""
  GATE:
    T metal below ~1600 K ........ vaporisation is disabled, so a runaway has
                                   no ceiling. Read the METAL column: T gas
                                   legitimately spikes near a droplet.
    pool volume > 0 AND CONTINUOUS  a pool that appears at each droplet and
                                   vanishes between them is not a weld pool,
                                   it is three separate solidification events.
                                   Every row must be non-zero, not just some.
    depth > 0.5 mm ............... more than one cell layer of the substrate
                                   is fused. At exactly 0.25 mm the pool is
                                   one cell deep and the number is the mesh
                                   talking, not the physics.
    depth < 6 mm ................. the pool has not reached the plate bottom
    bottom near 300 K ............ while true, the adiabatic bottom wall is
                                   free. Once it climbs the missing Robin
                                   condition (h = 80) is holding in heat that
                                   should be leaving.
    far met = 300.0 .............. the domain is still big enough. This is
                                   the gate; 'far gas' is not. Argon is
                                   advected, so it reports where the plume
                                   went, not how far heat has conducted.

  Al-5Mg solidus {TSOL:g} K, liquidus {TLIQ:g} K. 'Pool' is above the SOLIDUS,
  i.e. mushy plus fully liquid -- the region that is no longer solid.""")


if __name__ == "__main__":
    main()
