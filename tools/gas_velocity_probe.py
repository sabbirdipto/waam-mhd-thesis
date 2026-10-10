#!/usr/bin/env python3
"""Where is the fast argon, and is it a jet or an artefact?

THE QUESTION THIS SETTLES

    pool_stats_waam reports 'U gas' max: 4.2-5.1 m/s in M2, 7.0-13.1 m/s in G1.
    The shielding-gas force ledger priced the jet at its INLET peak of
    1.69 m/s, giving 0.5*rho*u^2 = 2.33 Pa against 228 Pa of arc pressure --
    1%, and dismissed it. At 13.1 m/s the same expression gives 139 Pa, 61% of
    arc pressure, which would make the shielding gas a first-order mechanical
    term and G1's flattened bead a real effect rather than an artefact.

    Which reading is right depends entirely on WHERE the fast cells sit.

      BULK GAS   alpha < 0.05. Up in the argon column, or in the constriction
                 between the nozzle and the bead. This is real flow: the cell
                 really is argon, it really moves that fast, and 0.5*rho*u^2
                 really is the pressure it can exert.

      INTERFACE  0.05 <= alpha <= 0.95. The two or three cells straddling the
                 metal surface. Spurious currents from the CSF surface-tension
                 model at a 1600:1 density ratio live here and essentially
                 nowhere else. 0.5*rho_gas*u^2 is meaningless for them: the
                 cell's density is the VOF blend, not argon's, and the momentum
                 is carried by the metal fraction that is barely moving.

    M2 showed 4-5 m/s with no jet imposed at all, so some of this is certainly
    artefact. The question is whether G1's extra 8 m/s joined the artefact or
    appeared somewhere new.

WHY A SEPARATE TOOL AND NOT A COLUMN

    A maximum over all gas cells is one number with no address. This reports
    the address, the alpha, and -- more useful than any maximum -- the 99th
    percentile of each class, so one pathological cell cannot set the headline.

USAGE
    python3 gas_velocity_probe.py <case> [tmin]
"""
import gzip, math, re, sys
from pathlib import Path


def _read(case, time, name):
    for c in (case / time / name, case / time / (name + ".gz")):
        if c.exists():
            if c.suffix == ".gz":
                with gzip.open(c, "rt", errors="replace") as fh:
                    return fh.read()
            return c.read_text(errors="replace")
    return None


def _parse(text, kind):
    """OpenFOAM writes the list length on the line AFTER the type, which is
    the detail that silently breaks naive single-line regexes."""
    if text is None:
        return None
    m = re.search(r"internalField\s+uniform\s+(.+?);", text, re.S)
    if m:
        v = m.group(1).strip()
        if kind == "vector":
            return [tuple(float(x) for x in v.strip("()").split())]
        return [float(v)]
    m = re.search(r"internalField\s+nonuniform\s+List<\w+>\s*(\d+)\s*\(", text, re.S)
    if not m:
        return None
    n, body = int(m.group(1)), text[m.end():]
    cut = body.find("\n)")
    chunk = body[:cut] if cut > 0 else body
    if kind == "vector":
        vals = [tuple(float(x) for x in t.split())
                for t in re.findall(r"\(([^)]+)\)", chunk)]
    else:
        vals = [float(x) for x in chunk.split()]
    return vals[:n] if len(vals) >= n else None


def p99(v):
    if not v:
        return float("nan")
    s = sorted(v)
    return s[min(len(s) - 1, int(0.99 * len(s)))]


case = Path(sys.argv[1]).expanduser().resolve()
tmin = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0

C = _parse(_read(case, "0", "C"), "vector")
if C is None:
    sys.exit(f"  no cell centres at {case}/0/C -- run "
             f"'postProcess -func writeCellCentres -time 0'")

times = sorted((d.name for d in case.iterdir()
                if d.is_dir() and re.fullmatch(r"\d+(\.\d+)?", d.name)
                and float(d.name) > tmin), key=float)

print(f"\n  case  {case.name}   {len(C)} cells   times > {tmin}")
print("  classes: metal alpha>0.95 | interface 0.05-0.95 | bulk gas alpha<0.05\n")
print(f"    {'time':>6} {'nBulk':>7} {'nIface':>7} | {'bulk max':>9} {'bulk p99':>9}"
      f" | {'ifc max':>8} {'ifc p99':>8} |  overall max is in...")
print("    " + "-" * 100)

for t in times:
    U = _parse(_read(case, t, "U"), "vector")
    A = _parse(_read(case, t, "alpha.metal"), "scalar")
    if U is None or A is None or len(U) != len(C):
        print(f"    {t:>6}  -- U or alpha.metal unreadable / wrong size")
        continue
    if len(A) == 1:
        A = A * len(C)

    bulk, ifc, best = [], [], (0.0, None, None)
    for i, a in enumerate(A):
        u = math.sqrt(U[i][0] ** 2 + U[i][1] ** 2 + U[i][2] ** 2)
        if a < 0.05:
            bulk.append(u); cls = "bulk-gas"
        elif a <= 0.95:
            ifc.append(u); cls = "INTERFACE"
        else:
            continue                      # metal: not a gas velocity
        if u > best[0]:
            best = (u, i, (cls, a))

    i = best[1]
    if i is None:
        print(f"    {t:>6}  -- no gas cells")
        continue
    x, y, z = (1000 * v for v in C[i])
    cls, a = best[2]
    print(f"    {t:>6} {len(bulk):>7} {len(ifc):>7} |"
          f" {max(bulk) if bulk else 0:>9.3f} {p99(bulk):>9.3f} |"
          f" {max(ifc) if ifc else 0:>8.3f} {p99(ifc):>8.3f} |"
          f"  {cls:<9} alpha {a:5.3f}  x {x:7.2f}  y {y:6.2f}  z {z:6.2f} mm")

print("""
  HOW TO READ THIS

  If the overall max is INTERFACE at every row, and 'bulk p99' stays near the
  1.69 m/s inlet peak, then the 13 m/s is a CSF artefact. 0.5*rho*u^2 = 139 Pa
  does not apply to it, the force ledger's 2.33 Pa stands, and G1's flattened
  bead needs a different explanation.

  If 'bulk max' is itself 7-13 m/s, the argon genuinely accelerates -- most
  likely through the gap between the 10 mm nozzle and the rising bead -- the
  ledger understated the load by using the inlet velocity, and the shielding
  gas is a first-order mechanical term.

  Compare G1 against M2. M2 had NO jet, so whatever M2 shows in the interface
  class is the artefact floor for this mesh and this surface-tension model.
""")
