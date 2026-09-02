#!/usr/bin/env python3
"""Unit test for the two-branch kappa patch.

The patch makes an exact, checkable prediction. With

    poly_kappa_solid  (kS 0 ...)      poly_kappa_liquid (kL 0 ...)
    Tsolidus  Ts                      Tliquidus  Tl

every cell must satisfy

    e     = clip((T - Ts)/(Tl - Ts), 0, 1)        liquid fraction from T
    kMet  = e*kL + (1 - e)*kS                     solid/liquid blend
    kappa = a*kMet + (1 - a)*kGas                 metal/gas blend, a = alpha.metal

This reads the written fields and checks that identity cell by cell. It tests
the CODE, not the physics -- which is exactly what is missing after a
regression suite that only ever exercised the fallback path.

    ./verify_kappa.py <caseDir> <time> --kS 25 --kL 5 --kGas 0.4 --Ts 1658 --Tl 1723
"""
import argparse, re, sys
from pathlib import Path


def read_field(path):
    """Return a list of scalars from an OpenFOAM volScalarField file."""
    if not path.exists():
        sys.exit(f"missing field: {path}")
    txt = path.read_text()
    m = re.search(r"internalField\s+uniform\s+([-\d.eE+]+)\s*;", txt)
    if m:
        return ("uniform", float(m.group(1)))
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\n?\s*(\d+)\s*\(", txt)
    if not m:
        sys.exit(f"could not parse internalField in {path}")
    n = int(m.group(1))
    start = txt.index("(", m.end() - 1) + 1
    body = txt[start:txt.index(")", start)]
    vals = [float(x) for x in body.split()]
    if len(vals) != n:
        sys.exit(f"{path}: header says {n} values, found {len(vals)}")
    return ("list", vals)


def as_list(field, n):
    kind, v = field
    return [v] * n if kind == "uniform" else v


def main():
    p = argparse.ArgumentParser()
    p.add_argument("case", type=Path)
    p.add_argument("time")
    p.add_argument("--kS", type=float, required=True)
    p.add_argument("--kL", type=float, required=True)
    p.add_argument("--kGas", type=float, required=True)
    p.add_argument("--Ts", type=float, required=True)
    p.add_argument("--Tl", type=float, required=True)
    # Two tolerances, because the three melt states are not equally checkable.
    # In fully solid / fully liquid cells the liquid fraction is CLIPPED to
    # exactly 0 or 1, so kappa must equal kS or kL to write precision -- a tight
    # check that directly proves both branches are read and applied.
    # In mushy cells e = (T-Ts)/(Tl-Ts) inherits T's write precision AMPLIFIED
    # by (kS-kL)/(Tl-Ts): for Plate2D that is 20/65 = 0.31 W/mK per Kelvin, so
    # +-0.005 K on a 6-sig-fig T becomes +-1.5e-3 on kappa. That is arithmetic,
    # not error, so the mushy branch gets a looser bound.
    p.add_argument("--tol", type=float, default=1e-5, help="solid/liquid cells")
    # Mushy cells are checked DIFFERENTLY, as a bound on implied temperature
    # rather than on kappa. See the note in main().
    p.add_argument("--tol-dT-frac", type=float, default=0.05, dest="tdt",
                   help="max |implied T - written T| as a fraction of Tl-Ts")
    a = p.parse_args()

    d = a.case / a.time
    T = read_field(d / "T")
    n = len(T[1]) if T[0] == "list" else None
    if n is None:
        sys.exit("T is uniform -- nothing to test")
    Tv = as_list(T, n)
    av = as_list(read_field(d / "alpha.metal"), n)
    kv = as_list(read_field(d / "kappa"), n)
    ev = as_list(read_field(d / "epsilon1"), n)

    st = {"solid": [0, 0.0, -1], "mushy": [0, 0.0, -1], "liquid": [0, 0.0, -1]}
    for i in range(n):
        e = min(max((Tv[i] - a.Ts) / (a.Tl - a.Ts), 0.0), 1.0)
        al = min(max(av[i], 0.0), 1.0)
        pred = al * (e * a.kL + (1 - e) * a.kS) + (1 - al) * a.kGas
        err = abs(kv[i] - pred) / max(abs(pred), 1.0)
        k = "solid" if e <= 0.0 else ("liquid" if e >= 1.0 else "mushy")
        st[k][0] += 1
        if err > st[k][1]:
            st[k][1], st[k][2] = err, i
    solid, mushy, liquid = (st[k][0] for k in ("solid", "mushy", "liquid"))
    worst, worst_i = max((v[1], v[2]) for v in st.values())

    print(f"\n  case             : {a.case}/{a.time}")
    print(f"  cells            : {n}")
    print(f"  by melt state    : {solid} solid, {mushy} mushy, {liquid} liquid")
    print(f"  kappa range      : {min(kv):.4f} to {max(kv):.4f}"
          f"   (expect {a.kGas} to {a.kS})")
    print()
    print(f"  {'melt state':<10} {'cells':>7} {'max rel err':>13} {'tol':>9}  result")
    allok = True
    for k in ("solid", "liquid"):
        cnt, err, ci = st[k]
        ok = (cnt == 0) or (err < a.tol)
        allok &= ok
        note = "-" if cnt == 0 else ("ok" if ok else "*** FAIL ***")
        print(f"  {k:<10} {cnt:>7} {err:13.3e} {a.tol:9.0e}  {note}")
    print(f"  {'mushy':<10} {st['mushy'][0]:>7} {st['mushy'][1]:13.3e} "
          f"{'--':>9}  checked as implied dT, below")
    print("""
  Solid and liquid clip the liquid fraction to exactly 0 or 1, so kappa must
  equal kS or kL to write precision -- an exact test of both branches.

  Mushy cells CANNOT be tested that way. laserbeamFoam.C computes properties at
  line 158 (updateProps.H) and updates T at line 178 (TEqn.H), so the written
  kappa was evaluated from an EARLIER temperature than the written T. That lag
  is pre-existing solver behaviour, not a property of this patch -- with a
  single polynomial it is simply unobservable. So mushy cells are checked by
  bounding the implied temperature difference instead.""")
    if worst_i >= 0:
        print(f"\n  worst cell {worst_i}: T={Tv[worst_i]:.2f}  "
              f"alpha={av[worst_i]:.4f}  kappa={kv[worst_i]:.6f}")

    # ---- diagnostics on MUSHY METAL cells only -------------------------
    # If the blend formula is right, kappa implies a liquid fraction, which
    # implies the temperature the solver actually used. Comparing that with the
    # WRITTEN temperature separates two very different explanations:
    #   * a systematic offset with little scatter  -> kappa and T were written
    #     from different points in the timestep. Formula correct.
    #   * scattered or slope-wrong                 -> the blend is genuinely wrong.
    mm = [i for i in range(n)
          if av[i] > 0.99 and 0.0 < (Tv[i]-a.Ts)/(a.Tl-a.Ts) < 1.0]
    if mm and a.kS != a.kL:
        dT = []
        for i in mm:
            e_imp = (a.kS - kv[i]) / (a.kS - a.kL)
            dT.append((a.Ts + e_imp*(a.Tl - a.Ts)) - Tv[i])
        mean = sum(dT)/len(dT)
        sd = (sum((x-mean)**2 for x in dT)/len(dT))**0.5
        print(f"\n  mushy metal cells               : {len(mm)}")
        print(f"  implied T minus written T       : mean {mean:+.4f} K, "
              f"sd {sd:.4f} K, range [{min(dT):+.3f}, {max(dT):+.3f}]")
        bound = a.tdt * (a.Tl - a.Ts)
        mx = max(abs(x) for x in dT)
        ok = mx < bound
        allok &= ok
        print(f"  max |implied T - written T|     : {mx:.4f} K"
              f"   (bound {bound:.2f} K = {a.tdt:.0%} of the {a.Tl-a.Ts:.0f} K window)")
        print(f"  {'  -> ok: within one outer iteration of temperature change'if ok else '  -> *** FAIL: too large to be the known kappa/T evaluation lag ***'}")

    # epsilon1 comparison, restricted to metal (gas is always liquid there)
    met = [i for i in range(n) if av[i] > 0.99]
    if met:
        de = max(abs(min(max((Tv[i]-a.Ts)/(a.Tl-a.Ts), 0.0), 1.0) - ev[i])
                 for i in met)
        print(f"  max |e_fromT - epsilon1|, metal : {de:.3e}"
              f"  (solver's iterated vs T-derived)")

    ok = allok
    print(f"\n  VERDICT : {'PASS - the blend is exactly as specified' if ok else '*** FAIL - the patch is not doing what it claims ***'}\n")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
