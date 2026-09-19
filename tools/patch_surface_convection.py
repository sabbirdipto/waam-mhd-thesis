#!/usr/bin/env python3
"""The convective half of Zhao Eq. 29 on the free surface. Radiation is already there.

WHAT IS MISSING AND WHAT IS NOT

    Zhao Eq. 29 is two terms, both CSF-converted to volumetric sources:

        qloss = -[ h(T - Ta) + eps sigma (T^4 - Ta^4) ] |grad F| 2rho/(rho_m+rho_g)
                  \_______/   \___________________/
                  CONVECTION      RADIATION

    TEqn.H already implements the RADIATION half, correctly and with a proper
    implicit/explicit split:

        Qrad_expl =  SB*em*(TRef^4 + 3*T^4)
        Qrad_inpl = -4*SB*em*T^3
        ... - fvm::Sp(Qrad_inpl*gradAlphaDamper, T) == ... + Qrad_expl*gradAlphaDamper

    which collapses to +SB*em*(TRef^4 - T^4)*gradAlphaDamper -- a loss above
    ambient, as it should be. gradAlphaDamper is mag(gradAlpha)*thermalDamper,
    i.e. Zhao's |grad F| times the density weighting. Nothing wrong with it.

    The CONVECTION half is simply absent. Nowhere in the solver is there an
    h(T - Ta)|grad F| term.

    THIS PATCH ADDS ONLY THAT. An earlier draft of a heat-loss patch in this
    project implemented the WHOLE of Eq. 29 as a new source; applying that one
    would double-count radiation. Do not use both.

HOW BIG, HONESTLY

    At the pool surface the two terms are not close:

        T = 1200 K, em = 0.15   radiation  17,600 W/m2   convection  72,000 W/m2
        T =  900 K, em = 0.20   radiation   7,300 W/m2   convection  48,000 W/m2
        T =  800 K, em = 0.30   radiation   6,800 W/m2   convection  40,000 W/m2

    So the model is currently capturing roughly a sixth of the surface loss.
    But the hot surface is only about 170 mm^2, which makes the missing term
    about 8 W against 1512 W of input -- 0.5% of the energy budget.

    This is therefore a CORRECTNESS fix, not a lever. It will not change the
    bead. It is worth doing because the model claims to implement Eq. 29 and
    at present implements half of it, and because that half-implementation
    means the model UNDER-cools -- which matters when arguing about why the
    pool is too small, since it removes "we are cooling too fast" from the
    list of explanations.

FULLY IMPLICIT, BECAUSE IT CAN BE

    h(T - Ta) is linear in T, so it splits exactly with no linearisation:

        LHS   + fvm::Sp(h*gradAlphaDamper, T)
        RHS   + h*Ta*gradAlphaDamper

    No Picard iteration, no timestep restriction, no error. The radiation term
    beside it needs the T^4 linearisation; this one does not, and writing it
    explicitly to match the neighbour's style would import a stability limit
    for no reason.

DIMENSIONS
    h [W/m2K] * gradAlphaDamper [1/m] = W/m3K, times T gives W/m3, which is
    what every other term in TEqn is. Checked rather than assumed: the
    dimensionedScalar below carries [1 0 -3 -1 0 0 0] and OpenFOAM will refuse
    the matrix if that is wrong.

SAFETY
    h defaults to ZERO, so the term vanishes and every existing case is
    bit-identical until h_conv_Wm2K is set. Gate A and Gate B check that.

USAGE
    python3 patch_surface_convection.py            # apply
    python3 patch_surface_convection.py --revert
"""
import sys, shutil
from pathlib import Path

SOLVER = Path.home() / "thesis/src/lbf-com/applications/solvers/laserbeamFoam"
CF, TE = SOLVER / "createFields.H", SOLVER / "TEqn.H"
SUFFIX = ".cbak"

# ---------------------------------------------------------------- createFields
# TRef is the ambient reference the radiation term already uses. Anchoring here
# keeps the two halves of Eq. 29 declared together, and avoids the emS anchor
# that patch_marangoni_mushy.py attaches to.
CF_ANCHOR = 'const dimensionedScalar TRef ("TRef",transportProperties);'

CF_ADD = r'''

// --- Convective half of Zhao Eq. 29, free surface (WAAM-MHD thesis) --------
// The radiation half is already in TEqn.H via Qrad_expl / Qrad_inpl. This is
// the h(T - Ta) term, which was absent.
//
// DEFAULT ZERO: absent from transportProperties, the term vanishes and every
// existing case is bit-identical. That is what Gate A and Gate B check.
//
// Ta deliberately reuses TRef rather than introducing a second ambient
// temperature. Two names for one physical constant is how they drift apart.
const dimensionedScalar hConvSurf
(
    "hConvSurf",
    dimensionSet(1, 0, -3, -1, 0, 0, 0),
    transportProperties.lookupOrDefault<scalar>("h_conv_Wm2K", 0.0)
);'''

# ---------------------------------------------------------------- TEqn
TE_OLD = '''          - fvm::Sp(Qrad_inpl*gradAlphaDamper, T)
          ==
            laser.deposition()
          + arcSource
          + Joulesource
          + (Qrad_expl - Qv)*gradAlphaDamper'''

TE_NEW = '''          - fvm::Sp(Qrad_inpl*gradAlphaDamper, T)
          // Convective half of Eq. 29. Linear in T, so it goes in fully
          // implicit and exactly -- unlike the T^4 radiation beside it, which
          // has to be linearised. Zero unless h_conv_Wm2K is set.
          + fvm::Sp(hConvSurf*gradAlphaDamper, T)
          ==
            laser.deposition()
          + arcSource
          + Joulesource
          + (Qrad_expl - Qv)*gradAlphaDamper
          + hConvSurf*TRef*gradAlphaDamper'''

EDITS = [(CF, CF_ANCHOR, CF_ANCHOR + CF_ADD),
         (TE, TE_OLD, TE_NEW)]


def revert():
    n = 0
    for path, _, _ in EDITS:
        bak = path.with_suffix(path.suffix + SUFFIX)
        if bak.exists():
            shutil.copy2(bak, path); bak.unlink(); n += 1
            print(f"  reverted {path.name}")
    print(f"  {n} file(s) restored" if n else "  nothing to revert")


def apply():
    for path, old, _ in EDITS:
        if not path.exists():
            sys.exit(f"ERROR: {path} not found")
        text = path.read_text()
        if "hConvSurf" in text:
            sys.exit(f"ERROR: {path.name} looks already patched. --revert first.")
        if text.count(old) != 1:
            hint = ""
            if path is TE and "arcSource" not in text:
                hint = ("\n  The anchor includes arcSource. Apply "
                        "patch_goldak.py first.")
            sys.exit(f"ERROR: {path.name}: expected one match for:\n"
                     f"-----\n{old}\n-----\nfound {text.count(old)}. "
                     f"Nothing modified." + hint)

    for path, old, new in EDITS:
        shutil.copy2(path, path.with_suffix(path.suffix + SUFFIX))
        path.write_text(path.read_text().replace(old, new, 1))
        print(f"  patched {path.name}   (backup: {path.name}{SUFFIX})")

    print("""
  BUILD
    of2506 && cd ~/thesis/src/lbf-com && ./Allwmake > /tmp/b.log 2>&1
    grep -c 'error:' /tmp/b.log

  REGRESSION -- h defaults to 0, so this MUST be bit-identical.

  ENABLE, in constant/transportProperties (top level, beside TRef):
    h_conv_Wm2K   80;      // Zhao Table 4

  DO NOT also apply the earlier whole-Eq.29 heat-loss patch. TEqn.H already
  carries the radiation half; that one would count it twice.

  EXPECT NO VISIBLE CHANGE. This is ~8 W against 1512 W. If the bead or the
  pool moves noticeably, something is wrong -- most likely the term has
  landed on the gas side, which would show up as a collapse in 'T gas'.""")


if __name__ == "__main__":
    revert() if "--revert" in sys.argv else apply()
