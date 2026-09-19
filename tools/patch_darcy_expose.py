#!/usr/bin/env python3
"""Make the mushy-zone Darcy constants case parameters instead of literals.

THE PROBLEM: THE BIGGEST LEVER ON MUSHY-ZONE FLOW IS NOT EXPOSED

    createFields.H hardcodes both constants of the enthalpy-porosity momentum
    sink, read from no dictionary anywhere in the solver:

        const dimensionedScalar DarcyConstantlarge
            ("DarcyConstantlarge", dimensionSet(1,-3,-1,0,0), 1.0e6);
        const dimensionedScalar DarcyConstantsmall
            ("DarcyConstantsmall", dimless, 1.0e-12);

    and UEqn.H uses them as

        DC = A (1 - epsilon1)^2 / (epsilon1^3 + q),    fvm::Sp(DC, U)

    Zhao Eq. 7 is the same expression. He does not state his value of A. He
    states only the range: "values between 10^4 and 10^7 are commonly
    recommended". That is THREE ORDERS OF MAGNITUDE, and A sets, linearly, how
    far metal in the mushy zone can travel before it freezes.

    Balancing the Marangoni body force against the Darcy sink and integrating
    over the ~0.73 s a cross-section stays above the solidus after the arc
    passes:

        eps1   T(K)     A = 1e5      A = 1e6 (current)    A = 1e7
        0.70   878.7    604 mm       60 mm                6.0 mm
        0.50   860.5     57 mm        5.7 mm              0.57 mm
        0.30   842.3    3.7 mm        0.37 mm             0.04 mm
        0.10   824.1   0.03 mm        0.003 mm            0.000 mm

    (The eps1 > 0.9 rows are not achievable -- there the liquid is inertia and
    viscosity limited, and the measured U metal is 0.2-0.5 m/s. The rows that
    decide bead width are 0.3 to 0.7.)

    So against a width deficit of 2.3 mm, a choice of A that Zhao never pins
    down moves the mush by anywhere from 0.04 mm to 3.7 mm at the same point
    in the freezing range. It is a free parameter of the same order as the
    thing being measured, and at present it cannot be changed without
    recompiling.

    THE POINT IS NOT THAT 1e6 IS WRONG. It is in Zhao's range and is a
    defensible default. The point is that a parameter this influential must be
    in the case file, recorded per run, and swept -- not a literal in a header.

A SECOND, SMALLER FINDING: q IS NOT ZHAO'S, AND IT DOES NOT MATTER

    Zhao Eq. 7 uses eps = 0.001 "to prevent division by 0". laserbeamFoam uses
    1e-12. That is a factor of 10^9, and materials/alMg5Properties even
    carries 'darcyEpsilon 1.0e-3;  // - [ZHAO] Eq. 7' -- which nothing reads.
    Dead configuration of exactly the kind Sec 7.3 warns about.

    Checked rather than assumed: q only enters where eps1^3 is comparable to
    it, i.e. eps1 below about 0.1, and at eps1 = 0.1 the metal travels 0.003 mm
    with q = 1e-12 and 0.006 mm with q = 1e-3. Both are frozen. The mismatch
    is immaterial and should be RECORDED, not corrected -- but it should stop
    being written into a file that implies it is in force.

SAFETY
    Defaults are exactly the current literals, so with no keys added every
    case is bit-identical. Gate A and Gate B check that.

USAGE
    python3 patch_darcy_expose.py            # apply
    python3 patch_darcy_expose.py --revert
"""
import sys, shutil
from pathlib import Path

SOLVER = Path.home() / "thesis/src/lbf-com/applications/solvers/laserbeamFoam"
CF = SOLVER / "createFields.H"
SUFFIX = ".dkbak"

OLD = '''const dimensionedScalar DarcyConstantlarge
(
    "DarcyConstantlarge", dimensionSet(1,-3,-1,0,0), 1.0e6
);

const dimensionedScalar DarcyConstantsmall
(
    "DarcyConstantsmall", dimless, 1.0e-12
);'''

NEW = '''// --- Mushy-zone Darcy constants (WAAM-MHD thesis) --------------------------
// Zhao Eq. 7 is this same sink. He never states A, only that "values between
// 10^4 and 10^7 are commonly recommended" -- three orders of magnitude, over
// which the distance metal at epsilon1 = 0.3 travels before freezing runs
// from 0.04 mm to 3.7 mm. Against a 2.3 mm width deficit that is a free
// parameter the size of the measurement, so it belongs in the case file and
// in the per-run record, not in a header.
//
// DEFAULTS ARE THE PREVIOUS LITERALS, so with no keys present every existing
// case is bit-identical. That is what Gate A and Gate B check.
//
// Note on darcy_q: Zhao Eq. 7 uses 0.001, laserbeamFoam 1e-12. It only enters
// where epsilon1^3 is comparable to it -- below epsilon1 ~ 0.1 -- and there
// the metal travels 0.003 mm either way. The difference is immaterial; the
// default stays at the solver's value so nothing moves, and the deviation is
// recorded rather than silently harmonised.
const dimensionedScalar DarcyConstantlarge
(
    "DarcyConstantlarge", dimensionSet(1,-3,-1,0,0),
    transportProperties.lookupOrDefault<scalar>("darcy_A", 1.0e6)
);

const dimensionedScalar DarcyConstantsmall
(
    "DarcyConstantsmall", dimless,
    transportProperties.lookupOrDefault<scalar>("darcy_q", 1.0e-12)
);'''


def revert():
    bak = CF.with_suffix(CF.suffix + SUFFIX)
    if bak.exists():
        shutil.copy2(bak, CF); bak.unlink()
        print(f"  reverted {CF.name}")
    else:
        print("  nothing to revert")


def apply():
    if not CF.exists():
        sys.exit(f"ERROR: {CF} not found")
    text = CF.read_text()
    if "darcy_A" in text:
        sys.exit("ERROR: createFields.H looks already patched. --revert first.")
    if text.count(OLD) != 1:
        sys.exit(f"ERROR: expected one match for the two Darcy declarations, "
                 f"found {text.count(OLD)}. Nothing modified.")

    shutil.copy2(CF, CF.with_suffix(CF.suffix + SUFFIX))
    CF.write_text(text.replace(OLD, NEW, 1))
    print(f"  patched {CF.name}   (backup: {CF.name}{SUFFIX})")

    print("""
  BUILD
    of2506 && cd ~/thesis/src/lbf-com && ./Allwmake > /tmp/b.log 2>&1
    grep -c 'error:' /tmp/b.log

  REGRESSION -- defaults are the old literals, so this MUST be bit-identical.

  SET, in constant/transportProperties (top level):
    darcy_A   1.0e6;      // Zhao Eq. 7 range is 1e4 to 1e7; he does not say
    darcy_q   1.0e-12;    // solver value. Zhao Eq. 7 says 1e-3; immaterial

  Write BOTH explicitly even though they are the defaults. A value that is
  swept has to be visible in the case that produced each run, or the sweep
  cannot be reconstructed afterwards.

  THE SWEEP THIS ENABLES
    darcy_A is now the one parameter that directly controls how far mushy
    metal travels before it freezes, and Zhao leaves it free. Run 1e5 / 1e6 /
    1e7 and read POOL WIDTH, not bead width:

      python3 ~/thesis/tools/pool_width.py | tee poolwidth.txt

    Do this only AFTER patch_marangoni_mushy.py is enabled -- with the
    liquidus gate still in place the Marangoni force is identically zero
    everywhere below 906 K, so darcy_A has nothing to act on and the sweep
    would return three identical answers and look like a null result.""")


if __name__ == "__main__":
    revert() if "--revert" in sys.argv else apply()
