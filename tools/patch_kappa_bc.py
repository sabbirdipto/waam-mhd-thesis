#!/usr/bin/env python3
"""Give kappa and cp real boundary conditions. Edit 8.

    NOT the planned Edit 3 (arc pressure). This is a NEW eighth change,
    found by testing rather than planned: an upstream defect that makes
    every thermal boundary condition inert. The guide's ledger of seven
    changes becomes eight.

THE DEFECT
    createFields.H declares both fields with the constructor

        volScalarField kappa(IOobject(...), mesh, dimensionedScalar(...,0.0));

    which gives every patch a 'calculated' boundary condition initialised to
    ZERO. updateKappaCp.H fills the interior and then calls
    kappa.correctBoundaryConditions() -- but 'calculated' has no evaluation
    rule, so that call does nothing and the boundary values stay at zero.

    Confirmed empirically in waam-arc-01: 0.04/kappa reports
        inlet / outlet / walls :  type calculated;  value uniform 0

WHY IT MATTERS
    TEqn.H transports heat by two routes, and zero coefficients kill both AT
    THE BOUNDARY ONLY:

        - fvm::laplacian(kappa, T)      conduction:  kappa_f = 0  -> no flux
        + fvm::div(rhophicp, T)         advection:   rhophicp uses
                                        fvc::interpolate(cp), also 0 -> no flux

    The consequence is that EVERY thermal boundary condition in EVERY case is
    inert. The domain is a perfectly insulated box whatever the dictionaries
    say. Measured in waam-arc-01: walls at a fixed 300 K, walls adiabatic, and
    walls held at 2000 K all produced identical output to nine significant
    figures, and the energy balance closed at 100% because nothing could ever
    leave.

    This also means the planned wall-thermal-BC sensitivity study would have
    measured nothing and reported "no sensitivity" -- the same species of
    silent false negative as the sigma_e_gas_ratio bug found in Phase 2.

THE FIX
    Pass zeroGradientFvPatchScalarField::typeName as the patch-field type, so
    correctBoundaryConditions() extrapolates the adjacent cell value to the
    face. That is the physically correct choice: the conductivity and specific
    heat AT a wall are those of the material touching it.

    zeroGradientFvPatchField is already used four times in createFields.H, so
    no new include is needed.

THIS PATCH CHANGES RESULTS. IT IS NOT BIT-IDENTICAL.
    Unlike Edits 1 and 2, there is no legacy path to preserve -- the whole
    point is that heat can now cross a boundary where before it could not.
    Any case in which heat reaches a boundary WILL give different numbers, and
    those numbers are the correct ones. Re-run both gates and record the
    change rather than expecting a match.

USAGE
    python3 patch_kappa_bc.py            # apply
    python3 patch_kappa_bc.py --revert
"""
import shutil, sys
from pathlib import Path

SOLVER = Path.home() / "thesis/src/lbf-com/applications/solvers/laserbeamFoam"
CF = SOLVER / "createFields.H"
BAK = CF.with_suffix(".H.kbak")

ZG = "    zeroGradientFvPatchScalarField::typeName"

EDITS = [
    # (anchor, replacement) -- anchors taken verbatim from the upstream file
    ('    mesh,\n'
     '    dimensionedScalar("cp", dimensionSet(0, 2, -2, -1, 0), 0.0)\n'
     ');',
     '    mesh,\n'
     '    dimensionedScalar("cp", dimensionSet(0, 2, -2, -1, 0), 0.0),\n'
     + ZG + '\n'
     ');'),

    ('    mesh,\n'
     '    dimensionedScalar("kappa",dimensionSet(1, 1, -3, -1, 0),0.0)\n'
     ');',
     '    mesh,\n'
     '    dimensionedScalar("kappa",dimensionSet(1, 1, -3, -1, 0),0.0),\n'
     + ZG + '\n'
     ');'),
]


def revert():
    if not BAK.exists():
        print("  nothing to revert"); return
    shutil.copy2(BAK, CF); BAK.unlink()
    print(f"  reverted {CF.name}")


def apply():
    if not CF.exists():
        sys.exit(f"ERROR: {CF} not found")
    text = CF.read_text()

    # Guard must be SPECIFIC to these two declarations. zeroGradientFvPatch-
    # ScalarField::typeName already appears three times in this file for other
    # fields, so testing for it alone reports every clean file as patched.
    # The trailing comma after the dimensionedScalar only exists once patched.
    if ('dimensionedScalar("kappa",dimensionSet(1, 1, -3, -1, 0),0.0),' in text
            or 'dimensionedScalar("cp", dimensionSet(0, 2, -2, -1, 0), 0.0),'
            in text):
        sys.exit("ERROR: looks already patched. Run --revert first to reapply.")

    for old, _ in EDITS:
        n = text.count(old)
        if n != 1:
            sys.exit(f"ERROR: expected exactly one match for:\n-----\n{old}\n"
                     f"-----\nfound {n}. The upstream file has changed. "
                     f"Nothing was modified.")

    # zeroGradientFvPatchField must already be available; it is used elsewhere
    # in this file, but check rather than assume.
    if "zeroGradientFvPatch" not in text:
        sys.exit("ERROR: createFields.H does not reference zeroGradientFvPatch "
                 "anywhere, so the type may not be declared. Nothing modified.")

    shutil.copy2(CF, BAK)
    for old, new in EDITS:
        text = text.replace(old, new, 1)
    CF.write_text(text)
    print(f"  patched {CF.name}   (backup: {BAK.name})")
    print("    cp    -> zeroGradient boundaries")
    print("    kappa -> zeroGradient boundaries")
    print("""
  Next:
    of2506
    cd ~/thesis/src/lbf-com && ./Allwmake > /tmp/b.log 2>&1
    grep -c 'error:' /tmp/b.log

  VERIFY the defect is gone -- the 2000 K wall test must now respond:
    cd ~/thesis/cases/00-tutorials/waam-arc-01
    sed -n '/boundaryField/,$p' 0.04/kappa | head -20
        walls should read 'zeroGradient', NOT 'calculated ... uniform 0'

    Then re-run the A/B pair at endTime 0.05. The 'bottom' column for a
    2000 K wall must now be far above the adiabatic one. If the two are
    still identical, the patch did not take and something else is wrong.

  EXPECT BOTH GATES TO MOVE.
    Gate A (hartmann1d) is isothermal, so it may still be bit-identical.
    Gate B (Plate2D) has real temperature gradients reaching its boundaries
    and WILL change. That is the patch working, not the patch failing.
    Record the new numbers as the post-Edit-8 baseline.""")


if __name__ == "__main__":
    revert() if "--revert" in sys.argv else apply()
