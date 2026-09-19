#!/usr/bin/env python3
"""Let Marangoni act in the mushy zone. A switch, defaulting to upstream.

THE DEFECT

    UEqn.H builds the Marangoni force correctly --

        Marangoni = epsilon1*Marangoni_Constant
                  * (gradT - nHatM*(nHatM & gradT))*mag(gradAlpha)

    where epsilon1 is the liquid fraction, so the force already ramps from
    zero at the solidus to full at the liquidus. That is the physically
    correct treatment and it is all that is needed.

    It then throws that away:

        forAll(MarangoniI, celli)
            if (TI[celli] < Tliquidus1Value)
                MarangoniI[celli] = vector::zero;

    Every cell below 906 K gets exactly zero, regardless of what epsilon1
    says. The ramp is replaced by a cliff.

WHY IT MATTERS HERE AND NOT UPSTREAM

    laserbeamFoam is a laser-welding solver. For steel, or for any alloy with
    a narrow freezing range, the mushy zone is a thin skin and zeroing it
    costs nothing. Al-5%Mg freezes over 815-906 K -- a 91 K range -- and in
    this case that band is most of the pool:

        measured at t = 0.9161 s:   462 cells above solidus
                                     73 cells above liquidus
                                -> 84% of the pool is mushy

    So the gate switches off Marangoni over 84% of the pool, and it does so
    precisely at the OUTER EDGE, where the surface temperature gradient is
    steepest and where outward flow is what sets the bead width.

WHY THAT IS THE SUSPECT FOR THE WIDTH DEFICIT

    A Rosenthal thin-plate estimate at 1512 W into a 4 mm plate at 20 mm/s
    gives a melt 2.9 mm wide on a cold plate and 4.3 mm once the plate has
    warmed to 450 K. Zhao's bead is 8.2 mm. Conduction alone cannot make that
    pool; convection has to roughly double it, and Zhao says twice which
    convection does it:

        "The liquid metal at the top of molten pool is driven and DOMINATED
         by Marangoni force to flow from the center to periphery."

    and, on the mushy zone specifically:

        "the solid-liquid coexistence region of aluminum alloy is large,
         where the temperature ranges from 815 to 906K, AND THE VELOCITY
         FIELD STILL EXISTS."

    The second sentence is a direct contradiction of the gate.

WHAT THIS PATCH DOES

    Wraps the gate in a switch read from transportProperties. The default is
    FALSE, which keeps the gate and therefore reproduces upstream behaviour
    bit for bit -- so Gate A and Gate B must be unchanged after this patch,
    and that is the regression test. Setting

        marangoni_mushy   true;

    in the metal sub-dictionary removes the cliff and leaves epsilon1 to do
    the ramping, which it was already doing.

    Nothing else changes. Darcy damping still resists motion in the mush
    through fvm::Sp(DC, U), with DC proportional to (1-epsilon1)^2, so the
    mushy zone is still far more viscous than the free liquid. This patch
    removes a third, redundant suppression -- not the first one.

USAGE
    python3 patch_marangoni_mushy.py            # apply
    python3 patch_marangoni_mushy.py --revert
"""
import sys, shutil
from pathlib import Path

SOLVER = Path.home() / "thesis/src/lbf-com/applications/solvers/laserbeamFoam"
CF, UE = SOLVER / "createFields.H", SOLVER / "UEqn.H"
SUFFIX = ".mbak"

# ---------------------------------------------------------------- createFields
CF_ANCHOR = ('[[maybe_unused]]const scalar '
             'emS(readScalar(transportProperties.lookup("emS")));')

CF_ADD = r'''

// --- Marangoni in the mushy zone (WAAM-MHD thesis) -------------------------
// UEqn.H hard-zeroes the Marangoni force below the liquidus, which discards
// the epsilon1 ramp already built into the expression. For Al-5Mg, whose
// 815-906 K freezing range covers ~84% of the pool, that switches off the
// dominant width mechanism over most of the molten region and all of its
// outer edge.
//
// DEFAULT FALSE: the gate stays, upstream behaviour is reproduced exactly,
// and every existing case is bit-identical. That is what Gate A and Gate B
// check.
//
// TO ENABLE, put  marangoni_mushy true;  at the TOP LEVEL of
// constant/transportProperties, beside sigma -- NOT inside the metal { }
// sub-dictionary. The lookup below is on the top-level IOdictionary, so a key
// nested inside metal { } is not found, the default false applies, the run
// completes normally and the result looks like "removing the gate changed
// nothing". A null input that reads as a null result is the worst kind of
// failure, so it is spelled out here and again at the end of the patch.
const bool marangoniMushy =
    transportProperties.lookupOrDefault<bool>("marangoni_mushy", false);'''

# ---------------------------------------------------------------- UEqn
UE_OLD = '''const scalar Tliquidus1Value = Tliquidus1.value();
const scalarField& TI = T;
vectorField& MarangoniI = Marangoni;
forAll(MarangoniI, celli)
{
    if(TI[celli] < Tliquidus1Value)
    {
        MarangoniI[celli] = vector::zero;
    }
}'''

UE_NEW = '''// The liquidus cliff. epsilon1 already ramps this force from zero at the
// solidus to full at the liquidus, and Darcy damping already resists motion
// in the mush, so this is a third and redundant suppression -- and the only
// one that is a step rather than a ramp. Zhao (2021) Sec 4.2 states the
// opposite for this alloy: "the solid-liquid coexistence region of aluminum
// alloy is large, where the temperature ranges from 815 to 906K, and the
// velocity field still exists."
if (!marangoniMushy)
{
    const scalar Tliquidus1Value = Tliquidus1.value();
    const scalarField& TI = T;
    vectorField& MarangoniI = Marangoni;
    forAll(MarangoniI, celli)
    {
        if(TI[celli] < Tliquidus1Value)
        {
            MarangoniI[celli] = vector::zero;
        }
    }
}'''

EDITS = [(CF, CF_ANCHOR, CF_ANCHOR + CF_ADD),
         (UE, UE_OLD, UE_NEW)]


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
        if "marangoniMushy" in text:
            sys.exit(f"ERROR: {path.name} looks already patched. --revert first.")
        if text.count(old) != 1:
            sys.exit(f"ERROR: {path.name}: expected one match for:\n"
                     f"-----\n{old}\n-----\nfound {text.count(old)}. "
                     f"Nothing modified.")

    for path, old, new in EDITS:
        shutil.copy2(path, path.with_suffix(path.suffix + SUFFIX))
        path.write_text(path.read_text().replace(old, new, 1))
        print(f"  patched {path.name}   (backup: {path.name}{SUFFIX})")

    print("""
  BUILD
    of2506 && cd ~/thesis/src/lbf-com && ./Allwmake > /tmp/b.log 2>&1
    grep -c 'error:' /tmp/b.log

  REGRESSION -- the default is false, so this MUST be bit-identical.
    Run Gate A and Gate B before enabling anything. If either moves, the
    switch is not defaulting the way it should and nothing downstream is
    trustworthy.

  ENABLE -- at the TOP LEVEL of constant/transportProperties, beside sigma.
  NOT inside the metal { } sub-dictionary. createFields.H looks this up on the
  top-level IOdictionary; a key nested inside metal { } is simply not found,
  the default false applies, and the run completes looking exactly like "the
  gate made no difference". Check it landed before launching:

      grep -n "marangoni_mushy" constant/transportProperties
      awk '/^metal/,/^}/' constant/transportProperties | grep marangoni_mushy
        -- the first must print a line, the SECOND must print nothing

    marangoni_mushy   true;

  THE TEST
    Re-run to the same time as A2 and compare POOL WIDTH, not bead width --
    this changes how far the melt spreads, and the bead follows the melt:

      python3 ~/thesis/tools/pool_width.py | tee poolwidth.txt

    'w fuse' is the number. If it widens toward 8 mm, this was the width
    mechanism and the force ablation was chasing the wrong thing. If it does
    not move, Marangoni in the mush is not the lever and the gate should be
    left in place as an honest upstream limitation.""")


if __name__ == "__main__":
    revert() if "--revert" in sys.argv else apply()
