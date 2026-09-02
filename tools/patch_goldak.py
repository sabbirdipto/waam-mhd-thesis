#!/usr/bin/env python3
"""Goldak double-ellipsoid arc heat source for laserbeamFoam.

WHY
    laserbeamFoam heats with a ray-traced laser. GMAW-WAAM heats with a moving
    arc, which Zhao et al. (2021) Eq. 27 models as Goldak's double ellipsoid --
    narrow in front, wide behind, because the torch is travelling.

    TEqn.H has no fvOptions, so a codedFvOption would be SILENTLY IGNORED: the
    case would run, the arc would do nothing, and nothing would say why. The
    source therefore has to enter the equation directly.

HOW
    Follows the pattern already in TEqn.H two lines away:
        Joulesource = mthd.valid() ? mthd->Joulesource() : zeroJoulesource
    arcSource is a peer of laser.deposition() -- another optional volumetric
    source on the RHS, identically zero when no arc is configured.

THREE THINGS THE NAIVE IMPLEMENTATION GETS WRONG
    1. Goldak is a HALF ellipsoid in depth. The 6*sqrt(3) normalisation with
       ff+fr=2 is correct only below the surface; applied symmetrically it
       delivers 2x eta*U*I. Enforced here by skipping dz > 0.
    2. In a VOF domain the ellipsoid overlaps argon (rho 1.6, cp 520). Gated
       on alpha.metal.
    3. A HALF domain contains half the arc. The target is scaled by sym.

    The source then renormalises so delivered power equals the target exactly,
    and PRINTS what fraction it had to move. A large correction means the
    ellipsoid is mostly outside the metal -- worth knowing, not hiding.

SAFETY
    Reads constant/caseParameters with READ_IF_PRESENT. No file, or no
    goldak_af_m key, means arcSource stays zero and every existing case is
    bit-identical.

USAGE
    python3 patch_goldak.py            # apply
    python3 patch_goldak.py --revert
"""
import sys, shutil
from pathlib import Path

SOLVER = Path.home() / "thesis/src/lbf-com/applications/solvers/laserbeamFoam"
CF, TE, MAIN = (SOLVER / "createFields.H", SOLVER / "TEqn.H",
                SOLVER / "laserbeamFoam.C")
NEW = SOLVER / "updateArcSource.H"

CF_ANCHOR = 'const Polynomial<8> polykappa_mL(readMetalKappa("poly_kappa_liquid"));'
CF_ADD = r'''

// --- Goldak arc heat source (WAAM-MHD thesis addition) ---------------------
// Read from constant/caseParameters, which generate_cases.py writes from
// caseMatrix.csv. Absent, or lacking goldak_af_m, the arc is simply off.
IOdictionary arcDict
(
    IOobject
    (
        "caseParameters",
        runTime.constant(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    )
);

const bool arcOn = arcDict.found("goldak_af_m");

scalar arcPower_ = 0.0, arc_af_ = 1.0, arc_ar_ = 1.0, arc_b_ = 1.0,
       arc_c_ = 1.0, arc_ff_ = 1.0, arc_fr_ = 1.0, arcV_ = 0.0,
       arc_x0_ = 0.0, arc_y0_ = 0.0, arc_z0_ = 0.0, arcSym_ = 1.0;

if (arcOn)
{
    arcPower_ = readScalar(arcDict.lookup("arc_power_W"));
    arc_af_   = readScalar(arcDict.lookup("goldak_af_m"));
    arc_ar_   = readScalar(arcDict.lookup("goldak_ar_m"));
    arc_b_    = readScalar(arcDict.lookup("goldak_b_m"));
    arc_c_    = readScalar(arcDict.lookup("goldak_c_m"));
    arcV_     = readScalar(arcDict.lookup("travel_speed_m_s"));
    arc_x0_   = readScalar(arcDict.lookup("arc_x0_m"));
    arc_y0_   = readScalar(arcDict.lookup("arc_y0_m"));
    arc_z0_   = readScalar(arcDict.lookup("arc_z0_m"));

    // Goldak's front/rear split. ff + fr = 2 by construction.
    arc_ff_ = 2.0*arc_af_/(arc_af_ + arc_ar_);
    arc_fr_ = 2.0*arc_ar_/(arc_af_ + arc_ar_);

    // A half domain holds half the arc, so it must receive half the power.
    word arcSymWord("full");
    if (arcDict.found("sym")) { arcDict.lookup("sym") >> arcSymWord; }
    arcSym_ = (arcSymWord == "half") ? 0.5 : 1.0;

    Info<< "Arc: Goldak double ellipsoid" << nl
        << "     af=" << arc_af_*1e3 << "mm  ar=" << arc_ar_*1e3
        << "mm  b=" << arc_b_*1e3 << "mm  c=" << arc_c_*1e3 << "mm" << nl
        << "     ff=" << arc_ff_ << "  fr=" << arc_fr_
        << "  travel=" << arcV_ << " m/s" << nl
        << "     power " << arcPower_ << " W, domain " << arcSymWord
        << " -> target " << arcPower_*arcSym_ << " W" << endl;
}
else
{
    Info<< "Arc: none (no goldak_af_m in caseParameters)" << endl;
}

volScalarField arcSource
(
    IOobject
    (
        "arcSource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimPower/dimVolume, 0.0)
);
'''

ARC_FILE = r'''// Goldak double-ellipsoid arc heat source -- Zhao et al. (2021) Eq. 27.
// Recomputed once per timestep from the torch position x0 + V*t.
//
// Three corrections the bare formula needs in a VOF free-surface domain:
//   * dz > 0 skipped -- Goldak is a HALF ellipsoid in depth. The 6*sqrt(3)
//     normalisation with ff+fr=2 integrates to eta*U*I only below the
//     surface; applied symmetrically it delivers exactly TWICE that.
//   * gated on alpha.metal -- otherwise power lands in argon (rho 1.6,
//     cp 520), which would reach absurd temperatures.
//   * renormalised to the target, with the raw fraction reported. A small
//     correction is the ellipsoid tail leaving the metal; a large one means
//     the arc is mostly off the workpiece.
{
    if (arcOn)
    {
        const scalar tNow = runTime.value();
        const scalar xc = arc_x0_ + arcV_*tNow;

        const vectorField& Cc = mesh.C().primitiveField();
        const scalarField& a1 = alpha1.primitiveField();
        scalarField& q = arcSource.primitiveFieldRef();

        const scalar kf = 6.0*Foam::sqrt(3.0)*arc_ff_*arcPower_
                        /(arc_af_*arc_b_*arc_c_*M_PI*Foam::sqrt(M_PI));
        const scalar kr = 6.0*Foam::sqrt(3.0)*arc_fr_*arcPower_
                        /(arc_ar_*arc_b_*arc_c_*M_PI*Foam::sqrt(M_PI));

        forAll(q, celli)
        {
            const scalar dz = Cc[celli].z() - arc_z0_;
            if (dz > 0.0) { q[celli] = 0.0; continue; }   // above the surface

            const scalar dx = Cc[celli].x() - xc;
            const scalar dy = Cc[celli].y() - arc_y0_;
            const scalar a  = (dx > 0.0) ? arc_af_ : arc_ar_;
            const scalar k  = (dx > 0.0) ? kf : kr;

            const scalar e = 3.0*( dx*dx/(a*a)
                                 + dy*dy/(arc_b_*arc_b_)
                                 + dz*dz/(arc_c_*arc_c_) );
            if (e > 30.0) { q[celli] = 0.0; continue; }   // exp(-30) ~ 1e-13

            q[celli] = Foam::max(Foam::min(a1[celli], 1.0), 0.0)
                     * k * Foam::exp(-e);
        }

        // What the mesh actually holds, before any correction.
        scalar raw = 0.0;
        forAll(q, celli) { raw += q[celli]*mesh.V()[celli]; }
        reduce(raw, sumOp<scalar>());

        const scalar target = arcPower_*arcSym_;
        if (raw > SMALL)
        {
            const scalar scale = target/raw;
            forAll(q, celli) { q[celli] *= scale; }

            if (runTime.outputTime() || tNow < SMALL)
            {
                Info<< "Arc: " << raw << " W in metal below the surface, "
                    << "renormalised x" << scale << " to " << target << " W"
                    << " (raw " << 100.0*raw/target << "% of target)" << endl;
            }
        }
        else
        {
            Info<< "Arc: WARNING no metal under the arc at t=" << tNow
                << " -- source is zero" << endl;
        }

        arcSource.correctBoundaryConditions();
    }
}
'''

MAIN_OLD = '            #include "TEqn.H"'
MAIN_NEW = ('            #include "updateArcSource.H"\n'
            '            #include "TEqn.H"')

TE_OLD = """            laser.deposition()
          + Joulesource"""
TE_NEW = """            laser.deposition()
          + arcSource
          + Joulesource"""

EDITS = [(CF, [(CF_ANCHOR, CF_ANCHOR + CF_ADD)]),
         (TE, [(TE_OLD, TE_NEW)]),
         (MAIN, [(MAIN_OLD, MAIN_NEW)])]


def revert():
    n = 0
    for path, _ in EDITS:
        bak = path.with_suffix(path.suffix + ".gbak")
        if bak.exists():
            shutil.copy2(bak, path); bak.unlink(); n += 1
            print(f"  reverted {path.name}")
    if NEW.exists():
        NEW.unlink(); print(f"  removed {NEW.name}")
    print(f"  {n} file(s) restored" if n else "  nothing to revert")


def apply():
    for path, subs in EDITS:
        if not path.exists():
            sys.exit(f"ERROR: {path} not found")
        text = path.read_text()
        if "arcSource" in text:
            sys.exit(f"ERROR: {path.name} looks already patched. --revert first.")
        for old, _ in subs:
            if text.count(old) != 1:
                sys.exit(f"ERROR: {path.name}: expected one match for:\n"
                         f"-----\n{old}\n-----\nfound {text.count(old)}. "
                         f"Nothing modified.")

    for path, subs in EDITS:
        shutil.copy2(path, path.with_suffix(path.suffix + ".gbak"))
        text = path.read_text()
        for old, new in subs:
            text = text.replace(old, new, 1)
        path.write_text(text)
        print(f"  patched {path.name}   (backup: {path.name}.gbak)")
    NEW.write_text(ARC_FILE)
    print(f"  created {NEW.name}")
    print("\n  Next:")
    print("    of2506 && cd ~/thesis/src/lbf-com && ./Allwmake > /tmp/b.log 2>&1")
    print("    grep -c 'error:' /tmp/b.log")
    print("  Then re-run Gate A and Gate B -- neither has caseParameters, so")
    print("  both must be BIT-IDENTICAL. That is the regression test.")


if __name__ == "__main__":
    revert() if "--revert" in sys.argv else apply()
