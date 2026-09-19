#!/usr/bin/env python3
"""Electromagnetic body force in the molten pool -- Zhao Eqs. 22-24, DERIVED.

WHY THIS IS NOT A TRANSCRIPTION OF THE PAPER

    Eq. 21 lists six momentum sources for the molten pool. Five are implemented
    (gravity, buoyancy, surface tension, Marangoni, arc pressure). The sixth,
    the electromagnetic force, is absent.

    It cannot be typed in from the paper, because Eqs. 22-24 as printed have
    the radial and axial expressions TRANSPOSED, plus a stray factor 2.

THE DERIVATION

    Gaussian surface current in Zhao's own convention -- note his prefactors
    carry 2*pi^2 rather than 4*pi^2, which fixes sigma_j as sqrt(2) times the
    standard Gaussian sigma:

        j_z(r,z) = (I / (pi sj^2)) exp(-r^2/sj^2) * f(z)      [downward]
        f(z)     = 1 - d/L,   d = z0 - z,   0 <= d <= L

    Enclosed current and Ampere's law give the azimuthal self-field:

        I_enc(r) = I [1 - exp(-r^2/sj^2)]
        B_theta  = -(mu0 I / (2 pi r)) [1 - exp(-r^2/sj^2)] f(z)

    Continuity, div J = 0, fixes the radial spreading current:

        (1/r) d(r j_r)/dr = -dj_z/dz = (I/(pi sj^2 L)) exp(-r^2/sj^2)
        =>  j_r = (I / (2 pi L r)) [1 - exp(-r^2/sj^2)]

    and F = J x B, using r_hat x th_hat = z_hat and z_hat x th_hat = -r_hat:

        F_r = -j_z B_th = -(mu0 I^2)/(2 pi^2 sj^2 r) E (1-E) f^2      INWARD
        F_z =  j_r B_th = -(mu0 I^2)/(4 pi^2 L r^2)  (1-E)^2 f        DOWNWARD

    with E = exp(-r^2/sj^2). Both negative: the arc current pinches the pool
    inward and drives it down. That is the whole point of the term -- it is a
    PENETRATION mechanism, not a widening one.

WHAT ZHAO PRINTED, AND THE SIX SIGNATURES

    Eq. 22/23 (his x,y):  mu0 I^2 x / (4 pi^2 r^3 L)   (1-E)^2  (1 - (z-z0)/L)
    Eq. 24    (his z)  :  mu0 I^2  / (2 pi^2 r sj^2)  E (1-E)  (1 - (z-z0)/L)^2

    feature                belongs to   Zhao puts it on
    ---------------------  -----------  ----------------
    E(1-E) envelope        RADIAL       Eq. 24, axial
    (1-E)^2 envelope       AXIAL        Eq. 22/23, radial
    (1-z/L)^2 depth        RADIAL       Eq. 24, axial
    (1-z/L)^1 depth        AXIAL        Eq. 22/23, radial
    1/sj^2 prefactor       RADIAL       Eq. 24, axial
    1/(L r^2) prefactor    AXIAL        Eq. 22/23, radial

    Six independent signatures, one transposition. His Eq. 24 IS the correct
    radial expression; his Eqs. 22-23 ARE the correct axial expression carrying
    a spurious x/r.

    SEPARATELY, Eq. 23's denominator reads 4 pi^2 r^3 *2* L against Eq. 22's
    4 pi^2 r^3 L. That halves F_y and breaks axisymmetry -- a pool cannot know
    the difference between +x and +y under an axisymmetric arc.

    AND THE SIGN TRAP: Zhao writes (1 - (z-z0)/L), which is correct only if z
    is DEPTH, positive downward. In this solver z is height, positive up, and
    the plate lies below z0 = 0. Transcribed literally the factor reaches 2 at
    the plate bottom instead of 0 -- the force would grow with depth and point
    the wrong way. Below it is written as (1 - (z0-z)/L) and clamped to the
    slab, which is the same physics in this solver's coordinates.

VALIDATION -- two independent routes agree

    F_r peaks at r = 0.60 sj = 1.20 mm  ->  5.10e4 N/m^3
    F_z peaks at r = 1.12 sj = 2.24 mm  ->  1.48e4 N/m^3

    at I = 135 A, sj = 2 mm, L = 4 mm. These reproduce the force table computed
    independently earlier in this project to three significant figures. The
    solver prints both at startup; if they do not appear, something is wrong
    before a single timestep has run.

WHY J AND B ARE STORED AS FIELDS, NOT JUST F

    Phase 3 applies an EXTERNAL magnetic field, and that force is
    J_arc x B_external -- the same current, a different B. Storing J lets one
    field serve both phases with no new physics.

    The numerical alternative already in the solver (mthd->JcrossB(), a
    potential solve) is not needed here: Rm = mu0 sigma U L = 9e-3, and the
    motional current sigma*U*B at 10 mT is ~1.8e4 A/m^2 against the arc's peak
    j_z of 1.07e7 A/m^2. The arc current dominates by three orders of
    magnitude. This term is ADDED BESIDE JcrossB, not in place of it, so the
    numerical path stays available and the two never both carry the arc.

SAFETY
    Off unless  arc_lorentz  is true in constant/caseParameters. Absent, the
    field is identically zero and every existing case is bit-identical. That is
    what Gate A and Gate B check.

USAGE
    python3 patch_lorentz.py            # apply
    python3 patch_lorentz.py --revert
"""
import sys, shutil
from pathlib import Path

SOLVER = Path.home() / "thesis/src/lbf-com/applications/solvers/laserbeamFoam"
CF = SOLVER / "createFields.H"
UE = SOLVER / "UEqn.H"
MAIN = SOLVER / "laserbeamFoam.C"
NEW = SOLVER / "updateArcLorentz.H"
SUFFIX = ".lzbak"

# ------------------------------------------------------------------ new header
NEW_SRC = r'''// Electromagnetic body force from the arc current -- Zhao Eqs. 22-24, in
// their DERIVED form. See patch_lorentz.py for the derivation and for the six
// signatures showing the printed radial and axial expressions are transposed.
//
//   F_r = -(mu0 I^2)/(2 pi^2 sj^2 r) E (1-E) f^2        inward  (pinch)
//   F_z = -(mu0 I^2)/(4 pi^2 L r^2)    (1-E)^2 f        downward
//   E   = exp(-r^2/sj^2),   f = 1 - d/L,   d = z0 - z clamped to [0, L]
//
// J and B are stored separately so Phase 3 can cross the same current with an
// applied field. Gated on alpha.metal: current flows in metal, not argon.
{
    if (lorentzOn)
    {
        const scalar tNow = runTime.value();
        const scalar xc   = arc_x0_ + arcV_*tNow;

        const vectorField& Cc = mesh.C().primitiveField();
        const scalarField& a1 = alpha1.primitiveField();
        vectorField& Jf = Jarc.primitiveFieldRef();
        vectorField& Bf = Bself.primitiveFieldRef();
        vectorField& Ff = arcLorentz.primitiveFieldRef();

        const scalar mu0 = 4.0*M_PI*1.0e-7;
        const scalar I   = lorentzI_;
        const scalar sj2 = lorentzSj_*lorentzSj_;
        const scalar L   = lorentzL_;

        forAll(Ff, celli)
        {
            Jf[celli] = vector::zero;
            Bf[celli] = vector::zero;
            Ff[celli] = vector::zero;

            const scalar d = arc_z0_ - Cc[celli].z();
            if (d < 0.0 || d > L) continue;          // outside the slab

            const scalar aM = max(min(a1[celli], 1.0), 0.0);
            if (aM < SMALL) continue;                // no current in argon

            const scalar dx = Cc[celli].x() - xc;
            const scalar dy = Cc[celli].y() - arc_y0_;
            const scalar r2 = dx*dx + dy*dy;
            const scalar r  = Foam::sqrt(r2);
            const scalar fz = 1.0 - d/L;

            const scalar E  = Foam::exp(-r2/sj2);
            const scalar Em = 1.0 - E;

            // axial current, downward (-z)
            const scalar jz = -(I/(M_PI*sj2))*E*fz;
            // radial spreading current, outward, finite at r -> 0
            const scalar jr = (r > SMALL) ? (I/(2.0*M_PI*L*r))*Em : 0.0;
            // self-field, azimuthal, negative sense because J points in -z
            const scalar Bt = (r > SMALL) ? -(mu0*I/(2.0*M_PI*r))*Em*fz : 0.0;

            const scalar ex = (r > SMALL) ? dx/r : 0.0;   // r_hat
            const scalar ey = (r > SMALL) ? dy/r : 0.0;

            Jf[celli] = aM*vector(jr*ex, jr*ey, jz);
            Bf[celli] = aM*vector(-Bt*ey, Bt*ex, 0.0);    // th_hat = (-y,x,0)/r
        }

        // F = J x (B_self + B_applied). One expression serves both phases; with
        // lorentzBext_ zero this reduces exactly to the boxed F_r, F_z above.
        arcLorentz = Jarc ^ (Bself + dimensionedVector("Bext", Bself.dimensions(), lorentzBext_));

        if (!lorentzReported_)
        {
            lorentzReported_ = true;
            scalar fmax = 0.0;
            forAll(Ff, celli) fmax = max(fmax, mag(Ff[celli]));
            reduce(fmax, maxOp<scalar>());
            const scalar pr = mu0*I*I/(2.0*M_PI*M_PI*sj2*lorentzSj_);
            const scalar pz = mu0*I*I/(4.0*M_PI*M_PI*L*sj2);
            Info<< "Arc Lorentz: Zhao Eqs. 22-24, derived (radial/axial"
                << " transposition and Eq.23 factor 2 corrected)" << nl
                << "     I=" << I << " A  sigma_j=" << lorentzSj_*1e3
                << " mm  L=" << L*1e3 << " mm" << nl
                << "     analytic peaks  F_r " << 0.3515*pr
                << "  F_z " << 0.4070*pz << " N/m^3"
                << "   (expect 5.10e4 and 1.48e4 at 135 A, 2 mm, 4 mm)" << nl
                << "     mesh max |F| " << fmax << " N/m^3" << endl;
        }
    }
}
'''

# ------------------------------------------------------------- createFields.H
CF_ANCHOR = "const volVectorField zeroJcrossB"

CF_ADD = r'''// --- Arc electromagnetic body force (WAAM-MHD thesis) ----------------------
// Zhao Eqs. 22-24 in their DERIVED form -- the printed radial and axial
// expressions are transposed and Eq. 23 carries a spurious factor 2. See
// patch_lorentz.py. OFF unless arc_lorentz is true, so every existing case is
// bit-identical; that is what Gate A and Gate B check.
//
// arc_lorentz      true;        // enable
// sigma_j_m        0.002;       // Zhao Table 4, "distribution radius of current"
// substrate_L_m    0.004;       // plate thickness. Zhao S2 says 4 mm and his
//                               // figures show 4 mm; S3's "5 mm" is a typo.
// B_applied_T      (0 0 0);     // Phase 3: external field, crossed with the
//                               // SAME arc current
const bool lorentzOn =
    arcDict.lookupOrDefault<Switch>("arc_lorentz", false);

scalar lorentzI_ = 0.0, lorentzSj_ = 1.0, lorentzL_ = 1.0;
vector lorentzBext_ = vector::zero;
bool   lorentzReported_ = false;

if (lorentzOn)
{
    lorentzI_  = readScalar(arcDict.lookup("arc_current_A"));
    lorentzSj_ = readScalar(arcDict.lookup("sigma_j_m"));
    lorentzL_  = readScalar(arcDict.lookup("substrate_L_m"));
    lorentzBext_ = arcDict.lookupOrDefault<vector>("B_applied_T", vector::zero);
}

volVectorField Jarc
(
    IOobject("Jarc", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedVector("zero", dimensionSet(0, -2, 0, 0, 0, 1, 0), Zero),
    zeroGradientFvPatchVectorField::typeName
);

volVectorField Bself
(
    IOobject("Bself", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedVector("zero", dimensionSet(1, 0, -2, 0, 0, -1, 0), Zero),
    zeroGradientFvPatchVectorField::typeName
);

volVectorField arcLorentz
(
    IOobject("arcLorentz", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::AUTO_WRITE),
    mesh,
    dimensionedVector("zero", dimensionSet(1, -2, -2, 0, 0, 0, 0), Zero),
    zeroGradientFvPatchVectorField::typeName
);

'''

# --------------------------------------------------------------------- UEqn.H
UE_OLD = "  - JcrossB*damper"
UE_NEW = ("  - JcrossB*damper\n"
          "    // Arc self-field, Zhao Eqs. 22-24 derived. Zero unless\n"
          "    // arc_lorentz is set; never double-counts with JcrossB,\n"
          "    // which is the numerical potential-solve path.\n"
          "  - arcLorentz*damper")

# ------------------------------------------------------------ laserbeamFoam.C
MAIN_OLD = '#include "updateArcSource.H"'
MAIN_NEW = ('#include "updateArcSource.H"\n'
            '            #include "updateArcLorentz.H"')

EDITS = [(CF, CF_ANCHOR, CF_ADD + CF_ANCHOR),
         (UE, UE_OLD, UE_NEW),
         (MAIN, MAIN_OLD, MAIN_NEW)]


def revert():
    n = 0
    for path, _, _ in EDITS:
        bak = path.with_suffix(path.suffix + SUFFIX)
        if bak.exists():
            shutil.copy2(bak, path); bak.unlink(); n += 1
            print(f"  reverted {path.name}")
    if NEW.exists():
        NEW.unlink(); print(f"  removed {NEW.name}")
    print(f"  {n} file(s) restored" if n else "  nothing to revert")


def apply():
    for path, old, _ in EDITS:
        if not path.exists():
            sys.exit(f"ERROR: {path} not found")
        text = path.read_text()
        if "arcLorentz" in text:
            sys.exit(f"ERROR: {path.name} looks already patched. --revert first.")
        if text.count(old) != 1:
            hint = ""
            if path is MAIN:
                incs = [l for l in text.splitlines() if "#include \"update" in l]
                hint = ("\n  includes found in laserbeamFoam.C:\n    "
                        + "\n    ".join(incs or ["(none)"])
                        + "\n  Add  #include \"updateArcLorentz.H\"  beside the"
                          " other per-timestep updates by hand, then re-run"
                          " with the other two edits only.")
            sys.exit(f"ERROR: {path.name}: expected one match for:\n"
                     f"-----\n{old}\n-----\nfound {text.count(old)}. "
                     f"Nothing modified." + hint)

    for path, old, new in EDITS:
        shutil.copy2(path, path.with_suffix(path.suffix + SUFFIX))
        path.write_text(path.read_text().replace(old, new, 1))
        print(f"  patched {path.name}   (backup: {path.name}{SUFFIX})")
    NEW.write_text(NEW_SRC)
    print(f"  created {NEW.name}")

    print("""
  BUILD
    of2506 && cd ~/thesis/src/lbf-com && ./Allwmake > /tmp/b.log 2>&1
    grep -c 'error:' /tmp/b.log

  REGRESSION -- off by default, so this MUST be bit-identical.
    Run Gate A and Gate B before enabling anything.

  ENABLE, in constant/caseParameters (NOT transportProperties -- this reads
  arcDict, the same dictionary arc_current_A and sigma_p_m already live in):

    arc_lorentz       true;
    sigma_j_m         0.002;     // Zhao Table 4
    substrate_L_m     0.004;     // 4 mm; S3's "5 mm" is a typo, figures show 4

  arc_current_A is already present for arc pressure and is reused.

  CHECK THE STARTUP LINE. It prints the analytic peaks and the mesh maximum:

      analytic peaks  F_r 5.1e+04  F_z 1.48e+04 N/m^3

  Those two numbers are the validation. If they differ, sigma_j or L is wrong.
  If 'mesh max |F|' is far below F_r, the 0.5 mm mesh is not resolving a force
  that peaks 1.2 mm from the axis -- which is itself a result worth recording.

  WHAT TO EXPECT
    This is a PENETRATION term. It pinches inward and drives down, so expect
    d_sol to increase and w_fuse to stay put or narrow slightly. If the bead
    gets WIDER, check the sign of the depth factor first -- see the SIGN TRAP
    note in this file's header.""")


if __name__ == "__main__":
    revert() if "--revert" in sys.argv else apply()
