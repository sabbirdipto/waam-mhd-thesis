#!/usr/bin/env python3
"""Edit 4 -- droplet deposition as injected spheres. Task A4.

WHAT THIS MODELS, AND WHAT IT DOES NOT
    GMAW does not melt the workpiece with the arc. Measured on this model, the
    arc alone saturates at ~832 K: past the 815 K solidus by 17 K and 74 K
    short of the 906 K liquidus. No pool forms. The melting energy arrives
    inside superheated filler droplets, which Zhao's assumption 1 models as
    liquid entering at 1506 K.

    Droplet FORMATION is not modelled. Spheres of a fixed radius are injected
    at a fixed rate, both taken from Zhao Eq. 30 and validated there against a
    simulated free-flying droplet (1.105 mm theoretical, 1.114 mm simulated,
    0.88% error). This follows Zhou [17] and Bai [18], whom Zhao cites and
    criticises -- their objection is that a fixed sphere at a fixed frequency
    does not PREDICT detachment, which is true and is the trade being made.
    Detachment is imposed here from a validated result rather than predicted.

THE BUDGET, which the solver prints so it can be checked rather than trusted
    wire 1.2 mm at 0.15 m/s        -> 1.696e-7 m3/s  ->  0.450 g/s
    droplet r 1.105 mm             -> 5.652e-9 m3    ->  30.0 Hz
    enthalpy above 300 K           -> 1.492 MJ/kg
      cp(815-300) + Lf + cp(906-815) + cp(1506-906)
    droplet power                  ->  671 W   against the arc's 1512 W: +44%
    momentum, m_dot * v_impact     -> 4.05e-4 N at v = 0.9 m/s

    v_impact is 0.9 m/s, not the 0.15 m/s wire feed. The paper measures it:
    "the maximum flow velocity within droplet reaches 0.9 m/s after droplet
    separation, which will have a significant effect on the fluid flow in the
    molten pool". The droplet is accelerated between necking and detachment by
    an electromagnetic pinch of order 1e6 N/m3. Using the feed speed would
    understate the impingement momentum by a factor of six.

THREE THINGS READING THE SOURCE CHANGED, none of which were guessable
    1. alphaEqnSubCycle.H does  alpha1 = alpha1.prevIter()  on every PIMPLE
       iteration after the first. Anything written into alpha1 inside the
       loop is therefore DISCARDED. The injection must happen once per
       timestep, before  while (pimple.loop()).
    2. UEqn computes Darcy damping from epsilon1 and runs BEFORE TEqn
       recomputes epsilon1 from T. Injecting a 1506 K sphere without also
       setting epsilon1 = 1 gives liquid metal a solid's Darcy coefficient,
       and the droplet is pinned where it lands.
    3. UEqn SOLVES for U, so a patched velocity is only an initial guess and
       is overwritten. Momentum has to enter as a body force, exactly as
       Marangoni and JcrossB already do -- and picks up the same
       damper = 2*rho/(rho1+rho2) weighting, which is Zhao's
       2*rho/(rho_metal+rho_gas) from Eqs 14-15.

SO: mass and enthalpy are placed geometrically as a sphere; momentum is a
    matched body force on the same metal. Impulse over one impact integrates
    to exactly m_drop * v_impact.

MASS IS EXACT BY CONSTRUCTION
    Firing on a clock would quantise badly: a 2.21 mm sphere spans 4.4 cells
    at 0.5 mm, so how much volume a sphere actually covers jitters with where
    it lands on the grid. Instead the solver accumulates owed volume at
    Vdot*dt and fires when it exceeds one droplet, then subtracts the volume
    ACTUALLY placed. Errors carry forward instead of accumulating, so the
    time-averaged mass flow is 0.450 g/s however the cells fall.

USAGE
    python3 patch_droplet.py            # apply
    python3 patch_droplet.py --revert
"""
import shutil, sys
from pathlib import Path

SOLVER = Path.home() / "thesis/src/lbf-com/applications/solvers/laserbeamFoam"

NEWFILE = SOLVER / "updateDropletSource.H"

# ---------------------------------------------------------------- the new file
DROPLET_H = r'''// Droplet deposition -- Edit 4.  Zhao et al. (2021), assumption 1 and Eq. 30.
//
// CALLED ONCE PER TIMESTEP, BEFORE the PIMPLE loop. That placement is not
// stylistic: alphaEqnSubCycle.H resets alpha1 to alpha1.prevIter() on every
// outer iteration after the first, so anything written to alpha1 inside the
// loop is thrown away without a word.
//
// Mass and enthalpy are placed as a sphere. Momentum is a body force on the
// same metal, because UEqn solves for U and a patched velocity would simply
// be overwritten. The two are matched: the impulse delivered over one impact
// integrates to exactly m_drop * v_impact.
{
    if (dropOn)
    {
        dropletForce.primitiveFieldRef() = vector::zero;

        const scalar tNow = runTime.value();
        const scalar dt   = runTime.deltaTValue();

        // The torch position, from THE SAME expression the arc uses. If these
        // two ever disagree the metal lands away from the hot spot and every
        // result downstream is quietly wrong.
        const scalar xc = arc_x0_ + arcV_*tNow;
        const scalar yc = arc_y0_;

        const vectorField& Cc = mesh.C().primitiveField();
        const scalarField& Vc = mesh.V();
        scalarField& a1 = alpha1.primitiveFieldRef();
        scalarField& Tf = T.primitiveFieldRef();
        scalarField& e1 = epsilon1.primitiveFieldRef();

        // THE OLD-TIME FIELDS TOO, and this is not optional.
        //
        // primitiveFieldRef() calls storeOldTimes(), which snapshots the
        // CURRENT values as old-time BEFORE returning write access. So the
        // three lines above have already recorded "this cell is argon at
        // 300 K" as the state at the start of the timestep. Every equation
        // then restarts from that:
        //
        //     advector.advect()  ->  alpha1 = alpha1.oldTime() - dt*fluxes
        //     fvm::ddt(rhoCp,T)  ->  anchored on T.oldTime()
        //
        // Writing only the current field is therefore erased within the same
        // timestep, silently. Measured before this fix: three droplets
        // injected, 18000.0 mm^3 of metal before and after, to the last digit,
        // and a temperature field bit-identical to a run with no droplets.
        //
        // Euler ddt uses exactly one old-time level, so one is enough. With
        // CrankNicolson a second level would also need setting.
        //
        // Note this works only because TEqn is written non-conservatively:
        //     fvm::ddt(rhoCp,T) - fvm::Sp(fvc::ddt(rhoCp), T)
        //       ==  rhoCp_old * (T - T_old) / dt
        // so the OLD heat capacity is the coefficient and T_old is the
        // anchor. In a conservative form the rhoCp jump from argon (832) to
        // aluminium (2.49e6) on the injection step would have driven T to
        // about 0.5 K instead.
        scalarField& a1o = alpha1.oldTime().primitiveFieldRef();
        scalarField& Tfo = T.oldTime().primitiveFieldRef();
        scalarField& e1o = epsilon1.oldTime().primitiveFieldRef();

        // AND the two DERIVED old-time fields, or the run aborts on the
        // injection step with a floating-point exception.
        //
        // rho and rhoCp are recomputed from alpha and T inside the PIMPLE
        // loop, so their CURRENT values look after themselves. Their old-time
        // values do not, and two terms in TEqn read them:
        //
        //   TRHS = LatentHeat*(fvc::ddt(rho, epsilon1) + ...)
        //        = 358000 * (rho*eps - rho_old*eps_old)/dt
        //        = 358000 * (2650*1 - 1.6*1)/1e-5   =  9.5e13 W/m^3
        //   the arc, for scale,                     =  7.6e9  W/m^3
        //
        // a heat SINK 12,000 times the arc, acting against a coefficient of
        // rhoCp_old = 832 (argon) instead of 2.49e6 (aluminium):
        //
        //   dT = -9.5e13 * 1e-5 / 832  =  -1.1e6 K
        //
        // T goes hugely negative, pow(2*pi*Mm*R*T, 0.5) in the vaporisation
        // term takes the root of a negative number, and trapFpe aborts.
        // Measured: signal 8 at the exact timestep droplet 1 was injected.
        //
        // Setting both to their metal values makes ddt(rho, eps) vanish in
        // those cells -- which is correct, because nothing about that metal
        // is changing between the two time levels; it simply is there.
        scalarField& rhoo   = rho.oldTime().primitiveFieldRef();
        scalarField& rhoCpo = rhoCp.oldTime().primitiveFieldRef();
        const scalar rhoDrop   = rho1.value();
        const scalar rhoCpDrop = rhoDrop*polycp_m.value(dropT_);

        // ---- is a droplet due? ------------------------------------------
        // Accumulate owed volume rather than firing on a clock, so that
        // quantisation of the sphere onto the mesh cannot bias the mass flow.
        dropOwed_ += dropVdot_*dt;

        if (dropOwed_ >= dropVdrop_)
        {
            // Free surface under the torch: the highest metal within one
            // droplet radius of the axis. Reduced across processors, because
            // the column can straddle a decomposition boundary.
            scalar zTop = -GREAT;
            forAll(Cc, celli)
            {
                const scalar dx = Cc[celli].x() - xc;
                const scalar dy = Cc[celli].y() - yc;
                if (dx*dx + dy*dy <= dropR_*dropR_ && a1[celli] > 0.5)
                {
                    zTop = Foam::max(zTop, Cc[celli].z());
                }
            }
            reduce(zTop, maxOp<scalar>());
            if (zTop < -GREAT/2) { zTop = arc_z0_; }   // nothing there yet

            // Clear of existing metal by a gap, so the sphere never creates
            // mass inside material that is already there -- that would be a
            // straight volume error, not an approximation.
            dropZc_ = zTop + dropGap_ + dropR_;
            dropXc_ = xc;
            dropYc_ = yc;

            scalar placed = 0.0;
            label  nCell  = 0;
            forAll(Cc, celli)
            {
                const scalar dx = Cc[celli].x() - dropXc_;
                const scalar dy = Cc[celli].y() - dropYc_;
                const scalar dz = Cc[celli].z() - dropZc_;
                if (dx*dx + dy*dy + dz*dz <= dropR_*dropR_)
                {
                    // Only the argon fraction is new volume. Counting the
                    // whole cell would over-report whenever a sphere clips
                    // metal that was already present.
                    placed += (1.0 - a1[celli])*Vc[celli];
                    a1[celli]  = 1.0;   a1o[celli]  = 1.0;
                    Tf[celli]  = dropT_; Tfo[celli] = dropT_;
                    e1[celli]  = 1.0;   e1o[celli]  = 1.0;
                    rhoo[celli]   = rhoDrop;
                    rhoCpo[celli] = rhoCpDrop;
                    // Both levels, every time. Setting only the current field
                    // is the bug this patch was rewritten to fix: epsilon1
                    // also matters because UEqn reads it for Darcy damping
                    // BEFORE TEqn recomputes it from T.
                    nCell++;
                }
            }
            reduce(placed, sumOp<scalar>());
            reduce(nCell,  sumOp<label>());

            dropOwed_ -= placed;
            dropImpactEnd_ = tNow + dropTau_;
            dropNum_++;

            alpha1.correctBoundaryConditions();
            T.correctBoundaryConditions();
            epsilon1.correctBoundaryConditions();

            Info<< "Droplet " << dropNum_ << " at t=" << tNow
                << "  centre (" << dropXc_*1e3 << " " << dropYc_*1e3
                << " " << dropZc_*1e3 << ") mm"
                << "  " << nCell << " cells"
                << "  placed " << placed*1e9 << " mm^3"
                << " of " << dropVdrop_*1e9 << " nominal"
                << "  (carried " << dropOwed_*1e9 << " mm^3)" << endl;

            if (nCell < 4)
            {
                Info<< "Droplet: WARNING only " << nCell << " cells in the"
                    << " sphere -- the mesh cannot resolve a "
                    << 2000.0*dropR_ << " mm droplet. Refine, or switch to"
                    << " the smooth-source form." << endl;
            }
        }

        // ---- momentum, during the impact window -------------------------
        // Spread over dropTau_ = d/v, the time the droplet takes to pass its
        // own diameter into the pool. The centre tracks downward at v_impact
        // so the force follows the droplet instead of staying where it was
        // born -- over 2.5 ms at 0.9 m/s it travels a full diameter.
        if (tNow < dropImpactEnd_ && dropTau_ > SMALL)
        {
            const scalar tSince = tNow - (dropImpactEnd_ - dropTau_);
            const scalar zc = dropZc_ - dropV_*tSince;

            // Weight by the metal actually present, so the force acts on
            // liquid and not on argon.
            scalar wsum = 0.0;
            scalarField w(Cc.size(), 0.0);
            forAll(Cc, celli)
            {
                const scalar dx = Cc[celli].x() - dropXc_;
                const scalar dy = Cc[celli].y() - dropYc_;
                const scalar dz = Cc[celli].z() - zc;
                if (dx*dx + dy*dy + dz*dz <= dropR_*dropR_)
                {
                    w[celli] = Foam::max(Foam::min(a1[celli], 1.0), 0.0)
                             * Vc[celli];
                    wsum += w[celli];
                }
            }
            reduce(wsum, sumOp<scalar>());

            if (wsum > SMALL)
            {
                // Total force over the window integrates to m_drop*v_impact.
                const scalar Ftot = dropMass_*dropV_/dropTau_;
                vectorField& fI = dropletForce.primitiveFieldRef();
                forAll(fI, celli)
                {
                    if (w[celli] > 0.0)
                    {
                        fI[celli] =
                            vector(0, 0, -1)*Ftot*w[celli]/(wsum*Vc[celli]);
                    }
                }
            }
        }

        dropletForce.correctBoundaryConditions();
    }
}
'''

# ---------------------------------------------------------------- createFields
CF_ANCHOR = '''volScalarField arcSource
(
    IOobject
    (
        "arcSource",'''

CF_NEW = '''// ---------------------------------------------------------------- Edit 4
// Droplet deposition. Enabled by the PRESENCE of droplet_radius_m, the same
// convention the arc uses -- a case without these keys behaves exactly as it
// did before this patch existed.
const bool dropOn = arcDict.found("droplet_radius_m");

scalar dropR_ = 0.0, dropV_ = 0.0, dropT_ = 0.0, dropGap_ = 0.0,
       dropVdot_ = 0.0, dropVdrop_ = 0.0, dropMass_ = 0.0, dropTau_ = 0.0,
       dropOwed_ = 0.0, dropImpactEnd_ = -1.0,
       dropXc_ = 0.0, dropYc_ = 0.0, dropZc_ = 0.0;
label dropNum_ = 0;

if (dropOn)
{
    dropR_    = readScalar(arcDict.lookup("droplet_radius_m"));
    dropV_    = readScalar(arcDict.lookup("droplet_v_impact_m_s"));
    dropT_    = readScalar(arcDict.lookup("droplet_T_K"));
    dropVdot_ = readScalar(arcDict.lookup("wire_volume_rate_m3_s"));
    dropGap_  = arcDict.lookupOrDefault<scalar>("droplet_standoff_m", 5e-4);

    // A half domain receives half the wire, exactly as it receives half the
    // arc. arcSym_ is already 0.5 or 1.0 from the block above.
    dropVdot_ *= arcSym_;

    dropVdrop_ = (4.0/3.0)*M_PI*dropR_*dropR_*dropR_*arcSym_;
    dropMass_  = dropVdrop_*rho1.value();

    // Impact duration: the time to pass one diameter at the impact speed.
    dropTau_ = (dropV_ > SMALL) ? 2.0*dropR_/dropV_ : 0.0;

    Info<< "Droplet: spheres, Zhao Eq. 30" << nl
        << "     r=" << dropR_*1e3 << "mm  v_impact=" << dropV_ << " m/s"
        << "  T=" << dropT_ << " K" << nl
        << "     volume rate " << dropVdot_ << " m3/s"
        << "  -> " << dropVdot_*rho1.value()*1e3 << " g/s"
        << "  -> " << dropVdot_/dropVdrop_ << " Hz" << nl
        << "     droplet " << dropVdrop_*1e9 << " mm^3, "
        << dropMass_*1e6 << " mg"
        << "  impact " << dropTau_*1e3 << " ms"
        << "  duty " << 100.0*dropTau_*dropVdot_/dropVdrop_ << "%" << nl
        // TWO forces, because quoting one invites the wrong comparison. The
        // time-averaged m_dot*v is what balances against gravity and the
        // Lorentz force; the peak is what acts during the 7.4% of the time a
        // droplet is actually landing, and it is 13.6x larger.
        << "     force  peak " << dropMass_*dropV_/dropTau_ << " N"
        << "   time-averaged " << dropVdot_*rho1.value()*dropV_ << " N" << nl
        << "     standoff " << dropGap_*1e3 << " mm above the free surface"
        << endl;
}
else
{
    Info<< "Droplet: none (no droplet_radius_m in caseParameters)" << endl;
}

// Momentum source, added to UEqn beside Marangoni and JcrossB. zeroGradient
// boundaries and not the default 'calculated', which Edit 8 established are
// initialised to zero and never evaluated.
volVectorField dropletForce
(
    IOobject
    (
        "dropletForce",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero", dimForce/dimVolume, Zero),
    zeroGradientFvPatchVectorField::typeName
);


''' + CF_ANCHOR

# ---------------------------------------------------------------- UEqn
UEQN_OLD = '''  - Marangoni*damper
  - JcrossB*damper
 ==
    fvOptions(rho, U)'''

UEQN_NEW = '''  - Marangoni*damper
  - JcrossB*damper
  - dropletForce*damper
 ==
    fvOptions(rho, U)'''

# ---------------------------------------------------------------- solver
MAIN_OLD = '''        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())'''

MAIN_NEW = '''        Info<< "Time = " << runTime.timeName() << nl << endl;

        // Droplet injection -- Edit 4. MUST be here, outside the PIMPLE loop:
        // alphaEqnSubCycle.H resets alpha1 to alpha1.prevIter() on every outer
        // iteration after the first, so an injection inside the loop is
        // silently discarded on iterations 2 and beyond.
        #include "updateDropletSource.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())'''

EDITS = [("createFields.H", CF_ANCHOR, CF_NEW),
         ("UEqn.H",         UEQN_OLD,  UEQN_NEW),
         ("laserbeamFoam.C", MAIN_OLD, MAIN_NEW)]


def revert():
    done = False
    for name, _, _ in EDITS:
        bak = SOLVER / (name + ".dbak")
        if bak.exists():
            shutil.copy2(bak, SOLVER / name); bak.unlink()
            print(f"  reverted {name}"); done = True
    if NEWFILE.exists():
        NEWFILE.unlink(); print(f"  removed {NEWFILE.name}"); done = True
    if not done:
        print("  nothing to revert")


def apply():
    if not SOLVER.exists():
        sys.exit(f"ERROR: {SOLVER} not found")
    if NEWFILE.exists():
        sys.exit("ERROR: updateDropletSource.H already exists. "
                 "Run --revert first to reapply.")

    # Check every anchor BEFORE writing anything. A half-applied patch to a
    # solver is worse than an unapplied one.
    texts = {}
    for name, old, _ in EDITS:
        p = SOLVER / name
        if not p.exists():
            sys.exit(f"ERROR: {p} not found")
        t = p.read_text()
        n = t.count(old)
        if n != 1:
            sys.exit(f"ERROR: expected exactly one anchor in {name}, found "
                     f"{n}. The upstream file has changed. Nothing modified.\n"
                     f"-----\n{old}\n-----")
        texts[name] = t

    for name, old, new in EDITS:
        p = SOLVER / name
        shutil.copy2(p, SOLVER / (name + ".dbak"))
        p.write_text(texts[name].replace(old, new, 1))
        print(f"  patched {name}")

    NEWFILE.write_text(DROPLET_H)
    print(f"  wrote   {NEWFILE.name}")
    print("""
  Build:
    of2506
    cd ~/thesis/src/lbf-com && ./Allwmake > /tmp/b.log 2>&1
    grep -c 'error:' /tmp/b.log        # must be 0

  A case with no droplet_radius_m in caseParameters is UNCHANGED by this
  patch -- the solver prints "Droplet: none" and behaves exactly as before.
  Re-run both gates to confirm that before trusting any droplet result.

  Then the startup banner should read, for the full domain:
    Droplet: spheres, Zhao Eq. 30
         r=1.105mm  v_impact=0.9 m/s  T=1506 K
         volume rate 1.696e-07 m3/s  -> 0.4495 g/s  -> 30 Hz
         droplet 5.652 mm^3, 14.98 mg  impact 2.456 ms  duty 7.37%
         force  peak 0.00549 N   time-averaged 0.000405 N

  Check those by hand before trusting anything downstream. The four the
  model rests on are 0.45 g/s, 30 Hz, 5.652 mm^3 and 4.05e-4 N, and each is
  derived from the wire feed and Zhao Eq. 30 alone -- no simulation needed.

  THE VALIDATION TARGET, also available before running: mass conservation
  fixes the bead cross-section at 84.8 mm^3 over 10 mm of travel = 8.5 mm^2.
  Zhao measures the SPSL deposit at ~8 mm wide and 1.80 mm high, a parabolic
  cap of ~9.6 mm^2. Ours coming out slightly smaller is correct, because some
  of that volume goes into remelted substrate rather than standing proud.""")


if __name__ == "__main__":
    revert() if "--revert" in sys.argv else apply()
