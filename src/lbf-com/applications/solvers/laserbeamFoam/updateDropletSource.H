// Droplet deposition -- Edit 4.  Zhao et al. (2021), assumption 1 and Eq. 30.
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
