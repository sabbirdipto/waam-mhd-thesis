// Goldak double-ellipsoid arc heat source -- Zhao et al. (2021) Eq. 27.
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
