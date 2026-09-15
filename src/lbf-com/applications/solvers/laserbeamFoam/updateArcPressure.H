// Arc pressure -- Zhao et al. (2021) Eqs. 25 and 26.
//
//     P_arc(r) = mu0*I^2/(8*pi*sigma_p^2) * exp(-r^2/(2*sigma_p^2))
//
// The arc column's azimuthal self-field pinches it non-uniformly (strong at the
// constricted cathode, weak at the broad anode), driving plasma toward the
// workpiece. The momentum flux arriving at the anode is this pressure.
// Peak 227.8 Pa at I = 135 A, sigma_p = 2 mm; surface integral mu0*I^2/(4*pi)
// = 1.82 mN. For scale, 2 mm of liquid aluminium weighs 52 Pa.
//
// Recomputed each outer iteration from the torch position x0 + V*t, exactly as
// updateArcSource.H does, so the two cannot drift apart.
//
// z GATE: pArc is masked by snGrad(alpha1) and so acts at ANY interface,
// including a free-flying droplet. At 228 Pa a droplet feels ~6x its own
// weight, which would double-count the plasma action on droplets this model
// already injects at a prescribed velocity. Applied only below
// arc_z0 + arc_p_zmax_m (default 3 mm: covers a ~2.4 mm bead, excludes
// droplets spawning near 3.5-4 mm).
{
    if (arcOn && arcI_ > SMALL)
    {
        const scalar tNow = runTime.value();
        const scalar xc   = arc_x0_ + arcV_*tNow;
        const scalar mu0  = 4.0*M_PI*1.0e-7;
        const scalar kp   = mu0*arcI_*arcI_/(8.0*M_PI*arcSigP_*arcSigP_);
        const scalar den  = 2.0*arcSigP_*arcSigP_;

        const vectorField& Cc = mesh.C().primitiveField();
        scalarField& pa = pArc.primitiveFieldRef();

        forAll(pa, celli)
        {
            const scalar dz = Cc[celli].z() - arc_z0_;
            if (dz > arcPzMax_) { pa[celli] = 0.0; continue; }

            const scalar dx = Cc[celli].x() - xc;
            const scalar dy = Cc[celli].y() - arc_y0_;
            const scalar e  = (dx*dx + dy*dy)/den;

            pa[celli] = (e > 30.0) ? 0.0 : kp*Foam::exp(-e);
        }

        pArc.correctBoundaryConditions();
    }
}
