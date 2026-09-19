// Electromagnetic body force from the arc current -- Zhao Eqs. 22-24, in
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
