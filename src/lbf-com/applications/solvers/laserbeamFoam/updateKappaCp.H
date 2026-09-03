{
    const scalarField& alpha1I = alpha1;
    const scalarField& TI = T;
    scalarField& kappaI = kappa;
    scalarField& cpI = cp;
    forAll(kappaI, celli)
    {
        // alpha1 is the METAL fraction (VOF). The LIQUID fraction is taken
        // from temperature via Zhao Eq. 8 -- the same linear relation the
        // solver itself uses at createFields.H:349 and inverts in TEqn.H:89.
        // The epsilon1 FIELD is not used here: this file is first included at
        // createFields.H:282, before epsilon1 exists. Deriving it from T also
        // keeps kappa a pure function of temperature, which is what the
        // polynomial fits assume.
        const scalar aM = max(min(alpha1I[celli], 1.0), 0.0);

        scalar kappaMetal;
        if (twoBranchKappa)
        {
            const scalar eL =
                max(min((TI[celli] - TsolidusMetal)
                       /(TliquidusMetal - TsolidusMetal), 1.0), 0.0);

            kappaMetal =
                eL*polykappa_mL.value(TI[celli])
              + (1.0 - eL)*polykappa_mS.value(TI[celli]);
        }
        else
        {
            // EXACTLY the original expression. Writing eL*k + (1-eL)*k with
            // identical branches is algebraically equal but differs by ~1 ULP
            // about 5% of the time, and a chaotic case amplifies that to tens
            // of percent. Branching keeps the legacy path BIT-identical, so a
            // regression test can distinguish "unchanged" from "barely changed".
            kappaMetal = polykappa_m.value(TI[celli]);
        }

        kappaI[celli] =
            (aM*kappaMetal)
          + ((1.0 - aM)*polykappa_g.value(TI[celli]));

        cpI[celli] =
            (
                max(min(alpha1I[celli], 1.0), 0.0)
               *polycp_m.value(TI[celli])
            )
          + (
                (1.0 - max(min(alpha1I[celli], 1.0), 0.0))
               *polycp_g.value(TI[celli])
            );
    }

    kappa.correctBoundaryConditions();
    cp.correctBoundaryConditions();
}
