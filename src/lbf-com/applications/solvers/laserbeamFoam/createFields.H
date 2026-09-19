#include "createRDeltaT.H"

const dimensionedScalar deltaN = 1e-8/pow(average(mesh.V()), 1.0/3.0);

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Temperature\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();


// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);


IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

const word interfaceTrackingScheme
(
    transportProperties.lookup("interfaceTrackingScheme")
);

if
(
    interfaceTrackingScheme != "MULES"
 && interfaceTrackingScheme != "isoAdvector"
)
{
    FatalErrorIn("laserbeamFoam.C")
        << "Unknown interface tracking scheme "
        << interfaceTrackingScheme
        << ". Available options are 'MULES' and 'isoAdvector'"
        << exit(FatalError);
}

dictionary& transportPropertiesMetal = transportProperties.subDict("metal");
dictionary& transportPropertiesGas = transportProperties.subDict("gas");


const dimensionedScalar Tsolidus1
(
    "Tsolidus",
    dimensionSet(0, 0, 0, 1, 0),
    transportPropertiesMetal//.lookup("Tsolidus")
);
// Info<<Tsolidus1<<endl;
const dimensionedScalar Tsolidus2
(
    "Tsolidus",
    dimensionSet(0, 0, 0, 1, 0),
    transportPropertiesGas//.lookup("Tsolidus")
);
// Info<<Tsolidus2<<endl;

const dimensionedScalar Tliquidus1
(
    "Tliquidus",
    dimensionSet(0, 0, 0, 1, 0),
    transportPropertiesMetal//.lookup("Tliquidus")
);

const dimensionedScalar Tliquidus2
(
    "Tliquidus",
    dimensionSet(0, 0, 0, 1, 0),
    transportPropertiesGas//.lookup("Tliquidus")
);

const dimensionedScalar LatentHeat1
(
    "LatentHeat",
    dimensionSet(0, 2, -2, 0, 0),
    transportPropertiesMetal//.lookup("LatentHeat")
);

const dimensionedScalar LatentHeat2
(
    "LatentHeat",
    dimensionSet(0, 2, -2, 0, 0),
    transportPropertiesGas//.lookup("LatentHeat")
);

const dimensionedScalar beta1
(
    "beta",
    dimensionSet(0, 0, 0, -1, 0),
    transportPropertiesMetal//.lookup("beta")
);

const dimensionedScalar beta2
(
    "beta",
    dimensionSet(0, 0, 0, -1, 0),
    transportPropertiesGas//.lookup("beta")
);

const dimensionedScalar Marangoni_Constant
(
    "Marangoni_Constant",
    dimensionSet(1, 0, -2, -1, 0),
    transportProperties//.lookup("dsigmadT")
);

const dimensionedScalar p0
(
    "p0",
    dimensionSet(1, -1, -2, 0, 0, 0, 0),
    transportProperties//.lookup("p0")
); // atmospheric pressure

const dimensionedScalar Tvap
(
    "Tvap",
    dimensionSet(0, 0, 0, 1, 0),
    transportProperties//.lookup("Tvap")
);

const dimensionedScalar Mm
(
    "Mm",
    dimensionSet(1, 0, 0, 0, -1),
    transportProperties//.lookup("Mm")
);

const dimensionedScalar LatentHeatVap
(
    "LatentHeatVap",
    dimensionSet(0, 2, -2, 0, 0),
    transportProperties//.lookup("LatentHeatVap")
);

const dimensionedScalar R
(
    "R",
    dimensionSet(1, 2, -2, -1, -1, 0, 0),
        8.314
);

volScalarField em
(
    IOobject
    (
        "em",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

const dimensionedScalar SB
(
    "SB",
    dimensionSet(1, 0, -3, -4, 0, 0, 0),
    scalar(5.67e-8)
); // Stefan-Boltzmann Constant

[[maybe_unused]]const scalar emL(readScalar(transportProperties.lookup("emL")));
[[maybe_unused]]const scalar emS(readScalar(transportProperties.lookup("emS")));

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
    transportProperties.lookupOrDefault<bool>("marangoni_mushy", false);
const dimensionedScalar TRef ("TRef",transportProperties);

// --- Convective half of Zhao Eq. 29, free surface (WAAM-MHD thesis) --------
// The radiation half is already in TEqn.H via Qrad_expl / Qrad_inpl. This is
// the h(T - Ta) term, which was absent.
//
// DEFAULT ZERO: absent from transportProperties, the term vanishes and every
// existing case is bit-identical. That is what Gate A and Gate B check.
//
// Ta deliberately reuses TRef rather than introducing a second ambient
// temperature. Two names for one physical constant is how they drift apart.
const dimensionedScalar hConvSurf
(
    "hConvSurf",
    dimensionSet(1, 0, -3, -1, 0, 0, 0),
    transportProperties.lookupOrDefault<scalar>("h_conv_Wm2K", 0.0)
);

const Polynomial<8> polykappa_m(transportPropertiesMetal.lookup("poly_kappa"));
const Polynomial<8> polykappa_g(transportPropertiesGas.lookup("poly_kappa"));

// --- Two-branch metal conductivity (WAAM-MHD thesis addition) ---------------
// Aluminium's thermal conductivity halves on melting (162 -> 81 W/mK across a
// 91 K window). One Polynomial<8> cannot represent that step. These optional
// entries give the solid and liquid branches separately; updateKappaCp.H
// blends them by liquid fraction, derived from T via Zhao Eq. 8.
// Absent from the dictionary, both fall back to poly_kappa, so every existing
// case runs unchanged.
auto readMetalKappa = [&](const word& key) -> Polynomial<8>
{
    return transportPropertiesMetal.found(key)
         ? Polynomial<8>(transportPropertiesMetal.lookup(key))
         : polykappa_m;
};
const Polynomial<8> polykappa_mS(readMetalKappa("poly_kappa_solid"));
const Polynomial<8> polykappa_mL(readMetalKappa("poly_kappa_liquid"));

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
       arc_x0_ = 0.0, arc_y0_ = 0.0, arc_z0_ = 0.0, arcSym_ = 1.0,
       // Zhao Eq. 25 arc pressure. arcI_ = 0 switches it off entirely,
       // mirroring the way goldak_af_m switches the heat source on.
       arcI_ = 0.0, arcSigP_ = 0.002, arcPzMax_ = 0.003;

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

    // Arc pressure, Zhao Eq. 25. OPTIONAL: absent or zero current means the
    // term is never applied, so older cases stay bit-identical.
    arcI_     = arcDict.lookupOrDefault<scalar>("arc_current_A", 0.0);
    arcSigP_  = arcDict.lookupOrDefault<scalar>("sigma_p_m", 0.002);
    arcPzMax_ = arcDict.lookupOrDefault<scalar>("arc_p_zmax_m", 0.003);

    if (arcI_ > SMALL)
    {
        const scalar mu0_ = 4.0*M_PI*1.0e-7;
        Info<< "Arc pressure: Zhao Eq. 25, I=" << arcI_ << " A, sigma_p="
            << arcSigP_*1e3 << " mm" << nl
            << "     peak " << mu0_*arcI_*arcI_/(8.0*M_PI*arcSigP_*arcSigP_)
            << " Pa, total " << mu0_*arcI_*arcI_/(4.0*M_PI)*1e3 << " mN" << nl
            << "     applied for z < arc_z0 + " << arcPzMax_*1e3 << " mm"
            << endl;
    }
    else
    {
        Info<< "Arc pressure: OFF (arc_current_A absent or zero)" << endl;
    }

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

// ---------------------------------------------------------------- Edit 4
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

// --- Droplet state across restarts (WAAM-MHD thesis) -----------------------
// The scalars above are plain C++ and do not survive a restart. Resuming the
// M1-mushy run at t = 1.16 s discarded 4.63 mm^3 (0.82 of a droplet) from
// dropOwed_, starved the pool for 0.04 s, collapsed it by 37%, and ruined
// ~4 mm of bead. See patch_droplet_restart.py for the measurements.
//
// READ_IF_PRESENT: a case with no <time>/uniform/dropletState gets exactly the
// defaults above, so a run from t = 0 is bit-identical. Gate A and Gate B
// check that.
//
// Only the three ACCUMULATED values are stored. Everything else -- dropR_,
// dropV_, dropT_, dropVdot_, dropVdrop_, dropMass_, dropTau_, dropGap_ -- is
// re-read from caseParameters at startup and is identical on a restart.
// Persisting a value that is also read from a dictionary is how the two drift
// apart later.
IOdictionary dropletState
(
    IOobject
    (
        "dropletState",
        runTime.timeName(),
        "uniform",
        runTime,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    )
);

if (dropOn)
{
    const scalar owedPrev = dropOwed_;
    dropOwed_      = dropletState.lookupOrDefault<scalar>("dropOwed", dropOwed_);
    dropImpactEnd_ = dropletState.lookupOrDefault<scalar>
                     ("dropImpactEnd", dropImpactEnd_);
    dropNum_       = dropletState.lookupOrDefault<label>("dropNum", dropNum_);

    if (dropletState.found("dropOwed"))
    {
        Info<< "Droplet state: restored from "
            << runTime.timeName() << "/uniform/dropletState" << nl
            << "     dropOwed " << dropOwed_*1e9 << " mm^3 ("
            << dropOwed_/max(dropVdrop_, VSMALL) << " droplet)"
            << "   dropNum " << dropNum_ << endl;
    }
    else if (owedPrev == 0.0)
    {
        Info<< "Droplet state: fresh start (no uniform/dropletState)" << endl;
    }
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


// The METAL's own solidus/liquidus, as scalars. NOT the TSolidus/TLiquidus
// volScalarFields declared further down -- those are blended between metal and
// gas by alpha1, and they do not exist yet at the point updateKappaCp.H is
// first included (line 282, before their declarations at 285 and 299).
const scalar TsolidusMetal(readScalar(transportPropertiesMetal.lookup("Tsolidus")));
const scalar TliquidusMetal(readScalar(transportPropertiesMetal.lookup("Tliquidus")));
const bool twoBranchKappa = transportPropertiesMetal.found("poly_kappa_solid");
Info<< "Metal kappa: "
    << (transportPropertiesMetal.found("poly_kappa_solid")
        ? "two-branch (solid/liquid blended by liquid fraction)"
        : "single poly_kappa (legacy)")
    << endl;


const Polynomial<8> polycp_m(transportPropertiesMetal.lookup("poly_cp"));
const Polynomial<8> polycp_g(transportPropertiesGas.lookup("poly_cp"));

volScalarField cp
(
    IOobject
    (
        "cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//change to no write when working
    ),
    mesh,
    dimensionedScalar("cp", dimensionSet(0, 2, -2, -1, 0), 0.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE//change to no write when working
    ),
    mesh,
    dimensionedScalar("kappa",dimensionSet(1, 1, -3, -1, 0),0.0),
    zeroGradientFvPatchScalarField::typeName
);

#include "updateKappaCp.H"


volScalarField TSolidus
(
    IOobject
    (
        "TSolidus",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1*Tsolidus1 + (1.0 - alpha1)*Tsolidus2
);


volScalarField TLiquidus
(
    IOobject
    (
        "TLiquidus",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    alpha1*Tliquidus1 + (1.0 - alpha1)*Tliquidus2
);

volScalarField DC
(
    IOobject
    (
        "DC",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("DC", dimensionSet(1,-3,-1,0,0), 1.0e14)
);

volScalarField metalaverage
(
    IOobject
    (
        "metalaverage",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    fvc::average(alpha1)
);

volScalarField epsilon1
(
    IOobject
    (
        "epsilon1",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    max(min((T - TSolidus)/(TLiquidus - TSolidus), 1.0), 0.0)
);

volVectorField n_filtered
(
    IOobject
    (
        "n_filtered",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("n_filtered", dimless, vector::zero)//,
    // zeroGradientFvPatchVectorField::typeName
);

volScalarField epsilon1mask
(
    IOobject
    (
        "epsilon1mask",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilon1mask",dimless,0.0)
);

volScalarField epsilon1mask2
(
    IOobject
    (
        "epsilon1mask2",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("epsilon1mask2", dimless, 0.0)
);

volVectorField gradT
(
    IOobject
    (
        "gradT",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("gradT", dimensionSet(0, -1, 0, 1, 0), vector::zero)
);

volScalarField LatentHeat
(
    IOobject
    (
        "LatentHeat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*LatentHeat1 + (1.0 - alpha1)*LatentHeat2
);

volScalarField TRHS
(
    IOobject
    (
        "TRHS",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("TRHS", dimensionSet(1, -1, -3, 0, 0), 0.0)
);

volScalarField Tcorr
(
    IOobject
    (
        "Tcorr",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Tcorr", dimensionSet(0, 0, 0, 1, 0), 0.0)
);

volScalarField beta
(
    IOobject
    (
        "beta",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    alpha1*beta1 + (1.0-alpha1)*beta2
);

volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhok", dimless, 0.0)
);

volScalarField Num_divU
(
    IOobject
    (
        "Num_divU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("Num_divU", dimensionSet(0,0,-1,0,0), 0.0)
);

volVectorField Marangoni
(
    IOobject
    (
        "Marangoni",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("Marangoni", dimensionSet(1, -2, -2, 0, 0), vector::zero)
);

// --- Arc electromagnetic body force (WAAM-MHD thesis) ----------------------
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

const volVectorField zeroJcrossB
(
    IOobject
    (
        "zeroJcrossB",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("zeroJcrossB", dimensionSet(1, -2, -2, 0, 0), vector::zero)
);

const volScalarField zeroJoulesource
(
    IOobject
    (
        "zeroJoulesource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("zeroJoulesource", dimensionSet(1, -1, -3, 0, 0), 0.0)
);

// Arc pressure (Zhao Eq. 25). AUTO_WRITE deliberately, like arcSource: a force
// you cannot look at is a force you cannot verify. Expect max 227.8 Pa at
// I = 135 A, sigma_p = 2 mm.
volScalarField pArc
(
    IOobject
    (
        "pArc",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
);

// Recoil Pressure
volScalarField pVap
(
    IOobject
    (
        "pVap",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(1,-1,-2,0,0,0,0), 0.0)
);

// Evaporation Cooling
volScalarField Qv
(
    IOobject
    (
        "Qv",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimensionSet(1,0,-3,0,0,0,0), 0.0)
);


volScalarField damper
(
    IOobject
    (
        "damper",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("damper", dimless, 1.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField thermalDamper
(
    IOobject
    (
        "thermalDamper",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("thermalDamper", dimless, 1.0),
    zeroGradientFvPatchScalarField::typeName
);

volScalarField alpha_filtered
(
    IOobject
    (
        "alpha_filtered",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("alpha_filtered", dimensionSet(0, 0, 0, -0, 0), 0.0)
);

volScalarField alpha_smoothed
(
    IOobject
    (
        "alpha_smoothed",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("alpha_smoothed", dimensionSet(0, 0, 0, -0, 0), 0.0)
);

volScalarField e1temp
(
    IOobject
    (
        "e1temp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("e1temp", dimensionSet(0, 0, 0, -0, 0), 0.0)
);
e1temp = fvc::average(epsilon1);

volScalarField rhoCp
(
    IOobject
    (
        "rhoCp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoCp", dimensionSet(1, -1, -2, -1, 0), 0.0)
);
rhoCp = (rho*cp);
// volScalarField rhoCp(rho*cp);
rhoCp.oldTime();


// surfaceScalarField rhophicp(fvc::interpolate(cp)*rhoPhi);

surfaceScalarField rhophicp
(
    IOobject
    (
        "rhophicp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhophicp", dimensionSet(1, 2, -3, -1, 0), 0.0)
);
rhophicp = fvc::interpolate(cp)*rhoPhi;











// --- Mushy-zone Darcy constants (WAAM-MHD thesis) --------------------------
// Zhao Eq. 7 is this same sink. He never states A, only that "values between
// 10^4 and 10^7 are commonly recommended" -- three orders of magnitude, over
// which the distance metal at epsilon1 = 0.3 travels before freezing runs
// from 0.04 mm to 3.7 mm. Against a 2.3 mm width deficit that is a free
// parameter the size of the measurement, so it belongs in the case file and
// in the per-run record, not in a header.
//
// DEFAULTS ARE THE PREVIOUS LITERALS, so with no keys present every existing
// case is bit-identical. That is what Gate A and Gate B check.
//
// Note on darcy_q: Zhao Eq. 7 uses 0.001, laserbeamFoam 1e-12. It only enters
// where epsilon1^3 is comparable to it -- below epsilon1 ~ 0.1 -- and there
// the metal travels 0.003 mm either way. The difference is immaterial; the
// default stays at the solver's value so nothing moves, and the deviation is
// recorded rather than silently harmonised.
const dimensionedScalar DarcyConstantlarge
(
    "DarcyConstantlarge", dimensionSet(1,-3,-1,0,0),
    transportProperties.lookupOrDefault<scalar>("darcy_A", 1.0e6)
);

const dimensionedScalar DarcyConstantsmall
(
    "DarcyConstantsmall", dimless,
    transportProperties.lookupOrDefault<scalar>("darcy_q", 1.0e-12)
);



// volVectorField gradAlpha(fvc::grad(alpha1));

volVectorField gradAlpha
(
    IOobject
    (
        "gradAlpha",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("gradAlpha", dimensionSet(0, -1, 0, -0, 0), vector::zero)
);
gradAlpha = fvc::grad(alpha1);

// volVectorField nHatM(gradAlpha/(mag(gradAlpha) + deltaN));


volVectorField nHatM
(
    IOobject
    (
        "nHatM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("nHatM", dimensionSet(0, 0, 0, -0, 0), vector::zero)
);
nHatM = gradAlpha/(mag(gradAlpha) + deltaN);


// This variable is 1 if (alpha.metal >= 0.5 and epsilon1 >=0.5), 0 otherwise.
// These are the cells that are liquid and metallic at the current time
volScalarField liquidMetalCells
(
    IOobject
    (
        "liquidMetalCells",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("liquidMetalCells", dimensionSet(0, 0, 0, -0, 0), 0.0)
);


typedef incompressibleInterPhaseTransportModel
    <
        immiscibleIncompressibleTwoPhaseMixture
    > transportModelType;

autoPtr<transportModelType> turbulence
(
    new transportModelType
    (
        rho, U, phi, rhoPhi, mixture
    )
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"





autoPtr<mthdModel> mthd(nullptr);

if (mthdModel::electromagneticsEnabled(mesh, "transportProperties"))
{
    mthd.reset(new mthdModel(mesh, "transportProperties"));
}

const dimensionedScalar eleccond1 =
    mthdModel::readElectricalConductivity(transportPropertiesMetal);

volScalarField electrical_resistivity
(
    IOobject
    (
        "electrical_resistivity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    1.0/eleccond1
);




// Create laser heat source
laserHeatSource laser(mesh);

#include "updateYDim.H"
#include "updateProps.H"

// volScalarField depositionTemp
// (
//     IOobject
//     (
//         "depositionTemp",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     mesh,
//     dimensionedScalar("depositionTemp", dimensionSet(1, -1, -2, -1, 0), 0.0),
//     zeroGradientFvPatchScalarField::typeName
// );


// depositionTemp=laser.deposition();
// depositionTemp.correctBoundaryConditions();

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh*rhok
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
    p_rgh = p - rho*gh*rhok;
}

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

#include "createMRF.H"
#include "createFvOptions.H"


#include "isoAdvector/createPorosity.H"
isoAdvection advector(alpha1, phi, U);
