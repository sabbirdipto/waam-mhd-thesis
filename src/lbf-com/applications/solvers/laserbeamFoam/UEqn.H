MRF.correctBoundaryVelocity(U);


if (laser.powderSim())
{
    DC =
        DarcyConstantlarge*Foam::pow((1.0-epsilon1mask),2)
       /(Foam::pow(epsilon1mask, 3) + DarcyConstantsmall);
}
else
{
    DC =
        DarcyConstantlarge*Foam::pow((1.0 - epsilon1), 2)
       /(Foam::pow(epsilon1, 3) + DarcyConstantsmall);
}

DC.correctBoundaryConditions();  // ADDED

gradAlpha = fvc::grad(alpha1);
gradAlpha.correctBoundaryConditions();  // ADDED - CRITICAL!

nHatM = gradAlpha/(mag(gradAlpha) + deltaN);
nHatM.correctBoundaryConditions();  // ADDED

Marangoni =
    epsilon1*Marangoni_Constant*(gradT - nHatM*(nHatM & gradT))*mag(gradAlpha);

// The liquidus cliff. epsilon1 already ramps this force from zero at the
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
}

Marangoni.correctBoundaryConditions();

if (damperSwitch)
{
    // To account for high density difference
    damper = 2.0*rho/(rho1 + rho2);
    damper.correctBoundaryConditions();  // ADDED
}

pVap = 0.54*p0*Foam::exp(LatentHeatVap*Mm*((T - Tvap)/(R*T*Tvap)));
pVap.correctBoundaryConditions();  // ADDED

// Computed HERE, not beside updateArcSource.H, because that include sits AFTER
// UEqn.H in laserbeamFoam.C while pEqn.H -- where the force actually takes
// effect, momentumPredictor being no -- runs later still. Here it is current
// for both.
#include "updateArcPressure.H"


if (mthd.valid())
{
    mthd->updateLorentzForce();
}

const volVectorField& JcrossB =
    mthd.valid()
  ? mthd->JcrossB()
  : zeroJcrossB;

fvVectorMatrix UEqn
(
    fvm::ddt(rho, U) + fvm::div(rhoPhi, U)
  + MRF.DDt(rho, U)
  + turbulence->divDevRhoReff(rho, U)
  + fvm::Sp(DC, U)
  - Marangoni*damper
  - JcrossB*damper
    // Arc self-field, Zhao Eqs. 22-24 derived. Zero unless
    // arc_lorentz is set; never double-counts with JcrossB,
    // which is the numerical potential-solve path.
  - arcLorentz*damper
  - dropletForce*damper
 ==
    fvOptions(rho, U)
);

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
                mixture.surfaceTensionForce()*fvc::interpolate(damper)
                + (fvc::interpolate(pVap)
                    *fvc::snGrad(alpha1)
                    *fvc::interpolate(damper))
              - ghf*fvc::snGrad(rho*rhok)
              - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        )
    );

    fvOptions.correct(U);
}
