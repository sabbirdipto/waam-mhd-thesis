/*--------------------------------*- C++ -*----------------------------------*\
  waam-arc-01 -- built by make_waam_case.py. Concrete values, no tokens.
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      alpha.metal;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    inlet
    {
        type            fixedValue;
        value           uniform 0;
    }
    outlet
    {
        type            inletOutlet;
        inletValue      uniform 0;
        value           uniform 0;
    }
    symmetry
    {
        type            symmetryPlane;
    }
    substrateBottom
    {
        type            zeroGradient;
    }
    farField
    {
        type            zeroGradient;
    }
}

// ************************************************************* //
