#!/usr/bin/env python3
"""Build ONE runnable WAAM arc case by hand. Task A1.

WHY THIS EXISTS
    templates/modelP holds only a mesh and a controlDict, so no generated case
    can start and the Goldak arc source -- written, compiled and verified
    arithmetically -- has never actually fired. This builds a single concrete
    case so that can be tested, before any of it is generalised back into the
    template (task A2).

WHAT IT DOES NOT DO
    No @tokens@. Every value is a real number. Debugging the physics setup and
    the token substitution at the same time is how you end up unable to tell
    which one is broken.

DONOR
    cases/00-tutorials/Plate2D-mhd. Proven on this build: Gate B ran it 6906
    steps with MHD, melting, latent heat, Darcy damping, VOF and Marangoni all
    active, holding div(muH) at 1e-12. fvSchemes and fvSolution are COPIED
    UNCHANGED from it -- those are validated settings and re-deriving them
    would be inventing risk.

THE PATCH PROBLEM THIS SOLVES
    The donor is 2D with patches  lowerWall / atmosphere / rightWall /
    leftWall / frontAndBack, two of them 'empty'. The modelP mesh is a 3D half
    domain with  inlet / outlet / symmetry / walls. Copying the donor's field
    files onto the modelP mesh fails on every single field. Every boundaryField
    below is therefore written fresh against the modelP patch names.

USAGE
    python3 make_waam_case.py            # build
    python3 make_waam_case.py --force    # overwrite an existing case
"""
import re, shutil, sys
from pathlib import Path

# ============================================================ parameters
HOME   = Path.home()
DONOR  = HOME / "thesis/cases/00-tutorials/Plate2D-mhd"
MATS   = HOME / "thesis/materials/transportProperties"
CASE   = HOME / "thesis/cases/00-tutorials/waam-arc-01"

# --- arc, from caseMatrix.csv row R00-SPSL-paper --------------------------
ARC_I, ARC_U, ETA = 135.0, 16.0, 0.70
ARC_POWER = ETA * ARC_U * ARC_I          # 1512.0 W
AF, AR, BB, CC    = 0.003, 0.006, 0.004, 0.004
TRAVEL            = 0.02                 # m/s  = 120 cm/min

# --- filler droplets, Edit 4 ----------------------------------------------
# The arc alone saturates at ~832 K: 17 K past the solidus, 74 K short of the
# liquidus. No pool. In GMAW the melting energy arrives inside superheated
# droplets, and these are the numbers that describe them.
#
# Radius and frequency are NOT free: Zhao Eq. 30 fixes them from mass
# conservation given the wire feed,
#     r = cbrt( 3*Vw*Rw^2 / (4*fd) ) = cbrt(3*0.15*(0.6e-3)^2/(4*30))
#       = 1.105 mm,  against 1.114 mm simulated in the paper, 0.88% error.
#
# v_impact is 0.9 m/s, NOT the 0.15 m/s wire feed. The paper measures it
# directly: the droplet is at ~0.3 m/s before necking and reaches 0.9 m/s
# after separation, accelerated by an electromagnetic pinch of order
# 1e6 N/m^3. Using the feed speed understates the impingement momentum 6x.
# It is the one number here imported rather than derived, so it deserves a
# sensitivity sweep: 0.15 / 0.9 / 1.5 m/s. Bead VOLUME cannot depend on it --
# mass conservation fixes that -- so the 8.5 mm^2 validation target is
# untouched and the sweep is free of side effects.
DROPLETS          = True
WIRE_D            = 1.2e-3               # m, ER5356
WIRE_V            = 0.15                 # m/s = 900 cm/min
DROP_R            = 1.105e-3             # Zhao Eq. 30
DROP_V            = 0.9                  # m/s, measured in the paper
DROP_T            = 1506.0               # K, Zhao Table 3 inlet
DROP_GAP          = 0.5e-3               # standoff above the free surface,
                                         # one cell: the sphere must not be
                                         # created inside existing metal
# --- half or full domain --------------------------------------------------
# "half"  y from 0 to +25 mm, a symmetryPlane at y = 0. Half the cells, so
#         half the cost. VALID ONLY for B = 0 or B along y.
# "full"  y from -25 to +25 mm, no symmetry plane, both y faces farField.
#         71,808 -> 143,616 cells. Required for any field direction that
#         breaks the mirror, which is longitudinal (x) or vertical (z).
#
# THE TRAP THIS EXISTS TO AVOID: a half domain with a symmetry-breaking field
# does not fail. It produces a smooth, plausible, completely wrong answer, and
# nothing in the output says so. The guard below refuses instead.
#
# The solver reads 'sym' from caseParameters and halves the arc power itself,
# so this one word also decides whether 756 W or 1512 W is delivered. Getting
# it wrong is a factor of two in the only energy input the model has.
SYM               = "full"
B_DIRECTION       = None                 # None, "x", "y" or "z" -- see guard

# Decomposition, pinned in the case rather than typed on a command line. It
# perturbs the arithmetic -- different summation order at every processor
# boundary -- and on a chaotic system that is a real perturbation, so an
# accidentally different -np is a silent extra variable in every comparison.
# Split along x: the torch travels in x over 60 mm, while y is where the pool
# sits, so splitting y would put a processor boundary through it.
# HPC: raise this. 143,616 cells wants ~20-50k per core, so 4 to 8 ranks on a
# laptop and (16 1 1) or (16 2 1) on a cluster node.
DECOMP            = (6, 1, 1)

# --- domain: sized from the physics, for a 0.5 s run ----------------------
# A truncated boundary must sit far enough out that the far field is
# undisturbed, or the answer is set by the boundary condition rather than by
# the arc. From the measured properties:
#
#   thermal diffusivity  a = k/(rho*cp) = 155/(2650*940) = 6.22e-5 m^2/s
#   diffusion length     L = sqrt(4*a*t) = 11.2 mm at 0.5 s   (22.3 mm at 2 s)
#
# The torch runs x = 0 -> 10 mm in 0.5 s and the Goldak ellipsoid spans
# [-ar, +af] = [-6, +3] mm about it, so the heated region is x in [-6, 13] mm.
# Boundaries are placed ~25 mm beyond that, a little over two diffusion
# lengths, and the zeroGradient / fixedValue-300 pair should then agree.
#
# z is DIFFERENT IN KIND. 6 mm is the substrate's real thickness, so the
# bottom is a physical surface, not a truncation -- it cannot be pushed
# further away, and its condition is a modelling choice rather than a
# convergence question. Adiabatic is defensible over 0.5 s; over a 2 s run
# the fixture contact matters and Zhao Table 3 should settle it.
CELL      = 0.5e-3          # FINE cell, through the pool. Production 0.25e-3.
X0, X1    = -0.015, 0.035   # 50 mm -- Zhao's coupon length
# y is mirrored for the full domain. The far boundary is Zhao's REAL plate
# edge: his coupon is 50 x 30 x 4 mm, so the long edges sit 15 mm from the
# bead. The previous 25 mm was chosen for domain independence, which is the
# right target when a boundary is an artificial truncation of a larger body
# and the wrong one when the boundary is the physical edge of the workpiece.
# Converging to the infinite-plate answer converges to a different experiment.
# The boundary IS thermally live at 15 mm, which is why Eq. 29 (h = 80) now
# sits on farField -- that condition is not a cost of trimming, it is what
# Zhao Table 3 specifies regardless of domain size.
# Consequence for reading pool_stats_waam.py: 'far met' is no longer a
# domain-size gate. It is a physical plate edge and is EXPECTED to warm. The
# gate is now that it stays well below the 815 K solidus.
Y1        =  0.015          # 30 mm full width -- Zhao's coupon
FY1       =  0.006          # fine to +-6 mm, one diffusion length sqrt(a*t)
Y0        = -Y1  if SYM == "full" else 0.0
FY0       = -FY1 if SYM == "full" else 0.0
# 4 mm substrate, 5 mm argon.
# SUBSTRATE: 4 mm is Zhao's actual baseplate. This is not a cost-saving. The
# conduction penetration depth over the ~0.7 s between the arc passing a
# cross-section and that metal freezing is sqrt(alpha*t) = 5.9 mm, so at 4 mm
# the bottom face is INSIDE the thermally active region and sets the local
# cooling rate during the window that decides bead height. At 6 mm the plate
# was 50% too heavy a heat sink. It also sets L in Zhao Eqs. 22-24.
# ARGON: multilayer is out of scope (SPSL only), so the 14 mm a five-layer
# wall would need is irrelevant. A ~2 mm bead plus a 1.1 mm droplet radius
# plus 0.5 mm standoff needs 3.6 mm of clearance; 5 mm leaves margin.
Z0, Z1    = -0.004, 0.005

# Region held at CELL. Everything outside grades away from it, so the mesh is
# fine where the pool is and coarse where only diffusion happens. Uniform
# 0.5 mm over this box would be 108,000 cells and uniform 0.25 mm 864,000;
# grading gives 53,856.
FX0, FX1  = -0.012, 0.032   # the WHOLE 40 mm bead, not just the pool.
                            # Was [-10,+20]: the arc reaches x = +30 mm over a
                            # 2.0 s pass, so the bead end and the arc-off
                            # collapse sat in 1.0-1.6 mm cells -- and that is
                            # exactly where 'deposit length' is measured
                            # (Zhao Table 5, 42.47 mm). Invisible until now
                            # because every run stopped at t <= 0.5 s, where
                            # the arc has only reached x = 0.
                            # FY0 / FY1 are set above, with the y extent, so
                            # the half/full mirror is decided in one place
EXPAND    = 2.0             # coarsest/finest cell within a graded section.
                            # Was 4.0. Trimming leaves only 3 mm of graded
                            # margin in x, and 4.0 over 3 cells means adjacent
                            # cells growing 2.0x, well past the ~1.2 guideline.
                            # 2.0 gives 1.19x adjacent for ~17k extra cells.

END_TIME, WRITE_EVERY, DT = 0.5, 0.02, 1e-5

SIGMA_E_GAS = 3.65             # metal 3.65e6 * the 1e-6 baseline ratio

# Solve the magnetic module, or not. WRITTEN EXPLICITLY, NEVER OMITTED.
#
# The materials file had no 'electromagnetics' entry at all, so the solver fell
# back on an internal default -- and that default is TRUE. Measured from the
# log of the run this comment was written during:
#
#     Solving for T        1 per timestep     the energy equation
#     Solving for p_rgh    3 per timestep     pressure
#     Solving for Hx/y/z   3 each             induction equation
#     Solving for pH       6 per timestep     magnetic divergence cleaning
#                         ---
#                          19 solves, 15 of them magnetic = 79%
#
# pH is elliptic, so two thirds of the EXPENSIVE solves were cleaning the
# divergence of a field that is identically zero. Phase 1 measured the cost of
# switching electromagnetics on at 2.7x, and that is what was being paid.
#
# It was left on originally for a good reason -- an inert MHD run that diverges
# is worth discovering before a field is ever applied. That question is now
# answered: it ran 12,531 steps at B = 0 alongside VOF, melting and Marangoni
# without diverging. So stop paying for the answer.
#
# With B = 0 the Lorentz force and Joule heating are identically zero, so this
# cannot change the physics. It will NOT be bit-identical -- removing 15 solves
# per step changes the arithmetic and this system amplifies that -- but it must
# agree within the noise floor, and that is checkable.
#
# Track B sets this true. So does Phase 7. It is a case-matrix column now.
ELECTROMAGNETICS = False

HEAD = """/*--------------------------------*- C++ -*----------------------------------*\\
  waam-arc-01 -- built by make_waam_case.py. Concrete values, no tokens.
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    location    "{loc}";
    object      {obj};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""

SYMM = {"type": "symmetryPlane"}


def symm():
    """The symmetry patch entry, or nothing at all on a full domain.

    Written as a splat (**symm()) into each field's patch dict rather than as
    a fixed key, so a full-domain case simply has no 'symmetry' patch instead
    of having one that names a face which no longer exists. blockMesh would
    reject the latter, but only after the case had been built.
    """
    return {} if SYM == "full" else {"symmetry": SYMM}


def walls(body):
    """The two halves of the old 'walls' patch, given the same condition.

    substrateBottom (z = z_min) is the plate's REAL underside: heat reaches it
    and its condition is a physical modelling choice. farField (y = y_max) is a
    TRUNCATION: nothing is there, the plate continues, and if the domain is big
    enough its condition should not matter.

    They shared one patch until now, which made the domain-independence check
    impossible -- varying the far side also varied the bottom, and the two
    effects landed in the same number with no way to separate them. Only T
    distinguishes them; every other field treats them alike.
    """
    return {"substrateBottom": dict(body), "farField": dict(body)}


def vol(obj, cls, dims, internal, patches, loc="0"):
    s = HEAD.format(cls=cls, loc=loc, obj=obj)
    s += f"dimensions      {dims};\n\ninternalField   uniform {internal};\n\n"
    s += "boundaryField\n{\n"
    for name, bodyd in patches.items():
        s += f"    {name}\n    {{\n"
        for k, v in bodyd.items():
            s += f"        {k:<16}{v};\n"
        s += "    }\n"
    s += "}\n\n// ************************************************************* //\n"
    return s


# ============================================================ fields
# Dimensions taken from the donor's own initial/ files, not guessed.
FIELDS = {
 "alpha.metal": vol("alpha.metal", "volScalarField", "[0 0 0 0 0 0 0]", "0", {
     # inletOutlet for the same reason as T. Academic (metal never reaches
     # z = +5 mm) but kept consistent so the patch is not two-way in one
     # field and clamped in another.
     "inlet":    {"type": "inletOutlet", "inletValue": "uniform 0",
                  "value": "uniform 0"},
     "outlet":   {"type": "inletOutlet", "inletValue": "uniform 0",
                  "value": "uniform 0"},
     **symm(),
     **walls({"type": "zeroGradient"})}),

 # inlet was fixedValue (0 0 0), which is hydrodynamically a LID: zero
 # velocity means zero through-flow, so the buoyant argon plume rises, cannot
 # leave, spreads along the roof and reaches farField -- a no-slip wall it
 # also cannot pass. That is why 'far gas' kept climbing even after gravity
 # was put on the correct axis. pressureInletOutletVelocity lets the plume
 # out while still admitting inflow where the pressure field wants it.
 # This is the MINIMUM fix. Zhao's actual top boundary (Table 3, Eq. 28) is a
 # 20 L/min annular argon jet directed downward between the wire and the
 # nozzle, which sweeps the plume out through the side pressure outlets.
 "U": vol("U", "volVectorField", "[0 1 -1 0 0 0 0]", "(0 0 0)", {
     "inlet":    {"type": "pressureInletOutletVelocity",
                  "value": "uniform (0 0 0)"},
     "outlet":   {"type": "pressureInletOutletVelocity",
                  "value": "uniform (0 0 0)"},
     **symm(),
     **walls({"type": "noSlip"})}),

 # Both substrateBottom and farField are now PHYSICAL surfaces of Zhao's
 # 50 x 30 x 4 mm coupon, exposed to air. Neither is a truncation any more,
 # so zeroGradient -- a perfect insulator -- is no longer defensible for
 # either. Zhao Table 3 puts Eq. 29 on every wall patch:
 #     -k dT/dn = h(T - Ta) + eps*sigma*(T^4 - Ta^4)
 # with h = 80 W/m2K. Note that 80 is high for natural convection off a small
 # plate (5-15 W/m2K would be typical); it is a lumped coefficient absorbing
 # the conduction path into the clamp and workbench. Say so in the writeup
 # rather than presenting it as a convection coefficient.
 # Magnitude: 21.4 cm2 of wall at dT = 300 K loses 51 W, which is 128 J over a
 # 2.5 s pass against 5402 J in -- about 2%. Locally under the arc it is
 # 0.032 W/mm2 against 30 W/mm2 arriving, a factor of 1000. So this will not
 # transform the SPSL answer; it is here because it is what the paper does and
 # because an adiabatic wall 15 mm from the bead reflects heat that should be
 # leaving.
 # kappaMethod lookup works only because patch_kappa_bc.py gave kappa real
 # boundary conditions; without that this fails at construction.
 "T": vol("T", "volScalarField", "[0 0 0 1 0 0 0]", "300", {
     # inletOutlet, NOT fixedValue. U here is now
     # pressureInletOutletVelocity, so the patch permits OUTFLOW and the
     # buoyant plume leaves through it. fixedValue would clamp that outgoing
     # gas to 300 K -- an artificial heat sink that also kills the buoyancy
     # driving the plume these runs are measuring. inletOutlet is zeroGradient
     # on outflow, 300 K on inflow: the correct partner for a two-way velocity
     # condition, and what 'outlet' has always used for that reason.
     "inlet":    {"type": "inletOutlet", "inletValue": "uniform 300",
                  "value": "uniform 300"},
     "outlet":   {"type": "inletOutlet", "inletValue": "uniform 300",
                  "value": "uniform 300"},
     **symm(),
     # written separately, not via walls(), because they remain independently
     # variable for sensitivity runs even though both now carry Eq. 29
     "substrateBottom": {"type": "externalWallHeatFluxTemperature",
                         "mode": "coefficient", "h": "uniform 80",
                         "Ta": "uniform 300", "kappaMethod": "lookup",
                         "kappa": "kappa", "value": "uniform 300"},
     "farField":        {"type": "externalWallHeatFluxTemperature",
                         "mode": "coefficient", "h": "uniform 80",
                         "Ta": "uniform 300", "kappaMethod": "lookup",
                         "kappa": "kappa", "value": "uniform 300"}}),

 "p_rgh": vol("p_rgh", "volScalarField", "[1 -1 -2 0 0 0 0]", "0", {
     "inlet":    {"type": "totalPressure", "p0": "uniform 0",
                  "value": "uniform 0"},
     "outlet":   {"type": "totalPressure", "p0": "uniform 0",
                  "value": "uniform 0"},
     **symm(),
     **walls({"type": "fixedFluxPressure", "value": "uniform 0"})}),

 # B = 0 for this test. The magnetic module stays active and must still behave
 # -- an inert MHD run that diverges is worth knowing about before a field is
 # ever applied.
 "H": vol("H", "volVectorField", "[0 -1 0 0 0 0 0]", "(0 0 0)", {
     "inlet":    {"type": "fixedValue", "value": "uniform (0 0 0)"},
     "outlet":   {"type": "fixedValue", "value": "uniform (0 0 0)"},
     **symm(),
     **walls({"type": "fixedValue", "value": "uniform (0 0 0)"})}),

 # Divergence cleaning. MUST be pinned somewhere (the fixedValue on outlet) or
 # the pressure-like system for pH is singular and the solve will not converge.
 "pH": vol("pH", "volScalarField", "[2 2 -5 0 0 0 0]", "0", {
     "inlet":    {"type": "zeroGradient"},
     "outlet":   {"type": "fixedValue", "value": "uniform 0"},
     **symm(),
     **walls({"type": "zeroGradient"})}),

 # Ray-tracer orientation marker: +1 atmosphere side, -1 metal side. The laser
 # object is constructed unconditionally even at zero power, so this must exist.
 "Laser_boundary": vol("Laser_boundary", "volScalarField", "[0 0 0 0 0 0 0]",
                       "0", {
     "inlet":    {"type": "fixedValue", "value": "uniform 1.0"},
     "outlet":   {"type": "zeroGradient"},
     **symm(),
     **walls({"type": "fixedValue", "value": "uniform -1.0"})}),
}

# ============================================================ constant/
import math as _math
WIRE_VDOT = _math.pi*(WIRE_D/2.0)**2*WIRE_V      # 1.696e-7 m3/s, full domain

DROPLET_BLOCK = f"""
// Droplets -- Edit 4. The PRESENCE of droplet_radius_m switches them on, the
// same convention the arc uses. Delete these five lines and the case runs
// arc-only, bit-comparably with everything before this edit.
droplet_radius_m       {DROP_R};        // Zhao Eq. 30 from the wire feed
droplet_v_impact_m_s   {DROP_V};        // measured after separation, not fed
droplet_T_K            {DROP_T};        // Zhao Table 3 velocity-inlet
wire_volume_rate_m3_s  {WIRE_VDOT:.6g};   // pi*(d/2)^2*Vw, halved by 'sym'
droplet_standoff_m     {DROP_GAP};        // clear of existing metal
""" if DROPLETS else ""

CASE_PARAMS = HEAD.format(cls="dictionary", loc="constant",
                          obj="caseParameters") + f"""// Read by createFields.H with READ_IF_PRESENT. The presence of goldak_af_m
// is what switches the arc on; without this file arcSource stays identically
// zero and the run is bit-identical to stock laserbeamFoam.

arc_power_W       {ARC_POWER};      // eta * U * I = {ETA} * {ARC_U} * {ARC_I}
goldak_af_m       {AF};             // forward semi-axis, narrow
goldak_ar_m       {AR};             // rear semi-axis, wide (torch is moving)
goldak_b_m        {BB};             // half width, y
goldak_c_m        {CC};             // depth, z -- HALF ellipsoid below surface
travel_speed_m_s  {TRAVEL};
arc_x0_m          -0.010;           // torch start. Domain is x[-15,+35]
                                      // (50 mm, Zhao coupon); a 40 mm bead at
                                      // 20 mm/s runs -10 -> +30, leaving 5 mm
                                      // lead-in and 5 mm run-out as Zhao has.
                                      // Was 0, which put arc-off at +40 mm --
                                      // 5 mm PAST the plate edge.
arc_y0_m          0;                // on the symmetry plane
arc_z0_m          0;                // substrate top surface
sym               {SYM};            // {SYM} domain -> target {ARC_POWER*(0.5 if SYM=="half" else 1.0):g} W
{DROPLET_BLOCK}"""

G_DICT = HEAD.format(cls="uniformDimensionedVectorField", loc="constant",
                     obj="g") + """// WRITTEN, NOT COPIED, and that distinction is the point of this file.
//
// This was copied verbatim from the donor, and the donor is Plate2D-mhd: a
// TWO-DIMENSIONAL x-y case in which y is vertical, so its g is (0 -9.81 0).
//
// This mesh is z-up and travel is +x. Everything else already agrees:
//     setFieldsDict          fills z < 0 with metal, argon above
//     caseParameters         arc_z0_m 0 is "substrate top surface"
//     blockMeshDict          inlet is the z-max face, substrateBottom z-min
//     updateArcSource.H      skips cells with (Cc.z() - arc_z0_) > 0
//     updateDropletSource.H  applies the impact force along (0 0 -1)
//
// Copying the donor's g pointed gravity along -y, HORIZONTAL here. It reaches
// the momentum equation through ghf*fvc::snGrad(rho*rhok) in UEqn.H, so the
// failure was not a crash but a plausible-looking wrong answer: the crown
// drifted off the centreline, Boussinesq buoyancy drove the argon plume 25 mm
// sideways into farField instead of upward, and the deposit never spread
// under its own weight.
//
// Re-meshed on a different vertical axis? Edit this. Do not make it a copy.

dimensions      [0 1 -2 0 0 0 0];
value           (0 0 -9.81);

// ************************************************************* //
"""

LASER_OFF = HEAD.format(cls="dictionary", loc="constant",
                        obj="LaserProperties") + """// WAAM is heated by the arc, not a laser. The laser object is constructed
// unconditionally by the solver, so this dictionary must exist -- but every
// power in timeVsLaserPower is zero.
radialPolarHeatSource no;

lasers
(
    laser0
    {
        timeVsLaserPosition
        {
            file            "$FOAM_CASE/constant/timeVsLaserPosition";
            outOfBounds     clamp;
        }
        timeVsLaserPower
        {
            file            "$FOAM_CASE/constant/timeVsLaserPower";
            outOfBounds     clamp;
        }
        V_incident      (0 0 -1);   // straight down; irrelevant at zero power
        laserRadius     0.0005;
        wavelength       1.064e-6;
        e_num_density   5.83e29;
    }
);
"""

# ============================================================ system/
def _graded_n(L, d0, R):
    """Cells to span L, starting at d0, with last/first cell ratio R."""
    n = 2
    while n < 1000:
        q = R**(1.0/(n - 1))
        span = d0*(q**n - 1)/(q - 1) if abs(q - 1) > 1e-12 else d0*n
        if span >= L:
            return n
        n += 1
    sys.exit(f"ERROR: cannot grade {L*1e3:g} mm from {d0*1e3:g} mm at ratio {R}")


def blockmesh():
    Lx, Ly, Lz = X1 - X0, Y1 - Y0, Z1 - Z0

    n_xl = _graded_n(FX0 - X0, CELL, EXPAND)
    n_xf = round((FX1 - FX0)/CELL)
    n_xr = _graded_n(X1 - FX1, CELL, EXPAND)
    # y: two sections on a half domain (fine at the symmetry plane, grading
    # out), three on a full one (grading in, fine, grading out) -- exactly the
    # x pattern mirrored, so the mesh is symmetric about y = 0 to the last bit.
    # That symmetry matters: on a full domain with B = 0 the SOLUTION must be
    # mirror-symmetric too, and any asymmetry is accumulated numerical error.
    # An asymmetric mesh would make that check meaningless.
    n_yl = _graded_n(FY0 - Y0, CELL, EXPAND) if SYM == "full" else 0
    n_yf = round((FY1 - FY0)/CELL)
    n_yc = _graded_n(Y1 - FY1, CELL, EXPAND)
    nx, ny, nz = n_xl + n_xf + n_xr, n_yl + n_yf + n_yc, round(Lz/CELL)

    # multiGrading on a single block: each triple is
    #   (fraction of edge length, fraction of cells, last/first cell ratio)
    # The leading x section runs from the COARSE end toward the fine region,
    # so its ratio is 1/EXPAND -- cells shrink along it.
    gx = (f"( ({(FX0-X0)/Lx:.6g} {n_xl/nx:.6g} {1/EXPAND:g}) "
          f"({(FX1-FX0)/Lx:.6g} {n_xf/nx:.6g} 1) "
          f"({(X1-FX1)/Lx:.6g} {n_xr/nx:.6g} {EXPAND:g}) )")
    if SYM == "full":
        gy = (f"( ({(FY0-Y0)/Ly:.6g} {n_yl/ny:.6g} {1/EXPAND:g}) "
              f"({(FY1-FY0)/Ly:.6g} {n_yf/ny:.6g} 1) "
              f"({(Y1-FY1)/Ly:.6g} {n_yc/ny:.6g} {EXPAND:g}) )")
    else:
        gy = (f"( ({(FY1-Y0)/Ly:.6g} {n_yf/ny:.6g} 1) "
              f"({(Y1-FY1)/Ly:.6g} {n_yc/ny:.6g} {EXPAND:g}) )")
    gz = "1"

    # On a full domain the y = y_min face is no longer a mirror -- it is the
    # same kind of truncation as y = y_max, so both join farField and the
    # symmetryPlane patch disappears entirely.
    if SYM == "full":
        ysym = ""
        yfar = "((3 2 6 7) (0 1 5 4))"
    else:
        ysym = "    symmetry        { type symmetryPlane; faces ((0 1 5 4)); }\n"
        yfar = "((3 2 6 7))"

    s = HEAD.format(cls="dictionary", loc="system", obj="blockMeshDict")
    s += f"""convertToMeters 1;

vertices
(
    ({X0} {Y0} {Z0})
    ({X1} {Y0} {Z0})
    ({X1} {Y1} {Z0})
    ({X0} {Y1} {Z0})
    ({X0} {Y0} {Z1})
    ({X1} {Y0} {Z1})
    ({X1} {Y1} {Z1})
    ({X0} {Y1} {Z1})
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({nx} {ny} {nz})
    simpleGrading
    (
        {gx}
        {gy}
        {gz}
    )
);

edges ();

boundary
(
    inlet           {{ type patch;         faces ((4 5 6 7)); }}
    outlet          {{ type patch;         faces ((0 3 7 4) (1 2 6 5)); }}
{ysym}    substrateBottom {{ type wall;          faces ((0 1 2 3)); }}
    farField        {{ type wall;          faces {yfar}; }}
);

mergePatchPairs ();

// ************************************************************* //
"""
    return s, nx*ny*nz, (nx, ny, nz)


SETFIELDS = HEAD.format(cls="dictionary", loc="system", obj="setFieldsDict") + """// Everything starts as argon; every cell below z = 0 becomes metal. The box is
// deliberately larger than the domain in x and y so the substrate is complete.
defaultFieldValues
(
    volScalarFieldValue alpha.metal 0
);

regions
(
    boxToCell
    {
        box (-1 -1 -1) (1 1 0);
        fieldValues
        (
            volScalarFieldValue alpha.metal 1
        );
    }
);

// ************************************************************* //
"""

# 4 subdomains split along x only. NOT a free choice: every comparison in this
# study is between runs that must differ in exactly one thing, and the
# decomposition perturbs the arithmetic (different summation order at every
# processor boundary). On a system where a 1-ULP difference amplifies to tens
# of percent, that is a real perturbation, so it is pinned here in the case
# rather than typed on a command line where it can silently change between
# runs. Split along x because the arc travels in x over a 60 mm domain; y is
# only 25 mm and all the physics sits against the symmetry plane, so splitting
# y would put a processor boundary through the pool.
# .format() is used on this, so every literal brace is doubled. The FoamFile
# header's braces are literal; only {NPROC} and {DX}/{DY}/{DZ} are fields.
DECOMPOSE_T = """FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      decomposeParDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

numberOfSubdomains  {NPROC};

method              simple;

simpleCoeffs
{{
    n               ({DX} {DY} {DZ});
    delta           0.001;
}}

// ************************************************************* //
"""

ALLRUN = """#!/bin/sh
cd "${0%/*}" || exit 1
. $WM_PROJECT_DIR/bin/tools/RunFunctions

# The fields are written to initial/, not 0/, so that a re-run always starts
# from a pristine copy instead of from whatever the last run left behind.
# THIS COPY IS NOT OPTIONAL AND NOT COSMETIC. Without it:
#   setFields  warns "Field alpha.metal not found" and does nothing, so the
#              whole domain stays argon and the arc finds no metal to heat
#   postProcess falls back to writing into constant/ instead of 0/, and
#              pool_stats_waam.py -- which only looks in time directories --
#              then reports no mesh at all
# Both failures are warnings, not errors, so the case runs and produces
# confident nonsense. Use this script rather than typing the steps by hand.
# REFUSE BEFORE DESTROYING. runApplication skips any stage whose log already
# exists, so on a second invocation the stages below all decline -- but the
# rm above them had already deleted 0/ and the processor directories. The
# case was left wiped AND unrun, which is the worst of both. Check first.
# -f FIRST. Previously the guard below ran first and rejected everything,
# including the -f that exists to get past it, so the flag was unreachable.
if [ "$1" = "-f" ]; then
    echo "Cleaning: time directories, processor*, 0/ and logs"
    foamListTimes -rm > /dev/null 2>&1
    rm -rf 0 processor* log.*
fi

if ls log.* >/dev/null 2>&1; then
    echo "ERROR: this case has already been run -- log.* files are present."
    echo "  NOTHING HAS BEEN TOUCHED."
    echo
    echo "  To re-run it from scratch (deletes 0/, processor*/, log.* and"
    echo "  every time directory):"
    echo "      ./Allrun -f"
    echo
    echo "  To keep this run, copy the case aside first."
    exit 1
fi

echo "Copying 'initial' to 0"
rm -rf 0 processor*
cp -r initial 0

runApplication blockMesh
runApplication setFields
runApplication decomposePar
runParallel laserbeamFoam
runApplication reconstructPar

# Cell centres and volumes LAST. writeCellCentres also writes a field called
# 'C', and C is the mesh cell-centre symbol in OpenFOAM; putting it into 0/
# before decomposePar means decomposing a field with a reserved name. The
# reference run generated these after the solve and is proven, so keep that
# order. They describe the mesh, so once per mesh is enough.
runApplication -s cellCentres postProcess -func writeCellCentres -time 0
runApplication -s cellVolumes postProcess -func writeCellVolumes -time 0
"""


# ============================================================ build
def main():
    force = "--force" in sys.argv

    # A half domain is a mirror about y = 0. A field along x or z is not
    # mirror-symmetric, so imposing one on a half domain silently solves a
    # different problem -- and the result is smooth and plausible, which is
    # what makes it dangerous. Refuse rather than assume.
    if SYM == "half" and B_DIRECTION in ("x", "z"):
        sys.exit(f"ERROR: B along {B_DIRECTION} breaks the mirror at y = 0, so "
                 f"a half domain cannot represent it.\n"
                 f"  Set SYM = \"full\". The answer would look completely "
                 f"reasonable and be wrong.")
    if SYM not in ("half", "full"):
        sys.exit(f"ERROR: SYM is {SYM!r}; it must be \"half\" or \"full\". "
                 f"The solver reads it to decide whether to deliver "
                 f"{ARC_POWER*0.5:g} W or {ARC_POWER:g} W.")
    for p, what in ((DONOR, "donor case"), (MATS, "property set")):
        if not p.exists():
            sys.exit(f"ERROR: {p} not found ({what})")
    if CASE.exists():
        if not force:
            sys.exit(f"ERROR: {CASE} exists. Re-run with --force to replace it.")
        shutil.rmtree(CASE)

    (CASE / "initial").mkdir(parents=True)
    (CASE / "constant").mkdir()
    (CASE / "system").mkdir()

    for name, text in FIELDS.items():
        (CASE / "initial" / name).write_text(text)

    # --- properties: substitute the one token the file still carries ---------
    props = MATS.read_text()
    if "@sigma_e_gas@" in props:
        props = props.replace("@sigma_e_gas@", str(SIGMA_E_GAS))
    elif "@" in props:
        sys.exit("ERROR: transportProperties still holds an unexpected @token@. "
                 "Substitute it before this case can run.")

    # --- the magnetic switch, written explicitly ----------------------------
    # Omitting it is not the same as setting it false: the solver defaults to
    # TRUE and silently spends 79% of its linear solves on a zero field. This
    # writes the entry whether or not the materials file carries one, and
    # verifies afterwards, because "I assumed it was off" is how 100 minutes
    # went into divergence-cleaning a field that does not exist.
    em = "true" if ELECTROMAGNETICS else "false"
    if re.search(r"^\s*electromagnetics\s+\S+;", props, re.M):
        props = re.sub(r"^\s*electromagnetics\s+\S+;",
                       f"electromagnetics {em};", props, count=1, flags=re.M)
        how = "replaced in the materials file"
    else:
        # Insert after the FoamFile header's separator rather than appending,
        # so it reads as a deliberate setting and not as an afterthought.
        m = re.search(r"^//\s*\*.*$", props, re.M)
        line = f"\nelectromagnetics {em};   // Track A: no field. Track B sets true.\n"
        if m:
            props = props[:m.end()] + "\n" + line + props[m.end():]
            how = "inserted after the header"
        else:
            props = props.rstrip() + "\n" + line
            how = "appended (no header separator found)"

    n = len(re.findall(r"^\s*electromagnetics\s+\S+;", props, re.M))
    if n != 1:
        sys.exit(f"ERROR: transportProperties ended up with {n} "
                 f"'electromagnetics' entries. Refusing to write a dictionary "
                 f"whose most expensive switch is ambiguous.")

    (CASE / "constant" / "transportProperties").write_text(props)

    (CASE / "constant" / "caseParameters").write_text(CASE_PARAMS)
    (CASE / "constant" / "LaserProperties").write_text(LASER_OFF)
    (CASE / "constant" / "timeVsLaserPower").write_text("(\n(0    0)\n(1000 0)\n)\n")
    (CASE / "constant" / "timeVsLaserPosition").write_text(
        "(\n(0    (0 0 1))\n(1000 (0 0 1))\n)\n")

    # constant/g is WRITTEN, never copied: the donor is 2D and y-up, this mesh
    # is z-up, and copying its g put gravity on a horizontal axis in every WAAM
    # case this script built. Only geometry-FREE settings may be copied below.
    (CASE / "constant" / "g").write_text(G_DICT)

    # --- proven settings: copied, never re-derived ---------------------------
    # NOTHING GEOMETRY-DEPENDENT IN THIS TUPLE. Re-adding "constant/g" here is
    # the specific mistake this comment exists to prevent.
    copied = []
    for rel in ("constant/turbulenceProperties",
                "system/fvSchemes", "system/fvSolution"):
        src = DONOR / rel
        if not src.exists():
            sys.exit(f"ERROR: {src} missing -- cannot proceed without it")
        shutil.copy2(src, CASE / rel); copied.append(rel)

    bm, ncells, dims = blockmesh()
    (CASE / "system" / "blockMeshDict").write_text(bm)

    # controlDict is COPIED from the donor and edited, never written fresh.
    # Rewriting it is how adjustTimeStep / maxCo / maxAlphaCo / maxDeltaT went
    # missing on the first attempt: they sit at the bottom of the donor's file,
    # they are easy to overlook, and omitting maxAlphaCo aborts the run on the
    # first timestep. Same principle as fvSchemes and fvSolution -- if the donor
    # has a working version, copy it.
    ctrl = (DONOR / "system/controlDict").read_text()
    for key, val in (("endTime", END_TIME), ("writeInterval", WRITE_EVERY),
                     ("deltaT", DT), ("maxDeltaT", DT)):
        pat = re.compile(rf"^{key}\s+\S+;", re.M)
        if len(pat.findall(ctrl)) != 1:
            sys.exit(f"ERROR: expected exactly one '{key}' line in the donor "
                     f"controlDict, found {len(pat.findall(ctrl))}")
        ctrl = pat.sub(f"{key:<15} {val};", ctrl, count=1)
    for must in ("adjustTimeStep", "maxCo", "maxAlphaCo", "maxDeltaT"):
        if must not in ctrl:
            sys.exit(f"ERROR: donor controlDict has no '{must}' -- the solver "
                     f"needs it and this script will not invent one")
    # libs: externalWallHeatFluxTemperature (Zhao Eq. 29) is NOT in
    # laserbeamFoam's default patchField table. VERIFIED on this machine:
    # log.laserbeamFoam gave "Unknown patchField type
    # externalWallHeatFluxTemperature for patch type wall", and the printed
    # table of valid types had zero entries matching 'external' or
    # 'wallHeatFlux'. The BC lives in
    #   src/thermoTools/derivedFvPatchFields/externalWallHeatFluxTemperature/
    # and ships as libthermoTools.so, which loads cleanly at run time.
    # MUST be written here, not added by hand: --force rewrites controlDict
    # from the donor, so a hand-added libs line vanishes on the next rebuild
    # and the run dies at startup. Exactly the constant/g failure mode.
    if "libs" not in ctrl:
        pat = re.compile(r"^application\s+\S+;", re.M)
        if len(pat.findall(ctrl)) != 1:
            sys.exit("ERROR: expected exactly one 'application' line in the "
                     "donor controlDict to anchor the libs entry")
        ctrl = pat.sub('libs            ("libthermoTools.so");\n\n'
                       + pat.findall(ctrl)[0], ctrl, count=1)

    (CASE / "system" / "controlDict").write_text(ctrl)
    (CASE / "system" / "setFieldsDict").write_text(SETFIELDS)
    (CASE / "system" / "decomposeParDict").write_text(
        DECOMPOSE_T.format(NPROC=DECOMP[0]*DECOMP[1]*DECOMP[2],
                           DX=DECOMP[0], DY=DECOMP[1], DZ=DECOMP[2]))
    run = CASE / "Allrun"; run.write_text(ALLRUN); run.chmod(0o755)

    # Read the patch names back out of a field that was actually written,
    # rather than restating them. The hardcoded version of this line still
    # said "walls" after the patch split and reported a case that no longer
    # existed -- a summary that cannot be wrong is worth the three lines.
    plist = "/".join(re.findall(r"^    (\w+)$", FIELDS["U"], re.M))

    print(f"  built {CASE}")
    print(f"    initial/   {len(FIELDS)} fields, patches {plist}")
    print(f"    magnetic   electromagnetics {em}   ({how})")
    print(f"    constant/  transportProperties (sigma_e_gas = {SIGMA_E_GAS}), "
          f"caseParameters, laser OFF")
    print(f"    gravity    (0 0 -9.81) written, NOT copied -- mesh is z-up")
    print(f"    copied     {', '.join(copied)}")
    print(f"    mesh       {ncells:,} cells  {dims[0]}x{dims[1]}x{dims[2]}  graded, {CELL*1e3:g} mm through the pool")
    print(f"    arc        {ARC_POWER:g} W, {SYM} domain "
          f"-> target {ARC_POWER*(0.5 if SYM=='half' else 1.0):g} W")
    print(f"    system/    decomposeParDict  {DECOMP[0]*DECOMP[1]*DECOMP[2]} subdomains, simple {DECOMP}")
    print(f"""
  Run it -- USE Allrun, do not type the steps by hand. It copies initial/ to
  0/ first; skipping that makes setFields a no-op and the whole domain stays
  argon, with nothing but a warning to say so.
      of2506
      cd {CASE}
      ./Allrun

  GATE -- what to check in log.laserbeamFoam:
    1.  "Arc: Goldak double ellipsoid"  with af=3mm ar=6mm b=4mm c=4mm
    2.  "Arc: ... W in metal below the surface, renormalised x..."
        raw should be NEAR 100% of target. Well below means the ellipsoid is
        hanging off the plate -- widen the domain, do not adjust the power.
    3.  "Metal kappa: two-branch"  -- confirms Edit 1 is live
    4.  max(T) climbing but staying under ~1600 K. Vaporisation is disabled
        (Tvap = 1e6) so runaway temperature has nowhere to go but up.
    5.  It reaches {END_TIME} s without the timestep collapsing.""")


if __name__ == "__main__":
    main()
