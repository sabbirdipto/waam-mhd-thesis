# Zhao 2021 vs our model: every difference, and which ones matter

**Date:** 2026-09-23. Consolidates comparisons scattered across `M1`, `M3`, `M4`,
`shielding-gas-eq28` and the 0.25 mm plan, plus what the paper itself says.

Reference: Zhao, Wei, Long, Chen, Liu, Ou, *Welding in the World* (2021)
65:1571-1590, "Modeling and simulation of heat transfer, fluid flow and
geometry morphology in GMAW-based wire arc additive manufacturing".

**Claims below are sourced.** Where the paper is silent this says so, because
an unpublished parameter is a different problem from a disagreement.

---

## 1. What matches, exactly

These were deliberately taken from the paper and are not in question.

| | value | source |
|---|---|---|
| alloy | 5A05 substrate, ER5356 wire, properties assumed identical | §2, §3.6 |
| current, voltage | 135 A, 16 V | Table 2 |
| travel speed | 20 mm/s (120 cm/min) | Table 2 |
| wire feed | 150 mm/s (900 cm/min), 1.2 mm diameter | Table 2 |
| shielding gas | 20 L/min argon | Table 2 |
| substrate | 50 x 30 x 4 mm | §2 (errata 2) |
| solidus / liquidus | 815 / 906 K | Table 4 |
| density | 2650 kg/m^3 | Table 4 |
| surface tension, d(sigma)/dT | 0.85 N/m, -1.55e-4 N/m.K | Table 4 |
| latent heat | 3.58e5 J/kg | Table 4 |
| thermal expansion | 3e-5 /K | Table 4 |
| heat transfer coefficient | 80 W/m^2.K | Table 4 |
| Goldak semi-axes af, ar, b, c | 3, 6, 4, 4 mm | §3.5 |
| arc-pressure / current radii | sigma_p = sigma_j = 2 mm | Table 4 |
| droplet temperature | 1506 K | §3.5, Table 3 |
| droplet radius | 1.105 mm (Eq. 30) | §4.1 |
| **pool momentum terms (Eq. 21)** | gravity, buoyancy, EM, surface tension, Marangoni, arc pressure | §3.4 |

All six of Eq. 21's momentum sources are implemented. That is the core of the
model and it is matched.

## 2. Differences forced by the solver

Not choices; consequences of OpenFOAM not being Fluent.

| | Zhao | ours |
|---|---|---|
| code | Fluent + self-coded UDFs (§3.7) | OpenFOAM v2506, laserbeamFoam vendored |
| interface | VOF, Torrey lineage [25] (§3.2) | **isoAdvector + MULES**, geometric reconstruction |
| pressure-velocity | PISO (§3.7) | PIMPLE, 3 correctors |
| phase change | enthalpy-porosity, Voller & Prakash [26] | same method |
| linear solvers | **not stated** | DICPCG (M2/G1); GAMG-PCG on R0 |

isoAdvector is a sharper interface scheme than Fluent's default VOF. It should
smear less -- which matters directly for the creep discussed in §6.

## 3. Differences by our choice, with reasons

### 3.1 Arc power: 1113.7 W against his 1512 W

Zhao §3.5 states eta = 0.70 is "the arc efficiency for heating molten pool,
**excluding heating droplets**" -- so pool-only. But Table 4 sources 0.70 to
ref [29], a welding textbook, where it is a **total** GMAW efficiency
(**errata 5**). Applied his way: pool 1512 W + droplet enthalpy 703-801 W =
**2215-2313 W against U.I = 2160 W, i.e. 103-107%**.

We split it: eta_total = 0.85, eta_2 = 0.334 to droplets (Zhu [32] Eq. 9),
leaving **eta_1 = 0.516 -> 1113.7 W** to the arc, plus 722 W as droplet mass.

**This was tested, not assumed.** M4 ran Zhao's literal 1512 W: the plate
underside reached 846-859 K against an 815 K solidus at every write and was
still climbing at ~94 K/s. His own 4 mm coupon would have melted through
before his 2.0 s track ended, and his macrograph shows no melt-through.

### 3.2 Half-ellipsoid with renormalisation

Eq. 27 as printed carries the 6sqrt(3)/(abc.pi.sqrt(pi)) prefactor with
f_f + f_r = 2, which integrates to eta.U.I **only if restricted to z < z0**.
The printed equation states only the x branches (**errata 3**). Applied
literally it delivers 2.eta.U.I, half of it into argon.

Ours gates on `if (dz > 0.0) continue` and renormalises to exactly eta.U.I in
metal. The measured factor is 1.016-1.031 and falling as the bead builds --
i.e. 97-98.4% of the source lands in metal before correction.

### 3.3 Full domain, no symmetry plane

Zhao uses a half-domain with XZ symmetry (§3, Fig. 3, boundaries CDHG/HGNM).
We run the full 50 x 30 mm (`sym full`). Costs 2x the cells; removes any
question about symmetry suppressing asymmetric flow -- relevant later, since a
transverse magnetic field breaks that symmetry by construction.

### 3.4 Single layer only

Zhao simulates SPSL, SPTL (two interlayer-idle cases) and SPFL (five layers),
using **dynamic meshing** to grow the air subdomain per layer (§3). We do
SPSL only, on a fixed mesh. The multilayer cases are not part of this thesis.

### 3.5 Timestep and Courant number

| | Zhao | ours |
|---|---|---|
| dt | **4e-5 s fixed** (§3.7) | adjustable, `maxDeltaT 1e-5`, `maxCo 0.1` |
| dt achieved at 0.25 mm | 4e-5 | **~5.2e-6** (R0, measured) |
| implied Courant at 0.9 m/s droplet | ~0.14 | 0.1 cap |
| implied Courant at 4-5 m/s gas | **0.5-1.3** | 0.1 cap |

His Courant numbers in the gas exceed unity. Advection schemes there are
strongly diffusive. *(Computed from his stated dt and our measured gas speeds;
he does not report Courant numbers.)* This is the most likely reason his
resolved gas transports less than ours, and it bears directly on §5.

### 3.6 Mesh

| | Zhao | ours |
|---|---|---|
| cells | 256,512 initial (§3.7) | 862,400 (R0) / 88,200 (M2, G1) |
| minimum size | **0.25 mm** | 0.25 mm (R0) / 0.5 mm (M2, G1) |
| domain | half | full |
| air subdomain height | **7 mm** (§3) | 7 mm (R0) / 5 mm (M2, G1) |

R0 matches his resolution. M2 and G1 are a factor of two coarser, which is why
`w fuse` from them carries a cell-scale uncertainty of about +/-0.3 mm per side.

## 4. The droplet model: the largest structural difference

**Zhao resolves droplet detachment. We prescribe it.**

He solves the full force balance at the wire tip -- gravity, buoyancy, surface
tension, Marangoni, electromagnetic (Eqs. 16-19) and plasma drag (Eq. 20) --
and lets necking and detachment emerge (§3.3, §4.1). Detachment occurs at
t = 24 ms and the free droplet reaches 0.9 m/s internally.

We inject a sphere of radius 1.105 mm at 1506 K, placed
`droplet_standoff_m = 0.5 mm` clear of existing metal, with `0.9 m/s`
imposed -- his *result* used as our *input*.

Consequences:

- **Not a defect for pool physics.** Everything downstream of impact -- mass,
  momentum, enthalpy -- is the same, and his 0.9 m/s already embeds the plasma
  drag and electromagnetic pinch we don't model.
- **But we cannot reproduce his droplet chapter**, and metal-transfer-mode
  questions (globular vs spray) are outside our model by construction.
- His plasma drag uses **rho_f = 6e-6 kg/m^3** and v_f = 100 m/s (§3.3). That
  density is 270,000x below argon at 300 K. It is stated for the plasma, not
  the VOF gas phase -- but it is the only gas density printed anywhere in the
  paper. See §5.

Note also **errata 4**: Eq. 30 with f_d = 30 Hz gives r = 1.105 mm (his printed
value), but his *simulated* droplet is 1.114 mm, implying f_d = 29.30 Hz,
-2.35%. Mass is conserved either way.

## 5. What we cannot match, because he did not publish it

### 5.1 The Darcy mushy-zone constant

Zhao gives only the range: *"values between 1e4 and 1e7 are commonly
recommended"* (§3.2, quoting the Fluent manual). **Three orders of magnitude,
no value.**

M3 measured what that costs. Dropping `darcy_A` from 1e6 to 1e5 -- one
coefficient, `diff -r` returns a single hunk:

```
w toe    5.942 -> 11.953 mm    +101%
h crown  2.770 ->  1.265 mm    -54%
w fuse   3.560 ->  3.650 mm     within scatter -- unchanged
```

**Bead geometry is controlled by a parameter he never states, while the fused
zone is not.** Matching his bead width therefore means tuning an unpublished
number until an artefact agrees. Matching `w fuse` does not.

### 5.2 The gas phase properties -- and this one decides the pool

Table 4 is titled "Thermo-physical properties of 5A05 and other parameters".
It lists metal density, solidus, liquidus, specific heat, conductivity,
viscosity, surface tension, latent heat, expansion coefficient, h, and
constants. **There is no argon density, specific heat or thermal conductivity
anywhere in the paper.**

His Eq. 4 is solved over the whole domain and Eq. 6 blends properties by VOF,
so his gas *should* advect enthalpy. How much depends entirely on
rho_gas . cp_gas, which is unpublished. We use argon at **1.6228 kg/m^3**,
correct at 300 K.

**This is a new erratum and arguably the most consequential one: the paper
specifies a two-phase energy equation without specifying the second phase.**

### 5.3 Where the gas leaves

| | Zhao | ours |
|---|---|---|
| top, r > R_nozzle | Vz = 0 (Table 3) | Vz = 0 -- matched |
| outlets | AEFB, ADHE, BCGF -- **three** pressure outlets (Table 3) | **two** x-faces only |
| y-faces | one is an outlet *(inferred from Fig. 3 lettering; not verified)* | `farField`, **noSlip walls** |

G1 sealed the ceiling outside the annulus, faithfully. But if one of his three
outlets is a y-face, his trapped gas has an escape ours lacks -- and
`gas_velocity_probe.py` shows G1's fast argon moved from *behind* the torch
(buoyant plume, M2) to *ahead of it hugging the plate* (forced lateral flow).
**Worth confirming against Fig. 3 before attributing anything to physics.**

### 5.4 Not published at all

- linear solver settings, convergence tolerances
- under-relaxation factors
- emissivity used in Eq. 29 (only "eps_b" symbolically)
- initial substrate temperature (300 K implied by Table 3)

## 6. The creep artefact -- shared, and he says so three times

Both models use enthalpy-porosity, which treats solid as a fluid of very high
viscosity. Liquid can therefore creep over substrate that never melted.
`pool_width.py` exists to detect it via `w toe / w fuse`:

```
M2   5.552 / 3.560  =  1.56x
G1   6.782 / 3.238  =  2.09x
M3   11.953 / 3.655 =  3.27x     (darcy_A 1e5)
```

**Zhao documents the same defect in his own results:**

> §4.2: "all the liquid metal has solidified when the time is equal to 1.48s;
> **however, the width of the deposition layer has slightly expanded**, which
> is caused by the inherent defect of enthalpy-porosity method that regards the
> solid as a fluid with weak momentum"

> §4.3: "the width of previous deposited layers **will increase slightly as the
> calculation proceeds and the deposit profile would be flatter**"

> §4.4: "the width difference of the 1st layer of SPFL between the simulated
> and experimental results **becomes larger** than that of SPSL; one of the
> reasonable explanations is that caused by the inherent shortcoming of the
> enthalpy-porosity method"

**And he never validates width.** Table 5 tabulates height and length only;
Conclusion (1) quotes errors "only 4.3% in height and 0.5% in length". Width is
never given an error bar, because width is what the method corrupts.

So the 7.76 mm target is his **simulation's** width, read from Fig. 16, from
the **SPFL** run (not SPSL), and subject to the creep he just described
(**errata 1**).

**Implication: calibrate b and c on `w fuse`, not `w toe`.** `w fuse` measures
melted substrate, is creep-immune by construction, is insensitive to the
unpublished `darcy_A`, and is the quantity the Lorentz force will act within.

## 7. Consolidated errata

1. **The 8.2 mm width is not in the paper.** Printed widths are 7.76 mm (layer 1)
   and 4.52 mm (layer 4), §4.3 from Fig. 16, both simulated, both SPFL. 8.2 is
   our parabolic inversion of Table 5. Quote 7.76; present 8.21 as a check.
2. **Substrate thickness is 4 mm.** §3's "50x15x5" half-domain is typographical;
   §2 and the figures give 4 mm.
3. **Eq. 27 lacks the z-restriction its normalisation requires.** Applied
   literally it delivers 2.eta.U.I.
4. **Droplet frequency inconsistent**: 30 Hz printed, 29.30 Hz implied by the
   simulated radius.
5. **eta = 0.70 applied as pool-only but sourced as total** -> 103-107% of U.I.
6. **Table 5 corroborates the over-delivery**: simulated height 1.55 mm against
   experimental 1.80, a 14% under-prediction. At conserved area a shorter bead
   is a wider one.
7. **NEW -- the gas phase is never specified.** No argon density, cp or k. The
   only gas density printed (6e-6 kg/m^3, §3.3) is for the plasma drag force.
8. **NEW -- width is known-corrupted and never validated.** §4.2, §4.3 and §4.4
   each attribute width error to enthalpy-porosity creep; Table 5 omits width.

## 8. What actually matters, ranked

1. **`darcy_A` unpublished (§5.1)** -- controls bead width 2:1, and bead width is
   the headline geometric comparison. Fatal to width-based calibration.
2. **Gas phase unpublished (§5.2)** -- controls how much heat the argon carries,
   and G1 shows the pool is sensitive to it at the 35% level.
3. **Arc power split (§3.1)** -- resolved in our favour by measurement: his
   literal reading melts his own coupon.
4. **Timestep, 8x (§3.5)** -- his gas Courant exceeds 1; ours is capped at 0.1.
   Different effective gas transport from identical boundary conditions.
5. **Droplet model (§4)** -- large structurally, small in consequence, because
   his result is our input.
6. **Outlet topology (§5.3)** -- possibly explains G1 entirely. Cheap to check.
7. **Mesh (§3.6)** -- eliminated by R0, which matches his 0.25 mm.

Items 1, 2 and 6 are all "the paper does not determine the answer". That is
worth stating plainly in the thesis: **a faithful replication of Zhao 2021 is
not fully specified by Zhao 2021.**
