# Open issues, deviations and deferred fixes

WAAM SPSL replication of Zhao et al., *Welding in the World* 65:1571–1590 (2021).
Status as of the A1′ reference run (88,200-cell mesh, gravity fixed, no arc
pressure, no Lorentz force).

Each entry records **what**, **evidence**, **magnitude**, and **what to do**.
Anything without measured evidence is marked as a hypothesis.

---

# 1. BLOCKING — must be fixed before the Table 5 comparison

## 1.1 Arc pressure not implemented (Zhao Eqs. 25, 26)

**What.** P_arc = μ₀I²/(8πσ_p²)·exp(−r²/2σ_p²), peak **227.8 Pa** at I = 135 A,
σ_p = 2 mm, surface integral 1.82 mN.

**Magnitude.** Third-largest term in the force ledger. For scale, 2 mm of liquid
aluminium weighs 52 Pa and the Laplace pressure of a 3 mm cap at γ = 0.85 N/m is
283 Pa. So arc pressure is ~4× the weight of the pool it presses on.

**Status.** Patch written and tested (`apply_stage2_arcpressure.py`), **not yet
compiled or run**. Adds `updateArcPressure.H`, a `pArc` AUTO_WRITE field, and
the term to `pEqn.H` (the live path — see 3.1).

**Do.** Apply after A1′ finishes. Verify max(pArc) = 227.8 Pa from a written
time directory before trusting it; a force that loads but computes zero is
worse than one that crashes.

## 1.2 Pool Lorentz force not implemented (Zhao Eqs. 22–24)

**What.** Kumar & DebRoy closed form for a Gaussian current density entering the
free surface and spreading to exit through a plate of thickness L. Radial
components inward, axial component downward, scaled by (1 − (z−z₀)/L).

**Magnitude.** Peak ~**435 Pa**, the largest body force in the problem. Opposes
the outward Marangoni surface flow and deepens penetration.

**Validity.** Rm = μ₀σUL ≈ 4×10⁻³ ≪ 1, so the quasi-static prescribed-current
form is valid and no induction equation is needed.

**Note.** L is the substrate thickness, now **4 mm** after the Stage 1 trim.
Using the old 6 mm would make the depth dependence wrong by 50%.

**Do.** Write `updateEMSource.H` in the style of `updateArcSource.H`, gate on
`alpha1`, clip r at half a cell (the closed form has 1/r and 1/r² factors),
add `- emForce*damper` to both `UEqn.H` and `pEqn.H`. Keep
`electromagnetics false`; the MTHD induction path is not what the paper uses.

## 1.3 Combined magnitude of what is missing

| Term | Magnitude | Status |
|---|---|---|
| Gravity, hydrostatic at h = 2 mm | 52 Pa | present |
| Marangoni | 16 Pa | present |
| Arc pressure | 228 Pa | **missing** |
| Pool Lorentz | 435 Pa | **missing** |

**91% of the spreading force is not in the model.** Any bead-height comparison
against Zhao before 1.1 and 1.2 land is precise and wrong in a known direction:
too tall, too narrow.

---

# 2. TO CHECK — measured but not yet explained

## 2.1 Gas velocity is ~5× the buoyancy scale, and it sets the time step

**Evidence.** Back-calculated from Courant, U = Co·Δx/Δt:

| | Co | dt | Implied U |
|---|---|---|---|
| t ≈ 0.04 | 0.048 | 4.0e-5 | 0.60 m/s |
| t = 0.24 | 0.09997 | 1.02e-5 | **4.91 m/s** |

`pool_stats_waam.py` agrees: `U gas` 4.84 m/s against `U metal` 0.46 m/s.
Courant is binding, not `maxDeltaT`.

**Why it looks wrong.** Buoyancy scaling with the argon β = 3.33×10⁻³ from
`transportProperties`:

| Path | Predicted v |
|---|---|
| Rise through 5 mm of argon | 0.63 m/s |
| Circulate 15 mm horizontally | 1.08 m/s |
| √(gβΔT·L) | 0.44 m/s |

Measured 4.9 m/s is ~5× the largest of these.

**Two hypotheses, neither tested:**

1. **Boussinesq outside its range in the gas.** The approximation needs
   βΔT ≪ 1. `T gas` reached 1488 K at t = 0.04 and 1264 K at 0.24, so
   ΔT ≈ 1000–1200 K and **βΔT ≈ 3–4**. The linearisation *overestimates*
   buoyancy, because real gas density falls as 1/T, not linearly.
2. **Parasitic VOF currents.** Known artifact, scales with γ/μ, always worst in
   the light low-viscosity phase. Gas at 10× the melt velocity is the classic
   signature.

**Cheap test.** Find *where* the fast gas is. A plume above the pool means
buoyancy; hugging the interface means parasitic currents. Adding the cell
location of the `U gas` maximum and its local α to `pool_stats_waam.py` would
settle it in one line.

**Impact.** Cost only. It sets dt, and dt drives every wall-clock estimate. The
melt side is unaffected (`U metal` 0.46 m/s, ample Courant margin). Belongs in
the limitations section either way.

## 2.2 Does the 4 mm plate stay solid?

**Evidence.** `bottom` reached **685 K at t = 0.24** against an 815 K solidus,
still climbing. The plate is now 6,000 mm³ (Zhao's coupon exactly), 3× lighter
than the old 18,000 mm³, and Eq. 29 removes only ~2% of the energy budget.

**Check on the production run.** Does `bottom` plateau once the moving source
reaches quasi-steady state, or keep climbing linearly? And does pool `depth`
stay well under 4 mm (gate text still says 6 mm — see 4.2).

**Watch the ends specifically**, not just the global maximum — see 3.2.

## 2.3 Crown symmetry — the still-open question

**History.**

| Run | h max / h crown | Decomposition |
|---|---|---|
| Old mesh, scotch | 1.062 | y-asymmetric partitions |
| Old mesh, simple (4 1 1) | 1.031 | y-symmetric |
| New mesh, transient proxy | ~1.020 | y-symmetric |
| **A1′** | pending | y-symmetric |

Under z-gravity the problem is exactly invariant under y → −y (symmetric mesh
grading, arc at y₀ = 0, Goldak even in y, droplets at y = 0), so **1.000 is
forced**. Anything above ~1% is a third source not yet found. Candidates:
droplet placement relative to cell centres (y = 0 falls between cells at
±0.25 mm), or an impact instability at 4.4 cells per droplet.

**Note.** `bead_profile.py` only prints the transverse section when a developed
window exists (`if win:`, line 448), which needs t > 0.428 s. A1′ at 0.46 s is
the first run that produces it.

---

# 3. KNOWN DEVIATIONS FROM ZHAO — record in the thesis

## 3.1 `momentumPredictor no`, so `UEqn.H` body-force block never executes

**Evidence.** `fvSolution` line 128. The `if (pimple.momentumPredictor())` block
in `UEqn.H` is dead; surface tension, `pVap` and buoyancy take effect through
`phig` in `pEqn.H` lines 31–33.

**Consequence.** Any new body force must go in **`pEqn.H`**. Patching only
`UEqn.H` compiles, runs, and does nothing. Both are kept in step so the code
stays correct if the predictor is ever switched on.

## 3.2 Plate end faces are on `outlet`, not Eq. 29

**Evidence.** `outlet` is the x = −15 and x = +35 faces, each spanning the full
z = −4 to +5, so 4 mm of metal plus 5 mm of argon on one patch with one
condition. Zhao's pressure outlets are in the **air subdomain only**; every
metal face in his Table 3 gets Eq. 29.

**What `inletOutlet` does there.** `valueFraction = 1 − pos0(φ)`. In solidified
metal the Darcy sink drives U → 0, so φ ≈ 0 and `pos0(0) = 1` gives
zeroGradient, i.e. **effectively adiabatic**. But the sign of a near-zero flux
is not robust, so individual faces could flip to `fixedValue 300`, a strong
local sink.

**Coverage.**

| Plate face | Patch | T condition | Area |
|---|---|---|---|
| bottom z = −4 | `substrateBottom` wall | Eq. 29 | 1500 mm² |
| long sides y = ±15 | `farField` wall | Eq. 29 | 400 mm² |
| end faces x = −15, +35 | `outlet` patch | `inletOutlet` 300 | 240 mm² |

**89% of the wetted plate boundary has Eq. 29.** The missing 11% is worth
**0.18%** of the energy budget (3.8 W at ΔT = 200 K).

**But.** The end faces sit only **5 mm** from the bead ends, and √(αt) = 5.9 mm
over the solidification window, so the thermal wave does reach them. Being
effectively adiabatic biases arc-on and arc-off hot — exactly where Zhao's
Fig. 8 shows the height dropping and where `deposit length` is measured.

**Fix if needed.** `topoSet` + `createPatch` to split `outlet` into a metal band
(z < 0) and a gas band, Eq. 29 on the metal part. Mesh-topology job, not a
dictionary edit. Low priority; monitor `bottom` near the ends first.

## 3.3 Droplet detachment is prescribed, not resolved

Zhao's headline novelty is resolving droplet formation, necking and detachment
from the electromagnetic pinch, surface tension and plasma drag. This model
injects pre-formed spheres at r = 1.105 mm, 30 Hz, v = 0.9 m/s, T = 1506 K —
all taken from Zhao's own results.

**Consequence.** Figs. 4, 5, 6 and the 0.88% droplet-radius validation cannot be
reproduced. Everything downstream of impact can. State this explicitly.

## 3.4 Zhao's plasma drag coefficient looks wrong

Eq. 20 states ρ_f ≈ 6×10⁻⁶ kg/m³ "argon density". Argon at 300 K is 1.62 kg/m³;
at 10,000 K arc temperature ~0.05 kg/m³. 6×10⁻⁶ is five orders below either —
probably a typo for 6×10⁻². Only matters if resolved droplet detachment is ever
attempted; treat as a free parameter and say so.

## 3.5 Surface tension sourced from pure aluminium

γ = 0.85 N/m and ∂γ/∂T = −1.55×10⁻⁴ N/m/K both come from Mills' **pure
aluminium** chapter. Mg is strongly surface-active in Al, so the true Al-5Mg γ
is lower. Marangoni drives the paper's headline result, making this the largest
physical uncertainty in the whole replication.

**Do.** Sensitivity runs at γ = 0.70 N/m and ∂γ/∂T ±30%.

## 3.6 No shielding gas inlet (Zhao Eq. 28)

Zhao imposes a 20 L/min annular argon jet, ~1–2 m/s downward, between wire and
nozzle. Not implemented. The top boundary is passive two-way venting instead.

**Consequence.** No forced convective cooling of the free surface, no jet
momentum, nothing sweeping the plume. Related to 2.1.

## 3.7 Vaporisation disabled

`Tvap = 1e6` as a switch. Matches Zhao's assumption 3, but removes the negative
feedback that would cap temperature. Monitor max(T metal); above ~1600 K nothing
in the model stops it.

## 3.8 `maxCo 0.1` is 3× more conservative than Zhao

Zhao ran 0.25 mm cells at a **fixed** dt = 4×10⁻⁵ s, which at ~2 m/s gas is
Co ≈ 0.3. The Brackbill capillary limit is 2.5×10⁻⁴ s, 25× above any step used
here, so surface tension is not the constraint, and isoAdvector is geometric and
tolerates higher CFL than the MULES compression 0.1 was inherited for.

**Do.** Back-to-back 0.46 s runs at maxCo 0.1 and 0.3; confirm h(x), A and pool
depth agree. Doubles as the temporal-convergence check the thesis needs, and
licenses ~3× speedup for the production run.

## 3.9 Mesh is 4.4 cells per droplet diameter

0.5 mm core against a 2.21 mm droplet. Below the 5–6 usually taken as the VOF
curvature floor, and `updateDropletSource.H` warns below 4.

**Acceptable** for the ablation study, which compares runs on identical grids so
curvature error largely cancels. **Not acceptable** for the Table 5 validation,
which needs 0.25 mm (684k cells full domain, 8.8 cells per droplet).

## 3.10 Other property and modelling notes

- ρ held constant at 2650 kg/m³; liquid Al-5Mg near the liquidus is ~2350, so
  the melt is ~13% heavy. Kept for fidelity with Zhao.
- ε_S = 0.3, ε_L = 0.15 are estimates for oxidised Al; radiation overtakes
  convection near 1300 K. One sensitivity run.
- h = 80 W/m²K is high for natural convection off a small plate (5–15 typical).
  It is a lumped coefficient absorbing the fixture conduction path. Say so
  rather than presenting it as a convection coefficient.
- Solver is OpenFOAM v2506 + laserbeamFoam against Zhao's Fluent + UDFs;
  different VOF scheme (isoAdvector vs geometric reconstruction).

---

# 4. TOOLING — stale text and missing scripts

## 4.1 `pool_stats_waam.py`: `far met` gate text is now wrong

It reads as a domain-size gate: "the domain is still big enough". After the trim
to Zhao's 50 × 30 mm coupon, `farField` is a **physical plate edge**, not a
truncation, and it is *expected* to warm.

**Replace with:** "plate edge temperature. This is a physical edge (Zhao's
coupon is 50 × 30 × 4 mm), not a truncation, so it will rise. Gate is that it
stays well below the 815 K solidus."

## 4.2 `pool_stats_waam.py`: two more stale gates

- `depth < 6 mm` → **`depth < 4 mm`** (plate is now 4 mm)
- `bottom near 300 K` → no longer true or expected; Eq. 29 is active and the
  plate is 3× lighter. Gate is < 815 K solidus.

## 4.3 Droplet-clearance note assumes 6 mm of argon

The `Z0, Z1` comment in `make_waam_case.py` computes clearance for 6 mm of
argon. It is now 5 mm. Recompute or the note misleads.

## 4.4 `bead_profile.py`: report `y` at `h max`

`h max > h crown` has two possible causes: an off-centre crown, or a symmetric
profile with a central dip and two shoulders. The second becomes physically real
once arc pressure creates a depression. Reporting the y-location of the maximum
makes the symmetry test unambiguous. One-line change, do it before A2.

## 4.5 Missing post-processing scripts

| Script | For |
|---|---|
| `height_profile.py` | Fig. 8, h(x) along the track plus deposit length from h > 0.1 mm |
| `cross_section.py` | Fig. 18a, α = 0.5 contour as a polyline via `sampledInterface` |
| `flow_field.py` | Figs. 10, 11, with scalar diagnostics (surface radial velocity sign, bottom circulation sign) rather than eyeballed streamlines |
| `energy_budget.py` | arc in, latent, sensible, radiative, convective, closure error |
| `compare_runs.py` extension | the ablation table and the SPSL row of Table 5 |

## 4.6 `Allrun` does not fail on solver crash

**Evidence.** `MPI_ABORT` fired at rank 0 and `Allrun` still ran
`reconstructPar` and both `postProcess` stages. A 30-second startup crash
produced output that looked like a completed run.

**Fix.**

```sh
runParallel laserbeamFoam
grep -q "^End" log.laserbeamFoam || { echo "ERROR: solver did not reach End -- see log.laserbeamFoam"; exit 1; }
```

Same family as the existing `Allrun -f` guard: fail loudly rather than produce
plausible-looking empty output.

## 4.7 `E kept` in `pool_stats_waam.py` reads 130–139%

Almost certainly the denominator. Arc energy to t = 0.9161 is 1385 J, but
droplets carry their own enthalpy in: 0.404 g from 300 K to 1506 K is ~536 J
sensible plus ~145 J latent ≈ 681 J. (1385 + 681)/1385 ≈ 149%, same order as
reported. **Check how `E kept` is defined** before treating it as a
conservation failure.

---

# 5. LESSONS THAT COST TIME — process notes

## 5.1 Anything hand-set with `foamDictionary` is lost on `--force`

Bitten three times: `constant/g` (copied from a 2-D y-up donor into five
z-up cases), `libs ("libthermoTools.so")`, and `arc_x0_m` (reset from −0.010 to
0, putting arc-off 5 mm past the plate edge). All three now live in
`make_waam_case.py`.

**Rule.** If it must survive a rebuild, it goes in the generator. `endTime` and
`maxDeltaT` are the remaining exceptions and must be re-set after every
`--force`.

## 5.2 A value stated in one place and derived in another will drift

The build summary hardcoded "4 subdomains, simple (4 1 1)" while the dictionary
was written from `DECOMP`. It reported the wrong thing for a full rebuild cycle.
Now derived.

## 5.3 Decomposition is not a free choice

`scotch` minimises cut edges but its subdomains are not mirror images, so it
seeds a +y/−y difference every timestep. It produced an **8.9% false asymmetry
signal** that was nearly attributed to physics. `simple (6 1 1)` slices x only,
so every processor holds complete y-columns and mirror symmetry is exact.

The generator's own comment already said this; it was overridden on my advice.

**Rule.** Freeze `DECOMP` across any set of runs being compared, and never use a
y-cutting decomposition while y-symmetry is a diagnostic.

## 5.4 Check the registered type table before assuming a BC exists

`externalWallHeatFluxTemperature` is not in laserbeamFoam's default patchField
table. It lives in `src/thermoTools/` and needed
`libs ("libthermoTools.so")`. Two wrong causes were proposed from memory before
the log was read. The log prints the full valid-type list on failure — that is
definitive evidence already on disk.

## 5.5 Bugs hide until the run is long enough

The fine region stopped at x = +20 mm while the arc reaches +30 mm, putting the
arc-off collapse in 1.0–1.6 mm cells — exactly where `deposit length` is
measured. Invisible because every run so far stopped at t ≤ 0.5 s, where the arc
has only reached x = 0.

---

# 6. FROZEN CONFIGURATION — ablation set A1′ to A4

Every case must match exactly. `--force` resets `endTime` and `maxDeltaT`.

```
mesh        98 x 50 x 18 = 88,200 cells
            x[-15,+35] y[-15,+15] z[-4,+5] mm  = Zhao's 50 x 30 x 4 mm coupon
            fine region x[-12,+32] y[-6,+6], CELL 0.5 mm, EXPAND 2.0
arc_x0_m    -0.010        bead -10 -> +30, 5 mm lead-in and run-out
DECOMP      (6, 1, 1)
endTime     0.46          first developed window is t > 0.428
deltaT      4e-05   maxDeltaT 4e-05   maxCo 0.1   maxAlphaCo 0.1
libs        ("libthermoTools.so")
gravity     (0 0 -9.81)
top         U pressureInletOutletVelocity; T/alpha inletOutlet
walls       externalWallHeatFluxTemperature, h = 80, Ta = 300
```

| Case | Arc pressure | Lorentz | Status |
|---|---|---|---|
| A1′ | off | off | running |
| A2 | **on** | off | patch ready, not compiled |
| A3 | off | **on** | not written |
| A4 | **on** | **on** | — |

**Record per case at t = 0.44 and 0.46:** A, h crown, h max, w toe, w half,
ripple, pool depth, U metal, U gas, bottom, far met.

---

# 7. VALIDATION TARGETS — SPSL only

| Quantity | Experiment | Zhao simulation |
|---|---|---|
| Height at mid-bead | 1.80 mm | 1.55 mm |
| Deposit length | 43.28 mm | 42.47 mm |
| Cross-section A | — | 8.48 mm² (mass-conserved, Q_v/v) |
| Width (toe) | — | ~8.2 mm |

Measure at **x = +10 mm** in this coordinate system (mid-bead), not x = +25 mm.
Zhao's X = 25 mm is mid-bead in a 0–50 mm plate; here the bead runs −10 to +30.

Qualitative targets: Fig. 8 height profile, Fig. 9 cross-section evolution,
Figs. 10–11 flow fields, Fig. 18a cross-section overlay.
