# Open issues, deviations and deferred fixes

WAAM SPSL replication of Zhao et al., *Welding in the World* 65:1571–1590 (2021).

Revision 2. Supersedes the version written at the A1′ reference run. Adds the
arc-efficiency finding, measured targets read from Fig. 9, and retractions.

Each entry records **what**, **evidence**, **magnitude**, **what to do**.
Anything without measured evidence is marked as a hypothesis. Anything I got
wrong earlier is marked RETRACTED, with what replaced it.

---

# 0. MEASURED TARGETS FROM THE PAPER

Read off Fig. 9 panels b1–b5 (PDF page 12), the cross-section at X = 20 mm with
a 300–906 K colour scale and z spanning the 4 mm plate. These are pixel-colour
estimates from a printed figure, so **±0.3 mm and ±70 K (one colour band)**.

| Quantity | Zhao | Source |
|---|---|---|
| Pool depth `d sol` at pool maximum | **1.0–1.3 mm** | b3, b4 (t = 0.84, 1.04 s), red ≥839 K reaching z ≈ 3.0 mm in a plate topped at z = 4 mm |
| Bottom face temperature at y = 0 | **435–570 K** | b2–b4, cyan-to-green at z = 0 |
| Far substrate, \|y\| > 7 mm, t = 1.48 s | **435–500 K** | b5 |
| Far substrate, t = 0.72 s | ~300 K | b1, dark blue |

From Fig. 18a (PDF page 18), the SPSL cross-section at X = 25 mm:

| Quantity | Zhao |
|---|---|
| Bead height above the z = 4 mm substrate | **1.5–1.6 mm** (Table 5 says 1.55) |
| Bead full width at the base | **~8 mm** |

From Table 5:

| Quantity | Experiment | Zhao simulation |
|---|---|---|
| Height at mid-bead | 1.80 mm | 1.55 mm |
| Deposit length | 43.28 mm | 42.47 mm |
| Cross-section A | — | 8.48 mm² (mass-conserved, Q_v/v) |

**Fig. 18 cannot show penetration.** The red half is metal (α > 0.5), and
remelted substrate is metal before and after fusion, so the fusion line is
invisible in the simulated half by construction. It shows the free-surface
profile only.

---

# 1. RETRACTIONS

Recording these because each nearly drove a parameter change on no evidence.

## 1.1 RETRACTED: "Zhao's penetration is ~1 mm, from Figs. 11 and 18"

Invented. Zhao never reports a penetration depth in the text or Table 5, and
Fig. 18 physically cannot show one. The figure was cited repeatedly as if it
were a measured target and it drove an entire "2.8× too deep" argument.

**Replaced by** §0: `d sol` 1.0–1.3 mm read from Fig. 9 b3/b4. The direction
turned out right — yours is 3.25 mm, about 2.7× deeper — but that was luck,
not method.

## 1.2 RETRACTED: "η_total = 0.85 is the most defensible next value"

0.85 is the top of Zhu's stated 0.60–0.85 band, but Zhao states 0.70 and the
goal is replication. 0.85 was reached for because it would make the numbers
work, which is fitting, not replicating.

**Replaced by** keeping η = 0.70 and testing whether the mesh, not the power,
was the reason 790 W produced no pool. See §3.

## 1.3 RETRACTED: "Goldak c = 4 mm is the suspect"

The power-fraction arithmetic is correct — at c = 4 mm, 21% of the arc source
lands below 2 mm depth, against 1.4% at c = 2 mm — but the hypothesis does not
survive Fig. 9: **Zhao uses c = 4 mm as well**, with the same plate and the same
z₀. A 2.7× difference in penetration cannot be explained by a parameter both
models share.

Kept here because the arithmetic stays useful if c is ever revisited for a
different reason.

## 1.4 RETRACTED: "externalWallHeatFluxTemperature fails because of kappa lookup ordering"

Asserted from reasoning, then a second wrong cause was asserted from memory,
before the log was read. The log prints the full valid-type table on failure —
definitive evidence already on disk. Actual cause: the type is not registered
in laserbeamFoam; it ships in `libthermoTools.so`, which loads fine.

---

# 2. RESOLVED

## 2.1 RESOLVED: gravity on the wrong axis

`constant/g` was copied from `Plate2D-mhd`, a 2-D x–y donor with y vertical,
into five z-up WAAM cases. Gravity acted horizontally.

**Fixed**, and the generator now writes it rather than copying. Confirmed by
the A1′ transverse section: h(−y) vs h(+y) summed to −1.7%, centroid
**−16.9 µm = 3.4% of one 500 µm cell**. Symmetry is forced to exactly 1.000 by
the problem's invariance under y → −y, so that residual is round-off.

History of h_max/h_crown: 1.062 (scotch, g=−y) → 1.031 (simple 4 1 1, g=−y)
→ **1.003** (simple 6 1 1, g=−z, solidified bead).

## 2.2 RESOLVED: scotch decomposition produced a false asymmetry

`scotch` minimises cut edges but its subdomains are not mirror images, so it
seeds a +y/−y difference every timestep. It produced an 8.9% asymmetry signal
that was nearly attributed to physics. `simple (6 1 1)` slices x only, so every
processor holds complete y-columns and mirror symmetry is exact.

The generator's own comment already argued this; it was overridden on my advice
and then restored.

## 2.3 RESOLVED: `E kept` reading 119–164%

The denominator counted arc energy only, while droplets enter as mass at
1506 K carrying their own enthalpy. Measured values tracked (arc+722)/arc
exactly: 1.478 predicted vs 1.19–1.33 reported at 1512 W, 1.914 vs 1.50–1.64 at
790 W. Not a conservation failure.

**Fixed** in `pool_stats_waam.py`: denominator is now arc + droplet enthalpy.

## 2.4 RESOLVED: the plume was trapped under a lid

`inlet` had `U = fixedValue (0 0 0)`, which is hydrodynamically a lid: zero
velocity is zero through-flow. The buoyant argon rose, could not leave, spread
along the roof and reached `farField`, a no-slip wall.

**Fixed**: `pressureInletOutletVelocity` on U, `inletOutlet` on T and
alpha.metal. Effect measured via Courant: implied gas velocity fell from
~4.05 m/s to ~0.60 m/s, and dt rose from 1.0e-5 to 4.0e-5.

---

# 3. THE OPEN QUESTION: heat input and mesh resolution

## 3.1 The arc-efficiency finding

Zhao Eq. 27 (PDF p. 7) calls η "the arc efficiency for heating molten pool,
**excluding heating droplets**, which is described in detail in Ref. [32]", and
Table 4 sets it to 0.7 citing Wu [29].

Ref. [32] is Zhu, Cheon, Tang, Na & Cui (2018), *Int J Heat Mass Transf*
126:1206–1221. Its Eqs. 7–9:

```
q_arc = η₁·U·I / (2πσ_xσ_y) · exp(...)                            (7)
η₁ = η − η₂                                                       (8)
η₂ = 4πr_d³ρ[C_s(T_s−T_0) + C_l(T_d−T_l) + h_sl]·f_d / (3UI)      (9)
```

verbatim: *"η₁ is the arc efficiency for heating the molten pool. The TOTAL arc
efficiency η ranges from 0.60 to 0.85 for GMAW process"* and *"η₂ represents
the percentage of arc heat absorbed by the droplets"*.

So η is the **total and includes the droplets**; η₁ is what remains for the arc
source. Zhao's 0.70 is a total efficiency placed in the η₁ slot.

η₂ with Zhao's numbers: mass rate 0.450 g/s, enthalpy 300 → 1506 K =
541 (solid) + 358 (latent) + 708 (liquid superheat) = 1607 kJ/kg, giving
**722 W = η₂ 0.334**.

| Reading | Arc source | Total to workpiece |
|---|---|---|
| Zhao as written (η in the η₁ slot) | 1512 W | **2234 W = 103% of U·I** |
| Zhu Eq. 8 applied (η₁ = 0.70 − 0.334) | **790 W** | 1512 W = 70% of U·I |

103% is impossible. This is an internal inconsistency in the paper and it
applies to Zhao's own model, not only to this replication.

## 3.2 What each power actually does, measured

Both on the frozen 88,200-cell configuration, 0.5 mm core cells.

| | 1512 W | 790 W | **Zhao (Fig. 9)** |
|---|---|---|---|
| Pool continuity | continuous from t = 0.08 | **0 cells at t = 0.06 AND 0.12** | continuous |
| `d sol` | 3.25 mm at t = 0.66 | **negative** (melt above z = 0, substrate never fuses) | **1.0–1.3 mm** |
| `bottom` | **810.1 K** at t = 0.66, 4.9 K under solidus | 465.8 K at t = 0.12 | **435–570 K** |
| T metal at t = 0.12 | 896.2 K | **790.7 K, below the 815 K solidus** | — |

So 1512 W is too hot at the bottom face by ~300 K, and 790 W does not sustain a
weld pool at all. Neither matches Fig. 9.

## 3.3 The bulk heating arithmetic, which is not in dispute

2234 W into 15.9 g of aluminium at c_p ≈ 1050 gives **134 K/s**. Measured
`bottom` rise from t = 0.46 to 0.66 is **124 K/s**. The local gradient has
saturated (the arc has run 5.3 thermal lengths, α/v = 2.5 mm), so what remains
is the whole plate warming together, and that does not saturate.

Fig. 9 b5 is consistent with this for the **bulk**: far substrate at ~435–500 K
at t = 1.48 s, against 300 + 134×1.48 = 498 K predicted.

**So Zhao's bulk heating matches. His local excess at the bottom face does
not.** Zhao sits ~100 K above bulk; this model sits ~420 K above it.

Loss paths are all too small to absorb the difference:

| Path | Share of input |
|---|---|
| Eq. 29 walls, h = 80 | 2.3% |
| Free-surface radiation | 0.1% |
| Shielding gas jet (not modelled) | 0.9% |

## 3.4 CURRENT HYPOTHESIS — untested

**790 W is correct, and the 0.5 mm mesh is too coarse to sustain a pool at that
power.**

At 0.5 mm a 2.21 mm droplet spans 4.4 cells, below the 5–6 usually taken as the
VOF floor and barely above the injector's own warning threshold. Its enthalpy is
averaged over a large volume on arrival — numerical diffusion. At Zhao's 0.25 mm
it spans 8.8 cells and stays concentrated, giving a higher local temperature for
the same energy.

This is self-consistent with everything measured: it explains why 790 W gave no
pool here but works for Zhao, and it explains the excessive bottom temperature
at 1512 W as the consequence of over-supplying heat to compensate.

It also predicts the failure mode: if resolution is the issue, refining should
*raise* local pool temperature at fixed power.

**Test (running):** 790 W at 0.35 mm core, 6.3 cells per droplet.

| | arc | cell | cells/droplet | Expected if hypothesis holds |
|---|---|---|---|---|
| have | 790 W | 0.5 mm | 4.4 | no pool ✓ |
| **test** | 790 W | 0.35 mm | 6.3 | pool continuous, `d sol` 1.0–1.3 mm, `bottom` < 700 K |
| have | 1512 W | 0.5 mm | 4.4 | pool ok, `d sol` 3.25 mm, `bottom` 810 K ✗ |

If it still gives no pool, Zhao's arc really is 1512 W and the difference lies
elsewhere.

---

# 4. BLOCKING — required before the Table 5 comparison

## 4.1 Arc pressure not implemented (Zhao Eqs. 25, 26)

P_arc = μ₀I²/(8πσ_p²)·exp(−r²/2σ_p²), peak **227.8 Pa** at I = 135 A,
σ_p = 2 mm, surface integral 1.82 mN.

Patch written and tested (`apply_stage2_arcpressure.py`), **not compiled or
run**. Adds `updateArcPressure.H`, a `pArc` AUTO_WRITE field, and the term to
`pEqn.H` — the live path, since `momentumPredictor no` makes the equivalent
block in `UEqn.H` dead code. Includes a z gate at `arc_z0 + 3 mm` so the term
does not act on free-flying droplets, which already carry a prescribed velocity
(228 Pa is ~6× a droplet's own weight).

## 4.2 Pool Lorentz force not implemented (Zhao Eqs. 22–24)

Kumar & DebRoy closed form, valid at Rm ≈ 4×10⁻³ ≪ 1. Radial inward, axial
downward, scaled by (1 − (z−z₀)/L) with **L = 4 mm** after the Stage 1 trim.
Peak ~**435 Pa**, the largest body force in the problem.

Note it **increases** penetration, so it must not go in until §3 is settled.

## 4.3 What is missing, quantified

Laplace pressure holding the current crown up: γκ = 2γ/R with R = 2.40 mm gives
**709 Pa**.

| Term | Pressure | Share of Laplace |
|---|---|---|
| Hydrostatic ρgh | 72.8 Pa | 10.3% |
| Marangoni | 15.5 Pa | 2.2% |
| **Arc pressure (missing)** | 227.8 Pa | 32.1% |
| **Pool Lorentz (missing)** | 435.0 Pa | 61.4% |

**91% of the spreading force is absent.** No bead-height comparison against
Zhao means anything until both land.

Why gravity was never the lever: the bead is capillary-dominated. Capillary
length √(γ/ρg) = 5.72 mm against a 2.4 mm half-width, so Bo = 0.176 and surface
tension beats gravity 5.7×. A/ℓ_c² = 0.26, i.e. a fixed-area ridge settles as a
circular arc with gravity as a perturbation. The Young–Laplace shape
perturbation from gravity is Δκ/κ = ρgh/(γκ) = **10.3%**, which matches the
measured h change of +9.8% on switching −y → −z.

---

# 5. KNOWN DEVIATIONS FROM ZHAO — record in the thesis

## 5.1 `momentumPredictor no`, so the `UEqn.H` body-force block never executes

`fvSolution` line 128. Surface tension, `pVap` and buoyancy take effect through
`phig` in `pEqn.H` lines 31–33. **Any new body force must go in `pEqn.H`.**
Patching only `UEqn.H` compiles, runs, and does nothing.

## 5.2 Plate end faces are on `outlet`, not Eq. 29

`outlet` is x = −15 and x = +35, each spanning z = −4 to +5, so 4 mm of metal
plus 5 mm of argon on one patch. Zhao's pressure outlets are in the air
subdomain only; every metal face in his Table 3 gets Eq. 29.

`inletOutlet` gives `valueFraction = 1 − pos0(φ)`, and in solidified metal
φ ≈ 0 so `pos0(0) = 1` yields zeroGradient — effectively adiabatic. But the
sign of a near-zero flux is not robust.

Coverage: `substrateBottom` 1500 mm² + `farField` 400 mm² have Eq. 29;
`outlet` 240 mm² does not. **89% covered**, missing 11% worth **0.18%** of the
energy budget. The end faces sit only 5 mm from the bead ends and
√(αt) = 5.9 mm, so the thermal wave does reach them.

Fix if needed: `topoSet` + `createPatch` to split `outlet` into metal and gas
bands. Mesh-topology job, low priority.

## 5.3 Droplet detachment is prescribed, not resolved

Zhao's headline novelty is resolving droplet formation, necking and detachment
from the electromagnetic pinch, surface tension and plasma drag. This model
injects pre-formed spheres at r = 1.105 mm, 30 Hz, v = 0.9 m/s, T = 1506 K —
all from Zhao's own results. Figs. 4, 5, 6 and the 0.88% radius validation
cannot be reproduced; everything downstream of impact can.

## 5.4 Zhao's plasma drag density looks wrong

Eq. 20 states ρ_f ≈ 6×10⁻⁶ kg/m³ "argon density". Argon is 1.62 kg/m³ at 300 K
and ~0.05 kg/m³ at arc temperature. Probably a typo for 6×10⁻². Only matters if
resolved detachment is attempted.

## 5.5 Surface tension sourced from pure aluminium

γ = 0.85 N/m and ∂γ/∂T = −1.55×10⁻⁴ N/m/K both come from Mills' **pure
aluminium** chapter. Mg is strongly surface-active in Al, so the true Al-5Mg γ
is lower. Marangoni drives the paper's headline result, making this the largest
physical uncertainty in the replication. Sensitivity runs: γ = 0.70 N/m,
∂γ/∂T ±30%.

## 5.6 No shielding gas inlet (Zhao Eq. 28)

Zhao imposes a 20 L/min annular argon jet, ~1–2 m/s downward. Not implemented;
the top boundary is passive two-way venting. No forced convective cooling of
the free surface, no jet momentum.

## 5.7 Vaporisation disabled

`Tvap = 1e6` as a switch. Matches Zhao's assumption 3, but removes the feedback
that caps temperature. Monitor max(T metal); above ~1600 K nothing stops it.

## 5.8 `maxCo 0.1` is 3× more conservative than Zhao

Zhao ran 0.25 mm cells at a fixed dt = 4×10⁻⁵ s, i.e. Co ≈ 0.3 at ~2 m/s. The
Brackbill capillary limit is 2.5×10⁻⁴ s, 25× above any step used here, and
isoAdvector tolerates higher CFL than the MULES compression 0.1 was inherited
for. Test 0.1 vs 0.3 back to back; doubles as the temporal-convergence check.

## 5.9 Other property notes

- ρ constant at 2650 kg/m³; liquid Al-5Mg near the liquidus is ~2350, so the
  melt is ~13% heavy. Kept for fidelity with Zhao.
- ε_S = 0.3, ε_L = 0.15 are estimates for oxidised Al; radiation overtakes
  convection near 1300 K.
- h = 80 W/m²K is high for natural convection off a small plate (5–15 typical).
  It is a lumped coefficient absorbing the fixture conduction path. Say so.
- C_s = 1050, C_l = 1180 J/kg·K in Zhu Eq. 9 are integral averages used **only**
  for η₂. They do not feed the solver, which uses c_p(T) from
  `transportProperties`. ±10% on C_s moves η₂ by ±0.016, ~35 W on the source.
- Solver is OpenFOAM v2506 + laserbeamFoam against Zhao's Fluent + UDFs;
  different VOF scheme (isoAdvector vs geometric reconstruction).

---

# 6. TOOLING

## 6.1 DONE: `d liq` reported beside `d sol`

`depth` was solidus-based, i.e. measured to the bottom of the mushy zone, while
Zhao's figures colour to 906 K. With a 91 K mushy interval the two differ
substantially. `pool_stats_waam.py` now reports both.

## 6.2 DONE: three stale gate texts corrected

`depth < 6 mm` → `d sol < 4 mm`; `bottom near 300 K` → `bottom < 815 K`, with
the note that it will climb; `far met = 300.0` → no longer a domain-size gate
since `farField` is now a physical coupon edge.

## 6.3 `bead_profile.py`: report `y` at `h max`

`h max > h crown` has two causes: an off-centre crown, or a symmetric profile
with a central dip and two shoulders. The second becomes real once arc pressure
creates a depression. Reporting the y-location makes the test unambiguous.

## 6.4 `Allrun` does not fail on solver crash

`MPI_ABORT` fired at rank 0 and `Allrun` still ran `reconstructPar` and both
`postProcess` stages. A 30-second startup crash looked like a completed run.

```sh
runParallel laserbeamFoam
grep -q "^End" log.laserbeamFoam || { echo "ERROR: solver did not reach End"; exit 1; }
```

## 6.5 `Allrun` should write cell geometry right after `decomposePar`

`pool_stats_waam.py` has hit the missing-geometry error three times mid-run.
The hazard it was deferred for — decomposing a field named `C`, colliding with
OpenFOAM's reserved cell-centre symbol — is gone once `decomposePar` has run.

```sh
runApplication decomposePar
runApplication -s cellCentres postProcess -func writeCellCentres -time 0
runApplication -s cellVolumes postProcess -func writeCellVolumes -time 0
runParallel laserbeamFoam
```

## 6.6 Missing post-processing scripts

| Script | For |
|---|---|
| `height_profile.py` | Fig. 8, h(x) along the track plus deposit length from h > 0.1 mm |
| `cross_section.py` | Fig. 18a, α = 0.5 contour as a polyline via `sampledInterface` |
| `flow_field.py` | Figs. 10, 11, with scalar diagnostics rather than eyeballed streamlines |
| `energy_budget.py` | arc in, latent, sensible, radiative, convective, closure error |
| `compare_runs.py` extension | ablation table and the SPSL row of Table 5 |

---

# 7. PROCESS NOTES — things that cost real time

## 7.1 Anything hand-set with `foamDictionary` is lost on `--force`

Bitten four times: `constant/g`, `libs ("libthermoTools.so")`, `arc_x0_m` (reset
from −0.010 to 0, putting arc-off 5 mm past the plate edge), and run settings.
All now live in `make_waam_case.py`. `endTime` and `maxDeltaT` remain exceptions
and must be re-set after every `--force`.

## 7.2 `--force` calls `shutil.rmtree(CASE)`

Far more destructive than `Allrun -f`, which only clears `0/`, `processor*/`,
`log.*` and time directories. It deleted a completed A1′ run whose archive was
made **after** the rebuild, so the copy held dictionaries only. A guard now
refuses when `log.laserbeamFoam` exists. **Archive before `--force`.**

## 7.3 A value stated in one place and derived in another will drift

The build summary hardcoded "4 subdomains, simple (4 1 1)" while the dictionary
was written from `DECOMP`. It reported the wrong thing for a full rebuild cycle.

## 7.4 Decomposition is not a free choice

See §2.2. Freeze `DECOMP` across any set of runs being compared, and never use a
y-cutting decomposition while y-symmetry is a diagnostic.

## 7.5 Read the log before proposing a cause

See §1.4. OpenFOAM prints the full valid-type table on a BC failure.

## 7.6 Bugs hide until the run is long enough

The fine region stopped at x = +20 mm while the arc reaches +30 mm, putting the
arc-off collapse in 1.0–1.6 mm cells — exactly where `deposit length` is
measured. Invisible because every run stopped at t ≤ 0.5 s, where the arc has
only reached x = 0.

## 7.7 Do not invent comparison targets

See §1.1. If a number is not in the text, a table, or measurable from a figure,
it does not exist. Mark inferences as inferences.

---

# 8. FROZEN CONFIGURATION

`endTime` and `maxDeltaT` reset on every `--force`.

```
mesh        98 x 50 x 18 = 88,200 cells   (0.5 mm core)
            x[-15,+35] y[-15,+15] z[-4,+5] mm = Zhao's 50 x 30 x 4 mm coupon
            fine region x[-12,+32] y[-6,+6], EXPAND 2.0
arc_x0_m    -0.010        bead -10 -> +30, 5 mm lead-in and run-out
DECOMP      (6, 1, 1)
endTime     0.66          developed window opens at t > 0.4645
deltaT      4e-05   maxDeltaT 4e-05   maxCo 0.1   maxAlphaCo 0.1
libs        ("libthermoTools.so")
gravity     (0 0 -9.81)
top         U pressureInletOutletVelocity; T/alpha inletOutlet
walls       externalWallHeatFluxTemperature, h = 80, Ta = 300
```

**Record per run:** A, h crown, h max, w toe, w half, ripple, `d sol`, `d liq`,
pool cells, U metal, U gas, bottom, far met, E kept.

## Planned ablation, after §3 is settled

| Case | Arc pressure | Lorentz |
|---|---|---|
| A1′ | off | off |
| A2 | **on** | off |
| A3 | off | **on** |
| A4 | **on** | **on** |

All four must share one frozen configuration, including whatever arc power and
mesh §3 resolves to.
