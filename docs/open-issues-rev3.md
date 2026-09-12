# Open issues, deviations and deferred fixes

WAAM SPSL replication of Zhao et al., *Welding in the World* 65:1571–1590 (2021).

**Revision 3.** Supersedes rev 2. The heat-input question is now settled: the
arc source is **790 W** (η₁ = η − η₂, Zhu Eq. 8), confirmed twice from
measurement. Rev 2's "790 W gives no pool" is retracted.

Each entry records **what**, **evidence**, **magnitude**, **what to do**.
Inferences are marked as inferences. Things previously stated wrongly are
marked RETRACTED with what replaced them.

---

# 0. TARGETS

## 0.1 From Table 5 and mass conservation — high confidence

| Quantity | Experiment | Zhao simulation |
|---|---|---|
| Height at mid-bead | 1.80 mm | 1.55 mm |
| Deposit length | 43.28 mm | 42.47 mm |
| Cross-section A | — | 8.48 mm² (= Q_v/v, mass-conserved) |

## 0.2 From Fig. 18a — high confidence

Substrate flat at z = 4 mm, crown at ~5.5–5.6 mm, so height **1.5–1.6 mm**,
consistent with Table 5. Base half-width ~4 mm, so full width **~8 mm**.

**Fig. 18 cannot show penetration.** The red half is metal (α > 0.5), and
remelted substrate is metal before and after fusion, so the fusion line is
invisible in the simulated half by construction.

## 0.3 Zhao's bead is parabolic — derived, high confidence

For a parabolic cap, h = 3A/(2w):

| | A | w | Predicted h | Reported h | |
|---|---|---|---|---|---|
| **Zhao** | 8.480 | 8.200 | 1.551 | 1.550 | **−0.1%** |

Zhao's own numbers are self-consistent with a parabolic cross-section to one
part in a thousand. **This is the sharpest available shape target**, and it
does not depend on reading any figure.

## 0.4 From Fig. 12 Case 2 — medium confidence, and it settled the heat input

Case 2 is 4.04 s of continuous deposition, no interlayer idle. The substrate
reads **green, ~570–680 K** on the 400–906 K scale.

| Budget | Bulk rise over 4.04 s | Plate at |
|---|---|---|
| 1512 W arc + 722 droplets (103% of U·I) | 486 K | 786 K |
| **790 W arc + 722 droplets (70%)** | **329 K** | **629 K** |

The figure supports the 70% budget; 786 K would be orange on that scale.

## 0.5 From Fig. 9 b1–b5 — LOW CONFIDENCE, do not compare against these

Pixel-colour estimates from a printed figure: `d sol` 1.0–1.3 mm at the pool
maximum, bottom face 435–570 K.

**Flagged as unreliable.** A Rosenthal thin-plate calculation for 1512 W into
4 mm aluminium gives 730 K at 2 mm from the arc and 587 K at 3 mm, which
contradicts the 435–570 K reading. One of the two is wrong and the colour
estimate is the weaker. Do not use these as targets.

---

# 1. RETRACTIONS

## 1.1 RETRACTED: "790 W gives no continuous weld pool"

Rev 2 recorded this as a measured finding. It was drawn from runs stopped at
**t ≤ 0.12 s**, where the arc has moved 2.4 mm into a 40 mm bead and the plate
is at 311 K — the cold-plate start transient. Zhao's earliest pool evidence is
Fig. 7a at **t = 0.52 s**; there was no evidence at all about his model at
t = 0.06.

**Replaced by §2.1**: the pool establishes at t = 0.2676 and holds for 20
consecutive writes.

Process failure: the failure mode had been predicted in advance, so the first
zeros were read as confirmation rather than as a question about whether the
test was valid yet. See §7.7.

## 1.2 RETRACTED: "Zhao's penetration is ~1 mm, from Figs. 11 and 18"

Invented. Zhao reports no penetration depth, and Fig. 18 cannot show one.
Replaced by §0.5, itself now marked low-confidence.

## 1.3 RETRACTED: "η_total = 0.85 is the most defensible next value"

Reached for because it would make the numbers work. Zhao states 0.70 and the
goal is replication. Replaced by §2.2: 0.70 as a **total** is correct, and the
earlier failure was a premature stop, not a wrong efficiency.

## 1.4 RETRACTED: "Goldak c = 4 mm is the suspect"

The power-fraction arithmetic holds (at c = 4 mm, 21% of the source lands below
2 mm depth, against 1.4% at c = 2 mm) but **Zhao uses c = 4 mm too**. A
difference cannot be explained by a shared parameter.

## 1.5 RETRACTED: two wrong causes for the `externalWallHeatFluxTemperature` failure

Both asserted before the log was read. Actual cause: the type is not registered
in laserbeamFoam; it ships in `libthermoTools.so`, which loads cleanly.

---

# 2. RESOLVED

## 2.1 RESOLVED: the arc source is 790 W, and it sustains a pool

`waam-arc-01-0p5mm-790W`, 88,200 cells, run 0 → 0.1476 then restarted to
0.6476.

**Pool continuity.** Last zero in `pool cells` at t = 0.2476. Every one of the
20 writes after it is non-zero, 16 to 156 cells.

**Power scaling, confirmed twice from measurement:**

| | Measured ratio vs the 1512 W case | Predicted |
|---|---|---|
| t = 0.02, before the first droplet (0.0333), arc only | **0.520** | arc ratio 790/1512 = **0.522** |
| t = 0.65, droplets included | **0.670** | total ratio 1512/2234 = **0.677** |

The first is the cleaner test: with no droplet yet landed, only the arc is
heating, and the rise ratio lands on the arc power ratio to 0.4%.

## 2.2 RESOLVED: the η₁/η₂ split

Zhao Eq. 27 (PDF p. 7) calls η "the arc efficiency for heating molten pool,
**excluding heating droplets** … described in detail in Ref. [32]", Table 4
sets it to 0.7 citing Wu [29].

Ref. [32] is Zhu, Cheon, Tang, Na & Cui (2018), *Int J Heat Mass Transf*
126:1206–1221, Eqs. 7–9:

```
q_arc = η₁·U·I / (2πσ_xσ_y) · exp(...)                            (7)
η₁ = η − η₂                                                       (8)
η₂ = 4πr_d³ρ[C_s(T_s−T_0) + C_l(T_d−T_l) + h_sl]·f_d / (3UI)      (9)
```

verbatim: *"η₁ is the arc efficiency for heating the molten pool. The TOTAL arc
efficiency η ranges from 0.60 to 0.85 for GMAW process"*; *"η₂ represents the
percentage of arc heat absorbed by the droplets"*.

So η is the **total and includes the droplets**. Zhao's 0.70 is a total placed
in the η₁ slot. Taken literally it gives η = 0.70 + 0.334 = **1.034**, above
unity and outside the range stated in the sentence Zhao cites.

η₂ from Eq. 9 with Zhao's numbers: 0.450 g/s, enthalpy 541 (solid) + 358
(latent) + 708 (liquid superheat) = 1607 kJ/kg → **722 W = η₂ 0.334**.

**η₁ = 0.366 → arc_power_W = 790 W.** Cross-check: 790 + 722 = 1512 W = 0.70·U·I.

Implemented in `make_waam_case.py`, computed from Zhu Eq. 9 rather than
hardcoded, with a guard that errors if η₁ ≤ 0 and a check that the implied
droplet frequency matches Zhao's 30 Hz.

## 2.3 RESOLVED: gravity on the wrong axis

`constant/g` was copied from `Plate2D-mhd`, a 2-D x–y donor with y vertical,
into five z-up WAAM cases. Confirmed fixed by the A1′ transverse section:
h(−y) vs h(+y) summed to −1.7%, centroid **−16.9 µm = 3.4% of one 500 µm cell**.
Symmetry is forced to exactly 1.000 by invariance under y → −y.

History of h_max/h_crown: 1.062 (scotch, g=−y) → 1.031 (simple 4 1 1, g=−y)
→ **1.003** (simple 6 1 1, g=−z, solidified bead).

## 2.4 RESOLVED: scotch decomposition produced a false asymmetry

Its subdomains are not mirror images, so it seeds a +y/−y difference every
timestep — an 8.9% signal nearly attributed to physics. `simple (6 1 1)` slices
x only, so mirror symmetry is exact.

## 2.5 RESOLVED: `E kept` above 100%

The denominator counted arc energy only while droplets enter as mass at 1506 K.
Fixed: denominator is now arc + droplet enthalpy. Reads 78–86%, and the
shortfall is **latent heat**, which `Ein = ρc_p(T−T₀)V` omits — latent is 22% of
the droplet enthalpy. See §6.2.

## 2.6 RESOLVED: the plume was trapped under a lid

`inlet` had `U = fixedValue (0 0 0)`, hydrodynamically a lid. Fixed with
`pressureInletOutletVelocity` on U and `inletOutlet` on T and alpha.metal.
Implied gas velocity fell from ~4.05 to ~0.60 m/s; dt rose from 1.0e-5 to 4.0e-5.

---

# 3. THE REMAINING ERROR IS SHAPE, AND IT IS QUANTIFIED

## 3.1 Geometry at t = 0.6476, developed window x = −7.75 to −0.75 mm

| | Zhao | 1512 W | **790 W** | toward Zhao? |
|---|---|---|---|---|
| A (mm²) | 8.480 | 10.178 | **8.062** | yes |
| h crown (mm) | 1.550 | 2.788 | **2.258** | yes |
| w toe (mm) | 8.200 | 5.315 | **6.537** | yes |

**A is within 4.9% of the mass-conserved 8.48**, against +20.0% at 1512 W.
Metal is no longer piling up behind the arc.

## 3.2 The parabolic test isolates what is left

| | A | w | Parabolic h | Measured h | |
|---|---|---|---|---|---|
| Zhao | 8.480 | 8.200 | 1.551 | 1.550 | −0.1% |
| 1512 W | 10.178 | 5.315 | 2.872 | 2.788 | −2.9% |
| **790 W** | 8.062 | 6.537 | 1.850 | 2.258 | **+22.1%** |

With the right amount of metal present, the bead distributes it **22% more
peaked than parabolic**: crown too sharp, toes too thin. If it were parabolic
at A = 8.06, h would be 1.85 mm against Zhao's 1.55.

**This is now the whole remaining discrepancy**, and it is what the two
unimplemented forces are for.

## 3.3 What is missing, quantified

Laplace pressure holding the crown up: γκ = 2γ/R with R = 2.40 mm = **709 Pa**.

| Term | Pressure | Share of Laplace |
|---|---|---|
| Hydrostatic ρgh | 72.8 Pa | 10.3% |
| Marangoni | 15.5 Pa | 2.2% |
| **Arc pressure (missing)** | 227.8 Pa | 32.1% |
| **Pool Lorentz (missing)** | 435.0 Pa | 61.4% |

Why gravity was never the lever: capillary length √(γ/ρg) = 5.72 mm against a
2.4 mm half-width, so Bo = 0.176 and surface tension beats gravity 5.7×.
A/ℓ_c² = 0.26. The Young–Laplace perturbation from gravity is
Δκ/κ = ρgh/(γκ) = **10.3%**, matching the measured +9.8% h change on switching
−y → −z.

---

# 4. BLOCKING — required before the Table 5 comparison

## 4.1 Arc pressure (Zhao Eqs. 25, 26)

P_arc = μ₀I²/(8πσ_p²)·exp(−r²/2σ_p²), peak **227.8 Pa**, surface integral
1.82 mN. Patch written and tested (`apply_stage2_arcpressure.py`), **not
compiled or run**.

Goes in `pEqn.H` — the live path, since `momentumPredictor no` makes the
equivalent block in `UEqn.H` dead code (§5.1). Includes a z gate at
`arc_z0 + 3 mm` so it does not act on free-flying droplets, which already carry
a prescribed velocity (228 Pa is ~6× a droplet's own weight).

## 4.2 Pool Lorentz force (Zhao Eqs. 22–24)

Kumar & DebRoy closed form, valid at Rm ≈ 4×10⁻³ ≪ 1. Radial inward, axial
downward, scaled by (1 − (z−z₀)/L) with L = 4 mm. Peak ~**435 Pa**.

Note it **increases** penetration — Zhao: "brings more heat to the bottom of
molten pool … causing an increase in penetration". So it works against §4.3.

## 4.3 Melt-through over 2.5 s is MARGINAL — cannot be called

`bottom` reached 641.8 K at t = 0.6476, decelerating: **376 → 262 → 171 →
146 K/s** across successive 0.1 s windows. Margin to the 815 K solidus 173 K.

| If the rate | Reaches 815 K at |
|---|---|
| holds at 146 K/s | t = 1.83 s |
| decays to ~115 | t = 2.15 s |
| reaches the bulk rate 89 | t = 2.57 s |

A 2.5 s run sits inside that spread. For contrast, the 1512 W case at 810 K and
124 K/s would have crossed at **t = 0.70 s**, so the η₁ correction transformed
this without guaranteeing a pass.

---

# 5. KNOWN DEVIATIONS FROM ZHAO

## 5.1 `momentumPredictor no`, so the `UEqn.H` body-force block never executes

`fvSolution` line 128. Surface tension, `pVap` and buoyancy take effect through
`phig` in `pEqn.H` lines 31–33. **Any new body force must go in `pEqn.H`.**

## 5.2 Plate end faces are on `outlet`, not Eq. 29

`outlet` spans z = −4 to +5, so 4 mm of metal plus 5 mm of argon on one patch.
`inletOutlet` with φ ≈ 0 in solid metal reduces to zeroGradient, i.e.
effectively adiabatic. Coverage: 1900 mm² of 2140 mm² has Eq. 29 — **89%**, and
the missing 11% is worth **0.18%** of the energy budget. The end faces sit 5 mm
from the bead ends and √(αt) = 5.9 mm, so the thermal wave reaches them.

## 5.3 Droplet detachment is prescribed, not resolved

Spheres injected at r = 1.105 mm, 30 Hz, v = 0.9 m/s, T = 1506 K, all from
Zhao's own results. Figs. 4, 5, 6 and the 0.88% radius validation cannot be
reproduced; everything downstream of impact can.

## 5.4 Zhao's plasma drag density looks wrong

Eq. 20 states ρ_f ≈ 6×10⁻⁶ kg/m³ "argon density". Argon is 1.62 kg/m³ at 300 K
and ~0.05 at arc temperature. Probably a typo for 6×10⁻².

## 5.5 Surface tension sourced from pure aluminium

γ = 0.85 N/m and ∂γ/∂T = −1.55×10⁻⁴ N/m/K from Mills' **pure aluminium**
chapter. Mg is strongly surface-active in Al. Marangoni drives the headline
result, so this is the largest physical uncertainty. Sensitivity runs:
γ = 0.70 N/m, ∂γ/∂T ±30%.

## 5.6 No shielding gas inlet (Zhao Eq. 28)

Zhao imposes a 20 L/min annular argon jet, ~1–2 m/s downward. Not implemented.

## 5.7 Vaporisation disabled

`Tvap = 1e6` as a switch, matching Zhao's assumption 3. Monitor max(T metal).

## 5.8 `maxCo 0.1` is inherited, not chosen

Traced to the `Plate2D` tutorial via `Plate2D-mhd`. Upstream laserbeamFoam
tutorials disagree: Plate2D 0.1, LPBF 0.2, PowderBed 0.25,
fiveLasersOnSubstrate 0.5, each with a `//0.1;` comment marking it as tuned.

**Do not assume Zhao licenses a higher value.** He used a *fixed* dt = 4e-5
with no adaptive control, and Fluent's VOF can run implicit interface
advection, which is unconditionally stable past Co = 1. isoAdvector is
explicit. Four things constrain Co here, only one of which is the capillary
limit: explicit interface advection, explicit body-force sources in `pEqn.H`,
the stiff Darcy sink, and the melting corrector iteration count.

Test 0.1 vs 0.25 back to back before changing it. `maxAlphaCo` limits the
interface only and `maxCo` everything; relaxing the latter while keeping the
former tight targets the cost at the gas, which is what sets dt
(`U gas` ~10 m/s against `U metal` ~0.4).

## 5.9 Other property notes

- ρ constant at 2650; liquid Al-5Mg near the liquidus is ~2350, so the melt is
  ~13% heavy.
- ε_S = 0.3, ε_L = 0.15 estimates for oxidised Al.
- h = 80 W/m²K is a lumped coefficient absorbing the fixture conduction path
  (5–15 would be typical for natural convection alone). Say so.
- C_s = 1050, C_l = 1180 J/kg·K in Zhu Eq. 9 are integral averages used **only**
  for η₂; the solver uses c_p(T). ±10% on C_s moves η₂ by ±0.016, ~35 W.
- OpenFOAM v2506 + laserbeamFoam against Zhao's Fluent + UDFs; isoAdvector vs
  geometric reconstruction.

---

# 6. TOOLING

## 6.1 `ripple` conflates a trend with an oscillation — FIX BEFORE QUOTING

At t = 0.6476 it reported **51.2%**, which looks like severe humping. But A(x)
climbs monotonically across the window, 5.713 → 9.838. `ripple` is
peak-to-peak, so it cannot distinguish a gradient from droplet periodicity.
Detrended, the residual is 28.6%; adjacent-station wobble is ~5%.

**Detrend before reporting.** The 17.8% figure from the 1512 W case is subject
to the same problem.

## 6.2 `E kept` numerator omits latent heat

`Ein = ρ·c_p·(T−T₀)·V`. Latent is 22% of the droplet enthalpy, so a reading
near 80% is consistent with conservation. Add the latent term.

## 6.3 `E kept` denominator counts droplet power from t = 0

The first droplet lands at t = 0.0333, so early rows read low (53.0% at
t = 0.02, where the ceiling is 52.3%). Ramp it with injected mass. Cosmetic.

## 6.4 `d sol` at 0.25 mm is the mesh floor, not a measurement

From t = 0.5076 it sits at exactly 0.25 mm — one 0.5 mm cell. `d liq` is 0.00 or
negative throughout, so there is no fully-liquid substrate. **Penetration
cannot be compared with Zhao at this resolution.** Needs 0.25 mm cells.

## 6.5 DONE: `d liq` reported beside `d sol`

## 6.6 DONE: three stale gate texts corrected

`depth < 6 mm` → `d sol < 4 mm`; `bottom near 300 K` → `bottom < 815 K`;
`far met = 300.0` → no longer a domain-size gate.

## 6.7 `bead_profile.py`: report `y` at `h max`

`h max > h crown` has two causes — an off-centre crown, or a symmetric profile
with a central dip. The second becomes real once arc pressure creates a
depression.

## 6.8 `Allrun` does not fail on solver crash

`MPI_ABORT` fired and `Allrun` still ran `reconstructPar` and both
`postProcess` stages. Add:

```sh
runParallel laserbeamFoam
grep -q "^End" log.laserbeamFoam || { echo "ERROR: solver did not reach End"; exit 1; }
```

## 6.9 `Allrun` should write cell geometry right after `decomposePar`

`pool_stats_waam.py` has hit the missing-geometry error four times mid-run. The
hazard it was deferred for (a field named `C` colliding with OpenFOAM's
reserved cell-centre symbol) is gone once `decomposePar` has run.

## 6.10 `writeControl runTime` counts from `startTime`, not from zero

After a restart at t = 0.1475696396 the writes land on 0.1675696, 0.1875696 …
0.6475696. **t = 0.66 was never a write time.** Use `adjustableRunTime` if a
clean grid matters.

## 6.11 Missing post-processing scripts

| Script | For |
|---|---|
| `height_profile.py` | Fig. 8, h(x) plus deposit length from h > 0.1 mm |
| `cross_section.py` | Fig. 18a, α = 0.5 contour via `sampledInterface` |
| `flow_field.py` | Figs. 10, 11, with scalar diagnostics |
| `energy_budget.py` | arc in, latent, sensible, radiative, convective, closure |
| `compare_runs.py` extension | ablation table and the SPSL row of Table 5 |

---

# 7. PROCESS NOTES

## 7.1 Anything hand-set with `foamDictionary` is lost on `--force`

Bitten four times: `constant/g`, `libs`, `arc_x0_m` (reset from −0.010 to 0,
putting arc-off 5 mm past the plate edge), and run settings. All now in the
generator. `endTime` and `maxDeltaT` remain exceptions.

## 7.2 `--force` calls `shutil.rmtree(CASE)`

Far more destructive than `Allrun -f`. It deleted a completed A1′ run whose
archive was made **after** the rebuild, so the copy held dictionaries only. A
guard now refuses when `log.laserbeamFoam` exists. **Archive before `--force`.**

## 7.3 A value stated in one place and derived in another will drift

The build summary hardcoded "4 subdomains, simple (4 1 1)" while the dictionary
was written from `DECOMP`.

## 7.4 Decomposition is not a free choice

Freeze `DECOMP` across compared runs; never use a y-cutting decomposition while
y-symmetry is a diagnostic.

## 7.5 Read the log before proposing a cause

OpenFOAM prints the full valid-type table on a BC failure — definitive evidence
already on disk. Two wrong causes were asserted before it was read.

## 7.6 Bugs hide until the run is long enough

The fine region stopped at x = +20 mm while the arc reaches +30 mm, putting the
arc-off collapse in 1.0–1.6 mm cells — where `deposit length` is measured.
Invisible because every run stopped at t ≤ 0.5 s.

## 7.7 Do not call a run early on a metric still in its transient

**The most expensive mistake so far.** A run was killed at t = 0.12 because
`pool cells` read zero, and that was taken as proof 790 W could not sustain a
pool. It could: the pool establishes at t = 0.2676. Two consequences — a
0.35 mm mesh run was launched to test a hypothesis that only existed because of
the premature stop, and the `bottom` scaling result (which never depended on
the pool) was delayed by hours.

Before reading a diagnostic, state what condition makes it valid, and check
that condition is met.

## 7.8 Do not invent comparison targets

If a number is not in the text, a table, or measurable from a figure, it does
not exist. Mark inferences as inferences, and mark figure readings with their
confidence.

---

# 8. FROZEN CONFIGURATION

`endTime` and `maxDeltaT` reset on every `--force`.

```
mesh        98 x 50 x 18 = 88,200 cells   (0.5 mm core)
            x[-15,+35] y[-15,+15] z[-4,+5] mm = Zhao's 50 x 30 x 4 mm coupon
            fine region x[-12,+32] y[-6,+6], EXPAND 2.0
arc_power_W 789.7          eta_1 = eta - eta_2, Zhu Eq. 8
arc_x0_m    -0.010         bead -10 -> +30, 5 mm lead-in and run-out
DECOMP      (6, 1, 1)
endTime     0.66           developed window opens at t ~ 0.4645
deltaT      4e-05   maxDeltaT 4e-05   maxCo 0.1   maxAlphaCo 0.1
libs        ("libthermoTools.so")
gravity     (0 0 -9.81)
top         U pressureInletOutletVelocity; T/alpha inletOutlet
walls       externalWallHeatFluxTemperature, h = 80, Ta = 300
```

**Record per run:** A, h crown, h max, w toe, w half, ripple (detrended),
`d sol`, `d liq`, pool cells, U metal, U gas, bottom, far met, E kept, and the
parabolic ratio h/(3A/2w).

## Ablation, now runnable

| Case | Arc pressure | Lorentz | Status |
|---|---|---|---|
| **A1′** | off | off | **done** — `waam-arc-01-0p5mm-790W` |
| A2 | **on** | off | patch ready, not compiled |
| A3 | off | **on** | not written |
| A4 | **on** | **on** | — |

All four share the frozen configuration above.

---

# 9. CASE INVENTORY

| Directory | Arc | Mesh | g | Provenance |
|---|---|---|---|---|
| `waam-arc-01-0p5mm-790W` | 790 W | 88,200 | −z | **yes — A1′ reference** |
| `waam-arc-01-gminusY-baseline` | 1512 W | old | **−y** | yes — pre-fix control |
| `waam-arc-01-A1ref-1512W` | 1512 W | 88,200 | −z | yes — config only, no fields |
| `waam-arc-01-A1ref` | 1512 W | ? | −z | **missing** |
| `waam-arc-01-A1pre-Tclamp` | 1512 W | ? | −z | **missing** |
| `waam-arc-01-A1-gz-oldmesh` | 1512 W | ? | −z | **missing** |
| `waam-arc-01` | 790 W | 251,160 | −z | **missing** — 0.35 mm mesh test, t = 0.06 |
| `-3mm-argon`, `-emON`, `-half-arconly`, `-ref` | 1512 W | old | −z | pre-existing variants |

The four marked **missing** need provenance written from their actual contents,
not from memory. A dictionary with a guessed label is worse than none.
