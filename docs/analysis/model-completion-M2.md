# Model completion — Lorentz, droplet restart, wall radiation

**Date:** 2026-09-17
**Status:** `waam-arc-01-M2-full` running to t = 1.2 s (PID 141535, ~6–8 h), all five banners confirmed live.

---

## What went in

| # | change | how verified |
|---|---|---|
| 1 | **Lorentz force, Zhao Eqs. 22–24** (derived forms) | analytic peaks match an independently computed force table to 4 s.f. |
| 2 | **Droplet state survives restart** | write, read and accumulation checked end to end |
| 3 | **Wall radiation** → effective h = 88 W/m²K | arithmetic, no code change |
| 4 | `h_conv_Wm2K 80` — free-surface half of Eq. 29 | was compiled but inert |
| 5 | `darcy_q 1.0e-3` — Zhao Eq. 7 as written | solver default was 1e-12 |

Both gates ran **bit-identical twice** (after patch 1, and again after patch 2).

---

## 1. Lorentz force — the derivation

Not transcribable from the paper: **Eqs. 22–24 as printed have the radial and axial
expressions transposed**, plus a stray factor 2.

Gaussian surface current in Zhao's own convention (his 2π² prefactors fix σ_j as √2 × the
standard Gaussian σ):

```
j_z = (I/(pi sj^2)) exp(-r^2/sj^2) f(z),     f = 1 - d/L,  d = z0 - z,  0 <= d <= L
I_enc = I [1 - exp(-r^2/sj^2)]
B_th  = -(mu0 I/(2 pi r)) [1 - exp(-r^2/sj^2)] f
```

Continuity (div J = 0) fixes the radial spreading current, and F = J × B gives:

$$F_r = -\frac{\mu_0 I^2}{2\pi^2\sigma_j^2 r}E(1-E)f^2 \qquad F_z = -\frac{\mu_0 I^2}{4\pi^2 L r^2}(1-E)^2 f$$

with E = exp(−r²/σ_j²). Both negative — **inward pinch, downward drive. This is a penetration
term, not a widening one.**

### Six signatures of the transposition

| feature | belongs to | Zhao puts it on |
|---|---|---|
| E(1−E) envelope | radial | Eq. 24, axial |
| (1−E)² envelope | axial | Eqs. 22–23, radial |
| (1−z/L)² depth | radial | Eq. 24, axial |
| (1−z/L)¹ depth | axial | Eqs. 22–23, radial |
| 1/σ_j² prefactor | radial | Eq. 24, axial |
| 1/(Lr²) prefactor | axial | Eqs. 22–23, radial |

His Eq. 24 **is** the correct radial expression; his Eqs. 22–23 **are** the correct axial
expression wearing a spurious x/r.

**Separately:** Eq. 23's denominator reads 4π²r³**2**L against Eq. 22's 4π²r³L — halving F_y and
breaking axisymmetry under an axisymmetric arc.

**Sign trap:** Zhao writes (1 − (z−z₀)/L), valid only if z is depth positive *downward*. This
solver has z positive *up* with the plate below z₀ = 0; transcribed literally the factor reaches
**2** at the plate bottom instead of 0. Implemented as (1 − (z₀−z)/L), clamped to the slab.

### Validation

Solver startup prints, and these reproduce the §3.3 force table computed independently months
apart:

```
analytic peaks  F_r 50978.02  F_z 14756.80 N/m^3     (expect 5.10e4, 1.48e4)
mesh max |F|    45444.70 N/m^3                        = 89.2% of peak
```

**The 11% shortfall is depth sampling, not radial.** F_r peaks at r = 0.6σ_j = 1.20 mm and the
nearest cell centres sit at r ≈ 1.25 mm, still **99.3–99.6%** of peak — the function is flat near
its maximum. But F_r carries (1−d/L)², and the topmost substrate cell centre is at d = 0.25 mm
where that factor is already **0.879**. 0.993 × 0.879 = 0.873 against 0.892 observed.

Refining to 0.25 mm would put that centre at d = 0.125 mm and recover the peak to ~94% — a 5%
gain. Useful for the mesh chapter, not a reason to refine.

### Architecture

`Jarc` and `Bself` are stored as fields and `arcLorentz = Jarc ^ (Bself + Bext)`, so **Phase 3
crosses the same current with an applied field** — no new physics needed.

Not using the existing numerical `mthd->JcrossB()` potential solve, because it isn't needed:
Rm = μ₀σUL = **9×10⁻³**, and the motional current σUB at 10 mT is ~1.8×10⁴ A/m² against the arc's
peak j_z of **1.07×10⁷ A/m²**. The arc current dominates by three orders. The term is added
*beside* JcrossB, never in place of it, so the two can never both carry the arc.

---

## 2. Droplet restart — measured cost of the bug

`dropOwed_`, `dropImpactEnd_`, `dropNum_` were plain C++ scalars, re-initialised on every restart.

Resuming M1-mushy at t = 1.16 s discarded **4.63 mm³ = 0.82 of a droplet**. Evidence:

```
time   deposited    drops          time   >sol   >liq
1.14    190.253     33.66          1.14    428     63
1.16    190.253     33.66          1.18    323      0     <- no liquid anywhere
1.18    190.253     33.66          1.22    270      1     <- -37% from 1.14
1.20    195.747     34.64          1.36    408     69     <- recovered
```

Three consecutive writes with no new metal, 0.04 s during which 1.2 droplets were due. The pool
lost 37% of its volume and took 0.2 s — **4 mm of track** — to rebuild. The bead shows it as a
trough to 5.2 mm² against the conserved 8.482, with a compensating overshoot to 12.1 mm². The run
finished 59 droplets instead of 60.

**Fix:** `IOdictionary` in `<time>/uniform/dropletState`, READ_IF_PRESENT / AUTO_WRITE. Only the
three *accumulated* values are stored; everything else is re-read from `caseParameters` and would
only drift if duplicated.

**Verified end to end** on a 5-step serial run then a restart:

| | value | expected |
|---|---|---|
| written at t = 5e-5 | dropOwed 8.4823e-12 m³ | 169.646 mm³/s × 5e-5 = 8.4823e-3 mm³ ✓ |
| restored on restart | "restored from 5e-05/uniform/dropletState" | ✓ |
| accumulated to t = 1e-4 | **1.69646e-11 m³** | 1.69646e-2 mm³ ✓ (a reset would give 8.48e-12) |

---

## 3. Wall radiation — why no patch

`externalWallHeatFluxTemperature` in `coefficient` mode gives h(T−Tₐ) and nothing else. Radiation
linearises exactly:

$$\varepsilon\sigma(T^4-T_a^4) = \underbrace{\varepsilon\sigma(T^2+T_a^2)(T+T_a)}_{h_{rad}}(T-T_a)$$

At measured wall temperatures with ε = 0.3: **12.0 W/m²K** at the 759 K bottom, **4.2** at the
478 K coupon edge. So `h uniform 88` (80 convective + ~8 radiative).

**Justification for not doing it properly:** ε for aluminium runs 0.05 (polished) to 0.4
(oxidised) — a factor of eight, against a linearisation error of ~1.5×. Implementing T⁴ exactly on
top of an emissivity uncertain by ×8 is false precision, and the whole term is **2–3 W against
1836 W in (0.15%)**.

---

## Force inventory — for the thesis

At measured conditions (|U| ≈ 0.34 m/s, pool half-width 1.85 mm, bead h = 2.7 mm).

### Body forces in UEqn (N/m³)

| force | magnitude | vs gravity |
|---|---|---|
| Darcy sink, ε₁ = 0.3 | **6.2×10⁶** | 240× |
| surface tension, droplet | 2.3×10⁶ | 89× |
| surface tension, bead crown | 1.0×10⁶ | 39× |
| droplet impact, peak (7.4% duty) | 9.7×10⁵ | 37× |
| Darcy sink, ε₁ = 0.5 | 6.8×10⁵ | 26× |
| arc pressure | 3.4×10⁵ | 13× |
| inertia ρU²/L | 1.7×10⁵ | 6.4× |
| Marangoni, at droplet impact | 8.6×10⁴ | 3.3× |
| droplet impact, time-averaged | 7.2×10⁴ | 2.8× |
| **Lorentz, radial** | **5.1×10⁴** | 2.0× |
| gravity ρg | 2.6×10⁴ | 1.0× |
| **Lorentz, axial** | **1.5×10⁴** | 0.57× |
| Marangoni, between droplets | 1.4×10⁴ | 0.55× |
| viscous μ∇²U | 1.8×10³ | 0.07× |
| buoyancy ρgβΔT | 90–530 | 0.003–0.02× |

**Read with care.** Surface tension, Marangoni and arc pressure are CSF-converted — their
body-force value is a stress × |∇α| ≈ 1500 m⁻¹ on a 0.5 mm mesh, and rises on refinement while
their integrated effect stays fixed. Comparing them to volumetric forces as body forces is apples
to oranges.

### Surface stresses — the fair comparison (Pa)

| stress | magnitude |
|---|---|
| Laplace, droplet (2σ/R) | **1538** |
| droplet impact ½ρv² | 1073 |
| Laplace, bead crown (σκ) | **680** |
| arc pressure peak | 228 |
| hydrostatic head ρgh | 70 |
| Marangoni shear, at impact | **57** |
| Marangoni shear, quiescent | **10** |

**This explains the M1 null result.** Marangoni shear is 10–57 Pa against 680 Pa of capillary
pressure holding the bead's shape — **1.5% to 8%**. It stirs the pool vigorously (Ma ≈ 610, well
above the ~80 onset, hence +34% velocity) but cannot reshape a surface governed by an
order-of-magnitude larger stress.

Bond ρgL²/σ = **0.21** · Weber ρv²D/σ = **5.6** · Marangoni/Laplace = **0.015–0.08**

### Expectation for the Lorentz term

2× gravity, comparable to Marangoni at impact, **two orders below both Darcy and capillarity**.
Expect `d sol` to move off 1.75 mm and `w fuse` to hold or narrow. It belongs in the model because
it is in Eq. 21; it will not close a 33% width gap.

---

## Gate procedure — recorded, because reconstructing it cost five round trips

The pieces were scattered across three places.

**Gate A** — Hartmann, MHD path. Checker lives **inside the case**, not `tools/`:

```bash
cd cases/00-tutorials/hartmann1d-postpatch
./Allclean && ./Allrun
postProcess -latestTime -func sample
python3 hartmann_check.py postProcessing/sample/10/line_centreProfile_H_U.csv \
        --sigma 1000 --mu 1
```

**The `--sigma 1000 --mu 1` is mandatory** on this copy. Without it the script assumes unity,
computes Ha = 0.500 instead of 15.810, compares a genuine plug profile against a parabola, and
reports **FAIL at 6.2903e-01 on a case that passes at 1.3869e-02**. `patch_hartmann_check.py`
exists to fix this and has been applied to `hartmann1d/` but not to `hartmann1d-postpatch/`.

Expected: `u(0)` **0.0158753** · `Ha` **15.810** · `max err` **1.3869e-02**

**Gate B** — melt regression. Tool is `cases/00-tutorials/pool_stats.py`, **not**
`tools/pool_stats_waam.py` (which is WAAM-geometry only and refuses Plate2D). Takes case names as
arguments, run from `cases/00-tutorials/`:

```bash
cd cases/00-tutorials/Plate2D-mhd-postpatch && ./Allclean && ./Allrun
cd .. && python3 pool_stats.py Plate2D-mhd-postpatch
```

Expected, from `results/gateB-edit8.txt`: **1869.0 / 0.0504 / 3.4135 / 0.0552 / 12.7787**
("melt volume" is **cells**, not mm³.)

**`max|U| whole domain` = 12.7787 is the strongest check available.** The tool's own docstring
gives it a 25% noise floor and calls it useless for comparison — so a measure that noisy landing
on the same four decimals means the solution is byte-for-byte identical, not merely inside
tolerance.

**Always `./Allclean` first.** A previous false pass came from `Allrun` printing "already run"
while the stats tool read stale fields and returned identical numbers — which is exactly what a
passing regression looks like.

**Worth building:** a single `tools/run_gates.sh` that runs both and diffs against the banked
files. Three patches landed today and each one cost the same reconstruction.

---

## M2-full case spec

`cases/00-tutorials/waam-arc-01-M2-full`, from M1-mushy, 1114 W, 88,200 cells, 6 subdomains,
`writeInterval 0.02`, `endTime 1.2`.

```
constant/transportProperties   (TOP LEVEL, not inside metal { })
    h_conv_Wm2K       80;
    darcy_A           1.0e6;      // explicit; swept separately in M3
    darcy_q           1.0e-3;     // Zhao Eq. 7 as written

constant/caseParameters        (arcDict)
    arc_lorentz       true;
    sigma_j_m         0.002;
    substrate_L_m     0.004;

0/T   substrateBottom, farField
    h   uniform 88;               // 80 convective + ~8 linearised radiative
```

Placement verified 3 / 0 / 3 / 2 — the second number being keys that landed inside `metal { }`,
which would have been silently ignored.

---

## Next

1. **M2-full lands (~6–8 h).** Compare against M1 at **t = 1.14**, where a settled profile and a
   verified mass control already exist. This isolates *what completing the model did*.
2. **M3-darcy, `darcy_A 1.0e5`**, everything else identical. Isolates the Darcy constant.
   10⁶ is inherited from laserbeamFoam — not Zhao's, who gives only the 10⁴–10⁷ range while
   paraphrasing the Fluent manual, whose default is 10⁵.
3. **1512 W replication**, single launch to 2.0 s, no restart.
4. **Mesh convergence** 0.5 / 0.35 / 0.25 mm. Zhao's full-domain equivalent is ~513,000 cells
   against 88,200.

**Prediction for M2, on the record:** Lorentz deepens (`d sol` off 1.75 mm), `h_conv` and the
wall-h change both remove heat so the pool shrinks slightly, `darcy_q` permits a little solid
creep. Net: marginally **narrower**, not wider. That is the expected cost of completing the model
honestly, and it should be stated as such rather than hidden.
