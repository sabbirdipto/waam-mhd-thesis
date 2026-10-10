# M4: literal η = 0.70 melts into the plate underside, and why Zhao's heat-source geometry doesn't transfer

> **CORRECTION (2026-09-19, after the plate-temperature check, §2).** An earlier
> version of this doc recommended adopting M4's settings as the replication
> baseline. **That recommendation is withdrawn.** At 1512 W the plate underside
> sits at 846–859 K, above the 815 K solidus, at every write and is still
> climbing. The solidus front reaches the bottom cell layer, and liquid
> penetration is 2.5× Zhao's. His experimental coupon did not melt through, so
> M4 cannot be the baseline. The corrected plan is in §7.

Run: `waam-arc-01-M4-1512W`. Single launch from 0 to 1.2 s, 120000 steps,
24093 s wall (6.7 h), 6 ranks. Finished cleanly. Reconstructed 1.06–1.20 for
analysis. Tool outputs banked in `~/thesis/results/M4-1512W.*` and
`M2-vs-M4-1.14.compare.txt`.

**One-line change from M2.** It was cloned from M2, so `darcy_A` comes back to
1.0e6 by inheritance rather than by an edit. A full recursive diff against M2
returns a single hunk:

```
< arc_power_W       1113.7;    // eta_1*U*I, NOT eta*U*I.
> arc_power_W       1512.0;    // eta_1 = 0.70, Zhao Eq.27 read literally.
```

The launch banner confirms the value took effect:
`Arc: 1462.48578 W in metal below the surface, renormalised x1.033856206 to 1512 W (raw 96.72524999% of target)`.

---

## 0. M4 is the faithful replication of what Zhao states

Zhao §3.5, directly under Eq. 27:

> "η is the arc efficiency for heating molten pool, **excluding heating
> droplets**, which is described in detail in Ref. [32]."

So Zhao states that η = 0.70 in Eq. 27 is **pool-only** (our η₁), with droplet
heating kept separate. M4 does exactly that. The ~3.4% energy overshoot
(1512 + 722 = 2234 W against U·I = 2160 W) belongs to Zhao, not to us. His
Table 4 cites 0.7 to ref [29], a welding textbook, where 0.7 is a total GMAW
arc efficiency. Applying a total efficiency as a pool-only one is errata 5.

The faithful replication is precisely what fails in §2. **Zhao's stated
parameter set, run in a correctly normalised solver, melts into the underside
of Zhao's own 4 mm coupon.** That is a replication result in its own right,
alongside errata 3 and 5.

---

## 1. Pool width: power moves it, but can't take it far enough

`pool_width`, eight writes from 1.06 to 1.20 s. The means are over those eight
writes.

| | w fuse (mm) | w liq (mm) | vs Zhao 8.2 |
|---|---|---|---|
| M2 (1113.7 W, A=1e6) | **3.560** (3.504–3.645) | 0 at 5 of 8 writes | −56.6% |
| M3 (1113.7 W, A=1e5) | 3.650 (3.626–3.671) | 0 at all 8 | −55.5% |
| **M4 (1512 W, A=1e6)** | **4.956** (4.927–5.044) | **2.367**, non-zero at all 8 | **−39.6%** |

**w fuse +39.2% for +35.8% power**, an exponent of n = 1.08. There is a liquid
core at every write (w liq 1.93–2.86 mm), which M2 had only intermittently and
M3 never had. `w melt = w fuse` at every write, so there is still no
undercutting.

**M3 creep verdict is airtight.** M2's late w fuse (3.560) and M3's (3.650)
agree within scatter. Changing the Darcy constant 10× had no effect on the pool.
M3's doubled bead was entirely creep over unmelted plate.

### Power cannot reach Zhao's width

With the Goldak half-width fixed at b = 4 mm, this is the arc-source power
needed to reach Zhao's width:

| fit through (1113.7, 3.560), (1512, 4.956) | to 7.76 mm (printed) | to 8.2 mm (inverted) |
|---|---|---|
| power law, n = 1.08 | 2289 W | 2408 W |
| linear | 2312 W | 2438 W |

Every reading exceeds U·I = 2160 W for the arc source alone. §2 then shows the
depth rises about 3× faster than the width, so the plate would be molten
through its whole thickness long before the width got near 8 mm. **Power is
ruled out twice over.**

---

## 2. Plate temperatures: M4 melts into the underside

`pool_stats_waam`, eight writes from 1.06 to 1.20 s. Saved as
`M4-1512W.pool_stats.txt` and `M2-full.pool_stats.txt`.

| | M2 (1113.7 W) | M4 (1512 W) | limit |
|---|---|---|---|
| **bottom face T** | 746.4 → 751.6 K | **846.2 → 859.3 K** | **< 815 K solidus** |
| bottom climb rate | ~37 K/s, linear | ~94 K/s | a plateau, not a climb |
| **d sol** (to 815 K) | 1.75 mm | **3.75 mm, saturated** | < 4 mm plate |
| **d liq** (to 906 K) | 0.75–1.25, mean **0.88 mm** | 2.25–2.75, mean **2.50 mm** | Zhao **~1 mm** |
| pool volume (> solidus) | 44–51 mm³, mean 47.9 | 128–142 mm³, mean 136.6 | |
| T metal between droplets | 922–992 K | 1008–1037 K | < ~1600 K |
| far-edge metal T | 396 → 412 K | 416 → 435 K | well below solidus |
| energy kept in domain | ~81% | ~79% | |

**M4 fails the bottom gate at every write.** The underside is 31–44 K above the
solidus, so the bottom of the plate is mushy under the pool. `d sol = 3.75 mm`
is the deepest cell centre on this mesh (8 layers of 0.5 mm, centres from −0.25
to −3.75). The value is saturated: the solidus front has reached the bottom cell
layer and cannot be measured any deeper. On the current trend (~94 K/s) the
underside crosses the **906 K liquidus near t ≈ 1.7 s** and would be fully
liquid well before Zhao's 2.0 s track ends. Zhao's macrograph shows no
melt-through.

**M2 passes, with margin.** Its underside is 63 K below the solidus and climbing
linearly at ~37 K/s. That reaches only ~781 K at t = 2.0 s, so it stays below
the solidus through Zhao's full 40 mm track.

**Penetration is the decisive comparison.** `pool_stats_waam` measures
`d liq` to 906 K, the same boundary Zhao colours in Figs. 11 and 18, where his
liquid penetration reads about 1 mm. **M2's 0.88 mm matches Zhao within a cell.
M4's 2.50 mm is 2.5× too deep.**

### Where the extra power went

| M2 → M4 (+36% power) | ratio |
|---|---|
| w fuse | ×1.39 |
| d liq | **×2.86** |
| pool volume | **×2.85** |

The extra power went mostly **down**, not out. The pool's shape
(depth : width) is wrong, not just its size.

### The mechanism: c = 4 mm puts half the arc power below the pool

The Goldak depth distribution is ∝ exp(−3z²/c²). The fraction of the source
deposited deeper than depth d is erfc(√3·d/c):

| c | below 1 mm | below 2 mm | beyond 4 mm |
|---|---|---|---|
| **4 mm (Zhao, current)** | **54%** | 22% | 1.4% |
| 2 mm | 22% | 1.4% | ≈0 |
| 1.5 mm | 10% | 0.03% | ≈0 |

With Zhao's calibrated c = 4 mm, the depth semi-axis equals the entire plate
thickness. **More than half the arc power is injected deeper than 1 mm, below
Zhao's whole liquid pool.** At 1512 W that is about 800 W heating the plate from
the inside. That heat mostly warms the plate body and drives the underside
upward rather than melting a wide pool at the surface. It explains why depth
responds nearly 3× as strongly as width.

---

## 3. Why Zhao's heat-source geometry doesn't transfer

Zhao §3.5, same paragraph as the η statement:

> "The geometric dimensions of heat source, af, ar, b, and c are **obtained by
> comparing the calculated and experimental morphologies of molten pool** of
> first deposition layer, which are 3mm, 6mm, 4mm, and 4mm."

The values 3/6/4/4 mm are **calibration parameters**, fitted so that *his*
implementation reproduced *his* experimental pool. His implementation differs
from ours in at least three ways, and a calibration would absorb each of them:

1. **Eq. 27 as printed has no z-restriction** (errata 3). The 6√3 normalisation
   integrates to ηUI only below the surface. Applied over all space it delivers
   up to 2ηUI, with the upper half landing in the bead and droplets above
   z = 0. Our source is correctly half-ellipsoidal and renormalised to exactly
   ηUI in metal. If his code followed the printed equation, his heat
   *distribution* differed from ours even at the same stated η.
   **Consistent, not proven.**
2. **Mesh.** Zhao's minimum cell is 0.25 mm (§3.7, 256,512 cells,
   Δt = 4e-5 s). Ours is 0.510 × **0.600** × 0.500 mm, 2.4× coarser in y.
3. **Solver and properties.** His is Fluent with UDFs; ours is laserbeamFoam
   with the two-branch kappa. Both take k and cp from Mills [33].

Zhao also concedes a width discrepancy himself (§4.4), attributing it to "the
inherent shortcoming of the enthalpy-porosity method", which is the creep M3
exposed. His Table 5 reports only **height and length**, with no width.

**Implication:** recalibrating the heat-source geometry against experiment is
**Zhao's own stated procedure**, so it counts as replication, not fitting.
Carrying his fitted a/b/c into a different implementation and expecting his
pool is the actual inconsistency. §2 shows the calibration has to cover
**shape**, meaning b and c, not only width.

---

## 4. Bead geometry

### Station-matched, full common window

21 stations from −7.25 to 2.75 mm.

| | M2 | M4 | Δ |
|---|---|---|---|
| A (mm²) | 10.756 | 10.150 | −5.6% |
| h (mm) | 2.815 | 2.403 | −14.6% |
| w toe (mm) | 6.052 | 6.728 | +11.2% |
| h/w | 0.4651 | 0.3572 | −23.2% |
| vs Zhao sim 0.189 | 2.46× | 1.89× | gap closed +39.1% |

### The start hump, and the settled region

M4's A(x) is 13.1–13.9 mm² at x = −7.25 to −5.75 and decays to about 8.0 by
x ≈ −0.75. At 1512 W the start transient runs about 9 mm, but `bead_profile`
allows only 5 mm. Its 58–60% ripple measures this hump, not droplet periodicity.

Settled stations only, −0.75 to 2.75 mm (8 stations):

| | M2 | M4 | Zhao |
|---|---|---|---|
| A (mm²) | 9.545 | 8.064 | 8.482 conserved |
| h (mm) | 2.679 | 2.020 | 1.55 sim / 1.80 exp |
| w toe (mm) | 5.552 | 6.317 | 7.76 printed / 8.2 inverted |
| h/w | 0.4825 | 0.3198 | 0.189 |

M4's bead is closer to Zhao's than M2's. §2 shows it gets there by melting the
plate through, so **its bead geometry is not usable as a result**.

### Why matching width fixes height for free

A is conserved at 8.4823 mm², and a parabolic cap has A = (2/3)·w·h. So
**w = 7.76 mm forces h = 3A/(2w) = 1.64 mm**, between Zhao's simulated 1.55 and
experimental 1.80. Calibrating on width brings the height along automatically.
(The same relation gives the 8.2 mm figure: 3A/(2 × 1.55) = 8.21 mm, which is
errata 1.)

---

## 5. Mass: injector correct, loss tracks interface smearing

Every droplet was placed at 5.500 or 6.000 mm³ (44 or 48 cells), with the
remainder carried forward (−0.347 / +0.152 mm³), so the running total stays
exact at 5.6517 mm³. **All 36 droplets were delivered by 1.2 s.**

| | deposited at 1.2 s | drops | lost | bead a>.5 / aV |
|---|---|---|---|---|
| M3 | 155.626 mm³ | 27.54 | 23.5% | 0.811 |
| M4 | 186.618 mm³ | 33.02 | 8.3% | 0.922 |

Loss tracks interface smearing, which points to loss during transport, such as
alpha bounding or clipping in thin films. It is mostly early: 4.6% in the
1.06–1.20 window against 8.3% cumulative. Still unexplained.

---

## 6. Runtime notes

- 0.208 s/step against M3's 0.165 s/step. Total 6.7 h against 5.9 h.
- `deltaT` never dropped below `maxDeltaT = 1e-5`. Courant (`maxCo 0.1`) did
  not bind even at 1512 W.
- T took 12.2 linear solves per step against M3's 9.0, a stiffer enthalpy
  equation with a mushy region about 2.8× larger.

---

## 7. Where this leaves the replication (corrected)

**Keep the arc at 1113.7 W.** It is energy-consistent (η_total = 0.85, with 324 W
booked as losses). It already gives Zhao's ~1 mm liquid penetration (0.88 mm),
and it keeps the underside below the solidus through the full 2.0 s track.
Power is not the missing ingredient.

**Recalibrate the heat-source SHAPE at that power, following Zhao §3.5.** The
width has to rise without the depth rising. Two levers, each tested alone
against M2 as the shared control:

| run | change from M2 | expected if the lever is right |
|---|---|---|
| **M5a** | **c = 2 mm** (from 4) | the ~54% of heat now below 1 mm returns to the surface: w fuse up, bottom T down, d liq held near 1 mm |
| **M5b** | **b = 6 mm** (from 4) | flux moves out toward the pool edge (peak flux at lateral position y occurs for b = y√6): w fuse up, d liq down slightly |

Read three numbers from each: **w fuse** against 7.76 mm, **d liq** against
~1 mm, and **bottom T** against 815 K. Then combine the two levers. w fuse and
d liq are flat across 1.06–1.20 s, so 1.2 s runs are long enough. Run the final
calibrated case to 2.0 s. If the machine has 12 or more physical cores, run
M5a and M5b in parallel.

**Confound to control:** a b/c calibrated on a 0.6 mm y-mesh is
mesh-dependent. Penetration is resolved only to one 0.5 mm cell layer. The
pending refinement should verify the final calibration once. A cell-scale
boundary error is about ±0.3 mm per side, so it cannot close a 40% width gap,
but it can shift a calibrated value.
