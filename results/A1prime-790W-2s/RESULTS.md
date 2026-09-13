# A1′ — 790 W, 0.5 mm, 2.0 s

`cases/00-tutorials/waam-arc-01-0p5mm-790W`, run 0 → 2.0 s in three segments
(0 → 0.1476, restart → 0.6476, restart → 1.999998). 88,200 cells, 6 ranks,
9.94 h wall clock for the final 1.35 s.

Reference case for the ablation study: gravity correct, energy budget correct,
**arc pressure and pool Lorentz force NOT implemented**.

---

## Table 5 comparison, mid-bead plane x = +10 mm

Zhao measures at X = 25 mm on a 0–50 mm plate, which is mid-bead. Here the bead
runs −10 → +30, so the analogue is **x = +10 mm**. Values are the mean over
x = 9–11 mm.

| | Experiment | Zhao sim | **Ours** | vs sim | vs exp |
|---|---|---|---|---|---|
| Height | 1.80 mm | 1.55 mm | **2.594 mm** | +67.4% | +44.1% |
| Deposit length | 43.28 mm | 42.47 mm | **44.41 mm** | **+4.6%** | **+2.6%** |
| Cross-section A | — | 8.48 mm² | 8.939 mm² | +5.4% | — |

**Deposit length is threshold-sensitive.** Zhao used Vernier callipers on a
coupon; here it depends on where the metal is judged to end:

| Threshold | Extent | Length |
|---|---|---|
| h > 0.1 mm | −12.66 → +33.19 | 45.85 mm |
| **h > 0.2 mm** | −12.66 → +31.75 | **44.41 mm** (quoted) |
| h > 0.5 mm | −12.21 → +30.25 | 42.46 mm |

3.4 mm of spread across that range, which is the whole discrepancy. The tail
beyond x = +29.75 (the arc position at the final write) has h = 0.12–0.39 mm,
under one 0.5 mm cell, so it is smeared interface rather than bead. The
h > 0.2 mm criterion is a stated choice, not a measurement.

---

## What is right

| | Evidence |
|---|---|
| **Energy budget η₁ = 0.366** | three independent checks, below |
| **Pool continuity** | last zero at t = 0.2476; non-zero on all 66 writes after |
| **No melt-through** | `bottom` plateaus at 678.9 K, 136 K under the 815 K solidus |
| **Symmetry** | transverse centroid +18 µm = 3.5% of one 500 µm cell |
| **Mass conservation** | window-mean A 8.223 mm² vs 8.48 expected, −3.0% |
| **Deposit length** | +4.6% vs Zhao sim |

### The energy budget, confirmed three ways

1. **t = 0.02, before the first droplet (0.0333).** Only the arc is heating. The
   `bottom` rise ratio against the 1512 W case is **0.520** against an arc power
   ratio of 790/1512 = **0.522**.
2. **t = 0.65, droplets included.** Ratio **0.670** against a total power ratio
   of 1512/2234 = **0.677**.
3. **`far met`, the plate edge, is the bulk.** 331.2 K (t = 0.6476) → 446.4 K
   (t = 1.9876) = **86 K/s**, against a predicted 1512 W into 16.8 g at
   c_p 1050 = **86 K/s**. Exact.

### No melt-through — and why it plateaus

| Window | `bottom` rate |
|---|---|
| 0.65 → 0.85 | 41.5 K/s |
| 1.05 → 1.25 | 40.5 K/s |
| 1.25 → 1.45 | 13.5 K/s |
| 1.45 → 1.65 | 9.0 K/s |
| 1.65 → 1.85 | 8.0 K/s |
| 1.85 → 1.99 | **−4.3 K/s** |

`bottom` is the maximum on the bottom face, which tracks the moving arc. Once
the source reaches quasi-steady state that local peak saturates: heat in
balances conduction into a widening heated volume. The **bulk** keeps rising —
that is `far met` at 86 K/s — but the local peak does not.

Peak 678.9 K, margin **136 K** to the solidus.

---

## What is wrong: a width deficit

| | Zhao | Ours | |
|---|---|---|---|
| A | 8.480 mm² | 8.223 | −3.0% |
| h crown | 1.550 mm | 2.360 | +52.3% |
| w toe | 8.200 mm | 5.932 | **−27.7%** |
| parabolic h/(3A/2w) | 0.999 | 1.135 | |

Zhao's bead is parabolic to 0.1%. Ours carries the right metal but distributes
it 13.5% more peaked than parabolic and 28% too narrow.

If our bead were parabolic at Zhao's h = 1.55 mm it would need
w = 3A/(2h) = **8.13 mm**, against Zhao's 8.20. So the height error is a
consequence of the width error, not independent of it.

The parabolic excess falls as the bead develops: **22.1% (t = 0.65) → 15.4%
(t = 1.51) → 13.5% (t = 1.99)**.

### What is missing, quantified

Laplace pressure holding the crown up, γκ = 2γ/R with R = 2.40 mm: **709 Pa**.

| Term | Pressure | Share |
|---|---|---|
| Hydrostatic ρgh | 72.8 Pa | 10.3% |
| Marangoni | 15.5 Pa | 2.2% |
| **Arc pressure — NOT implemented** | 227.8 Pa | 32.1% |
| **Pool Lorentz — NOT implemented** | 435.0 Pa | 61.4% |

**663 Pa of 736 Pa is absent.** That is the whole remaining discrepancy.

---

## Caveats

**`ripple` 49.4% is not droplet periodicity.** Detrending does not reduce it
(50.1%), because it is not a linear trend: A(x) climbs from 5.6 at x = −7.75 to
9.66 at x = +7.75 then falls to 8.4 — a long-wavelength hump over ~30 mm. Mean
adjacent-station change is **4.0%**, and droplets land every 0.665 mm. The
metric measures bead-scale non-uniformity, and the tool's explanation
("droplet periodicity … humping, not flow") is misleading.

The hump itself is physical: the bead is thickest around x = +5 to +9 and
thinner at both ends, which is heat accumulation changing how far metal spreads
as the plate warms — the same effect Zhao describes for his Fig. 8.

**A depends on where the window ends.** The window-mean A runs 8.062 (t = 0.65,
4 stations) → 8.399 (t = 1.51, 46) → 8.223 (t = 1.99, 66). It fell as the window
extended toward the arc, where metal is newer and still spreading. The
single-station value at x = +10 mm is 8.939 mm². Quote the plane, not the mean.

**`d sol` 0.75 mm is resolved; `d liq` is not.** `d sol` sits at 0.25 mm (one
cell, the mesh floor) until t ≈ 1.0, then 0.75 mm (three cells) through
t = 1.85. `d liq` stays 0.00–0.25 mm, so there is little fully-liquid substrate.
Zhao's Fig. 9 colour estimate was 1.0–1.3 mm, low confidence. Same order.

**`E kept` ~79–80%** because the numerator omits latent heat, which is 22% of
the droplet enthalpy. Consistent with conservation.

**The last write is 1.9875696396, not 2.0.** `writeControl runTime` counts from
`startTime`, and this run restarted from 0.6475696396. The arc is at
x = +29.75 mm there, 0.25 mm short of the 40 mm bead. There is no arc-off in
the model — `endTime` **is** the arc-off, so the pool at the leading end is
still liquid, as it is in Zhao's own Fig. 9 (a section takes ~0.73 s to solidify
after the arc passes).

---

## Files

| | |
|---|---|
| `stats.txt` | `pool_stats_waam.py`, all 67 writes |
| `bead-summary.txt` | `bead_profile.py`, developed-window history |
| `bead-final.txt` | per-x profile and transverse section at t = 1.9875696396 |
| `stats-to0.66.txt` | the earlier 0 → 0.66 s table, kept for the restart record |

Fields are not in git (gitignored). Rebuild from the dictionaries plus `Allrun`
in ~10 h if they are needed.
