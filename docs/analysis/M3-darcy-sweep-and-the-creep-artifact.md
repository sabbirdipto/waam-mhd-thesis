# M3 — the Darcy sweep, and why the width problem is a melting problem

Run: `waam-arc-01-M3-darcyA1e5`, 0 → 1.2 s, single launch, 120000 steps,
21181 s wall (5.9 h), 6 ranks. Finished cleanly.

**One-line change from M2.** `diff -r` over the whole case tree returns exactly
one hunk, in `constant/transportProperties:188`:

```
< darcy_A           1.0e6;     // explicit; swept separately in M3
> darcy_A           1.0e5;     // Fluent default; Zhao gives only the 1e4-1e7 range
```

`darcy_A` is read at runtime — `createFields.H:1059`,
`transportProperties.lookupOrDefault<scalar>("darcy_A", 1.0e6)` — so this
required no rebuild and no gate re-run. M2 and M3 ran the same binary. The
comparison is therefore attributable to that single coefficient with no
confounds.

---

## 1. The geometry result

Station-matched, 25 common stations, −7.25 to 4.75 mm, 0.5 mm steps.

| | M2 (A = 1e6) | M3 (A = 1e5) | Δ | Zhao |
|---|---|---|---|---|
| A (mm²) | 10.444 | 8.993 | −13.9% | 8.482 conserved |
| h crown (mm) | 2.770 | 1.265 | −54.3% | 1.55 sim / 1.80 exp |
| w toe (mm) | 5.942 | 11.953 | **+101.2%** | 7.76 printed / 8.2 inverted |
| h/w | 0.4661 | 0.1058 | −77.3% | 0.1890 |
| vs Zhao h/w | 2.47× | 0.56× | | 1.00× |
| h/(3A/2w) | 1.051 | 1.121 | | 1.000 |

Against the targets individually:

| | M2 | M3 |
|---|---|---|
| A vs 8.4823 | +23.1% | **+6.0%** |
| h vs 1.55 (sim) | +78.7% | −18.4% |
| h vs 1.80 (exp) | +53.9% | −29.7% |
| w vs 8.2 | −27.5% | +45.8% |
| ripple | 34.7–43.3% | 24.8–32.2% |

`compare_bead` reports "gap to Zhao closed +130.0%" — it closed the aspect-ratio
gap and overshot by 30%. In log-error terms the miss shrank from 0.393 to
0.252, so it is a genuinely smaller error, on the opposite side.

`h/(3A/2w)` stays near 1 in both, so both sections remain near-parabolic and
h and w are consistent measures of the same shape rather than artefacts of the
width definition.

**Dropping one coefficient by 10× doubled the bead width and more than halved
its height.** A parameter Zhao never reports is the single largest controller
of bead geometry in this model. That sensitivity is a result in its own right,
and arguably a more defensible one than any value of `darcy_A` tuned to match
his macrograph.

---

## 2. Why M3 is not physical

`pool_width` on M3, all eight writes from 1.06 to 1.20:

```
  w fuse max  3.626 – 3.671 mm      (dead flat across the window)
  w liq  max  0.000 mm              (every single time)
  w melt max  = w fuse exactly      (no undercutting anywhere)
```

Set against the bead:

```
  w fuse =  3.655 mm       w toe = 11.953 mm       ratio 3.27
```

**The bead is 3.27× wider than the fused zone.** 4.15 mm of metal on each side
sits on substrate that never reached the 815 K solidus. Liquid aluminium
cannot wet and spread across cold solid plate 4 mm beyond the fusion boundary —
the contact line would freeze. And `w liq = 0.000` at every write means no
substrate cell at the surface ever crossed the liquidus at all: there is no
molten pool, only a mushy smear.

So the width gain is **not** the bead spreading into a wider puddle. There is
no wider puddle. It is the metal creeping out over unmelted plate because the
Darcy sink,

$$S = -A\frac{(1-f_l)^2}{f_l^3 + q}\,\mathbf{u}$$

no longer rises fast enough at low liquid fraction to arrest it. This is the
solid-creep artifact Zhao himself acknowledges. M2's +18.1% width gain was
noted at the time as *partly* reproducing it; M3 reproduces it outright.

**The kinematic signature confirms it.** Lateral spreading over the last four
writes:

- M2: 5.939 → 5.969 mm, +0.5%, ≈ 0.38 mm/s
- M3: 11.650 → 11.894 mm, +2.1%, ≈ **3.05 mm/s — eight times faster**

A bead 1.2 s old, well behind the pool trailing edge, should have solidified
and stopped moving. M3's is still widening when the run ends.

---

## 3. What this reframes

**The width deficit was never a Darcy problem.** M3's fused zone is 3.655 mm
against Zhao's 8.2 — **−55%** — and it does not move across the window. A 10×
change in `darcy_A` doubled the bead and left the pool untouched. Tuning `A`
does not widen the weld; it lets metal walk out over plate that never melted,
buying a number closer to Zhao while being less like a weld.

The problem is upstream: **the arc is not melting the plate wide enough.**

Also worth recording: `w melt = w fuse` exactly at every time, so the pool shows
**zero undercutting** — precisely the effect the Lorentz force was added to
produce. Consistent with the force inventory, where Lorentz radial
(5.1e4 N/m³) sits two orders of magnitude below Darcy (6.2e6). The Lorentz
implementation is correct and currently irrelevant; it will matter when the
applied field goes in, not before.

---

## 4. The energy budget, settled

`make_waam_case.py` prints the split explicitly, and it closes exactly:

```
U*I 2160 W  ->  droplets   (eta_2 0.334)   722 W as mass at 1506 K
                arc source (eta_1 0.516)  1114 W
                unmodelled losses          324 W
                                          ----
                                          2160 W
```

So the model is already energy-consistent at **η_total = 0.85**, with 15%
booked as losses, and `arc_power_W 1113.7; // eta_1*U*I, NOT eta*U*I` records
that as deliberate.

This brackets Zhao's η = 0.70 cleanly, because the ambiguity is *what η
multiplies*:

| reading | arc source | total vs U·I | status |
|---|---|---|---|
| η on arc source alone (Eq. 27 literal) | **1512 W** | 103.4% | M4, running |
| current split | 1113.7 W | 85% | M2, M3 |
| η as a total | **790 W** | 70% | already run — produced no pool |
| losses set to zero | **1438 W** | 100% | the physical ceiling |

1512 W spends the 324 W of losses and then 74 W that does not exist. It is
worth running anyway — it answers "does Zhao's stated parameter set reproduce
Zhao's reported bead?", which is the question a replication chapter exists to
ask — but it must be labelled a **literal replication of a published
inconsistency**, not a physical model. 1438 W is the honest upper bound and is
only 5% below it, which is nothing against the ~2× pool widening needed.

---

## 5. Open: the mass deficit

`check_droplet` on M3 at t = 1.2 s reports **27.54 nominal droplets against 36
due** at 30 Hz — a **23.5% shortfall**. M2 ran at ~87% of nominal over its
early window (2.58/3 at 0.1 s, 5.22/6 at 0.2 s, 12.51/14.4 at 0.48 s); M3 sits
at ~76%.

The interface-smearing gap widens alongside it — `bead a>.5 / bead aV`:

```
M2  t=0.10   0.942
M2  t=0.48   0.883
M3  t=1.20   0.811
```

The tool's own legend attributes a steady fractional shortfall to cells being
clipped at placement, and says the injector reports it as "placed X of Y
nominal". That line has not yet been pulled from the log. **A width result
built on a bead missing a quarter of its metal needs this explained before it
is quoted.** Note it does not invalidate the M2↔M3 comparison, since both
carry the same defect, but it does bound what either can be compared to Zhao.

---

## 6. Corrections made during this analysis

Recorded because both were mine and both were stated before being checked.

**Claimed M3 restarted twice.** `runstat.sh` showed three `Arc:` renormalisation
lines under "startup gates" and I read repeated lines as repeated launches.
`grep -c 'Arc: Goldak double ellipsoid'` returns **1**. Single launch. The
`Droplet state: fresh start` banner is therefore correct behaviour, not a
failure — and the restart patch remains unexercised by any long run.

**Claimed the mesh had a 5 mm substrate against `substrate_L_m 0.004`.**
Inferred from 88,200 cells / 4,900 per layer = 18 z-layers, assuming isotropic
0.5 mm cells. The cells are not isotropic. `blockMeshDict` vertices give
50 × 30 × 9 mm on 98 × 50 × 18, so **dx 0.510, dy 0.600, dz 0.500 mm**, and the
substrate is z = −0.004 → 0, eight cells, **4 mm**. Confirmed independently:
50 × 30 × 4 = 6000 mm³, exactly the substrate volume `check_droplet` measures.
`substrate_L_m` is right and the Lorentz current-return depth is right.

**Worth keeping from that mistake:** dy = 0.600 mm is the coarsest axis, and it
is the axis bead width is measured along. `w fuse = 3.655 mm` is about **six
cells**. Every width number in this dispute is resolved by six cells in the
coarsest direction. It does not invalidate M2↔M3 — same mesh, 2× effect — but
the pending 0.5/0.35/0.25 mm convergence study is now aimed squarely at the
quantity under argument, and should refine **y first**.

---

## 7. M4 — the run this produced

`waam-arc-01-M4-1512W`, cloned from **M2** (not M3), so `darcy_A` returns to
1.0e6 by inheritance rather than by editing. Full recursive diff against M2
returns one hunk, `arc_power_W 1113.7 → 1512.0`. One variable.

Launch confirmed the target took:

```
Arc: 1462.48578 W in metal below the surface,
     renormalised x1.033856206 to 1512 W (raw 96.72524999% of target)
```

The raw fraction, 96.73%, matches M3's 96.32% at the same stage, as it should —
same geometry, same ellipsoid, 1.36× the power.

**What M4 tests:** whether a 36% increase in arc power widens the fused zone
from 3.655 mm toward Zhao's 8.2 mm. If it does, the width deficit was a heat
input problem all along and `darcy_A` should be left at 1e6. If `w fuse` stays
near 3.6 mm, the problem is in the heat *distribution* — Goldak `b` (4 mm
half-width), the two-branch kappa, or the η₁/η₂ split — not its magnitude, and
no power setting will fix it.

**The number to read first is `w fuse`, not the bead.**
