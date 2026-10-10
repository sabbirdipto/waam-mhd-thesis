# Marangoni mushy-zone gate — result, and what it rules out

**Updated:** 2026-09-17, after the M1 run completed to 1.16 s.
**Status:** M1-mushy extended and running to 2.0 s (PID 132865, ~7–8 h). Restart verified resumed
from 1.160007678, not from zero.

---

## THE RESULT

Removing the liquidus gate on the Marangoni force:

| | effect |
|---|---|
| mean \|U\| in the metal | **+34%** (0.2556 → 0.3416 m/s) |
| bead width `w toe` | **−1.30% ± 1.66%** — not significantly different from zero |
| cross-sectional shape `h/(3A/2w)` | **1.009 → 1.009** — identical to three decimals |
| longitudinal uniformity, CV of A(x) | **19.6% → 11.9%**, a 39% reduction |
| fully-liquid width `w liq` | −13.4% |
| pool volume at t = 1.14 | −1.2% |
| `E kept` | −1.8 points (80.0–80.6 vs 81.8–82.6) |

**Pool velocity rose by a third and bead width did not move.**

That is the finding. It is not "the change did nothing" — the flow field responded strongly and
measurably, and the bead became markedly more even along its length. It is that **pool velocity and
bead width are decoupled in this regime**: the melt boundary is set by where conduction puts the
815 K isotherm, and convection inside that boundary redistributes metal and heat without moving it.

Zhao attributes his bead width to Marangoni-driven outward flow (§4.2, and Conclusion 3). That
attribution does not survive a controlled test.

### Why the headline `h/w` number is misleading

`compare_bead` reports `h/w +4.8%` and `gap to Zhao closed −7.5%`, which reads as a worsening. It
isn't. Per-station over 29 points: `A` +2.54%, `h` +4.35% ± 1.74%, `w` −1.30% ± 1.66%. For a
parabolic cap, Δh ≈ ΔA − Δw = +3.84% expected against +4.35% observed — residual 0.5%, inside the
error bar. Height tracks mass station for station (at x = 3.25, A is +21.8% and h is +21.8%).

**Do not quote `h/w` as a worsening without this decomposition.** The mass-normalised shape
parameter is the honest number and it did not move at all.

### Validity of the comparison

- **Between-run mass control: PASSES.** Mean `A` +0.7%; `deposited` agreed to 0.061% at every early
  write; `drops` identical to two decimals.
- **The switch is provably active.** `T metal` identical at t = 0.02 (no melt yet) and at the
  droplet-placement instants (imposed 1506 K), different at every other write from t = 0.04 — the
  first write with any cell above solidus. Both runs `numberOfSubdomains 6`, so this is not
  summation order.
- **Steady-state validity: FAILS** — see the duration problem below.

---

## The duration problem (why the extension is running)

`A(x)` should be flat at **8.482 mm²** — wire feed over travel speed, conserved. Over the 29-station
developed window it arches:

```
x:  -9.25   -6.25   -2.25    1.25    4.75
A:   7.13    9.66   11.10    9.93    8.93     (gate OFF)
     6.41    9.34   12.55    9.57    7.80     (gate ON)
```

Both runs, same shape. The bead is still redistributing — metal has flowed backward out of the
near-arc region and piled up behind it, which is the mechanism Zhao describes ("the metal flowing
backward promotes the continuous increase in the height of deposition layer"), caught mid-flight.

**Zhao runs SPSL to 2.0 s. We ran to 1.16 s.** His X = 25 mm station is reached at t ≈ 1.0 s and has
**1.0 s to settle**; our window was laid between t = 0.145 and 0.845 s, so most of it had 0.3–0.8 s.

And `w fuse` has never stopped growing: block means 2.07 → 2.76 → 3.08 → 3.23 mm, reaching ~3.73 mm
by 1.156 s. **Every width number so far is from a still-expanding pool compared against Zhao's
settled one.** If it grows another 15–20% by 2.0 s the deficit shrinks with no parameter change —
which is why the extension runs before any power decision.

The A/B conclusion above survives this, because both runs carry the identical transient. The
absolute comparison to Zhao does not.

---

## Thermal state at 1114 W (measured, t = 0.98–1.16)

**No melt-through.** `d sol` = **1.75 mm** on every row of both runs (44% of a 4 mm plate),
unmoved over 0.18 s. `bottom` at 748 K with **67 K of margin** to the solidus, and its rate has
collapsed from ~1015 K/s early (331 → 656 K over t = 0.02–0.34) to **38 K/s** — a local quasi-steady
state under a moving source.

**But the coupon is a calorimeter, and it is now measured.** `far met` — the physical plate edge,
15 mm out — rises dead linear:

```
381.0  383.6  386.1  388.6  391.1  393.6  396.1  398.5  400.9  403.3
increments:  2.6  2.5  2.5  2.5  2.5  2.5  2.4  2.4  2.4     →  124 K/s, no saturation
```

An energy balance on 15.9 J/K of aluminium gives 93–113 K/s depending on how the ~19% leaving with
the argon plume is counted. Same order — confirmed.

After the domain trim, `farField` **is** the real coupon edge, so there is no larger body to conduct
into and the only sinks are h = 80 on the bottom and edges. That is faithful to Zhao's Table 3,
which means **his plate does this too** — his X = 25 mm measurement sits on a similarly preheated
plate. Preheat is a shared feature, not something separating our result from his.

---

## Errata — Zhao et al. (2021), consolidated

1. **The 8.2 mm bead width is not in the paper.** Only printed widths are **7.76 mm (layer 1)** and
   4.52 mm (layer 4), §4.3 p.1585, from Fig. 16, both simulated, both SPFL. Table 5 tabulates height
   and *length* only. 8.2 mm is our parabolic inversion of Table 5: 3(8.482)/2(1.55) = 8.21 mm —
   which is why `compare_bead` reports his parabolic ratio as 0.999, by construction.
   **Quote 7.76 mm; present 8.21 as a consistency check.**

2. **Substrate thickness: SETTLED at 4 mm.** §2 says 50×30×**4**; §3 says the metal subdomain is
   50×15×**5** (half-domain). The figures show 4 mm, so §3 is typographical. Consequences:
   our 6,000 mm³ coupon is correct and needs no rebuild; **`L` = 4 mm in Eqs. 22–24**, which
   unblocks A3; and Zhao has **no extra thermal margin** over us, so a power level that cooks our
   plate would have cooked his.

3. **Eq. 27 lacks the z-restriction its own normalisation requires.** With the standard
   6√3/(abcπ√π) prefactor and f_f + f_r = 2, the integral equals ηUI only when restricted to each
   lobe's x half-range **and** to z < z₀. Verified: full-space single-lobe integral is 2f·ηUI;
   halving twice and summing gives ηUI. Eq. 27 as printed states only the x branches. Applied
   literally it delivers **2ηUI = 3024 W**, half of it into argon (ρ 1.6, c_p 520) which then
   conducts back into the plate. Unverifiable from the paper — he may have gated on metal fraction
   and not written it — but the published equation does not integrate to ηUI.
   *(Our `updateArcSource.H` handles this correctly: `if (dz > 0.0) continue;`.)*

4. **Droplet frequency is internally inconsistent.** Eq. 30 with f_d = 30 Hz gives r = 1.1052 mm
   (matches his printed 1.105). His *simulated* droplet is r = 1.114 mm, implying **f_d = 29.30 Hz,
   −2.35%**. Mass is conserved either way. Also: Eq. 30 *is* mass conservation rearranged, so
   agreement with it is a consistency check, not independent validation.

5. **Energy.** η = 0.70 is applied as pool-only per his own text ("excluding heating droplets") but
   sourced from [29], a textbook *total* GMAW efficiency. Pool 1512 W + droplet enthalpy 703–801 W
   = **2215–2313 W against U·I = 2160 W, i.e. 103–107%**. Subtracting ~24 W of electrode-extension
   I²R (small in aluminium; ~170 W in steel) leaves a true overshoot of **1.4–6.0% of U·I** — sloppy
   rather than absurd, and at the low end inside voltage-measurement error.

6. **Table 5 corroborates the over-delivery.** Simulated height 1.55 mm against experimental 1.80 —
   his model under-predicts crown height by 14%. At fixed mass-conserved area a shorter bead is a
   wider one, which is what excess energy produces. And he never tabulates width.

---

## Droplet momentum — hypothesis raised and closed

Suspected the droplets were placed without momentum coupling. **Wrong.** `updateDropletSource.H`
applies `dropletForce` as a body force in `UEqn.H`, impulse-matched to m_drop·v_impact, over a
window τ = 2R/v = **2.46 ms** at 7.4% duty. Peak force 5.49 mN = **9.7×10⁵ N/m³ = 37× gravity**.
The author anticipated the failure mode explicitly ("a patched velocity would simply be
overwritten").

The earlier "no crater" note was a **sampling artifact**. Weber number

$$We = \rho v^2 D/\sigma = 2650(0.81)(2.21{\times}10^{-3})/0.85 = \mathbf{5.6}$$

so impact dynamic pressure (1073 Pa) is *below* the Laplace pressure (1538 Pa) — the crater is
shallow and springs back. Zhao's own crater appears at t = 0.820 and is gone by 0.832 s: a **12 ms**
feature, which he resolves with 4 ms frames. Our write interval is 20 ms.

**Both mechanisms Zhao names for pool expansion are now accounted for**, and neither widens the
pool. Marangoni stirs it (measured, +34% velocity, 0% width). Droplet impingement is present and
strong in bursts but recovers within ~10 ms.

---

## Model gaps still open

**Lorentz force (Eqs. 22–24) absent.** Five of Eq. 21's six momentum sources are implemented.
Build from the *derived* forms — Eq. 23 as printed has a spurious factor 2 breaking axisymmetry and
the radial/axial expressions are transposed. `L` = 4 mm (item 2 above). Expect it to deepen rather
than widen: it drives flow down and inward.

**Eq. 29 is half-implemented in each location:**

| | convection h(T−Tₐ) | radiation εσ(T⁴−Tₐ⁴) |
|---|---|---|
| free surface (`TEqn.H`) | compiled as `hConvSurf`, **inert** — `h_conv_Wm2K` unset | ✓ |
| walls (`0/T` BC) | ✓ h = 80 on `substrateBottom`, `farField` | **absent** |

Free-surface convection ≈ 8 W; wall radiation ≈ 2–3 W at measured temperatures (5190 W/m² under the
arc, 311 W/m² at the edge). **Together ~0.6% of 1836 W.** Both make the model *retain* heat, so
adding them shrinks the pool — the width deficit gets marginally worse, not better.

**h = 80 on a clamped coupon is unphysical.** Contact conductance to a metal bed is 500–5000 W/m²K.
It is Zhao's value, so keep it for the replication and record it as inherited; flag as a sensitivity
candidate if the thermal branch turns out to matter.

**Deviations to document:** droplets placed rather than grown/detached (so `0/T` inlet is plain
300 K, not Zhao's Table 3 T = 1506 K inside the wire radius); η_total = 0.85 vs his 0.70; mesh
0.5 mm vs his 0.25 mm; full-width domain (no symmetry patch — **required**, a magnetic field breaks
his XZ symmetry); and the arc source anchored at z₀ = 0 means **the deposited bead receives no arc
heat**, only droplet enthalpy and conduction (faithful to Zhao, questionable physics, becomes a real
problem multi-layer).

---

## The power question

No physically admissible power reproduces Zhao's bead:

| pool W | η_total implied | admissible? | result |
|---|---|---|---|
| 790 | 0.70 | yes | **no pool at all** |
| 1114 | 0.85 — top of Zhu's 0.60–0.85 band | yes | 5.17 mm, −33% |
| 1438 | 1.00 — the energy ceiling | no | untested |
| 1512 | 1.034 | no | 6.5 mm, −16% (old config) |

**Prediction for 1512 W, on the record before running.** I previously predicted melt-through; on the
measured thermal data I now think that is **probably wrong**. Expect: no melt-through; `d sol`
~2.2 mm; `bottom` ~850 K; `far met` ~576 K at 2.0 s; bead near 6.5 mm — closer to Zhao, still short.
*Caveat:* scaling `bottom`'s 345 K local excess linearly gives ~890 K, which would be melt-through,
while the empirical sub-linear penetration scaling (1.23× depth for 1.41× power) says no. Latent-heat
buffering decides. That makes the run informative rather than a foregone conclusion.

Run it **to 2.0 s, not 1.16 s** — stopping early produces another unsettled bead.

---

## Method rules (§7)

1. **Never read a result without first proving the run that produced it happened.** Caught three
   times in one day.

2. **The 30 Hz / 50 Hz beat.** Beat period exactly 0.1 s, 3 droplets per 5 writes. Write *k* lands at
   droplet phase (0.6k mod 1), cycling {0.6, 0.2, 0.8, 0.4, 0.0} — five evenly spaced phases. So
   **window means over 5 writes are phase-complete; single writes never are.** Writes at t = 0.1,
   0.2, 0.3 are triple-poisoned: droplet-impact temperature, a local minimum in `w fuse` against
   ±40% scatter, and a **±1 droplet (16.67%)** ambiguity in deposited volume because the accumulator
   sits at exactly zero. Depends on f_d being exactly 30.000 Hz; recompute as LCM(1/f_d,
   writeInterval) if wire feed or droplet volume changes.

3. **Fields are the asset; summary text is a derivative.** Archive fields before cleaning anything.
   *(Note: a false alarm that the reference had been deleted came from checking one guessed path.
   `waam-arc-01-P1114W` was intact all along — 59 times, 3.1 GB. An empty result from a guessed path
   is not evidence of deletion.)*

4. **One variable per run.** `h_conv_Wm2K 80` was briefly appended to the running case and removed.
   It would not have affected the run in flight (`hConvSurf` is a `const dimensionedScalar` built
   once before the time loop) but it would have been picked up on any restart — and the extension to
   2.0 s *is* a restart, so removing it is what made that safe.

5. **Report mass-normalised shape, not raw ratios.** `h/w` over a window where `A` has shifted
   attributes mass redistribution to shape. `h/(3A/2w)` divides it out.

---

## Next

1. **Extension to 2.0 s (running).** Three questions: does `A(x)` flatten toward 8.482; does
   `w toe` keep growing past 5.17 mm; does `bottom` stay short of 815 K.
2. **Then 1512 W to 2.0 s**, power only — no `h_conv`, no Lorentz. The power comparison needs one
   variable.
3. **Then the model-completion case**: `h_conv_Wm2K 80`, wall radiation, Lorentz from the derived
   forms, mesh convergence 0.5 / 0.35 / 0.25 mm (Zhao's full-domain-equivalent is ~513,000 cells
   against our 88,200 — roughly 1.8× finer linearly).
4. Target: Zhao's own claimed bound, *"maximum dimensional error does not exceed 15%."* Currently
   ~33% out on width.
5. Only then does the magnetic field go on.

**Microstructure note:** the target EMF paper (`s00170021084665.pdf`) builds a 60-layer wall. Not
simulable — Zhao needed 100 h for 2 s of one layer. The link must come from solidification
parameters (G, R) extracted from single-layer runs. Decide how before Phase 3, because it determines
what each run must save.
