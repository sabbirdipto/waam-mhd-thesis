# Zhao Fig. 11 measured — the creep-immune target we never had

**Date:** 2026-09-23. Source: Zhao et al. (2021) *Weld World* 65:1571, Fig. 11
(temperature, velocity and stream traces, t = 0.832 s), measured off the
published raster by pixel analysis rather than by eye.

**Why this matters.** `zhao-2021-vs-our-model-differences.md` §5.1 and §6 argued
that `w toe` is the wrong calibration target because Zhao's enthalpy-porosity
creep corrupts it and because the controlling coefficient, `darcy_A`, is
unpublished; the right target is `w fuse`, which is creep-immune. But the paper
never prints `w fuse` — Table 5 gives height and length only. Fig. 11 contains
it as a contour. This document extracts it.

---

## 1. Method, and the four checks that say it worked

The figure is a flat-colour contour plot, so the palette is exactly recoverable:
11 patches, 10 labelled boundaries (300 … 906 K), uniform 67.3 K bands.

**The end patches are half-width (21 and 20 px against ~47 for the interior
ones), and every label sits on a patch *boundary*, not a patch centre.** That
settles the one thing the whole measurement rests on: the red band is
open-ended, so **red = T ≥ 906 K = at or above liquidus**, and orange is
839–906 K. Read the other way — red as the 839–906 K band — every number below
would be wrong.

Axis calibration comes from the tick marks, found outside the frame, not from
the labels (the outer labels are pulled inward to stay inside the box and are
off by up to 7 px). Panels b/c/d: 43.3–44.2 px/mm, ticks evenly spaced to
better than 0.5 px. Panel a: 23.2 px/mm in X and 23.1 in Z.

Checks:

1. **Plate geometry, three ways.** In all three cross-sections the coloured body
   bottoms out at Z = 0.06 mm and the surface reaches the outer plate edges at
   Z = 3.95–4.05. A 4 mm substrate with its underside at Z = 0 — confirming
   errata 2 independently of the text.
2. **Cross-panel agreement.** Panel (a) is a longitudinal section through Y = 0
   and crosses panels b, c, d. Crown heights from (a) at X = 20/18/15 are
   +0.76/+1.28/+1.71 mm; from (b)/(c)/(d) they are 0.88/1.27/1.78. Agreement
   within 0.12 mm across two independently calibrated panels.
3. **Mass.** The integrated deposit cross-section at X = 15 mm is 8.02 mm²
   against the mass-conserved 8.4823 mm² — −5.4%. The parabolic estimate
   (2/3)·w·h on the same station gives 7.78 mm².
4. **Arc kinematics close.** The in-flight droplet in panel (a) sits at
   X = 21.2–22.2, centre 21.7, against the 21.64 mm the case matrix gives.
   Back-projecting at 20 mm/s puts arc-on at X = 5.00 mm. Clean number, not
   fitted.

**Uncertainty.** Not the pixel pitch (±0.023 mm) — the contour banding. A
contour position is resolved to about half a band in temperature, which in
space is of order one of Zhao's 0.25 mm cells. Treat widths as ±0.3 mm and
depths as ±0.2 mm. The 815 K solidus falls *inside* the 771–839 K band, so it
is reported as a bracket with a linear-in-T interpolation, and the bracket is
the hard bound.

---

## 2. The measurements

t = 0.832 s, arc centre X = 21.64 mm, arc-on X = 5.0 mm, travel 20 mm/s, plate
top Z = 4.00 mm. Lengths in mm, areas in mm².

| station | behind arc | h crown | w toe | A bead | w liq | d liq | w sol(815) | d sol(815) |
|---|---|---|---|---|---|---|---|---|
| b  X = 20 mm | 1.64 mm | 0.88 | 6.40 | 3.99 | 3.42 | 1.22 | **5.68** | **1.60** |
| c  X = 18 mm | 3.64 mm | 1.27 | 6.19 | 5.56 | 2.70 | 1.00 | 5.48 | 1.72 |
| d  X = 15 mm | 6.64 mm | 1.78 | 6.54 | 8.02 | 1.79 | 0.30 | 3.97 | 1.15 |

Solidus brackets (the interpolation above sits inside these):

| station | w at 839 K | w at 771 K | d at 839 K | d at 771 K |
|---|---|---|---|---|
| X = 20 | 5.26 | 6.45 | 1.47 | 1.84 |
| X = 18 | 5.00 | 6.37 | 1.48 | 2.15 |
| X = 15 | 3.28 | 5.23 | 0.84 | 1.72 |

Longitudinal section (panel a), melt pool in the Y = 0 plane: length 10.07 mm
overall, 7.1 mm behind the arc centre and 1.5 mm ahead at most depths, deepest
point 1.49 mm below the plate top. The rear boundary sits at X = 14.4–14.7 mm
over the whole depth range −0.6 to +0.8 mm relative to the plate top — a
near-vertical solidification front 1.4 mm tall, transitioning through one
0.25 mm cell of orange.

---

## 3. What the numbers say

### 3.1 The pool is the right depth and the wrong width

Set against our banked values:

| | Zhao Fig. 11 (peak, X = 20) | M2 (`darcy_A` 1e6) | ratio |
|---|---|---|---|
| w fuse, solidus | **5.68** | 3.560 | **1.60×** |
| d fuse, solidus | **1.60** | 1.75 (M1, 1114 W) | 0.91× |

**Depth already matches to 9%. Width is short by a factor of 1.6.** That is a
much sharper diagnosis than "the pool is too small", and it points at a
different knob: total power and the depth semi-axis `c` are doing the right
thing, and the transverse semi-axis `b` is not. Raising power to widen the pool
would overshoot depth, and M1 already measured what that costs — at 1512 W (M4)
the plate melts through.

**Caveat before acting on it.** Our `w fuse` is measured on dy = 0.600 mm cells
and 3.655 mm is about six of them; the 5.68 mm target is nine and a half. R0
(0.25 mm, matching Zhao) is the run that decides whether the 1.6× is physics or
discretisation. **Do not change `b` on the 0.5 mm evidence.**

### 3.2 Zhao's own creep ratio, and what it implies for `darcy_A`

The creep diagnostic is `w toe / w fuse`, and the fusion width it should be
divided by is the *maximum* that station ever reached, not the instantaneous
one. At X = 15 the solidus zone has already contracted from ~5.7 to 3.97 mm as
the substrate cooled back, so the instantaneous 1.65 there is an artefact of
timing. The stations still near the arc give the honest value:

```
Zhao, X = 20      6.40 / 5.68 = 1.13
Zhao, X = 18      6.19 / 5.48 = 1.13
Zhao, settled toe against peak fused width   6.54 / 5.7 = 1.15
```

Against the M3 sweep:

| | `darcy_A` | w toe / w fuse |
|---|---|---|
| **Zhao Fig. 11** | **unpublished** | **1.15** |
| M2 | 1e6 | 1.56 |
| G1 | 1e6 + gas | 2.09 |
| M3 | 1e5 | 3.27 |

One decade of `darcy_A` took (ratio − 1) from 2.27 to 0.56, a factor of 4.
Another decade predicts 0.14, i.e. a ratio of **1.14** — on top of Zhao's 1.15.
**Zhao's `darcy_A` is most consistent with 1e7, the top of the 1e4–1e7 range he
quotes from the Fluent manual.**

This is a prediction, not a result, because the ratio has two inputs and our
denominator is independently 1.6× too small — a narrower fused zone leaves more
cold plate to creep over at the same `darcy_A`. **Testing it costs one run**:
`darcy_A = 1e7` as a single-hunk change to `constant/transportProperties`,
runtime-looked-up exactly as M3 was, same binary, no gate re-run. If the ratio
lands near 1.15 with `w fuse` unmoved, the coefficient is recovered.

### 3.3 The 7.76 mm width target should be retired

Fig. 11 puts Zhao's *simulated* bead toe at 6.19–6.54 mm at t = 0.832 s. §4.3's
7.76 mm comes from Fig. 16, layer 1, after four layers were deposited on top of
it, and §4.3 says in as many words that the width "will increase slightly as
the calculation proceeds and the deposit profile would be flatter". The two
numbers are the same bead 1.2 mm apart, and the difference is the creep he
describes.

Area confirms they are the same bead: (2/3)(6.5)(1.8) = 7.8 mm² early,
(2/3)(7.76)(1.55) = 8.0 mm² late. Both near the conserved 8.48 — creep
redistributes the section, it does not create metal.

Consequences for the comparison table:

| | vs 7.76 (Fig. 16, retired) | vs 6.4 (Fig. 11, single layer) |
|---|---|---|
| M2  w toe 5.552 | −28.5% | −13.3% |
| M4  w toe 6.317 | −18.6% | −1.3% |
| G1  w toe 6.782 | −12.6% | +6.0% |
| M3  w toe 11.953 | +54.0% | +86.8% |

Note what this does *not* rescue: G1 lands within 6% of Zhao's Fig. 11 toe
width while its fused zone is 3.238 mm against his 5.68. It reaches the right
outline by creeping twice as far, which is why `w toe` was disqualified as a
calibration target in the first place.

Also worth recording: Fig. 11's settled crown of 1.78 mm is closer to Zhao's
**experimental** 1.80 mm than to his own tabulated simulated 1.55 mm. Table 5's
1.55 is the creep-flattened late-time value. His model is *better* at 0.832 s
than the number he published.

### 3.4 Melted metal below the bead — answered

Yes, and substantially. Liquidus penetration is 1.22 mm at X = 20 and 1.49 mm
at Y = 0 in the longitudinal section; solidus penetration is 1.60–1.72 mm, or
40–43% of the 4 mm plate. At X = 20 the ≥906 K region is 6.07 mm² in total, of
which 2.77 mm² — **46%** — is remelted substrate below the original plate top.

Our M1 value of `d sol` = 1.75 mm (44% of the plate) is the same number. The
plate is supposed to melt this much.

### 3.5 Zhao's bead section does not settle for ~6 mm behind the arc

| station | behind arc | A bead |
|---|---|---|
| X = 20 | 1.64 mm | 3.99 |
| X = 18 | 3.64 mm | 5.56 |
| X = 15 | 6.64 mm | 8.02 |

**A factor of two between 1.6 and 6.6 mm behind the arc**, converging on the
conserved 8.48 only at the far station. Panel (a) shows the crown height
plateauing at +1.75 mm around X = 12.5–15, the same 6–9 mm behind the arc.

This is the independent confirmation that the `bead_profile.py` window matters
and is not a nitpick — it is worth a factor of 2 in A. It also supplies the
number: **a trailing margin of 6–7 mm is the minimum**, and the ~11 mm implied
by M2's own `x_tail` is conservative but defensible. Our settled M2 window
(x = −0.75 to +2.75, which is 11.25–14.75 mm behind the arc at t = 1.2 s) is
comfortably past it. **G1's reported window starting at x = −9.2 mm is not**,
and this is the second, independent reason to reject it.

---

## 4. Open

- Test `darcy_A = 1e7` (§3.2). Single-hunk change, no rebuild.
- Re-measure `w fuse` on R0 at 0.25 mm before touching `b` (§3.1).
- Rebuild the bead geometry table on a station-matched window once
  `bead_profile.py`'s `developed()` is fixed, using 6–7 mm as the floor for the
  trailing margin (§3.5).
- Fig. 9 panels b4 (t = 1.04 s, liquid) and b5 (t = 1.48 s, solid) would give
  the creep *increment* directly, as a second and independent estimate of
  `darcy_A`. Same method as here; the figure has not been measured.
