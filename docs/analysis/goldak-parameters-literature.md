# Is c = 4 mm a realistic Goldak depth for our case? — the literature answer

**Date:** 2026-09-24. Triggered by M1's converged 2.0 s result: right depth
(1.75 mm vs Zhao's 1.60), right solidification timing (8.8 mm behind the arc vs
8.82), and a fused zone **35% too narrow** with essentially **no sustained
liquid** (`w liq` 0.99 mm, zero on alternate writes, against Zhao's 3.42 mm).

**Answer: no.** By Goldak's own stated rule, `c` should equal the measured
penetration. Ours is 1.60–1.75 mm. `c = 4 mm` is **2.3–2.5× too deep**, and it
sits numerically among values fitted to 15–20 mm thick-section steel at 4–27×
our heat input.

---

## 1. Notation — read this before quoting any rule

Three conventions are in circulation and almost every published rule of thumb
is stated without saying which one it uses.

| convention | longitudinal (travel) | transverse (width) | depth |
|---|---|---|---|
| **Goldak 1984 original** | `c_f`, `c_r` | `a` | **`b`** |
| **Zhao 2021 / our Eq. 27** | `a_f`, `a_r` | **`b`** | `c` |

So **`b` means penetration depth in Goldak and half bead width in Zhao.** Every
value below has been mapped into our (Zhao) convention. The Frontiers in
Materials DED review (8:747389) states our convention explicitly and is a
citable anchor for it.

---

## 2. What Goldak actually prescribes

From §III-B of Goldak, Chakravarti & Bibby (1984), *Metall. Trans. B*
15B:299–305, verbatim:

> **"If the cross-section of the molten zone is known from experiment, these
> data may be used to fix the heat source dimensions."**

> **"For example, the width and depth are taken directly from a cross-section
> of the weld."**

> **"Physically these parameters are the radial dimensions of the molten zone
> in front, behind, to the side, and underneath the arc."**

> **"In the absence of better data, the experience of these authors suggests it
> is reasonable to take the distance in front of the heat source equal to
> one-half the weld width and the distance behind the heat source equal to
> twice the width."**

And Goldak & Akhlaghi, *Computational Welding Mechanics* (Springer, 2005):

> **"The size and shape of the 'double ellipsoid' … can be fixed by recognizing
> that the solid-liquid interface is the melting point isotherm."**

> **"In the author's investigation it was found that the most accuracy was
> obtained when the ellipsoid size and shape were equal to that of the weld
> pool."**

Three things follow.

1. **`c` = measured penetration.** "Radial dimensions … underneath the arc" and
   "depth … taken directly from a cross-section" leave no room for reading `c`
   as anything else.
2. **The `a_f = w/2`, `a_r = 2w` rule is a fallback**, explicitly "in the
   absence of better data" — not the primary method. Note it is stated against
   the **full** width, so `a_f = b` and `a_r = 4b`, i.e. **`a_r/a_f = 4`**. A
   later and widely used modification instead sets the front length equal to
   the penetration and the rear to 4× the front — same 4:1 ratio, different
   anchor (MDPI *Metals* 14(11):1213).
3. **Nothing in the 1984 paper relates any parameter to plate thickness.**
   Verified absent. `c ≤ thickness` is not a published rule either.

### The semi-axis is a contour, not a container

§II-D: the semi-axes are defined such that **"the power density falls to
0.05 q(0) at the surface of the ellipsoid"** — which is where the constant 3
comes from, since e⁻³ = 0.0498.

So the Gaussian always deposits energy outside the nominal ellipsoid. Goldak
simply asserts the 5% contour coincides with the fusion boundary. That is the
honest technical basis for the erfc argument below: an inflated `c` does not
"spread the source", it **moves power into metal that never melts**.

Worth keeping straight: the 5% figure is an *intensity* ratio. The fraction of
**total power** beyond the semi-axis is **1.43%**, not 5%. The integral form
does not appear to be quoted anywhere in the literature.

---

## 3. Where c = 4 mm sits in the published record

| c (mm) | reference | material / process | plate | heat input |
|---|---|---|---|---|
| 5.0 | Hammad et al. (via PMC10707621) | EH36 steel WAAM | thick | very high |
| 4.72 | Pyo, Kim & Kim (2020), *Metals* 10(7):885 | SS400 GMAW butt | 15.6 mm | **2880 J/mm** |
| **4.0** | **Zhao et al. (2021)** | **5A05 GMAW-WAAM** | **4 mm** | **108 J/mm** |
| 3.0 | Wang, Zhang & Zhang (2024), *Micromachines* 15(12):1489 | 7075 Al CMT-PADV WAAM | 15 mm | ~165 J/mm |
| 2.96 | Liu et al., ICAS 2020, paper 0774 | 316L arc-fuse AM | 20 mm | 480 J/mm |

**Zhao's `c` is statistically indistinguishable from values fitted to
thick-section steel at up to 27× his own heat input, and it equals his entire
substrate thickness.**

### The one clean numerical confirmation of the rule

Liu et al. (ICAS 2020) is the only source found that reports both the
semi-axes and the measured pool:

```
depth semi-axis 2.96 / measured penetration 3.157 = 0.94
width semi-axis 4.46 / measured half-width 4.758 = 0.94
```

Both axes at 0.94× the measured pool — as clean a confirmation of
"semi-axis = measured dimension" as the literature offers.

### What low-current arc deposition actually penetrates

| source | process | current | penetration |
|---|---|---|---|
| Cambon et al. (2020), *Weld. World* | 316L WAAM-CMT | 124 A | **0.9 mm** measured, 0.85 computed |
| Tongov & Petkov (2023) | steel WAAM, **4 mm plate** | 76 A | **"slightly greater than 1 mm"** |
| Kim et al. (2023), *Appl. Sci.* 13(11):6653 | **5083-O / ER5356** tandem GMAW, 20.5–26 mm/s | pulsed | **0.00–3.0 mm** across the DOE |
| Zhao Fig. 11 (our measurement) | 5A05 GMAW-WAAM | 135 A | **1.60 mm** (solidus) |

Penetration in low-current aluminium arc deposition is an O(1 mm) quantity.
Our own 1.75 mm sits comfortably inside that band. **A depth semi-axis of 4 mm
describes a pool nobody in this process window has ever measured.**

---

## 4. The erfc criterion — and why it appears to be ours

Because the Goldak exponential factorises, the z-marginal of deposited power is
exactly ∝ exp(−3z²/c²), so the fraction deposited deeper than `d` for a
half-ellipsoid is **erfc(√3·d/c)**:

| c | below 1 mm | below 1.6 mm | below 2 mm |
|---|---|---|---|
| **4.0 mm (Zhao, ours)** | **54.0%** | **32.7%** | 22.1% |
| 3.0 mm | 41.4% | 19.1% | 10.2% |
| 2.0 mm | 22.1% | 5.0% | 1.4% |
| **1.6 mm** | 12.6% | **1.4%** | 0.2% |

**At c = 4 mm, a third of the arc power is deposited below the entire measured
pool.** At c = 1.6 mm it is 1.4%.

An independent search of the welding-simulation literature found **no paper
that states this criterion in any quantitative form**, and none that warns
about power deposited below the fusion zone. The qualitative precedent is
Goldak & Akhlaghi's own "the ellipsoid should equal the pool"; the only paper
that worries about misplaced power at all is Mokrov & Simon's MR10
modification, which blocks the rear ellipsoid from heating the plate beneath
the weld reinforcement.

So the erfc form is best framed as **a quantitative restatement of Goldak's own
rule**, not a citation. That makes it a small original contribution rather than
a gap in our reading.

---

## 5. The parameters are coupled — you cannot move depth without moving width

The prefactor is `6√3 f η U I / (a_f b c π√π)`. It contains the **product**
`a·b·c`, so shrinking `c` raises the peak intensity *everywhere*, including at
the surface. There is no parameter that moves depth without moving width.

Published support, assembled from fragments because nobody states it in one
place:

- **Mokrov & Simon (MR10, TU Graz/RWTH):** with standard Goldak the melt depth
  **"increases only slightly"** when the depth parameter is increased. `c` is a
  *weak* lever on penetration.
- **Wang, Zhao, Wu & Zhang (2003)**, an explicit inversion from measured width
  and depth to semi-axes: their lateral parameter depends only on width, but
  their depth parameter cannot be written without the width parameter.
- **Monash / *Additive Manufacturing* (2024)**, "Application of melt pool
  profiles for parameter calibration of Goldak's heat source model": calibrating
  on a transverse cross-section alone gives **12.97% error** in the longitudinal
  semi-axis, falling to 1.46% once a top view is added. *"Using the CS profile
  only cannot uniquely determine b value."*

That last result has a direct consequence for the replication: **Zhao's
`a_f = 3`, `a_r = 6` are essentially unconstrained by the single macrograph he
fitted them to.** A transverse section cannot identify them.

### What this predicts for reducing c

The surface intensity scales as 1/c, so c = 4 → 1.6 raises it **2.5×**, while
the density at 1.6 mm depth *falls* by a factor of 5:

```
q(0,0,z) ~ (1/c) exp(-3 z^2 / c^2)
c = 4.0 :  z=0 -> 0.250     z=1.6 mm -> 0.155
c = 1.6 :  z=0 -> 0.625     z=1.6 mm -> 0.031
```

Hotter and wider at the surface, cooler underneath. Both directions are the
ones we need: `w fuse` up, superheat up (so a *sustained* pool instead of
droplet pulses), penetration down from our +9%.

A crude static-Gaussian estimate — melt half-width `y_m = b√(ln(I₀/I_th)/3)`,
with `I₀/I_th = 1.89` inferred from our current 3.69 mm fused width — gives
**5.76 mm at c = 1.6 mm**, against Zhao's 5.68. That agreement is almost
certainly luckier than the model deserves: it ignores conduction, convection,
latent heat and the moving source. Take the direction, not the number.

The same estimate puts penetration at ~1.15 mm, i.e. it may *overshoot*
downward. If so, c ≈ 2.0–2.5 mm is the compromise. **The run decides; the
estimate only says which way to move.**

---

## 6. Why Zhao's calibration should not have transferred

§3.5: *"The geometric dimensions of heat source, a_f, a_r, b, and c are
obtained by comparing the calculated and experimental morphologies of molten
pool of first deposition layer, which are 3mm, 6mm, 4mm, and 4mm."*

No objective function, no algorithm, no uncertainty, no identifiability check —
visual matching of one macrograph. Measured against the calibration literature
(contour-distance genetic algorithms, multi-objective Pareto searches,
thermocouple-based NSGA-II) that is the weakest reporting tier, and it is worth
saying plainly in the thesis.

More importantly, **errata 3 gives a mechanism for why `c = 4 mm` was needed**.
Eq. 27 as printed carries no z-restriction, so the 6√3 normalisation — which
integrates to ηUI only over the half-space below the surface — delivers up to
2ηUI if applied over all space, with half landing in the bead and droplets
above z = 0. If his code followed the printed equation, **his `c = 4 mm` is
partly compensating for a source dumping half its energy above the plate.**

Our source is correctly half-ellipsoidal (`if (dz > 0.0) continue;`) and
renormalised to exactly 1113.7 W in metal. A calibration that absorbed an
implementation artefact cannot transfer to an implementation without it.

Supporting circumstantial evidence: **no other group reuses Zhao's 3/6/4/4.**
No citing paper was found that adopts those values.

### f_f and f_r — do not cite Goldak for 0.667/1.333

Goldak 1984 gives only `f_f + f_r = 2` and the fitted **`f_f = 0.6, f_r = 1.4`**,
which was a thermal-history fit, not a continuity condition. The relation

```
f_f = 2 a_f/(a_f + a_r),   f_r = 2 a_r/(a_f + a_r)
```

is a later addition. Zhao's 0.667/1.333 is exactly that relation at
`a_r/a_f = 2`, so our banner values are internally consistent with `a_f = 3`,
`a_r = 6`. **If `a_f` or `a_r` ever change, `f_f` and `f_r` must be recomputed.**

---

## 7. Recommended values

Target pool = Zhao Fig. 11 as measured in `zhao-fig11-measured.md`:
`w fuse` 5.68 mm, `d fuse` 1.60 mm.

| parameter | current | Goldak rule | with Liu's 0.94 | recommendation |
|---|---|---|---|---|
| `b` (half-width) | 4.0 | 2.84 | 2.67 | **hold at 4.0 for now** |
| `c` (depth) | 4.0 | **1.60** | 1.50 | **1.6 mm** |
| `a_f` | 3.0 | — | — | hold |
| `a_r` | 6.0 | — | — | hold |
| `f_f` / `f_r` | 0.667 / 1.333 | — | — | hold (consistent with a_r/a_f = 2) |
| power | 1113.7 W | — | — | hold |

**Change `c` alone.** Reasons: it is the largest single lever (2.5× on surface
intensity), it has the clearest primary-source justification, and one variable
per run is the rule that has made every comparison in this project
interpretable. `b = 4` against a target half-width of 2.84 is *not* the binding
constraint — intensity is — and narrowing `b` would confine the pool laterally
while raising intensity, which confounds the test.

**Priority.** This now outranks the `darcy_A = 1e7` test. `darcy_A` tunes how a
bead creeps *after* it freezes; `c` decides whether there is a pool at all, and
M1's converged result says there isn't one.

---

## 8. Open, and worth having in the thesis

Two questions the literature search could not answer, both of which are
genuinely open rather than merely unfound:

1. **The erfc energy-confinement criterion** (§4) is not published in any
   quantitative form.
2. **Should z₀ track the growing bead?** Our source origin sits at the original
   plate top while the bead builds 2.6 mm above it. No paper was found that
   addresses whether the Goldak origin should move up with the deposit, or what
   error fixing it introduces for low-penetration bead-on-plate. Mokrov &
   Simon's reinforcement modification is the nearest precedent. With `c = 4 mm`
   the source is not merely too deep — it is mis-centred by roughly the bead
   height.

Two papers worth pulling through a library, both inaccessible to this search:

- **Yi & Park (2022), *Metals* 12(6):891, Table 9** — Goldak parameters for
  Al 5083 H321, 6061 T6 and 6082 T6 GMAW. The single most on-point source found.
- **Li Peilin & Lu Hao (2011)**, *Trans. China Welding Institution* 32(11):89–91
  — the only paper found with "sensitivity analysis of double-ellipsoid heat
  source parameters" in its title; reports regression formulae from measured
  width and depth to the semi-axes.

## Sources

- Goldak, Chakravarti & Bibby (1984), *Metall. Trans. B* 15B:299–305 — https://link.springer.com/article/10.1007/BF02667333
- Goldak & Akhlaghi, *Computational Welding Mechanics* (2005) — https://link.springer.com/book/10.1007/b101137
- Kik (2024), "Calibration of Heat Source Models", *Metals* 14(11):1213 — https://www.mdpi.com/2075-4701/14/11/1213
- Monash (2024), "Application of melt pool profiles for parameter calibration of Goldak's heat source model", *Additive Manufacturing* — https://www.sciencedirect.com/science/article/pii/S2214860424004251
- Mokrov & Simon, "A fine modification of the double ellipsoid heat source" — https://openlib.tugraz.at/download.php?id=5d440d3492815&location=browse
- Bjelić et al. (2022), *Thermal Science* 26(3):2081 — https://thermalscience.rs/pdfs/papers-2021/TSCI210131181B.pdf
- Kumar & DebRoy (2004), *IJHMT* 47:5793 — https://modeling.matse.psu.edu/research_files/papers/2004IJHMT_Kumar.pdf
- Liu et al., ICAS 2020 paper 0774 — https://www.icas.org/icas_archive/ICAS2020/data/papers/ICAS2020_0774_paper.pdf
- Pyo, Kim & Kim (2020), *Metals* 10(7):885 — https://www.mdpi.com/2075-4701/10/7/885
- Wang, Zhang & Zhang (2024), *Micromachines* 15(12):1489 — https://www.mdpi.com/2072-666X/15/12/1489
- Cambon et al. (2020), *Welding in the World* — https://link.springer.com/article/10.1007/s40194-020-00951-x
- Tongov & Petkov (2023), *Environment. Technology. Resources.* — https://archive-journals.rtu.lv/etr/article/download/4658/4598/5854
- Kim et al. (2023), *Applied Sciences* 13(11):6653 — https://www.mdpi.com/2076-3417/13/11/6653
- Yi & Park (2022), *Metals* 12(6):891 — https://www.mdpi.com/2075-4701/12/6/891 (tables inaccessible)
- Li Peilin & Lu Hao (2011), *Trans. China Welding Institution* 32(11) — http://hjxb.hwi.com.cn/hjxb/article/id/20111123 (full text not retrieved)
- Frontiers in Materials 8:747389 — https://www.frontiersin.org/journals/materials/articles/10.3389/fmats.2021.747389/full
- *Materials* 15(7):2545, DED heat-source study — https://www.mdpi.com/1996-1944/15/7/2545
