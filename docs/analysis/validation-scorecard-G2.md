# Validation scorecard: Zhao experiment vs Zhao simulation vs our G2

**Date:** 2026-09-24. G2 = `waam-arc-01-G2-yopen`, converged at t = 2.0 s,
0.5 mm mesh, arc 1113.7 W, `darcy_A` 1e6, Lorentz on, Eq. 28 shielding gas,
air y-faces open.

**Bar used:** Zhao's own claimed accuracy, §5 Conclusion (1) — *"maximum
dimensional error does not exceed 15%"*, with 4.3% on height and 0.5% on
length. If we cannot meet the bar he sets for himself, we have not replicated
him.

---

## 0. The inputs match by construction

| parameter | Zhao Table 2 | ours | |
|---|---|---|---|
| current | 135 A | 135 A | ✔ |
| voltage | 16 V | 16 V | ✔ |
| travel speed | 120 cm/min = 20 mm/s | 20 mm/s | ✔ |
| wire feed | 900 cm/min = 150 mm/s | 150 mm/s | ✔ |
| wire diameter | 1.2 mm | 1.2 mm | ✔ |
| shielding gas | 20 L/min argon | 19.99998 L/min measured | ✔ |
| substrate | 5A05, 50 × 30 × 4 mm | same | ✔ |
| deposited area (feed/travel) | 8.4823 mm² | 8.4823 mm² | ✔ |

Nothing here is in dispute. Everything below is output.

---

## 1. What Zhao's EXPERIMENT actually gives us

This is the first thing to state, because it constrains what "validation"
can even mean.

| quantity | value | source |
|---|---|---|
| bead height | **1.80 mm** | Table 5 |
| bead length | tabulated, value not extracted | Table 5 |
| bead width | **never published** | errata 8 |
| penetration depth | **never published** | — |
| fusion zone width | **never published** | — |
| pool dimensions | **never published** | — |

**Zhao's experimental validation is two numbers, and one of them is height.**
Table 5 tabulates height and length only; §4.2, §4.3 and §4.4 each attribute
width error to enthalpy-porosity creep, and width is given no error bar
anywhere. So every other row in this scorecard is measured against his
*simulation*, not against reality.

---

## 2. The scorecard

Zhao's simulated column is measured from Fig. 11 (t = 0.832 s) by pixel
analysis — see `zhao-fig11-measured.md`. Widths are stated at **10% of peak
crown height**, matching our `bead_profile` definition, so the two are
comparable.

### Tier 1 — creep-immune. These are the real test.

| quantity | Zhao exp | Zhao sim | **G2** | error | verdict |
|---|---|---|---|---|---|
| `w fuse` — melted substrate width | — | **5.68 mm** | **3.531 mm** | **−37.8%** | ✘ FAIL |
| `d fuse` — melted substrate depth | — | **1.60 mm** | **1.75 mm** | +9.4% | ✔ PASS |
| `w liq` — liquid width | — | **3.42 mm** | **0.098 mm** | **−97.1%** | ✘ FAIL |
| `d liq` — liquid depth | — | **1.22 mm** | **0.00 mm** (7 of 11 writes) | **−100%** | ✘ FAIL |
| pool length | — | **10.07 mm** | not yet measured | — | open |
| solidification front behind arc | — | 8.82 mm (solidus) / 6.95 (liquidus) | **5.8 mm** | −34% / −17% | ✘ FAIL |

These cannot be faked by the enthalpy-porosity creep, because unmelted metal
cannot be spread by a numerical artefact. They are the quantities the
calibration must hit.

### Tier 2 — mass-conserved. Should pass if the injector is right.

| quantity | conserved | Zhao sim | **G2** | error vs conserved | verdict |
|---|---|---|---|---|---|
| bead cross-section `A` | 8.4823 mm² | 8.02 mm² | **7.722 mm²** | −9.0% | ✔ PASS |
| deposited mass loss | 0% | — | ~10% (G1 measure) | — | ✔ PASS |

### Tier 3 — creep-corrupted. Weak evidence, and Zhao says so himself.

| quantity | Zhao exp | Zhao sim | **G2** | error vs sim | verdict |
|---|---|---|---|---|---|
| `h crown` | **1.80 mm** | 1.78 (Fig. 11) / 1.55 (Table 5) | **1.781 mm** | **+0.1%** | ✔ PASS |
| `w toe` | — | 6.18 (Fig. 11) / 7.76 (Fig. 16, late) | **7.472 mm** | **+20.9%** | ✘ FAIL |
| `h/w` aspect ratio | — | 0.288 | **0.238** | −17.4% | ✘ FAIL |
| shape factor `h/(3A/2w)` | — | 0.91 | **1.149** | +26% | ✘ FAIL |
| **creep per side** (`w toe − w fuse`)/2 | — | **0.25 mm** | **1.97 mm** | **7.9×** | ✘ FAIL |

### Tier 4 — thermal gates. No Zhao data; these are sanity, not validation.

| quantity | requirement | **G2** | verdict |
|---|---|---|---|
| `bottom` T | < 815 K, no melt-through | 753.2 K, 62 K margin, rate 7 K/s | ✔ PASS |
| `far met` T | well below solidus | 478.6 K | ✔ PASS |
| `E kept` | stable, no drift | 82.9%, flat | ✔ PASS |
| bead uniformity (ripple) | Zhao's crown is flat to ±5% | **54.7%** | ✘ FAIL |

---

## 3. Tally

```
Tier 1  creep-immune        1 pass,  4 fail,  1 open
Tier 2  mass                2 pass
Tier 3  creep-corrupted     1 pass,  4 fail
Tier 4  thermal gates       3 pass,  1 fail
```

**The single quantity that matches Zhao's experiment is the bead height, to
0.1%. The single quantity that matches his simulation cleanly is the same
height. Everything creep-immune fails.**

Against Zhao's own 15% bar, we are at **−37.8% on fused width, −97% on liquid
width, and 7.9× on creep**.

---

## 4. Why the height match is not the result it looks like

Three reasons, all of which have to be stated together.

1. **Height and width are one measurement, not two.** The deposited area is
   fixed by wire feed over travel speed at 8.4823 mm². A near-parabolic cap has
   `A = (2/3)wh`, so at fixed mass matching the height *forces* the width. They
   are not independent evidence.

2. **They are the pair Zhao says his own method corrupts.** §4.2: *"the width
   of the deposition layer has slightly expanded, which is caused by the
   inherent defect of enthalpy-porosity method"*. §4.3: *"the width of previous
   deposited layers will increase slightly as the calculation proceeds and the
   deposit profile would be flatter"*. §4.4 repeats it. Table 5 omits width
   entirely.

3. **The interiors are different.** Our bead overhangs its own melted zone by
   **1.97 mm per side** against Zhao's 0.25 mm — eight times as much. Over
   roughly a third of its width our deposit sits *on* the plate rather than
   fused *into* it, which in welding has a name: lack of fusion.

---

## 5. What has to change to call it a success

In priority order, with the target and the current value:

| # | criterion | target | now | status |
|---|---|---|---|---|
| 1 | `w liq` — a **sustained** liquid pool, non-zero at every write | 3.42 mm | 0.098 mm | `c` bracket running |
| 2 | `w fuse` | 5.68 mm ± 15% | 3.531 mm | `c` bracket + R0 mesh |
| 3 | creep per side | ≤ 0.5 mm | 1.97 mm | follows from 1 and 2 |
| 4 | `d fuse` held | 1.60 mm ± 15% | 1.75 mm | must not degrade |
| 5 | pool length | 10.07 mm | unmeasured | **add to the tool suite** |
| 6 | bead uniformity | flat crown | 54.7% ripple | expected to follow from 1 |
| 7 | `h`, `w toe`, `h/w` | Fig. 11 values | mixed | **validation set — must not be calibrated on** |

**Rows 1–5 are the calibration set. Rows 6–7 are the validation set and must
stay out of the fit**, or the eventual claim means nothing. That separation is
what makes "our model reproduces Zhao" a statement rather than a tautology.

---

## 6. Open measurement gaps

- **Pool length** is in Zhao's Fig. 11a (10.07 mm, rear edge 6.95 mm behind the
  arc, leading edge 1.66 mm ahead) and we have no equivalent tool output. It is
  a Tier 1 quantity and it should be added.
- **`x_tail` definition** — `bead_profile` calls it the "pool trailing edge",
  and whether it tracks the liquidus or the solidus changes the row 6 error
  from −17% to −34%. Worth pinning down in the source.
- **G2's bead window** is `developed()`'s, which anchors its lower bound to a
  0.05 mm VOF smear detected 4.6 mm upstream of arc-on. The numbers above
  should be re-stated on the agreed convention window
  (`x_arc_on + 5 mm` to `x_tail`) before they go into the thesis.
- **Zhao's Table 5 bead length** was never extracted; it is one of only two
  experimental numbers he publishes and should be read off the paper.
