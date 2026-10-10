# Phase 2 — Materials and scaling (complete)

Deliverables: `alMg5Properties` (OpenFOAM dictionary), `alMg5_properties.csv` (property tables 300–2000 K), `make_props.py`, `scaling.py`, `scaling_report.txt`.

## The headline: laserbeamFoam's MTHD is in the wrong physical regime

Two different couplings share the name "Lorentz force", and the distinction decides the thesis.

| | Driver | Scaling | Effect | Needs |
|---|---|---|---|---|
| **Stirring** | J from the arc (135 A, injected) | linear in B | **drives** flow | a current source |
| **Braking** | J from motion, σ(u×B) | quadratic in B | **damps** flow | nothing |

At the melt pool, |J| = I/(πσ_j²) = **1.07e7 A/m²**. Comparing the two as stresses over the pool scale:

| B (mT) | stirring (Pa) | braking (Pa) | ratio |
|---|---|---|---|
| 10 | 322 | 0.33 | 981× |
| 20 | 645 | 1.31 | **491×** |
| 50 | 1611 | 8.21 | **196×** |
| 100 | 3223 | 32.9 | 98× |

They cross only at **B = 9.8 T**. At every field strength any experiment applies, current-driven stirring is 100–1000× the motional braking.

**MTHD as shipped has only the braking term.** laserbeamFoam is a laser code — there is no welding current, so `J = curl(H)` is purely motional. It is therefore two to three orders of magnitude off the dominant physics of magnetised GMAW, at every field strength.

This upgrades Phase 6 from "the novel technical result" to "the load-bearing step". Without current injection the model cannot show the effect at all — not weakly, not at higher field, not at finer mesh. It is in the wrong regime.

## This also explains Phase 1's null

| Configuration | Ha | N |
|---|---|---|
| Plate2D as run (σ = 1e6) | 0.28 | 7.6e-4 |
| Plate2D at true σ = 3.65e6 | 0.54 | 2.8e-3 |

N ~ 1e-3. With no current injection, Plate2D at 20 mT sits **four orders of magnitude** below the braking threshold. The Phase 1 "no resolved velocity effect" was not a resolution failure — **it is the physically correct answer for that configuration.** Re-running at the true conductivity would not change it; σ buys 3.65× against a 1e4 shortfall.

Gate B proved the machinery works. It could never have shown a magnetic effect on the melt, and no field strength would have made it. Revise the Phase 1 write-up accordingly: the honest sentence is "expected null, correctly reproduced", which is stronger than "unresolved".

## σ_e(T) for Al-5%Mg — the new number

Nobody has measured it. Constructed and cross-validated:

```
SOLID   300-815 K   rho = 1e-8(0.1743 + 7.502e-3 T + 4.028e-6 T^2) + 3.30e-8  Ohm.m
MUSHY   815-906 K   linear in liquid fraction
LIQUID  906-2000 K  rho = 1e-8(8.2019 + 1.9799e-2 T - 2.1227e-6 T^2) + 3.0e-8
```

σ = 1.64e7 (300 K) → 8.15e6 (solidus) → **3.65e6 S/m (liquidus)** → 2.36e6 (2000 K).

- Solid backbone: Desai, James & Ho (1984) JPCRD 13(4):1131, pure-Al recommended values (±2–5%), plus a temperature-independent Matthiessen offset from the measured 29% IACS of 5083/5356. Cross-checked against Brandt & Neuer's measured AlSiCu-19 (matching RT resistivity): **agrees to 1%.**
- Melting jump ×2.19–2.42: Brandt & Neuer (2007) Int J Thermophys 28(5):1429 measured this for 11 Al alloys, band ±5%, no compositional trend. **Transferred from Al-Si to Al-Mg — nobody has measured Al-Mg.** State as a limitation.
- Liquid branch: two independent routes (jump-ratio; pure-liquid-plus-solute) agree to 8%. Uncertainty ±15% at the liquidus, ±25% by 2000 K.

**The tutorial's σ = 1e6 S/m understates the real melt by 3.65×** (Ha by 1.91×, N by 3.65×).

## k(T) came free, and validates σ_e(T)

Wiedemann–Franz, k = L₀σT with L₀ = 2.44e-8. Three independent checks:

| Check | WF | Measured | Error |
|---|---|---|---|
| Pure Al, solid at T_m | 213 | 211 | **+1.0%** |
| Pure Al, **liquid** at T_m | 91.7 | 91 | **+0.8%** |
| Al-5Mg at 300 K | 120.2 | 123 (AZoM/eFunda 5182) | **−2.2%** |

k rises 120 → 162 W/mK through the solid (impurity-scattering regime — the opposite of pure Al), then **halves to 81 W/mK on melting**. Caveat: Int J Thermophys 40 (2019) reports that the ideal Lorenz number *overestimates* k for liquid Al-Si, so treat the liquid branch as an upper bound.

Pure Al at 300 K is the one case that misses (+11%) — expected, because inelastic phonon scattering breaks WF in pure metals at intermediate T. In a concentrated alloy elastic impurity scattering dominates and WF holds. The pattern is the physics working, not luck.

## Three problems found in Zhao's Table 4

1. **Surface tension is pure aluminium.** 0.85 N/m and dγ/dT = −1.55e-4 are footnoted to Mills ref **[34], the pure-Al chapter**, while every thermal property comes from ref [33], the LM5/5182 **alloy** chapter. Mg is strongly surface-active in Al, so a true Al-5Mg value is lower. 0.85 coincides with the *oxidized* pure-Al figure (0.865, Garcia-Cordovilla 1986); clean pure Al is 1.10. No measured Al-5Mg value was located. **Marangoni drives the paper's headline surface flow, so this is the largest single uncertainty in the replication.** Carry as a sensitivity case.
2. **Argon density 6e-6 kg/m³** (Sec 3.3, plasma drag). Argon is 1.62 kg/m³ at 300 K and ~0.05 even at 10,000 K. Off by four orders. Probable typo; plasma drag is then understated. Flagged, not silently corrected.
3. **The mushy-zone constant is never stated.** Zhao gives only the range 1e4–1e7. A_Darcy is therefore a free parameter of the replication, not a value inherited from the paper. Set to 1e6 and test.

Also: density is one constant 2650 kg/m³ for both phases. Real liquid Al-5Mg near the liquidus is ~2350, so the melt is modelled ~13% too heavy, biasing every inertial term.

## Scaling results that set the case matrix

- **Ha = 1 at 6.2 mT** — Hartmann layers begin to form
- **N = 1 at 269 mT** — order-unity bulk braking. No experiment reaches it.
- **Rm = 4.1e-3** ≪ 1 — quasi-static approximation safe, induced field negligible
- J×B_ext beats gravity at **2.4 mT** and the arc's own self-field at **13.5 mT**
- Joule heating **3.6 W against 1512 W** of arc = **0.24%** — negligible, now quantified with a real conductivity rather than asserted
- Force ledger as stresses (not CSF volumetrics, which are mesh-dependent): self-field 435 Pa > surface tension 283 > arc pressure 228 > gravity 78 > Marangoni 15.5 > buoyancy 0.5

The earlier "20–30 mT crossover" claim was in the right neighbourhood but imprecise. The defensible statement is **13.5 mT against the arc's own self-field**.

## Open

- Mills (2002) pp. 32–36 is paywalled; c_p and k are constructions, not Mills' tables. Worth one interlibrary-loan request — it would also settle whether Mills lists electrical resistivity for 5182, which would collapse most of the uncertainty above.
- Joshi et al. (2013) Adv Mater Res 665:76, "Electrical Resistivity of Liquid Al-Mg Binary Alloys" — theory only, paywalled, closest primary source to the exact system.
- Metals 10(11):1447 Table 1 — EN AW 5754 with an external field, names its property sources.
- Whether Zhao 2022's field is static or oscillating. If oscillating, a **third** mechanism (induction) applies and the scaling above needs a shielding-parameter analysis.

## Sources

- Desai, James & Ho (1984), JPCRD 13(4):1131 — https://srd.nist.gov/jpcrdreprint/1.555725.pdf
- Brandt & Neuer (2007), Int J Thermophys 28(5):1429 — https://link.springer.com/article/10.1007/s10765-006-0144-0
- Garcia-Cordovilla, Louis & Pamies (1986), J Mater Sci 21:2787 — https://link.springer.com/article/10.1007/BF00551490
- Joshi et al. (2013), Adv Mater Res 665:76 — https://www.scientific.net/AMR.665.76
- Thermal & Electrical Conductivity of Liquid Al–Si Alloys, Int J Thermophys 40 (2019) — https://link.springer.com/article/10.1007/s10765-019-2497-1
- Dinsdale & Quested (2004), J Mater Sci 39 — https://link.springer.com/article/10.1023/B:JMSC.0000048735.50256.96
- NDT aluminium conductivity chart — https://content.ndtsupply.com/media/Conductivity_Al%20Reference%20Chart.pdf
- eFunda AA5182 — https://www.efunda.com/materials/alloys/aluminum/show_aluminum.cfm?ID=AA_5182&show_prop=all
- Zhao et al. (2021), Weld World 65:1571 — https://doi.org/10.1007/s40194-021-01123-1

---

# Addendum — the solver modification (Phase 2b)

Phase 2 uncovered a representational limit in laserbeamFoam that had to be fixed before the property set was usable. Branch `waam-mhd` in the `lbf-com` clone.

## The problem

`transportProperties` reads **one `Polynomial<8>` per phase** for thermal conductivity. Aluminium's kappa **halves on melting** — 162 → 81 W/mK across a 91 K window. No smooth degree-7 polynomial can represent that step: the best full-range fit is 18% wrong in the solid and **84% wrong at the liquidus**.

That error sits at the **solidification front**, where Darcy damping drives velocity to zero so conduction locally dominates — and where ExaCA reads the cooling rate that the entire microstructure prediction depends on. The Pe = 35 argument that protects the melt interior does not protect the front.

## The fix

Two optional entries, `poly_kappa_solid` and `poly_kappa_liquid`, blended by liquid fraction inside the metal term, leaving the existing metal/gas blend by `alpha1` untouched. **Worst error falls from 84% to 2.6%**, and the residual 2.6% is a difference of mixing rule in the mushy window, not a fitting failure.

Liquid fraction is derived from T via Zhao Eq. 8 rather than read from the `epsilon1` field, because `updateKappaCp.H` is first included from `createFields.H:282` — before `epsilon1` (339), `TSolidus` (285) and `TLiquidus` (299) exist. Measured cost of that choice: `max |e_fromT − epsilon1| = 8.2e-04` in metal cells. Negligible.

Also uses the metal's own solidus/liquidus **scalars**, not the `TSolidus`/`TLiquidus` fields — those are blended between metal (815 K) and gas (1.0 K) by `alpha1` and would be wrong inside the metal branch near the free surface.

## Verification

| test | property | result |
|---|---|---|
| Gate A, Hartmann | non-chaotic — **decisive** | 1.3869e-02, bit-identical |
| Gate B, Plate2D | chaotic — **sensitive** | all five pool statistics bit-identical |
| Two-branch, solid cells | 2469 cells | 3.7e-10 |
| Two-branch, liquid cells | 675 cells | 5.4e-10 |
| Two-branch, mushy cells | 56 cells | implied ΔT 0.63 K of a 65 K window |

The fallback path **branches to the original expression verbatim** rather than computing `e·k + (1−e)·k` with identical branches. Those are algebraically equal but differ by ~1 ULP about 5% of the time, and the first attempt showed why that matters: Plate2D amplified it to a **10% velocity difference**. Making the legacy path bit-identical is what lets a regression test distinguish "unchanged" from "barely changed".

## Three methodological findings

**1. The chaotic and non-chaotic cases do different jobs, and you need both.** The chaotic case is *sensitive* — it amplifies any real change to visibility. The non-chaotic case is *decisive* — it can attribute. Neither alone is sufficient: a difference in Plate2D alone cannot be interpreted, and Hartmann alone cannot detect. This is a sharper statement than §05 of the roadmap currently makes.

**2. The noise floor did its job on a real question.** The first regression attempt produced a 10% velocity difference after a code change. Without the measured floor that is indistinguishable from a bug — you either chase it for a day or wave it through. With it, the answer took one line: 9.96% against a floor peaking at 21%, with the non-chaotic case bit-identical alongside.

**3. `kappa` lags `T` by construction, and always has.** `laserbeamFoam.C` calls `updateProps.H` at line 158 and `TEqn.H` at line 178, so written kappa was evaluated from an earlier temperature than written T. Pre-existing behaviour — invisible with a single polynomial because you cannot back out the temperature. It is one *outer iteration*, not one timestep, so it shrinks with `nOuterCorrectors`. Worth knowing before any future test assumes fields in a time directory are simultaneous.

## Also settled

- **`emS`/`emL` are emissivity**, blended by liquid fraction at `TEqn.H:15`. The tutorial's 0.4/0.1 are steel; oxidised aluminium is ~0.3/0.15. Not negligible — radiation overtakes convection at ~1300 K and is 2× it by 2000 K, and this pool is above 906 K throughout.
- **`hartmann1d-fine` is not fine.** Identical mesh to `hartmann1d` (both `half-span 0.4950 + half-cell 0.0050`). Almost certainly the `foamDictionary`/`#eval` trap recorded in Phase 0 — the lesson was written down but the case never fixed. Must be genuinely refined or renamed before it can be cited as a mesh-independence check.
- **Gate A's σ = 1000 is a dimensionless tuning knob**, not a material property — the tutorial pins L and μ_f to 1 so Ha is set directly. Gate A validates the MHD *coupling*; it says nothing about whether 3.65e6 S/m for Al-5Mg is right. Separate claims, keep them separate.

## Tooling added

- `tools/patch_kappa.py` — applies/reverts the change, verifies both files match before touching either, and checks the anchor precedes the include site
- `tools/verify_kappa.py` — unit test of the blend, checking each melt state against the exact prediction
