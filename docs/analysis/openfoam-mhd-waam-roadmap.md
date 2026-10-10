# MHD-WAAM thesis — plan rev. 4 (1 Sep 2026)

Full formatted plan: https://claude.ai/code/artifact/8c61870c-696c-4ccc-b5be-46d395b6096d
Companion PDFs in this project: *Chaotic Sensitivity in Melt-Pool CFD*, *The Case Matrix, Explained*.

**Status: Phases 0 and 1 complete. Phase 2 is next.**

## Phase 1 results (new in rev. 4)

Both gates passed. The magnetic half of the project is de-risked; everything hard from here is arc physics.

| Gate | Result | Number |
|---|---|---|
| A — Hartmann analytical benchmark | PASS, Ha **imposed** not fitted | 1.39% max deviation |
| B — MHD + melting + VOF + Marangoni | PASS, 6,906 steps, no collapse | div(μH) 1e-12 |
| Melt volume, field vs no field | **Resolved null** against measured floor | 0.14% vs 0.4-0.7% floor |
| Melt velocity, field vs no field | **Unresolved** — ordering flips | floor peaks at 21% |
| Cost of switching the field on | multiply every compute estimate | **2.7×** |

Build: ESI v2506 installed alongside OF10, `of10`/`of2506` shell functions. `OpenFoam_com_main` cloned and compiled at commit **7578023c** (2026-07-30); `libMTHD.so` present alongside `laserbeamFoam`, `compressibleLaserbeamFoam`, `laserMeltFoam`, `setSolidFraction`.

**Gate A detail.** `tutorials/.../hartmann1d` reproduces `u(y)/u(0) = [cosh(Ha) − cosh(Ha·y/L)]/[cosh(Ha) − 1]` to 1.39%. Two bugs in my own checker had to be fixed before the number meant anything: (1) Ha was being *fitted*, which slides along the Hartmann family and absorbs a wrong profile — it must be imposed as `B·L·√(σ/μ)` from the sampled H column; (2) `L` was taken from the cell-centre span (0.495) rather than the wall (0.5), which manufactured a 16% error at the first cell. Correct `L = half-span + half-cell`. Core agreement is 1e-4 to 1e-7; the residual is confined to an under-resolved wall layer (~3 cells across L/Ha). U_1 ≈ 1e-5 and U_2 = 0 confirm fully-developed unidirectional flow, which is the precondition for comparing to the analytical solution at all.

**Gate B detail.** `Plate2D` with `electromagnetics true` — the one claim in rev. 3 not supported by evidence. MHD coexists with VOF, melting, latent heat, Darcy damping and Marangoni. Divergence held at 1e-12 across a 1e6 conductivity jump. Timestep fell 6.5× but did not collapse.

**The tutorial's own conductivity floor is 1e-6** (metal 1e6 S/m, gas 1.0) — far more aggressive than the 1e-3 rev. 3 recommended, and it still holds divergence at 1e-12. Useful prior evidence for the S00/S01 sensitivity runs.

## Three refinements to the noise-floor protocol (new in rev. 4)

The Phase 0 protocol survived contact with a real MHD comparison, but not unchanged. Each of these caught a result I had already written down.

1. **One twin pair is not a noise floor.** It is one sample of a distribution. The Phase 1 pair gave 0.6% at one instant and 21% at another. Phase 7 headline numbers need **three to five** perturbed runs and a quoted standard deviation.
2. **Report the floor as a time series, not a number.** Across a single 10 ms run it grew from 0.6% to 21%. Plot floor against signal; if the curves cross anywhere in the window, the effect is not resolved there.
3. **Time-average only over a genuinely quasi-steady window.** A 5-write average produced a clean "15% braking at S/N ≈ 15" that the full series destroyed — the ordering of field and no-field runs *flips* at t = 2, 5 and 8 ms. Averaging over a window that still contains transients manufactures a signal. Check the series, then choose the window, then average.

What survived: **melt volume is unaffected — 0.14% against a 0.4-0.7% floor.** That is a resolved null and a publishable result. **Melt velocity is unresolved at this field strength and mesh** — not "no effect" but "below the resolution of this experiment", which is a different sentence and the only honest one. This is a statement about Plate2D at 20 mT, not about WAAM.

## The headline change since rev. 2

**laserbeamFoam's ESI mainline already contains the magnetic module.** `src/MTHD` — induction equation for H with divergence cleaning, `J = curl(H)`, `JcrossB` already in every solver's UEqn, Joule heating, and σ/μ blended across the VOF interface by `alpha1`. It ships a Hartmann validation tutorial. That was rev. 2's Phase 5 and its highest-risk item — **and as of Phase 1 it is verified running on this machine.**

**But it does not solve the problem.** There is no current injection anywhere in MTHD — laserbeamFoam is laser-based, so the only current is motional, `J = σ(u×B)`.

| | Laser (MTHD as shipped) | GMAW-WAAM |
|---|---|---|
| Current source | motional EMF | welding arc, 135 A |
| \|J\| | ~2e5 A/m² | ~1.1e7 A/m² |
| Lorentz at 50 mT | ~1e4 N/m³ | ~5.4e5 N/m³ |
| vs gravity (2.6e4) | below — negligible | 20× above — dominant |

~50× and it crosses a qualitative threshold. **Injecting the welding current is the core contribution.**

## Toolchain (settled in Phase 1)

| | Keep | Added |
|---|---|---|
| Fork | OpenFOAM.org (Foundation) | **OpenFOAM.com (ESI)** |
| Version | 10, `/opt/openfoam10` | **v2506**, `/usr/lib/openfoam/openfoam2506` |
| Branch | `OpenFoam_Org_main` (no MHD) | **`OpenFoam_com_main`** @ 7578023c |
| Role | fallback | thesis solver |

**Not `DEDFoamv2`** — I recommended it, then checked the other branches. `OpenFoam_com_main` is the repo's default branch, has the **newest MTHD (484 lines vs 439)**, ships the same Hartmann tutorials. DEDFoamv2's only unique feature is a Lagrangian *powder* cloud; you deposit wire.

One OpenFOAM per terminal — `.bashrc` auto-source commented out, `of10` / `of2506` functions instead.

## Physics ledger

- **10 items ship free**: VOF, enthalpy-porosity, latent heat, CSF surface tension, Marangoni, Boussinesq, Lorentz, Joule, σ/μ blending, applied field via H BCs
- **6 are codedFvOption / codedFixedValue snippets** — no solver surgery: Goldak, arc pressure, plasma drag, surface heat loss, moving wire inlet, annular gas inlet
- **1 is genuinely new code**: welding-current injection

The paper's **assumption 1** models wire feeding as hot liquid from a virtual tube — that is a *boundary condition*, not a source term. The hardest-sounding item needs no solver edits.

## The two hard problems

**1. There is no arc anywhere.** Every branch heats with a ray-traced laser. Set laser power to zero, add Goldak (Eq 27) as a codedFvOption with a translating centre, plus arc pressure (Eqs 25-26) and plasma drag. This is the bulk of the replication work.

**2. Injecting 135 A.** Two routes, in order:
- **Route A** — Ampère's law: `H_θ(r) = I_enc(r)/2πr`, `I_enc(r) = I[1 − exp(−r²/2σ_j²)]`. Impose as an azimuthal H boundary profile; the induction equation propagates it inward. **Boundary conditions only, no code.** Try first.
- **Route B** — add a φ sub-solve, `∇·(σ∇φ)=0`, superpose conduction and motional current. Requires editing MTHD. Cadiou et al. (2020) is the template.

**Validation**: set B_ext = 0 and compare computed J×B against the paper's analytic force fields (Eqs 22-24). Quantitative, against published numbers, before any external field is added.

## Nine phases (20 weeks)

- **P0 — COMPLETE.** Toolchain verified, project structured, 31-case matrix with self-auditing generator, chaotic-sensitivity constraint found and documented.
- **P1 — COMPLETE.** ESI v2506 installed, `OpenFoam_com_main` built at 7578023c. Gate A 1.39%. Gate B passed. MHD costs 2.7×. Noise protocol revised. *Still open from P1: HPC access application, GitHub remote (work exists in one place on a virtual disk), ExaCA Al-Mg material file check.*
- **P2 (wk 2) — NEXT.** Zhao Table 4 → dictionary via Mills (2002). **σ_e(T) for Al-5Mg** — the most consequential new number. Unit checks (120 cm/min = 0.02 m/s, 900 cm/min = 0.15 m/s). Force-scaling analysis by hand.
- **P3 (wk 3-5)** Goldak + arc pressure + surface heat loss as codedFvOptions. Build the non-chaotic parallel regression case.
- **P4 (wk 5-7)** Wire feed as a moving `codedFixedValue` inlet. Droplet EM force + plasma drag. Gate: radius within 5% of 1.105 mm (Eq 30), f≈30 Hz, v_max 0.9 m/s.
- **P5 (wk 7-10)** Replication: SPSL at paper fidelity, mesh study, SPTL ×2, SPFL. Gate: within 15% of Table 5. **Safety point — submittable thesis from here.**
- **P6 (wk 10-12)** Current injection. Gate: J×B matches Eqs 22-24 at B=0. *The novel technical result.*
- **P7 (wk 12-14)** 13-case magnetic matrix on HPC. M00 + its noise twin FIRST. Now budgeted at **2.7× the field-off estimates**, and with 3-5 twins per headline number.
- **P8 (wk 14-17)** Microstructure — ExaCA or microstructureFoam, decide on evidence.
- **P9 (wk 17-20)** Verification, sensitivity, writing.

## Timelines — honest

Full plan = **20 weeks working time, 5-6 months with slack**. Fits the 6-8 month runway; does **not** fit 3 months.

Minimum viable (12 weeks): skip resolved droplet (equivalent mass source from the start), SPSL only, 6 magnetic cases, G-R map instead of grain simulation. Phases 1-7 are shared, so **build in MVT order and decide at week 12**.

## Microstructure: two routes, decide at P8

| | ExaCA | microstructureFoam |
|---|---|---|
| How | extract (x,y,z,tm,tl,cr), separate CA code | phase-field grains in the same run, coupled to flow |
| Extra software | ExaCA + Kokkos + a thermal solver | none — a branch of what you have |
| Risk | Al-5Mg response function may not exist | MTHD 45 lines behind mainline |

Neither models **dendrite fragmentation**, the dominant experimental refinement mechanism. Expect under-prediction; treat as the finding and add shear-rate-dependent nucleation as the follow-up.

## Carried forward, unchanged

- **Chaotic sensitivity** — perturbed *serial* run reached 3.54e-1 relative vs parallel 4.00e-1 by t=0.01 s. Parallel is one perturbation among many, not a source of error. Mean temperature agreed to 0.44%: energy conserved, local distribution differs. No cell-by-cell verification; integral/time-averaged quantities only. Enforced mechanically by `generate_cases.py`.
- **Symmetry trap** — half domain valid only for B=0 and B along y. Longitudinal/transverse-X break it and render plausibly wrong.
- **Force scaling** — external field first-order above ~20-30 mT. Joule heating ~3 W vs 1512 W arc: negligible, report it.
- **Domain cuts** — 40→20 mm, decouple droplet after one proof run, 0.35 mm + AMR, no 2D axisymmetric.
- **Case matrix** — 31 rows, 17 enabled, 4 noise twins. R05-SPFL is the paper's best-validated case (4.3%/0.5% error). *The `solver` column still reads `waamMhdFoam`/`additiveFoam` and needs revising now that the toolchain is ESI.*
- **Zhao 2022 (IJAMT 119)** as the magnetic validation target — same lab, same alloy, same rig, no simulation.

## Immediate next action

**Phase 2.** Zhao Table 4 into an OpenFOAM dictionary, temperature-dependent c_p, k, μ chased back to Mills (2002) rather than guessed, and σ_e(T) for solid and liquid Al-5%Mg established — the number that sets both the Lorentz force and the Joule heating, and the one the source paper never needed.
