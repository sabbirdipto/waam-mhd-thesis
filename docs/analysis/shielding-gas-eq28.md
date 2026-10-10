# Shielding gas (Zhao Eq. 28): implemented, and running as G1

**Date:** 2026-09-22. **Supersedes the status line in
`plan-remote-pc-0p25mm-exaca-ready.md` §4** ("NOT implemented"). It is now
implemented as boundary conditions only — no solver change, so Gates A and B
are untouched — and is running on the 0.5 mm mesh as
`waam-arc-01-G1-gas` = **M2 plus this one term**.

This was the last unmatched item in Zhao's Table 3.

---

## 1. The physics, and why nothing in it is tuned

20 L/min of argon leaves the annulus between the wire (Rw = 0.6 mm, from the
1.2 mm wire) and the nozzle (Rn = 10 mm). For steady, fully developed,
axisymmetric laminar flow in an annulus the momentum equation collapses to

```
(mu/r) d/dr( r du/dr ) = dp/dz
```

Integrating twice and applying no-slip at both walls, u(Rw) = u(Rn) = 0:

```
u(r) = K [ Rn^2 - r^2 + (Rn^2 - Rw^2) ln(r/Rn) / ln(Rn/Rw) ]
```

K is **not** a free parameter. Requiring the profile to carry the measured flow
rate, `integral u 2 pi r dr = Q`, fixes

```
K = (2Q/pi) / [ Rn^4 - Rw^4 - (Rn^2 - Rw^2)^2 / ln(Rn/Rw) ]
```

so the only inputs are the two radii and Zhao's Table 2 flow rate. Evaluated:

| quantity | value |
|---|---|
| annulus area | 313.03 mm² |
| mean velocity Q/A | **1.0649 m/s** |
| peak | **1.6933 m/s at r = 4.208 mm** (r_peak = sqrt((Rn²−Rw²)/(2 ln(Rn/Rw)))) |
| profile integral | 20.000000 L/min (numerical check) |
| Re = rho·U·Dh/mu on Dh = 2(Rn−Rw) | **1431** — laminar, so Zhao's assumption is self-consistent |
| mass flow | 0.541 g/s |

The solver computes K itself from the three dictionary keys and printed
**19.999998 L/min, peak 1.693342306 m/s at r = 4.208092647 mm, mean
1.064866367 m/s** at start-up, agreeing with the independent Python evaluation
to 7 significant figures.

## 2. What changed in the case

`tools/patch_shieldgas.py` (md5 `3e6c52b917be20fede5ae408adef8e12`), tested on
a synthetic copy of the case before touching anything: clone, patch,
double-apply guard, revert to byte-identical.

| file | change |
|---|---|
| `0.orig/U`, patch `inlet` | `pressureInletOutletVelocity` -> `codedFixedValue` with the profile above, centred on the moving torch, **Vz = 0 outside the nozzle and inside the wire**, as the paper states |
| `0.orig/p_rgh`, patch `inlet` | `totalPressure` -> `fixedFluxPressure`. Mandatory: once velocity is prescribed on a patch, the pressure equation needs the flux-consistent gradient, not a pressure value |
| `constant/caseParameters` | `gas_Q_m3_s 3.333333e-4`, `gas_Rn_m 0.010`, `gas_Rw_m 0.0006`. The BC reads `arc_x0_m`, `arc_y0_m`, `travel_speed_m_s` from the same file, so nozzle and arc cannot disagree |
| `system/controlDict` | `gasFlux`: a `surfaceFieldValue` summing `phi` over `inlet` at every write, so the delivered flow is recorded, not assumed |

`0.orig/T` and `0.orig/alpha.metal` were **not** touched: `inletOutlet` with
300 K and 0 already acts as a Dirichlet condition wherever gas enters and
zero-gradient where it does not.

**Decisions taken** (2026-09-22): faithful reading — the whole top patch is a
velocity inlet with Vz = 0 outside the nozzle, rather than a hybrid that leaves
the rest of the ceiling open; domain top kept at 5 mm rather than Zhao's 7 mm,
so G1's mesh is identical to M2's and the comparison stays one-variable. The
7 mm top belongs with the 0.25 mm case, which is built from scratch anyway.

**Consequence of the faithful reading:** all 0.541 g/s must leave through the
two x-face outlets, 300 mm² of gas cross-section, about 1.11 m/s mean exit.
The y-faces are no-slip walls and the ceiling is now closed except for the
annulus.

## 3. Nozzle clipping — a boundary effect to state, not to hide

The nozzle is 20 mm across and the torch travels from x = −10 mm; the domain
runs x = −15..35 mm. So part of the annulus lies outside the domain until the
torch reaches x = −5 mm, and again once it passes x = +25 mm:

| t (s) | torch x (mm) | delivered | % of 20 L/min |
|---|---|---|---|
| 0.00 | −10 | 17.05 L/min | 85.2 |
| 0.02 | −9.6 | 17.52 | 87.6 |
| 0.10 | −8.0 | 19.03 | 95.2 |
| 0.20 | −6.0 | 19.93 | 99.6 |
| **0.25 – 1.75** | **−5 → 25** | **20.000 L/min (−3.333e-4 m³/s)** | **100** |
| 1.90 | 28 | 19.03 | 95.2 |
| 2.00 | 30 | 17.05 | 85.2 |

G1 runs to 1.2 s, so only its start transient is affected — the region already
excluded from every comparison. The 2.0 s production runs lose up to 15% of the
jet in the last 0.25 s; the `gasFlux` monitor records it per write.

## 4. Prediction, on the record before the run finished

From the force ledger: jet stagnation pressure ½rho·u_peak² = **2.33 Pa**
against 228 Pa of arc pressure (1%); jet surface shear 0.04–0.17 Pa against
10–57 Pa of Marangoni (≤2%); forced convection h ≈ 15 W/m²K against the 80
already applied over the whole free surface. Argon density in
`transportProperties` is 1.6228 kg/m³, argon at 300 K.

So `w fuse` and `d liq` are expected to move by **less than the write-to-write
scatter** (±0.07 mm on M2's w fuse). The visible effect should be in the argon
field: the 3–8 m/s spurious VOF interface currents now compete with a real
~1 m/s downflow, and the bead surface may cool slightly faster.

If that is what G1 shows, the thesis sentence changes from "expected to be
negligible" to "measured and negligible", with Zhao's Table 3 completely
matched.

## 5. What to read when G1 lands

Compare against M2 over **1.06–1.20 s**, the window used for M2/M3/M4:

| tool | numbers | M2 reference |
|---|---|---|
| `pool_width.py` | w fuse, w liq | 3.560 mm (3.504–3.645), w liq 0 at 5 of 8 writes |
| `pool_stats_waam.py` | d liq, d sol, bottom T, pool cells, E kept | 0.88 mm, 1.75 mm, 746–752 K, 44–51 mm³, ~81% |
| `bead_profile.py` | A, h, w toe, h/w | 9.545 mm², 2.679 mm, 5.552 mm, 0.4825 (settled stations) |
| `check_droplet.py` | deposited mass, loss | 8–23% loss band |
| log | `gasFlux` sum(phi) per write | −3.333e-4 m³/s from t = 0.25 s |

Any difference smaller than the twin-run noise floor must be reported as
"below resolution", not as "no effect" — the standing rule from the
chaotic-sensitivity protocol.

## 6. Two gotchas found while doing this

1. **Every case's `Allrun` is stale.** It copies `initial/` into `0/`, but the
   cases ship their fields in `0.orig/` and no `initial/` exists. M2 and E1 were
   both started by hand, which is why nobody noticed. G1's copy was fixed with
   `sed -i 's/\binitial\b/0.orig/g' Allrun`. **`tools/make_waam_case.py`
   (md5 `4ade4cc742c57c005c549b5a41447fd6`) still writes the stale version**, so
   the 0.25 mm case it generates will inherit it. Fix it there once.
2. **`0.orig` carries `C`, `Cx`, `Cy`, `Cz` and `V`.** The `Allrun` comment warns
   that `C` must not be in `0/` before `decomposePar`, but in practice E1 and G1
   both decomposed cleanly with it present. Worth knowing before anyone "fixes"
   it.

## 7. Launch recipe that works (both M2/E1 and G1 used it)

```bash
cd ~/thesis/cases/00-tutorials/<case>
of2506
rm -rf 0 processor* && cp -r 0.orig 0
blockMesh    > log.blockMesh    2>&1
setFields    > log.setFields    2>&1
decomposePar > log.decomposePar 2>&1
git -C ~/thesis rev-parse HEAD > SOLVER_COMMIT
nohup mpirun -np 6 laserbeamFoam -parallel > log.laserbeamFoam.<tag> 2>&1 &
```

G1 start-up, 2026-09-22: BC compiled to `dynamicCode/shieldGasEq28`, banner as
in §1, `Arc: Goldak double ellipsoid` at 1113.7 W, no errors, 756 steps in
120 s (~0.16 s/step, so ~5.5 h for 1.2 s).
