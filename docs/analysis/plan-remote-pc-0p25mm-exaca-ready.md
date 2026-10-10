# Plan: remote PC, 0.25 mm mesh, ExaCA-ready runs

Decided 2026-09-19. From here on every production run happens on a remote PC,
on a 0.25 mm mesh, with settings chosen so its output feeds ExaCA directly.
This **supersedes the mesh choice in the M4 doc §7**: the M5a/M5b shape probes
run at 0.25 mm, not 0.5 mm. The heat-source calibration depends on the mesh, so
calibrating at 0.5 mm and then switching would have discarded it.

The 0.5 mm runs M1–M4 remain the development history. M2 becomes the coarse
level of a mesh study.

> **Updated 2026-09-20** with the E1 results (§7) and two corrections found while
> writing `claude/E1-ExaCA-Data-Explained.pdf`: current ExaCA needs the
> temperature data at its own cell size, so a resampling step is required (§3),
> and Zhao 2022's "80 to 31 µm" is a **pore** diameter, not a grain size (§3).

---

## 1. Version control: DONE (2026-09-19)

Repository: **github.com/sabbirdipto/waam-mhd-thesis** (**public**, branch
`main`). It is `~/thesis` itself, and the solver is **vendored** inside it as
plain files at `src/lbf-com/`. It is not a separate git repo; the upstream
origin is recorded in `src/lbf-com-PROVENANCE.txt`. One repo holds everything.

**Tag `0p5mm-M1-M4` = commit `dd395e2`** is the exact solver and case recipes
that produced P1114W and M1–M4. Solver changes since the previous push (A2,
`a6fb6ee`): 6 files, +285/−8. They are the Lorentz force (`updateArcLorentz.H`,
new), restart-safe droplet state (`syncDropletState.H`, new), surface convection
in `TEqn.H`, and mushy Marangoni plus `darcy_A`/`darcy_q` in `createFields.H` and
`UEqn.H`. Droplet transfer (`updateDropletSource.H`) was already in A2's tree, so
the commit message's droplet bullet overstates what that commit added.

To build the remote at exactly this state:
`git clone --branch 0p5mm-M1-M4 https://github.com/sabbirdipto/waam-mhd-thesis.git thesis`

**`.gitignore` fix made in that commit.** A `*.orig` rule meant for patch
backups also matched every case's `0.orig/` folder (the initial conditions). In
`.gitignore` the last matching rule wins, so the earlier `!**/0.orig/` was
overridden. It is re-included at the end of the file (36 `0.orig` files had
never been pushed). Also added: `*bak` (covers `.drbak`), `dynamicCode/`,
`nohup.out`, and `cases/**/zz-*/`. The `zz-` prefix now marks scratch cases,
which are never committed.

**Add before the first 0.25 mm case:** most of the 33 MB in that commit is the
cell-centre field `C` (3 MB each in `0.orig/`), which is derived from the mesh
and regenerable. At 864k cells each `C` is about 30 MB, so every new case would
add tens of MB. Append `**/0.orig/C`, `**/0.orig/V`, `**/0/C`, `**/0/V` to the
end of `.gitignore`, where they override the `0.orig` re-include. Make sure
`make_waam_case.py` / `Allrun` regenerate them
(`postProcess -func writeCellCentres -time 0`), because the analysis tools read
geometry from `0/`.

At every launch, write `git rev-parse HEAD` into the case directory as
`SOLVER_COMMIT`, so the question "which build produced this run" always has an
answer. E1 does this.

M2's full tables (`pool_width`, `bead_profile`, `check_droplet`,
`pool_stats_waam`) are banked in `results/` as the reference for checking the
remote build.

**Not yet committed:** `tools/export_exaca.py`, `tools/exaca_coverage.py` and
the E1 recipe (`system/exacaRecorder` + controlDict `functions` block).

## 2. Prove the remote runs the same solver

1. **Gate A** (`hartmann_check.py --sigma 1000 --mu 1`): u(0) 0.0158753,
   Ha 15.810, err 1.3869e-02.
2. **Gate B** (`pool_stats.py Plate2D-mhd-postpatch`):
   1869.0 / 0.0504 / 3.4135 / 0.0552 / 12.7787. Use the same 6-rank
   decomposition. `max|U|` has a 25% noise floor, so judge it within noise.
3. **WAAM smoke test:** run M2's recipe on the **old 0.5 mm mesh** to 0.2 s on
   the remote, then compare its tables against M2's. Compare **integrals
   only**: different hardware and core count put the chaotic flow on a
   different trajectory, so bit-for-bit agreement is not expected (this is a
   standing project rule).

## 3. What "ExaCA-ready" means

| setting | status and why |
|---|---|
| **Recorder**: tm, tl, cr, and G at tl, per cell, every timestep, with multiple rows per cell for remelting | **Built (2026-09-19)** as a coded function object, `system/exacaRecorder`, included from `controlDict`. **No solver change**, so Gates A/B are unaffected. It records exact liquidus crossings at every timestep instead of reconstructing them from 0.02 s snapshots, so there is no aliasing of the 30 Hz droplets (sampling at 50 Hz is below the 60 Hz Nyquist rate). cr is taken over the first 10 K below the liquidus, which also rejects dips that re-melt. Output is MBs, not GBs. First use: E1, §7. **Version 2 proposed** after E1's 19% coverage (§7). |
| **Uniform grid where the metal melts** | The recorder output has to sit on a regular lattice, and the 0.5 mm mesh gives one: its multi-section `simpleGrading` makes a uniform 0.5 mm core over x −12..32 mm, y −6..6 mm and all z, graded only in the outer bands, which is where the cell-volume ratio of 4 comes from. Build the 0.25 mm mesh the same way: uniform core, grading outside. |
| **Resampling onto the ExaCA grid** (added 2026-09-20) | **ExaCA does not interpolate.** From release 1.3.0 (Feb 2024) through the current 2.x, `HeatTransferCellSize` is ignored with "temperature data must be at the same resolution as the CA cell size" (`src/CAinputs.hpp`), and releases ≤ 1.2 could not use a coarser thermal grid together with remelting. Our 0.5/0.25 mm data must therefore be resampled to the CA cell size (5–10 µm) before any grain-resolving run. Read as-is it only works with CellSize = 500 µm (a smoke test of the input chain). The earlier statement here, "its own cell size must divide the thermal cell", described ExaCA ≤ 1.2 without remelting. Tool to write: `resample_exaca.py` (§7). |
| **Runs to 2.0 s** | That is Zhao's full 40 mm track. Take the ExaCA sub-volume from mid-track, where everything has solidified well before `endTime`. |
| **`writeInterval 0.04`, binary, compressed** | Each write is about 10× larger than now. With 0.04 s writes, 0.04 × 30 = 1.2 droplets per interval, so successive writes land at phases 0.2/0.4/0.6/0.8/0.0 and five consecutive writes cover every phase. **Never use 0.05 s** (only 2 phases) **or 0.1 s**: 0.1 × 30 = 3 exactly, so every write lands on the same droplet phase. The analysis tools parse ASCII, so check them before switching to binary. |
| **Al–5Mg interfacial response function** | ExaCA cannot run without it. It has been open since P1 and sits on the ExaCA side (thesis partner). ExaCA ships material files for Inconel 625, SS316 and SS316L (austenite/ferrite), plus AlSi10Mg on the development branch only (power law, A 7.3e-8, B 3.34, freezing range 31 K; a stand-in for smoke tests, not for results). Nothing for Al–Mg. Derive it from the liquidus slope, partition coefficient, solute diffusivity and Gibbs–Thomson coefficient using dendrite-growth theory (KGT), then fit ExaCA's polynomial form (cubic, quadratic, power or exponential). `freezing_range` 91 K to match the flow model. |

**Target microstructure (corrected 2026-09-20).** The "80 to 31 µm" in Zhao
2022 (IJAMT 119:2571) is the **maximum pore diameter** at the bottom of the
wall (their §3.3 and conclusion 4), **not a grain size**, as this doc
previously said. Their quantitative microstructure data is **maximum dendrite
arm spacing**: about 13 → 9 µm at the bottom of the wall with the field, with
more equiaxed crystals; 33 vs 32 µm in the middle, where reheating erased the
difference; the top is "more uniform and smaller". They report no grain sizes,
so a grain-size target must come from their micrographs or another Al–5Mg WAAM
study. Dendrite arm spacing can be compared with our cooling rates directly,
without ExaCA. Their process also differs from Zhao 2021: Panasonic "active"
(CMT-like) mode, 135 A, 18.5 V, 10 mm/s, wire 7.9 m/min, a 60-layer wall on
5A06. Neither ExaCA nor microstructureFoam models dendrite fragmentation, the
main refinement mechanism under electromagnetic stirring. Expect
under-prediction and report it as a finding (roadmap P8).

**Domain size:** 5 µm CA cells across the whole bead would be billions of
cells. A desktop can handle a 2 mm slab at 10 µm (~50 M cells) or a 2 mm cube
at 5 µm (64 M), taken at a steady-state x, mirroring Zhao's section at
X = 25 mm. `export_exaca.py --xmin/--xmax` cuts the slab and prints the cell
count, and the slab is what gets resampled. ExaCA's `XRegion`/`YRegion` inputs
can also cut a region from a larger file. Kokkos on an NVIDIA GPU makes larger
domains practical.

## 4. Shielding gas (Zhao Eq. 28): NOT implemented. Add it to the 0.25 mm case

Checked 2026-09-19. The top patch `inlet` is `pressureInletOutletVelocity`, an
open atmospheric boundary. **No argon jet is imposed**, and `caseParameters`
has no gas entries. The side outlets (x-faces) are also
`pressureInletOutletVelocity`, and `farField` (y-faces) and `substrateBottom`
are `noSlip` walls.

Zhao (Table 2, Table 3, Eq. 28): **20 L/min** of argon enters through the
annulus between the wire (Rw = 0.6 mm) and the nozzle (Rn = 10 mm) with the
laminar annular-Poiseuille profile, at T = 300 K. Vz = 0 outside the nozzle.
The nozzle moves with the torch. His air subdomain starts **7 mm** tall; ours is
5 mm. Evaluated: **mean 1.07 m/s, peak 1.69 m/s at r = 4.2 mm**, and the
profile integrates back to exactly 20 L/min.

**Expected effect on the pool is small on every channel:**

| channel | jet | pool-side comparison |
|---|---|---|
| stagnation pressure ½ρU² | ≤ 2.3 Pa | arc pressure 228 Pa (~1%) |
| surface shear μU/δ | 0.04–0.17 Pa | Marangoni 10–57 Pa (≤ 2%) |
| forced-convection cooling | h ≈ 15 W/m²K (impinging-jet estimate) | Eq. 29 already applies h = 80 everywhere, and the whole surface loss is ~2% of the budget |
| droplets | — | impact velocity 0.9 m/s is Zhao's simulated value, so gas and plasma effects are already in it |

So M1–M4 and E1 are not compromised. The shielding gas is still the one item in
Zhao's Table 3 we have not matched, and the 0.25 mm case is being built from
scratch anyway. Add it there as **boundary conditions only**:

- `U` on the top patch: `codedFixedValue` with the Eq. 28 profile, centred on
  the torch (`arc_x0_m` + `travel_speed_m_s`·t, y = 0), zero outside Rn
- `p_rgh` on the top patch: `fixedFluxPressure`
- `T` on the top patch: 300 K
- `alpha.metal` on the top patch: 0
- raise the domain top to 7 mm above the plate

Also noted: `pool_stats_waam` reports gas velocities of 3–8 m/s. With no jet
imposed, those cannot be shielding gas. They are most likely **spurious
interface currents**, a known VOF artefact at a metal/gas density ratio of
~1600 with σ = 0.85 N/m. They carry negligible momentum, but they are worth
characterising once on the 0.25 mm mesh.

## 5. Budget

- Mesh: 50 × 30 × 9 mm at 0.25 mm is 200 × 120 × 36 = **864,000 cells**,
  9.8× the current 88,200 (fewer with the outer bands graded).
- Timestep: at 0.5 mm the Courant limit (`maxCo 0.1`) never bound. Halving
  the cell size will probably make it bind, so budget for Δt of 3–5e-6.
- Cost: **20–35× per simulated second.** M4 used about 40 core-hours per
  1.2 s. At 0.25 mm that becomes ~800–1300 core-hours per 1.2 s, and
  **~1300–2200 core-hours per 2.0 s**. Zhao reported 1600 core-hours
  (100 h × 16 cores) for 2 s, which is a good sanity check.
- On 32 cores: roughly 1–2 days per 1.2 s probe and 2–3 days per 2.0 s
  production run. **The 13-case magnetic matrix is about a month of
  continuous machine time** on one such box.
- Disk: about 300 MB per write before compression, so ~15 GB per 2.0 s run at
  0.04 s writes, and double that once reconstructed.
- ExaCA input after resampling (new): a 2 mm slab at 10 µm is ~8 M cubes with
  data for E1's 19% coverage, ~20 M for a fully covered bead; with ~2.9 melt
  events per cell, write binary `.catemp` (48 B per event) or keep only the
  last melting.

## 6. Run order on the remote

Case family `waam-arc-02-*` = 0.25 mm; `waam-arc-01-*` stays = 0.5 mm.

1. **R0: M2's exact settings at 0.25 mm, 1.2 s**, plus the recorder and the
   shielding gas (§4). It is the control for everything that follows, and
   together with M2/E1 it gives a two-level mesh study (adding 0.35 mm later
   would allow a Richardson error estimate). It is also the first test of
   whether the finer mesh reduces the 8–23% mass loss, which tracked how thin
   the bead is relative to the cells. A droplet goes from ~45 cells to ~370.
2. **R1a: c = 2 mm. R1b: b = 6 mm.** Each changes a single parameter from R0.
   Read w fuse (target 7.76 mm), d liq (target ~1 mm) and bottom T (must stay
   below 815 K). Background on why c and b are the levers is in the M4 doc §2
   and §7. **Also read ExaCA coverage** (new, §7): a correctly calibrated pool
   should score close to 100%.
3. Combine the two into the **final baseline** and run it to 2.0 s. Then do
   the ExaCA export, resampling and G–R map.
4. **Magnetic matrix.**

## 7. Local ExaCA run E1 (while the server is offline)

`waam-arc-01-E1-exaca`, launched 2026-09-19 on the local PC: **M2 exactly**
(same 0.5 mm mesh and physics) run to **2.0 s** with the recorder. Only
`endTime` and the `functions` block differ from M2. The recorder file was
verified by checksum before launch (180 lines, md5
`c9575b2c6540a69e07b354c2ba9262c8`). It compiled cleanly, started at
t = 1e-05, and the arc is renormalised to 1113.7 W. The run takes about 10 h.

Purpose: data for the thesis partner to build and debug the ExaCA side. It
carries M2's known limitation, a fused zone ~3.6 mm wide against Zhao's
~7.8 mm. The calibrated 0.25 mm baseline will replace it in the same format.
Export with `python3 ~/thesis/tools/export_exaca.py` from the case directory
(md5 `f7f9c810ad1078fd2868b337a9a4fa82`). The package goes to
`results/exaca_waam-arc-01-E1-exaca/`, because `postProcessing/` is
git-ignored.

### Result (2026-09-20)

- **Run:** finished, 200,000 steps, 38,038 s (10.6 h) on 6 ranks, Δt 1e-5
  throughout. Recorder: 3,610 events written, 68 cells still open at 2.0 s,
  **8,064 events dropped** (metal fraction fell below 0.5 mid-event).
- **Export:** 0 malformed, 0 off-lattice, tl ≤ tm 0. 1,231 cells, 916 melted
  more than once (~2.9 each: the torch moves 0.67 mm per droplet, the droplet
  is 2.2 mm across). cr median 2,099 K/s; G median 77.7 K/mm; R = cr/G median
  33 mm/s, above the 20 mm/s travel speed, so events are dominated by droplet
  freezing, and 0.5 mm cells smear G low. Steady window **x 6–26 mm**
  (cr 1,830–2,120 K/s, G 75–85 K/mm). Plate underside 757–758 K at 1.8–2.0 s,
  57 K below the solidus.
- **Coverage 19%** (`tools/exaca_coverage.py`): of 1,234 bead cells above z = 0
  in x 6–26 mm at t = 2.0 s, 234 have an event, plus 406 melted-plate cells.
  The events are the droplet-footprint column (~2 mm wide). Two causes, split
  not yet measured: (b) metal arriving already mushy never opens an event (the
  recorder's melt definition is ExaCA's own: right for the plate, wrong for
  deposited metal); (c) the 8,064 surface drops. Underlying both: M2's pool is
  mostly mush (≤ 68 of 436 pool cells at or above the liquidus at 2.0 s), the
  same heat shortfall at the sides that makes w fuse 3.6 mm instead of 7.8 mm.
  **Coverage is now a calibration metric.**
- **Verdict:** E1 is setup/pipeline data for the ExaCA side, **not
  grain-structure results**. Explained in full in
  `claude/E1-ExaCA-Data-Explained.pdf`.

### Next

1. **`resample_exaca.py`** (new, needed for every dataset): 0.5/0.25 mm events
   → CA grid for a slab, binary `.catemp` output. *Copy* mode (each fine cell
   takes its coarse cell's events; keeps remelting; blocky) for testing;
   *smooth* mode (trilinear tm/tl/cr on the last melting, melt boundary placed
   where the interpolated peak temperature = 906 K) for results. E1 slab
   x 14–16 mm at 10 µm: ~8 M cubes with data, ~23 M events, ~1.1 GB binary.
2. **Recorder v2 (awaiting go-ahead):** alpha hysteresis (open > 0.5, drop
   < 0.1); "arrived mushy" events (tm = tl = arrival, cr over the next 10 K,
   flagged, so the export writes the ExaCA file with and without them; the
   approximation is that ExaCA counts undercooling from zero on arrival); peak
   temperature per event; coverage report in the export. Then **E2** overnight
   with M2 physics.
3. **Partner package** (commands given 2026-09-20): `results/for_partner_E1/`
   holding `exaca_waam-arc-01-E1-exaca.csv`, `exaca_E1_slab_x14-16mm.csv`,
   `gr_waam-arc-01-E1-exaca.csv`, `README.txt` (COVERAGE WARNING, RESOLUTION
   NOTE and SLAB FILE appended) and `E1-ExaCA-Data-Explained.pdf`, zipped to
   `D:\Downloads\for_partner_E1.zip`.

## 8. Open: remote PC specs needed

Physical core count, RAM, GPU (model, CUDA), native Linux or WSL, free disk.
