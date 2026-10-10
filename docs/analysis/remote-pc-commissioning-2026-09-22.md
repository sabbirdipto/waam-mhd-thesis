# Remote PC commissioned: build, gates, benchmark (DESKTOP-3NBGHNH)

**Date:** 2026-09-22. **Supersedes** the environment and timing sections of
`plan-remote-pc-0p25mm-exaca-ready.md`, which predates every number here.

The lab machine is validated end to end: OpenFOAM v2506 built from the
vendored solver on gcc 13, both verification gates reproduced **to every
printed digit**, and the full WAAM physics chain run in parallel for the
first time. The machine is 2.5x faster than the laptop. The 0.25 mm
production run is a ~28 hour job, and disk is the binding constraint.

---

## 1. The machine, and one thing it cannot do

| item | value |
|---|---|
| hostname | DESKTOP-3NBGHNH |
| CPU | 13th Gen Intel Core i7-13700 - **8 P-cores + 8 E-cores**, 16 cores / 24 threads |
| RAM visible to WSL | 15 GB (host likely 32 GB; WSL2 takes half by default) |
| distro | Ubuntu 24.04.3 noble under WSL2 |
| vhdx location | `C:\Users\PC\AppData\Local\wsl\{1b3f834c-5891-491f-8bff-321d7efde427}` |
| root device | `/dev/sdd`, 1007 GB **virtual**, 8.4 GB used after commissioning |
| OpenFOAM | v2506 ESI at `/usr/lib/openfoam/openfoam2506`, `WM_MPLIB SYSTEMOPENMPI`, gcc 13 |
| account | `sabbirdipto` (uid 1001, sudo). `thesis20` and its `qe-7.5` data untouched |
| `/etc/wsl.conf` | `[boot] systemd = true` + `[user] default = sabbirdipto` |

**WSL presents a fake CPU topology.** `lscpu -e` reports a perfectly
uniform **12 cores x 2 threads**, every core with identical cache IDs. The
real die is 8 P + 8 E. Hyper-V has synthesised a homogeneous layout that
does not correspond to the physical chip.

The consequence is permanent and worth writing down: **P-core affinity
cannot be set from inside WSL.** The guest's 24 vCPUs are placed on host
logical processors by the Windows hypervisor, and `--bind-to core` would
bind to a map that does not exist. Any `taskset` or OpenMPI affinity flag
is applied to fiction.

The working answer is to **use 8 ranks and pass no affinity flags at all**.
With only 8 busy threads on a 16-core chip, Windows' Thread Director
promotes them to P-cores by itself. Asking for 16 would force half onto
E-cores; because MPI is bulk-synchronous and `decomposePar` splits the mesh
into equal *cell counts* on the assumption of equal *time*, the E-core
ranks become stragglers and the whole run proceeds at their pace.

## 2. Build

```bash
cd ~/thesis/src/lbf-com && of2506
export WM_NCOMPPROCS=8
./Allwmake 2>&1 | tee ~/thesis/log.build.0922
```

Clone at `6612ce9`. Zero errors on gcc 13 - a generation newer than the
laptop's - so both machines can stay on one commit with no fork.

`$FOAM_USER_LIBBIN`: `libMTHD.so`, `libgeometricVoFLMFOAM.so`,
`libincompressibleTransportModelsLMFOAM.so`,
`libincompressibleTurbulenceModelLMFOAM.so`, `liblaserHeatSource.so`,
`liblaserbeamFoamFunctionObjects.so`, `libmultiphaseVapMixtureThermo.so`.
`$FOAM_USER_APPBIN`: `laserbeamFoam`, `laserMeltFoam`,
`compressibleLaserbeamFoam`, `setSolidFraction`.

`ldd $(which laserbeamFoam)` resolves `libMTHD.so` to the user lib path, not
a system copy. Worth keeping in the recipe: a solver can build perfectly and
still die at startup with `error while loading shared libraries` if
`FOAM_USER_LIBBIN` is missing from the launching shell, which looks like a
solver bug and is purely environmental.

## 3. Gate A - Hartmann. The README was wrong twice

**The README's Gate A line is stale on both counts.** It names
`hartmann1d-postpatch` and gives only the check invocation. Correct on both:

1. **The case is `cases/00-tutorials/hartmann1d`.** The patched
   `hartmann_check.py` - the one carrying the "63% FAIL" comment that
   `patch_hartmann_check.py` inserts - exists only there; a repo-wide grep
   finds no copy under `-postpatch`. The patch's own verification block says
   `cd ~/thesis/cases/00-tutorials/hartmann1d`.
2. **There is a required step nobody wrote down.** `system/sample` is *not*
   a runtime function object; `controlDict` has no `functions` block at all.
   The dict's own header says `Recommended call as: postProcess -time <final
   time> -func sample`. Without it `postProcessing/` is never created and the
   check dies on `FileNotFoundError`. Gate A had **two** remembered command
   lines; the patch fixed one and left the other implicit.

Working recipe:

```bash
cd ~/thesis/cases/00-tutorials/hartmann1d && of2506
rm -rf 0 processor* postProcessing log.*
cp -r initial 0
find 0 -name "*.orig" -exec bash -c 'cp "$0" "${0%.orig}"' {} \;
blockMesh > log.blockMesh 2>&1
laserbeamFoam > log.laserbeamFoam 2>&1          # serial, ~1 s, endTime 10
postProcess -time 10 -func sample > log.postProcess 2>&1
python3 hartmann_check.py postProcessing/sample/10/line_centreProfile_H_U.csv
```

**Result, identical with and without `--sigma 1000 --mu 1`:**

```
samples 100   u(0) 0.0158753   Ha 15.810   max|u/u0 - theory| 1.3869e-02 @ y = -0.4950
VERDICT PASS
```

Every digit matches the banked `results/gateA-edit8.txt`. That is stronger
than expected: the case runs **serial**, so there are no MPI reductions to
reorder the floating-point sums, and OpenFOAM compiles without
`-ffast-math`, so IEEE semantics survive the compiler change. Identical
source plus identical input gave an identical arithmetic sequence.

**Independent cross-check against theory, not just against the bank.**
Hartmann flow flattens Poiseuille into a plug core because fluid crossing
**B** induces **J** = sigma(**u** x **B**), and **J** x **B** opposes the
motion most where the fluid is fastest:

```
u(y)/u(0) = [cosh(Ha) - cosh(Ha y/L)] / [cosh(Ha) - 1],   Ha = B L sqrt(sigma/mu)
```

With B = 0.99993, L = 0.5, sigma = 1000, mu = 1: Ha = 0.99993 x 0.5 x
sqrt(1000) = **15.81**, and u/u0 at y = -0.4950 evaluates to **0.1462**,
against the banked theory column's 0.146238.

**The residual is physically coherent, not noise.** The Hartmann layer is
delta = L/Ha = 0.5/15.81 = 0.0316 thick and the cell size is 0.01, so only
~3 cells span it. That is exactly why the 1.39% maximum error sits in the
wall-adjacent cell while the core is accurate to 1e-6. The gate is not
passing by luck; the error is bounded under-resolution of a layer whose
thickness theory predicts.

## 4. Gate B - Plate2D

`diff -rq Plate2D-mhd Plate2D-mhd-postpatch` returns **exit 0, zero lines**:
byte-identical trees, one case under two names. The README/results-file
naming discrepancy is cosmetic. Run `Plate2D-mhd` so `pool_stats.py`'s
printed header matches `results/gateB-edit8.txt` verbatim.

`pool_stats.py` lives in `cases/00-tutorials/`, not `tools/`.

```
melt volume 1869.0 cells | mean|U| metal 0.0504 | max|U| metal 3.4135
alpha-weighted mean 0.0552 | max|U| domain 12.7787
```

Exact. Gate B matters more than Gate A because it exercises code Gate A
never touches - isoAdvector interface reconstruction, MULES bounding,
enthalpy-porosity solidification with Darcy damping, buoyancy, Marangoni -
all coupled, with no closed-form answer. Five independent statistics
reproducing to the last digit means the whole coupled path compiled to the
same arithmetic.

## 5. Two traps in the case tree, one of them committed

**`cp -r initial 0` against an existing `0/`.** POSIX `cp -r SRC DST`
creates DST as a copy of SRC when DST is absent, but copies SRC *into* DST
when it exists. Every tutorial `Allrun` runs it unguarded, and
`hartmann1d/0/` is committed - so the nesting happened, and **eleven files
under `cases/00-tutorials/hartmann1d/0/initial/` are tracked in git**. The
failure is silent: OpenFOAM resolves fields by name, so a stray `0/initial`
subdirectory is simply never read, and the case runs using whatever the
committed `0/` holds.

**It was benign, and the diff proves it.** `diff -rq initial 0` produced no
`Files ... differ` lines - every shared file is byte-identical. The only
extras are `0/alpha.metal` (generated from `alpha.metal.orig` by the
`Allrun` `find` step), `0/B` (magnetic flux density, computed by the solver
as mu_mag x H and therefore a post-run artifact; the run completed fine
without it), and `0/initial` itself. No banked number was ever produced from
contaminated initial conditions.

`Plate2D-mhd` does show `initial/alpha.metal` differing from `0/alpha.metal`
- because that case's `Allrun` runs `setFields` between the copy and the
solver. Confirmed directly: after `setFields`, `0/alpha.metal` is
**byte-identical** to the committed file. Also benign.

**Correction to `shielding-gas-eq28.md` section 6, gotcha 1.** The claim
"every case's `Allrun` is stale" is wrong. Tutorial cases (`hartmann1d`,
`Plate2D-mhd`) ship a real `initial/` and their `Allrun` is correct as
written. Only the WAAM cases generated by `tools/make_waam_case.py` ship
`0.orig/` and inherit the stale `initial/` reference.

## 6. Benchmark and smoke test - `waam-M2-bench`

A copy of `waam-arc-01-M2-full` with `endTime` cut to 0.2 s and
`numberOfSubdomains 8`, `simpleCoeffs n (4 2 1)` - four cuts along the 50 mm
travel direction, two across the bead, which keeps subdomain interfaces
small on a domain long in x.

**Mesh: 88,200 cells**, volume 0.0985 to 0.394 mm^3 (grading ratio 4),
geometry check clean on all of them. At 8 ranks that is **11,025 cells per
rank** - right at the threshold where halo exchange begins to dominate
useful work, since communication scales with subdomain *surface* while work
scales with *volume*. Twelve ranks would drop to 7,350 and sixteen to 5,500,
both into the losing regime. Another reason 8 is correct here.

| metric | value |
|---|---|
| steps to 0.2 s | 20,072 |
| ExecutionTime | 1,313.79 s (21.9 min) |
| **per step** | **0.0655 s** |
| laptop, 6 ranks | 0.16 s/step -> **2.5x faster here** |
| disk, decomposed | 386 MB |
| disk, + `reconstructPar` | 655 MB (10 writes) |
| peak memory | ~1 GB of 15 |

**Per-step cost rises then plateaus** - 0.0455, 0.0553, 0.0700, 0.0695 s
across the run. That is the pool growing: more cells crossing the solidus
means more Darcy source evaluation and a harder pressure solve, settling
once the pool reaches quasi-steady size. (The 3.5x speedup quoted early in
the session came from the cheap opening steps and is wrong; 2.5x is the
honest figure.)

**Timestep and Courant - the numbers that size the 0.25 mm run:**

```
maxDeltaT 1e-05, maxCo 0.1, maxAlphaCo 0.1, adjustTimeStep yes
peak Courant     0.1030      (breached the limit)
steps below cap  2,243 / 20,072  = 11.2%
min dt           6.452e-06   (65% of cap)
mean dt          9.964e-06   (99.6% of cap)
max dt           1.001e-05
```

So Courant pushed the solver down **repeatedly but gently**: the total
deficit against a uniform 1e-05 is 7.2e-04 s spread over 2,243 steps, an
average trim of 3.2%, with the 6.45e-06 floor a rare excursion.

The `max` sitting fractionally *above* `maxDeltaT` is `adjustableRunTime`
working correctly - to land exactly on a write time it computes
`nSteps = timeToNextWrite/deltaT` and uses `round()`, not `ceil()`, so
rounding down by a fraction of a step lengthens the remainder by ~0.1%.
Not a failure of the cap.

### Smoke test: the WAAM chain in parallel

Energy in matches M2's specification exactly - **arc 1,113.7 W + droplets
722 W (0.450 g/s at 1506 K) = 1,836 W**.

**Mass conservation is inside the banked band.** Zhao's wire - 1.2 mm
diameter at 150 mm/s - delivers pi(0.6)^2 x 150 = 169.65 mm^3/s, which at a
nominal droplet of 5.652 mm^3 is **30.0 droplets per second**: a clean
number confirming the injector is tuned to the feed rate rather than fitted.
By 0.2 s that is 6.00 droplets nominal against **5.22** measured
(29.485 mm^3), a **13% shortfall** - mid-band against M2's banked 8-23%.
Three falls in `deposited` (0.04->0.06, 0.10->0.12, 0.14->0.16) of ~0.4 mm^3
each make up the difference. Known M2 behaviour, not introduced by the port.

**`>liq` reads zero on 7 of 10 rows. This is a cold start, not a failure.**
`check_droplet.py`'s guidance - "for a real weld pool these must be non-zero
at EVERY row... deposition without a pool" - is written for the steady
state. At 0.2 s the substrate is still charging: `bottom T` climbs
monotonically and smoothly 331.6, 370.4, 413.4, 451.9, 488.6, 516.9, 540.7,
559.8, 576.6, 591.3 K, and pool cells trend upward (0, 10, 0, 31, 44, 17,
69, 24, 67, 71) to 8.88 mm^3 against M2's banked 44-51 mm^3 at 1.06-1.20 s.
The negative `d sol`/`d liq` say the same from another angle: the melted
region sits *above* z = 0, inside the bead, because nothing has yet
penetrated the plate. **The same row would be a genuine red flag at 1.1 s.**

Minor: `T bead = nan` at t = 0.02 is an empty bead. The interface-smearing
gap between `bead a>.5` and `bead aV` grows absolutely but shrinks as a
fraction (100% -> 15%), which is the bounded behaviour the tool requires.
`E kept` 84-91%, above the banked ~81% at 1.06-1.20 s, as expected while the
plate is still cool.

## 7. Disk - the binding constraint, and the format decision

| volume | free | total | note |
|---|---|---|---|
| C: | **51 GB** | 312.5 GB | holds the vhdx - **this is the budget** |
| D: | 1,476 GB | 1,863 GB | **off limits by instruction** |
| E: | 16.2 GB | 163.5 GB | 90% full - unusable, and *less* free than C: |

Linux reports 948 GB available. **Ignore it.** WSL2's ext4 lives in a
dynamically-expanding vhdx whose 1007 GB is a *virtual* maximum; real space
comes from the Windows volume holding the file. And the vhdx **grows but
never shrinks** - deleting 40 GB inside Linux frees ext4 but leaves the
Windows file 40 GB larger until explicitly compacted (see
`disk-cleanup-runbook`; never the `--set-sparse --allow-unsafe` shortcut).

Allowing Windows its 10-15% headroom leaves roughly **30-35 GB usable**.

**Use `writeCompression on`. Do not use `writeFormat binary`.**
`pool_stats_waam.py` carries a hand-written parser for OpenFOAM's *ASCII*
field syntax - it reads the declared list length rather than walking
parentheses - so binary would break the analysis pipeline outright. Gzip
gets 4-6x on OpenFOAM field data (enormously repetitive text), OpenFOAM
reads `.gz` transparently, and **all four analysis tools already handle it**:
`bead_profile.py`, `check_droplet.py`, `pool_stats_waam.py`, `pool_width.py`.

`export_exaca.py` and `exaca_coverage.py` are not in that list but do not
read fields at all - `export_exaca.py` has no `internalField` or
`nonuniform` handling, only four plain `open()` calls consuming the
recorder's own event output. **One thing to verify at the first write of the
0.25 mm run:** if the recorder emits through OpenFOAM's `OFstream` rather
than a plain `std::ofstream`, its CSV will be gzipped too and
`export_exaca.py` line 69's bare `open()` will choke.

Do **not** lower `writePrecision` from 10 to save space. ASCII truncates
doubles either way, and restarting a chaotic run from truncated fields is
exactly what `patch_droplet_restart.py` exists to manage.

Other disciplines: delete `processor*` immediately after `reconstructPar`
(a complete second copy), and treat this box as compute, not archive -
reduce each finished run to extracted quantities and the few time
directories actually cited.

## 8. The 0.25 mm run, sized

| quantity | value |
|---|---|
| cells | 705,600 (8x) |
| cells/rank at 8 ranks | 88,200 - comfortably efficient, unlike 0.5 mm |
| Courant | doubles at fixed dt, so dt falls toward ~5.5e-06 |
| steps for 1.2 s | ~218,000 |
| per step | ~0.46 s (8x cells at ~7x cost, better parallel efficiency) |
| **wall clock** | **~28 hours** |
| memory | ~8 GB of 15 |
| disk, uncompressed | **~31 GB - exceeds budget. Compression is mandatory** |

This is the first case where going above 8 ranks would genuinely help
(58,800 cells/rank at 12), but only if Windows keeps the extra threads off
the E-cores. Measure it on the real case rather than assuming.

Earlier estimates in this session of 9-12 h and 32 h are both superseded.

## 9. Open items

1. **Housekeeping commit, on the laptop not here.** Putting GitHub
   credentials in a home directory on a shared lab PC is not worth it; the
   remote stays read-only toward GitHub and pulls. Four fixes belong in one
   commit: defensive `rm -rf 0` in every `Allrun`; delete the eleven tracked
   `0/initial/*` files; add `postProcess -time 10 -func sample` to the
   README's Gate A line and correct the case name to `hartmann1d`; fix
   `tools/make_waam_case.py` to emit `0.orig` rather than `initial`.
2. Correct `shielding-gas-eq28.md` section 6 gotcha 1 (section 5 above).
3. G1 on the laptop - `gasFlux` sum(phi) and the M2 comparison over
   1.06-1.20 s.
4. Generate the 0.25 mm case with `writeCompression on`.
5. Laptop disk: `~/thesis` is **35 GB**, root filesystem at 56 GB used. Most
   of it is M1-M4, E1 and G1 field data that need not stay once extracted.
6. Set `reconstructionScheme` explicitly in `fvSolution` - it currently
   warns and falls back to a default, so the choice is implicit rather than
   on the record. Identical on both machines, so it does not affect
   comparability.
