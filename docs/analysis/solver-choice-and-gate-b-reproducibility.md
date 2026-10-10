# GAMG for the 0.25 mm mesh, and what it revealed about Gate B

**Date:** 2026-09-23, during R0 commissioning on DESKTOP-3NBGHNH.
**Two results.** One is a 2x speedup. The other is more important: **Gate B's
velocity statistics are not a valid regression metric**, and never were.

---

## 1. The problem: R0 was a five-day run

R0 (`waam-arc-02-R0`, 862,400 cells, 8 ranks) started at **1.814 s/step**,
measured over two intervals so the dynamic-code compile could be excluded.
At ~240,000 steps for 1.2 s that is **five days**.

The log said exactly where it went:

```
DICPCG: p_rgh ... No Iterations 121     (relTol 0.05)
DICPCG: p_rgh ... No Iterations 124     (relTol 0.05)
DICPCG: p_rgh ... No Iterations 167     (p_rghFinal, relTol 0, tol 1e-07)
smoothSolver: T ... No Iterations 3
```

**412 conjugate-gradient iterations per timestep on pressure, 3 on
temperature.** That asymmetry is structural, not a bug. The pressure equation
is elliptic - an incompressible fluid has no finite signal speed, so every
cell couples to every other cell and there is no local shortcut. Temperature
is parabolic: information travels a few cells per timestep and a local
smoother resolves it in three sweeps.

Worse, the first corrector needed 121 iterations to reduce the residual by a
factor of only 21. CG's convergence goes as sqrt(kappa), and for a Poisson
problem kappa ~ 1/h^2, so iterations ~ 1/h. The 1600:1 metal/argon density
ratio inflates kappa further, because `1/rho` jumps by that factor across the
interface and `DIC` - diagonal incomplete Cholesky - is a purely *local*
preconditioner that cannot see the long-wavelength pressure modes the jump
creates.

## 2. Why this never appeared at 0.5 mm

**Refinement created the problem; it did not reveal it.** Halving the cell
size doubles the CG iteration count and multiplies the cells by 9.78, so the
pressure solve costs roughly **20x more** at 0.25 mm than at 0.5 mm. On the
0.5 mm mesh it was one cost among several and runs finished in 5.5 h on the
laptop. Nobody had a reason to look.

Plate2D demonstrated this directly and is worth keeping as the counterexample:
at that size `DICPCG` at `relTol 0` needed only **21/9/1 iterations**, because
a small 2D mesh is well conditioned. There was no waste to remove, and GAMG
produced **no speedup at all there** (202.62 s vs 204.93 s). Any solver
benchmark run on Plate2D would have concluded, wrongly, that GAMG is useless.

## 3. The four-run experiment

The intent was to check that GAMG does not change the answer, using Gate B
because it has banked numbers. The result was something else.

| run | solver | relTol | melt vol | mean\|U\| | max\|U\| metal | alpha-wtd | max\|U\| domain |
|---|---|---|---|---|---|---|---|
| banked | DICPCG | 0.05 | 1869.0 | 0.0504 | 3.4135 | 0.0552 | 12.7787 |
| GAMG loose | GAMG | 0.05 | 1858.0 | 0.0427 | 3.4351 | 0.0475 | 10.5105 |
| **A** | GAMG | **0** | 1866.0 | 0.0445 | 3.3266 | 0.0437 | 7.2920 |
| **B** | DICPCG | **0** | 1859.0 | 0.0396 | 3.8742 | 0.0438 | 15.5120 |

**Compare B with the bank. Same solver. Only `relTol` changed, 0.05 to 0.
`mean|U|` moves 21.4%** - *larger* than the 11.0% between A and B, which is
the actual solver change.

Across all four configurations:

```
melt volume          1858 .. 1869      range 0.59%     ROBUST
mean|U| in metal     0.0396 .. 0.0504  range 27%
alpha-weighted mean  0.0437 .. 0.0552  range 26%
max|U| whole domain  7.29 .. 15.51     range 113%      MEANINGLESS
```

The mechanism is visible in the residuals. `relTol 0.05` says "reduce the
residual 20-fold and stop". DICPCG grinds out exactly that (1.64e-04 ->
7.67e-06, ratio 21) and quits. One GAMG iteration overshoots to 170-fold
(2.15e-03 -> 1.26e-05) and cannot be asked for a fraction of an iteration. Both
satisfy the stopping criterion; they stop at **different points inside the
tolerance ball**, and 7,000 timesteps compound the difference.

## 4. Determinism is not convergence

Gate B reproduced **bit-for-bit** across two machines and two gcc versions.
That looked like reproducibility. It was **determinism**: identical source,
identical inputs, identical arithmetic sequence, identical answer. It proved
the build was correct, which is what a gate is for. It said nothing about
whether those numbers are converged with respect to the numerics.

The only way to tell the two apart is to run the same case twice with a
deliberately perturbed numerical setting. Nobody had, because there was no
reason to suspect a problem. It surfaced by accident while chasing speed.

**Consequences:**

1. **Melt volume is the gate.** It is the thermally-driven integral, it agrees
   to 0.4% between A and B and to 0.59% across all four, and it is the number
   to test a build against.
2. **`mean|U|`, `max|U|` and the alpha-weighted mean are indicative at best.**
   The "~2.6%" printed by `pool_stats.py` is a write-to-write floor, not a
   numerical-path floor. The real path floor is ~20%. The README currently
   presents all five statistics as equivalent; it should not.
3. **Re-examine any banked conclusion that rests on a velocity statistic.**
   The Marangoni gate (M1) and the Darcy sweep (M3) both involve flow. A ~20%
   numerical uncertainty may not overturn either, but that has to be checked
   rather than assumed.
4. This is further evidence for `Chaotic-Sensitivity-in-WAAM-CFD.pdf`, now in
   a case that was believed reproducible.

## 5. Decision and measured result

**Adopted for R0: GAMG-preconditioned PCG on `p_rgh`, `relTol 0`.**

GAMG is used as a *preconditioner* inside CG rather than as the solver. With a
1600:1 density jump the agglomeration straddles the discontinuity and pure
GAMG can stall; wrapping it in CG keeps Krylov robustness while GAMG removes
the smooth error CG is bad at. `cacheAgglomeration true` matters on a static
mesh - the hierarchy is built once, not 240,000 times.

`relTol 0` rather than 0.05: with GAMG it costs 3 iterations instead of 1,
essentially free, and it removes the 21% sensitivity documented above.

**Measured on R0, 473 steps:**

| | DICPCG | GAMG-PCG |
|---|---|---|
| iterations per step | 121/124/167 = **412** | 3/3/2 = **8** |
| s/step | 1.814 | **0.9207** |
| RSS, 8 ranks | 3.09 GB | 3.14 GB |
| 1.2 s projection | ~5.0 days | **~2.5-3 days** |

**A correction worth recording: the pressure solve was about half the step,
not 85%.** If GAMG's pressure cost is near zero, the residual 0.92 s/step is
everything else - momentum, energy, isoAdvector/MULES interface advection, the
MHD and droplet sources, the recorder. So 51x fewer iterations bought 1.97x
overall, and further pressure tuning has little left to give. The remaining
cost is the physics, not the linear algebra.

## 6. Open

- Verify at the first write (t = 0.04) that the ExaCA recorder emits plain
  text rather than being gzipped by `writeCompression on` - `export_exaca.py`
  line 69 uses a bare `open()` and would choke on a `.gz`.
- Amend the README so Gate B is stated as a melt-volume test with velocity
  reported as indicative.
- `fvSolution.dicpcg.bak` is kept in both `waam-arc-02-R0` and `Plate2D-mhd`;
  Plate2D has been restored to DICPCG.
- The R0 recipe is still untracked and exists only on the lab PC.
