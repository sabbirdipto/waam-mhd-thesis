# WAAM-MHD thesis

Undergraduate thesis: the effect of magnetic fields on melt-pool flow and
microstructure in GMAW-based wire-arc additive manufacturing of Al-5Mg
(5A05 substrate, ER5356 wire). Replicates Zhao et al. (2021) Weld World
65:1571 in OpenFOAM v2506 (ESI), then adds external fields and couples the
thermal history to ExaCA.

- `src/lbf-com/`  laserbeamFoam, vendored (upstream in `src/lbf-com-PROVENANCE.txt`),
  WAAM additions in `applications/solvers/laserbeamFoam/`: Goldak arc source,
  arc pressure, droplet transfer, arc Lorentz force
- `cases/`    case recipes only; regenerate fields with
  `cp -r 0.orig 0 && blockMesh && setFields && decomposePar`
- `tools/`    solver patch scripts and analysis tools
- `results/`  banked text output for every reference run

Verification gates, after every build:
- A: `cases/00-tutorials/hartmann1d-postpatch`, `python3 hartmann_check.py --sigma 1000 --mu 1`
  -> u(0) 0.0158753, Ha 15.810, err 1.3869e-02
- B: `cd cases/00-tutorials && python3 pool_stats.py Plate2D-mhd-postpatch`
  -> 1869.0 / 0.0504 / 3.4135 / 0.0552 / 12.7787
