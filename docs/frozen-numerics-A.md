# Frozen configuration for the ablation set A1'-A4

Every case in the ablation study MUST match this exactly. `--force` resets
endTime and maxDeltaT to generator defaults, so re-set them after every rebuild.

mesh        98 x 50 x 18 = 88,200 cells
            domain x[-15,+35] y[-15,+15] z[-4,+5] mm  (50 x 30 mm plate,
            4 mm substrate, 5 mm argon) = Zhao's coupon
            fine region x[-12,+32] y[-6,+6], CELL 0.5 mm, EXPAND 2.0
arc_x0_m    -0.010          bead -10 -> +30, 5 mm lead-in and run-out
DECOMP      (6, 1, 1)       x-only: preserves exact mirror symmetry in y,
                            which scotch does not (cost us an 8.9% false
                            asymmetry signal)
endTime     0.46            first developed window is t > 0.428
deltaT      4e-05   maxDeltaT 4e-05   maxCo 0.1   maxAlphaCo 0.1
                            Courant binds late in the run as gas reaches
                            ~4.9 m/s; dt falls to ~1e-5
libs        ("libthermoTools.so")   for externalWallHeatFluxTemperature
gravity     (0 0 -9.81)
top         U pressureInletOutletVelocity; T/alpha inletOutlet
walls       externalWallHeatFluxTemperature, h = 80, Ta = 300  (Zhao Eq. 29)

## Cases
A1'  gravity only                        <- reference (this run)
A2   + arc pressure (Eq. 25, 228 Pa)
A3   + pool Lorentz (Eqs. 22-24, 435 Pa)
A4   + both                              <- full model

## Record per case
A, h crown, h max, w toe, w half, ripple at t = 0.44 and 0.46
pool depth, U metal, U gas, bottom, far met

## Open issues
- gas velocity ~4.9 m/s is ~5x the buoyancy scale sqrt(g*beta*dT*L); either
  parasitic VOF currents or Boussinesq applied at beta*dT ~ 4 in the argon.
  Sets the time step. Does not affect melt-side physics (U metal ~0.46 m/s).
- maxCo 0.1 is 3x more conservative than Zhao, who ran 0.25 mm cells at a
  fixed dt 4e-5, i.e. Co ~ 0.3. Test 0.1 vs 0.3 back to back before the
  production run.
- 0.5 mm gives 4.4 cells per droplet diameter. Acceptable for RELATIVE
  comparison on identical grids; the Table 5 validation needs 0.25 mm
  (684k cells, 8.8 cells per droplet).
