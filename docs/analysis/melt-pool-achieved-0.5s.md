# The melt pool exists: 0.5 s full-domain droplet run

Case `waam-arc-01`, full domain, 143,616 cells, 4 ranks, Edit 4 active,
`electromagnetics false`, adiabatic walls, 15 droplets at 30 Hz, ~3 h wall
clock. This run answers the question that has blocked the project since
`Arc-Alone-Cannot-Melt.pdf`.

## Result: continuous, liquid, convecting

**Mushy continuity from t = 0.07 s (droplet 2).** The last row with zero cells
above the 815 K solidus is t = 0.06. Every row after it is non-zero, rising to
390 cells / 48.75 mm^3 at t = 0.5 and still growing.

**Liquid continuity from t = 0.20 s (droplet 6).** The last row with zero cells
above the 906 K liquidus is t = 0.19. From droplet 6 onward the pool never
fully freezes between injections.

**Sustained convection.** `U metal` last read exactly 0.0000 at t = 0.06.
Thereafter it never returns to zero: 0.2-0.6 m/s between droplets, 1.2-1.6 m/s
at impact. This is the precondition for the entire magnetic study -- J x B
changes microstructure by changing flow, and there was nothing to act on until
now.

**Real fusion depth.** 0.25 mm (one cell, a mesh artefact) at droplet 1, then
0.75 mm at t = 0.14, 1.25 mm at t = 0.24, 1.75 mm at t = 0.44 -- 3.5 cells, and
still stepping upward at the end.

## The saturation prediction, checked

The mechanism argued beforehand was that droplet lifetime diverges as the
surrounding metal approaches solidus, because the time to cool from 1506 K to
815 K goes as -tau ln[(T_sol - T_inf)/(T_drop - T_inf)]. Fitting the
between-droplet maxima of `T metal` (13 points, t = 0.09 to 0.49) to
T = T_sat - A exp(-t/tau):

    T_sat = 958 K    A = 165 K    tau = 192 ms    rms 6.1 K

- **T_sat = 958 K clears the 906 K liquidus by 52 K.** The pool is
  self-sustaining, not a sequence of freezing events.
- The fitted tau = 192 ms is close to the predicted arc saturation timescale
  b^2/alpha = 257 ms for the Goldak b = c = 4 mm source. The 0.1 s run stopped
  at 39% of that, which was the whole reason it showed no pool.
- The linear power-scaling estimate made in advance was 1009 K (droplet
  sensible enthalpy only) to 1071 K (with latent). Actual 958 K, so the
  extrapolation **overestimated by 5-12%** -- respectable for a linear
  extrapolation of a nonlinear system, and the direction and the crossing of
  liquidus were both right. The overestimate is expected: droplet heat is
  delivered concentrated at the surface and a large share goes into heating the
  bead itself rather than raising the substrate peak.
- Timing was predicted at droplet 4-5. Mushy continuity actually arrived at
  droplet **2** (early) and liquid continuity at droplet **6** (close).

## Conservation, independently verified

- **Mass.** 84.499 mm^3 deposited against 84.823 mm^3 of wire delivered; the
  0.324 mm^3 difference is 0.057 of a droplet still accumulating in
  `dropOwed_`. Bead cross-section 8.450 mm^2 against the 8.482 mm^2 fixed by
  A = Vdot/v_travel before the run -- **-0.38%**, better than the 0.1 s run's
  -1.3%. 14.95 droplets placed.
- **Energy.** `E kept` 137.8%, between the 133.6% predicted from arc plus
  droplet sensible enthalpy alone and the 145.3% that would follow if every
  droplet's latent heat had been released. It should sit between: a 48.75 mm^3
  pool is still liquid and has not given its latent heat back.

## Open items this run created

1. **Edit 10 is now required.** `bottom` reached 557.7 K (arc-alone at 0.5 s was
   517.2 K). The adiabatic wall is holding in 45-89 W depending on the mean wall
   temperature, 2-4% of the 2,191 W input. It biases toward melting, so it did
   not manufacture this result, but it will distort any quantitative pool
   dimension compared against Zhao.

2. **The far-field gate was misread, and `pool_stats` rev 4 fixes it.** The
   single `far` column jumped 300.0 -> 321.4 at t = 0.12 and then wandered
   308 -> 311 -> 347 -> 374 -> 426 -> 415. Conduction into a far boundary is
   smooth and monotonic; that is advected argon, blown 25 mm sideways at the
   4-10 m/s shown in `U gas`. `far` is now split into `far met` (the actual
   domain-independence gate) and `far gas`. Hot argon piling against a wall it
   cannot leave through is a real issue for spurious buoyancy, but it is a
   gas-domain question, not a substrate one.

3. **Nothing has saturated geometrically.** Pool cells, pool volume and depth
   are all still climbing at t = 0.5. The peak temperature is ~99% of its fitted
   asymptote but the pool *shape* is not steady, so 0.5 s is not yet a state to
   compare against Zhao's steady bead. Longer runs are needed before the
   magnetic sweep has a baseline.

4. **No crater.** `bead aV` equals `deposited` to the digit at every row, so no
   metal has moved below z = 0 and the free surface is undeformed despite
   0.9 m/s impingement. Worth inspecting in ParaView before Edit 3 (arc
   pressure) goes in.
