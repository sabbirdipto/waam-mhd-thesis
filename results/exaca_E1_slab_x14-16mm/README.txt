ExaCA input from OpenFOAM case waam-arc-01-E1-exaca
==============================================

FILES
  exaca_waam-arc-01-E1-exaca.csv   ExaCA temperature file: x,y,z,tm,tl,cr
  gr_waam-arc-01-E1-exaca.csv      solidification map:    x,y,z,tl,cr,G,R

UNITS            x,y,z m | tm,tl s | cr K/s | G K/m | R = cr/G m/s

GRID             regular, spacing 500 um (HeatTransferCellSize = 500 um).
                 Points are cell centres at (k + 1/2) * 0.5 mm; every point in
                 the file was checked to lie on that lattice.
                 Rows sorted by position, then time: a cell that melted more
                 than once has consecutive rows in chronological order.

COORDINATES      z = 0 is the top of the 4 mm plate, +z is the build direction.
                 y = 0 is the bead centreline. The torch starts at x = -10 mm at
                 t = 0 and travels +x at 20 mm/s.

DEFINITIONS      liquid metal : metal fraction > 0.5 and T >= 906 K (liquidus)
                 tm : time the cell became liquid metal. Interpolated liquidus
                      crossing if it melted in place; arrival time if metal
                      arrived already molten (droplet)
                 tl : time it cooled back through the liquidus, interpolated
                 cr : 10 K / (time taken to cool from 906 K to 896 K),
                      i.e. the mean cooling rate over the first 10 K of
                      undercooling -- the quantity ExaCA's linear undercooling
                      uses. A dip below the liquidus that re-melts before
                      896 K does not count as solidification.
                 G  : |grad T| at tl (unreliable at the free surface, where
                      the neighbouring cell is gas)
                 Recorded every solver time step (~1e-5 s), not from snapshots.

MATERIAL         Al-5Mg (5A05 plate, ER5356 wire). Solidus 815 K, liquidus
                 906 K, freezing range 91 K.
                 ExaCA needs an interfacial response function for Al-5Mg; it is
                 not part of this package.

PROCESS          GMAW WAAM, 135 A, 16 V, travel 20 mm/s, wire 1.2 mm at 150 mm/s,
                 single bead on a 50 x 30 x 4 mm plate (Zhao et al. 2021).

CAVEATS          Model state "M2": energy-consistent arc (1113.7 W into the pool,
                 722 W carried by droplets) with Zhao's heat-source dimensions.
                 Penetration matches Zhao (~1 mm liquid) but the fused zone is
                 ~3.6 mm wide against his ~7.8 mm bead, so solidification
                 conditions are those of a narrower pool than the experiment.
                 Treat this as pipeline data. The calibrated 0.25 mm runs will
                 replace it with the same format.
                 The first ~5 mm of track is a start transient; use the per-x
                 table below to choose a steady region.
                 Cells inside the bounding box that never appear here never
                 melted: plate below z = 0, gas above the bead. Only cells
                 listed here have a solidification history.

SUMMARY
  events 203 in 63 cells, 54 remelted
  domain (cell centres) mm: x 14.25..15.75  y -0.75..0.75  z -1.25..1.25
  cr K/s median 1826   G K/m median 7.861e+04
