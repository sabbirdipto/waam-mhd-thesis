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
  events 3610 in 1231 cells, 916 remelted
  domain (cell centres) mm: x -11.25..29.75  y -1.75..2.25  z -1.25..2.75
  cr K/s median 2099   G K/m median 7.773e+04

COVERAGE WARNING (added 20 Sep 2026, checked against the bead at t = 2.0 s)
  In the steady window (x 6..26 mm) only 19% of the bead cells above the
  plate have a solidification event (234 of 1234). The events form a
  column about 2 mm wide on the centreline, where the droplets land, plus
  the melted plate below it (406 cells). The rest of the bead is missing
  for two reasons: metal that reached a cell already below the 906 K
  liquidus (mushy) never opened an event, and 8064 open events were
  discarded when the metal/gas surface moved through a cell and its metal
  fraction fell below 0.5. So the CAVEATS line "Cells inside the bounding
  box that never appear here never melted" is NOT true inside the bead
  for this dataset: ExaCA will treat 81% of the bead as never melted.
  USE THIS DATASET TO SET UP AND TEST THE ExaCA PIPELINE, NOT FOR
  GRAIN-STRUCTURE RESULTS.

RESOLUTION NOTE (added 20 Sep 2026)
  ExaCA 1.3 and later ignore HeatTransferCellSize: temperature data must
  be given at the CA cell size, and earlier versions could not combine a
  coarser thermal grid with remelting. This file is on a 500 um grid, so
  as it stands it can only be read with CellSize = 500 (a test of the
  input chain). Grain-resolving runs need the data resampled onto the CA
  grid first. See E1-ExaCA-Data-Explained.pdf, sections 3.3 and 9.1.

SLAB FILE
  exaca_E1_slab_x14-16mm.csv holds the same events cut to 14 <= x <= 16
  mm, the middle of the steady window: 203 events in 63 cells, 54 of them
  melted more than once. Cell centres x 14.25..15.75, y -0.75..0.75,
  z -1.25..1.25 mm, i.e. a box of 2.0 x 2.0 x 3.0 mm = 12.1 M ExaCA cells
  at 10 um or 96.4 M at 5 um. Melting 1.201..1.574 s, freezing
  1.204..1.575 s, so every event completed before the run ended. Median
  cr 1826 K/s, median G 78.6 K/mm. The torch passes x = 15 mm at
  t = 1.25 s.
