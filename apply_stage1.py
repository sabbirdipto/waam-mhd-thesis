#!/usr/bin/env python3
"""Stage 1: trim the domain to Zhao's coupon, thin the plate to 4 mm, extend
the fine region past the bead end, open the top boundary, add the Eq. 29
convective wall condition.

Run from the repository root:   python3 apply_stage1.py

REFUSES RATHER THAN GUESSES. Every edit is an exact-string replacement
asserted to match exactly once. If make_waam_case.py has moved on, the script
names the failing edit and writes nothing.

Idempotent. Does NOT touch constant/g, decomposeParDict (already simple 4 1 1
in the generator), or the solver source.

WHAT CHANGES AND WHY
--------------------
1. Domain 60x50x12 -> 50x30x9 mm.
   The domain is not a truncation of an infinite plate, it IS the workpiece.
   Zhao welded a 50 x 30 x 4 mm 5A05 coupon. Its edges are 15 mm from the bead
   and they genuinely warm up; converging to an infinite-plate answer converges
   to the wrong experiment. The old rationale (25 mm verified domain
   independent) optimised the wrong quantity.

2. Plate 6 -> 4 mm.
   sqrt(alpha*t) over the ~0.7 s between arc passage and solidification is
   5.9 mm, so plate thickness is INSIDE the thermally active depth and sets
   the local cooling rate during exactly the window that decides bead height.
   Also sets L in the Kumar-DebRoy Lorentz force (Zhao Eqs. 22-24).
   Thermal mass falls 3.0x overall, so the plate now heats as Zhao's does.

3. Argon 6 -> 5 mm.
   Multilayer is out of scope, so 14 mm of headroom is no longer needed. A
   ~2 mm bead plus a droplet radius plus standoff needs 3.6 mm; 5 mm is ample.

4. Fine region x [-10,+20] -> [-12,+32].     <-- BUG, independent of the trim
   The arc travels to x = +30 mm over a 2.0 s pass, but the fine region
   stopped at +20 mm, putting the arc-off collapse in 1.6 mm cells. That is
   exactly where 'deposit length' is measured (Zhao Table 5: 42.47 mm).
   Never showed up because every run so far stopped at t <= 0.5 s, where the
   arc has only reached x = 0.

5. EXPAND 4.0 -> 2.0.
   Trimming leaves 3 mm of graded margin in x. At EXPAND 4 that is 3 cells
   growing 2.0x each, too abrupt. 2.0 gives 1.19x adjacent, inside the usual
   guideline, for 17k extra cells.

6. inlet velocity fixedValue (0 0 0) -> pressureInletOutletVelocity.
   A zero-velocity inlet is hydrodynamically a LID. The buoyant argon plume
   rises, cannot leave, spreads along the roof and reaches farField, which is
   a no-slip wall it also cannot pass. This is why 'far gas' kept climbing
   after the gravity fix. Proper fix is Zhao Eq. 28 (20 L/min annular argon
   jet); this is the minimum that lets the plume escape.

7. substrateBottom + farField: zeroGradient -> Eq. 29 convective condition.
   zeroGradient is a perfect insulator. Zhao Table 3 puts Eq. 29 (h = 80,
   plus radiation) on every wall patch. It is ~1-2% of the energy budget over
   2.5 s, so it will not transform SPSL, but the walls are now physical coupon
   surfaces exposed to air and adiabatic is no longer defensible at 15 mm.
   h = 80 W/m2K is high for natural convection off a small plate (5-15 would
   be typical) because it is absorbing the fixture conduction path too. State
   that in the writeup rather than presenting 80 as a convection coefficient.

RESULT: 98 x 50 x 18 = 88,200 cells, 39% fewer than 143,616, with a LONGER
fine region than before.
"""
import sys
from pathlib import Path

GEN = Path.cwd() / "tools" / "make_waam_case.py"

EDITS = []

# ---------------------------------------------------------------- 1,2,3,4,5
EDITS.append(("domain extents + rationale", '''X0, X1    = -0.025, 0.035   # 60 mm''',
              '''X0, X1    = -0.015, 0.035   # 50 mm -- Zhao's coupon length'''))

EDITS.append(("y extent + rationale", '''# y is mirrored for the full domain. The far boundary stays at 25 mm on each
# side, because 25 mm is what was VERIFIED domain-independent: at 0.5 s the
# outermost cell layer sat at 300.2 K, 0.04% of the peak rise. Mirroring
# carries that result over unchanged. Narrowing to Zhao's real 30 mm plate
# would save 24% of the cells but put the boundary 1.3 diffusion lengths out,
# where it starts to matter and Edit 10's convective condition becomes
# necessary -- a poor trade for a quarter of the cells.
Y1        =  0.025
FY1       =  0.010          # fine region edge, +y''',
'''# y is mirrored for the full domain. The far boundary is Zhao's REAL plate
# edge: his coupon is 50 x 30 x 4 mm, so the long edges sit 15 mm from the
# bead. The previous 25 mm was chosen for domain independence, which is the
# right target when a boundary is an artificial truncation of a larger body
# and the wrong one when the boundary is the physical edge of the workpiece.
# Converging to the infinite-plate answer converges to a different experiment.
# The boundary IS thermally live at 15 mm, which is why Eq. 29 (h = 80) now
# sits on farField -- that condition is not a cost of trimming, it is what
# Zhao Table 3 specifies regardless of domain size.
# Consequence for reading pool_stats_waam.py: 'far met' is no longer a
# domain-size gate. It is a physical plate edge and is EXPECTED to warm. The
# gate is now that it stays well below the 815 K solidus.
Y1        =  0.015          # 30 mm full width -- Zhao's coupon
FY1       =  0.006          # fine to +-6 mm, one diffusion length sqrt(a*t)'''))

EDITS.append(("z extent + rationale", '''# 6 mm substrate, 6 mm argon. The argon was 3 mm and that was too short once
# a bead exists. Mass conservation fixes the bead cross-section at 8.5 mm^2,
# which over an ~8 mm width is a cap about 1.6 mm tall, leaving 1.4 mm of
# clearance -- less than the 2.21 mm diameter of a detached droplet (Zhao
# Eq. 30). Nothing warned about this; the bead would simply have grown into
# the inlet patch. 6 mm leaves 4.4 mm of clearance above a finished bead.
# Cost: nz goes 18 -> 24, so 53,856 -> 71,808 cells, +33%. The substrate and
# both graded directions are untouched, so the domain-independence result
# still holds -- argon carries a negligible share of the thermal energy.
Z0, Z1    = -0.006, 0.006''',
'''# 4 mm substrate, 5 mm argon.
# SUBSTRATE: 4 mm is Zhao's actual baseplate. This is not a cost-saving. The
# conduction penetration depth over the ~0.7 s between the arc passing a
# cross-section and that metal freezing is sqrt(alpha*t) = 5.9 mm, so at 4 mm
# the bottom face is INSIDE the thermally active region and sets the local
# cooling rate during the window that decides bead height. At 6 mm the plate
# was 50% too heavy a heat sink. It also sets L in Zhao Eqs. 22-24.
# ARGON: multilayer is out of scope (SPSL only), so the 14 mm a five-layer
# wall would need is irrelevant. A ~2 mm bead plus a 1.1 mm droplet radius
# plus 0.5 mm standoff needs 3.6 mm of clearance; 5 mm leaves margin.
Z0, Z1    = -0.004, 0.005'''))

EDITS.append(("fine region in x", '''FX0, FX1  = -0.010, 0.020   # the pool and its trailing track''',
'''FX0, FX1  = -0.012, 0.032   # the WHOLE 40 mm bead, not just the pool.
                            # Was [-10,+20]: the arc reaches x = +30 mm over a
                            # 2.0 s pass, so the bead end and the arc-off
                            # collapse sat in 1.0-1.6 mm cells -- and that is
                            # exactly where 'deposit length' is measured
                            # (Zhao Table 5, 42.47 mm). Invisible until now
                            # because every run stopped at t <= 0.5 s, where
                            # the arc has only reached x = 0.'''))

EDITS.append(("expansion ratio", '''EXPAND    = 4.0             # coarsest/finest cell within a graded section''',
'''EXPAND    = 2.0             # coarsest/finest cell within a graded section.
                            # Was 4.0. Trimming leaves only 3 mm of graded
                            # margin in x, and 4.0 over 3 cells means adjacent
                            # cells growing 2.0x, well past the ~1.2 guideline.
                            # 2.0 gives 1.19x adjacent for ~17k extra cells.'''))

# ------------------------------------------------------------------------ 6
EDITS.append(("open the top boundary", ''' "U": vol("U", "volVectorField", "[0 1 -1 0 0 0 0]", "(0 0 0)", {
     "inlet":    {"type": "fixedValue", "value": "uniform (0 0 0)"},''',
''' # inlet was fixedValue (0 0 0), which is hydrodynamically a LID: zero
 # velocity means zero through-flow, so the buoyant argon plume rises, cannot
 # leave, spreads along the roof and reaches farField -- a no-slip wall it
 # also cannot pass. That is why 'far gas' kept climbing even after gravity
 # was put on the correct axis. pressureInletOutletVelocity lets the plume
 # out while still admitting inflow where the pressure field wants it.
 # This is the MINIMUM fix. Zhao's actual top boundary (Table 3, Eq. 28) is a
 # 20 L/min annular argon jet directed downward between the wire and the
 # nozzle, which sweeps the plume out through the side pressure outlets.
 "U": vol("U", "volVectorField", "[0 1 -1 0 0 0 0]", "(0 0 0)", {
     "inlet":    {"type": "pressureInletOutletVelocity",
                  "value": "uniform (0 0 0)"},'''))

# ------------------------------------------------------------------------ 7
EDITS.append(("Eq. 29 wall condition", ''' # 'walls' is the bottom AND the far side, which are different in kind: the
 # far side is a truncation of a larger plate, the bottom is the substrate's
 # real underside. zeroGradient is the least-wrong single choice for both over
 # 0.5 s -- a truncation should not act as a sink, and fixture contact is slow
 # on this timescale. Verify by running this against fixedValue 300: with the
 # boundaries two diffusion lengths out the two must AGREE, and if they do not
 # the domain is still too small. Note that until Edit 8 this test was
 # meaningless, because kappa was zero on every boundary and no condition did
 # anything at all.''',
''' # Both substrateBottom and farField are now PHYSICAL surfaces of Zhao's
 # 50 x 30 x 4 mm coupon, exposed to air. Neither is a truncation any more,
 # so zeroGradient -- a perfect insulator -- is no longer defensible for
 # either. Zhao Table 3 puts Eq. 29 on every wall patch:
 #     -k dT/dn = h(T - Ta) + eps*sigma*(T^4 - Ta^4)
 # with h = 80 W/m2K. Note that 80 is high for natural convection off a small
 # plate (5-15 W/m2K would be typical); it is a lumped coefficient absorbing
 # the conduction path into the clamp and workbench. Say so in the writeup
 # rather than presenting it as a convection coefficient.
 # Magnitude: 21.4 cm2 of wall at dT = 300 K loses 51 W, which is 128 J over a
 # 2.5 s pass against 5402 J in -- about 2%. Locally under the arc it is
 # 0.032 W/mm2 against 30 W/mm2 arriving, a factor of 1000. So this will not
 # transform the SPSL answer; it is here because it is what the paper does and
 # because an adiabatic wall 15 mm from the bead reflects heat that should be
 # leaving.
 # kappaMethod lookup works only because patch_kappa_bc.py gave kappa real
 # boundary conditions; without that this fails at construction.'''))

EDITS.append(("T wall patches", '''     # written separately, not via walls(), because the domain check
     # varies farField alone while substrateBottom is held fixed
     "substrateBottom": {"type": "zeroGradient"},
     "farField":        {"type": "zeroGradient"}}),''',
'''     # written separately, not via walls(), because they remain independently
     # variable for sensitivity runs even though both now carry Eq. 29
     "substrateBottom": {"type": "externalWallHeatFluxTemperature",
                         "mode": "coefficient", "h": "uniform 80",
                         "Ta": "uniform 300", "kappaMethod": "lookup",
                         "kappa": "kappa", "value": "uniform 300"},
     "farField":        {"type": "externalWallHeatFluxTemperature",
                         "mode": "coefficient", "h": "uniform 80",
                         "Ta": "uniform 300", "kappaMethod": "lookup",
                         "kappa": "kappa", "value": "uniform 300"}}),'''))


def main():
    if not GEN.exists():
        sys.exit(f"ERROR: {GEN} not found. Run from the repository root.")
    src = GEN.read_text()

    if "Zhao's coupon length" in src:
        print("already patched -- nothing to do"); return

    for label, old, new in EDITS:
        n = src.count(old)
        if n != 1:
            sys.exit(f"REFUSING: edit '{label}' matched {n} times, expected 1.\n"
                     f"  NOTHING WRITTEN. Apply that edit by hand.")
        src = src.replace(old, new, 1)

    compile(src, str(GEN), "exec")
    GEN.write_text(src)
    print(f"patched {GEN}  ({len(EDITS)} edits)")
    print("\nnext:")
    print("  python3 tools/make_waam_case.py --force")
    print("  cd cases/00-tutorials/waam-arc-01 && ./Allrun -f")


if __name__ == "__main__":
    main()
