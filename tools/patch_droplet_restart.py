#!/usr/bin/env python3
"""Make the droplet accumulator survive a restart. It currently does not.

THE DEFECT

    createFields.H declares the droplet state as plain C++ scalars:

        scalar ... dropOwed_ = 0.0, dropImpactEnd_ = -1.0, ...
        label  dropNum_ = 0;

    These are NOT registered IOobjects. Nothing writes them, nothing reads them
    back. Every restart re-initialises them, which throws away the partial
    droplet in the accumulator and cancels any impact window in progress.

MEASURED COST, from the M1-mushy run stopped and resumed at t = 1.16 s

    The accumulator held  169.646 mm^3/s * 1.16 s - 34 * 5.6517 mm^3
                        = 4.63 mm^3 = 0.82 of a droplet.  Discarded.

    check_droplet across the seam -- three consecutive writes, no new metal:

        time   deposited    drops
        1.14    190.253     33.66
        1.16    190.253     33.66     <-- identical
        1.18    190.253     33.66     <-- identical
        1.20    195.747     34.64

    0.04 s during which 1.2 droplets were due. Starved of 37 g impacts and of
    1506 K enthalpy, the pool collapsed:

        time   >sol   >liq
        1.14    428     63
        1.18    323      0      <-- no liquid anywhere
        1.22    270      1      <-- -37% from 1.14
        1.36    408     69      <-- recovered, 0.2 s later

    0.2 s at 20 mm/s is 4 mm of track. The bead shows it as a trough to
    5.2 mm^2 against the conserved 8.482, with a compensating overshoot to
    12.1 mm^2 where the recovery caught up. The run finished 59 droplets
    instead of 60.

    Total: 4.63 mm^3 lost outright, ~11 mm^3 redistributed over 10 mm of bead,
    and roughly 4 mm of geometry that cannot be used.

WHY IT MATTERS MORE FROM HERE

    Runs are getting longer -- 2.0 s is 15 h at the current timestep -- and a
    long run that cannot be resumed safely is a run that must be repeated in
    full if anything interrupts it. That is the real cost; the 4.63 mm^3 is
    just the receipt.

THE FIX

    Standard OpenFOAM idiom for solver state that must survive a restart: an
    IOdictionary in <time>/uniform, READ_IF_PRESENT so a fresh case starts from
    the same defaults, AUTO_WRITE so it is written with every field write.

    Three values are enough. dropR_, dropV_, dropT_, dropVdot_, dropVdrop_,
    dropMass_, dropTau_ and dropGap_ are all re-derived from caseParameters at
    startup and are identical on a restart, so they are deliberately NOT stored
    -- persisting a value that is also read from a dictionary is how the two
    drift apart later.

    dropXc_, dropYc_, dropZc_ are recomputed from the torch position each step
    and likewise do not need storing.

SAFETY
    READ_IF_PRESENT means a case with no <time>/uniform/dropletState gets
    exactly the current defaults, so a run from t = 0 is bit-identical. That is
    what Gate A and Gate B check.

USAGE
    python3 patch_droplet_restart.py            # apply
    python3 patch_droplet_restart.py --revert
"""
import sys, shutil
from pathlib import Path

SOLVER = Path.home() / "thesis/src/lbf-com/applications/solvers/laserbeamFoam"
CF = SOLVER / "createFields.H"
MAIN = SOLVER / "laserbeamFoam.C"
NEW = SOLVER / "syncDropletState.H"
SUFFIX = ".drbak"

# ------------------------------------------------------------------ new header
NEW_SRC = r'''// Persist the droplet accumulator so a restart does not discard it.
// Included immediately after updateDropletSource.H, so the values written are
// the ones that step just produced. See patch_droplet_restart.py.
{
    if (dropOn)
    {
        dropletState.set("dropOwed",      dropOwed_);
        dropletState.set("dropImpactEnd", dropImpactEnd_);
        dropletState.set("dropNum",       dropNum_);
    }
}
'''

# ------------------------------------------------------------- createFields.H
# anchored AFTER the whole droplet setup block, so dropVdrop_ is already set
CF_ANCHOR = ("// Momentum source, added to UEqn beside Marangoni and JcrossB."
             " zeroGradient")

CF_ADD = r'''// --- Droplet state across restarts (WAAM-MHD thesis) -----------------------
// The scalars above are plain C++ and do not survive a restart. Resuming the
// M1-mushy run at t = 1.16 s discarded 4.63 mm^3 (0.82 of a droplet) from
// dropOwed_, starved the pool for 0.04 s, collapsed it by 37%, and ruined
// ~4 mm of bead. See patch_droplet_restart.py for the measurements.
//
// READ_IF_PRESENT: a case with no <time>/uniform/dropletState gets exactly the
// defaults above, so a run from t = 0 is bit-identical. Gate A and Gate B
// check that.
//
// Only the three ACCUMULATED values are stored. Everything else -- dropR_,
// dropV_, dropT_, dropVdot_, dropVdrop_, dropMass_, dropTau_, dropGap_ -- is
// re-read from caseParameters at startup and is identical on a restart.
// Persisting a value that is also read from a dictionary is how the two drift
// apart later.
IOdictionary dropletState
(
    IOobject
    (
        "dropletState",
        runTime.timeName(),
        "uniform",
        runTime,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    )
);

if (dropOn)
{
    const scalar owedPrev = dropOwed_;
    dropOwed_      = dropletState.lookupOrDefault<scalar>("dropOwed", dropOwed_);
    dropImpactEnd_ = dropletState.lookupOrDefault<scalar>
                     ("dropImpactEnd", dropImpactEnd_);
    dropNum_       = dropletState.lookupOrDefault<label>("dropNum", dropNum_);

    if (dropletState.found("dropOwed"))
    {
        Info<< "Droplet state: restored from "
            << runTime.timeName() << "/uniform/dropletState" << nl
            << "     dropOwed " << dropOwed_*1e9 << " mm^3 ("
            << dropOwed_/max(dropVdrop_, VSMALL) << " droplet)"
            << "   dropNum " << dropNum_ << endl;
    }
    else if (owedPrev == 0.0)
    {
        Info<< "Droplet state: fresh start (no uniform/dropletState)" << endl;
    }
}


'''

# ------------------------------------------------------------ laserbeamFoam.C
MAIN_OLD = '#include "updateDropletSource.H"'
MAIN_NEW = ('#include "updateDropletSource.H"\n'
            '        #include "syncDropletState.H"')

EDITS = [(CF, CF_ANCHOR, CF_ADD + CF_ANCHOR),
         (MAIN, MAIN_OLD, MAIN_NEW)]


def revert():
    n = 0
    for path, _, _ in EDITS:
        bak = path.with_suffix(path.suffix + SUFFIX)
        if bak.exists():
            shutil.copy2(bak, path); bak.unlink(); n += 1
            print(f"  reverted {path.name}")
    if NEW.exists():
        NEW.unlink(); print(f"  removed {NEW.name}")
    print(f"  {n} file(s) restored" if n else "  nothing to revert")


def apply():
    for path, old, _ in EDITS:
        if not path.exists():
            sys.exit(f"ERROR: {path} not found")
        text = path.read_text()
        if "dropletState" in text:
            sys.exit(f"ERROR: {path.name} looks already patched. --revert first.")
        if text.count(old) != 1:
            sys.exit(f"ERROR: {path.name}: expected one match for:\n"
                     f"-----\n{old}\n-----\nfound {text.count(old)}. "
                     f"Nothing modified.")

    for path, old, new in EDITS:
        shutil.copy2(path, path.with_suffix(path.suffix + SUFFIX))
        path.write_text(path.read_text().replace(old, new, 1))
        print(f"  patched {path.name}   (backup: {path.name}{SUFFIX})")
    NEW.write_text(NEW_SRC)
    print(f"  created {NEW.name}")

    print("""
  BUILD
    of2506 && cd ~/thesis/src/lbf-com && ./Allwmake > /tmp/b.log 2>&1
    grep -c 'error:' /tmp/b.log

  REGRESSION -- READ_IF_PRESENT, so a fresh case is bit-identical.
    Gate A  : centreline u(0) 0.0158753, Ha 15.810, err 1.3869e-02
              (run hartmann_check.py with  --sigma 1000 --mu 1)
    Gate B  : 1869.0 / 0.0504 / 3.4135 / 0.0552 / 12.7787
              cd cases/00-tutorials && python3 pool_stats.py <case>

  THE TEST THAT MATTERS -- restart across a seam and look for the gap.

    Take any WAAM case with droplets, run to ~0.30 s, stop, set
    startFrom latestTime, run to ~0.40 s, then:

        python3 ~/thesis/tools/check_droplet.py | grep -E '^ +[0-9]'

    BEFORE this patch the seam shows three consecutive writes with an
    identical 'deposited' value and a frozen 'drops' column. AFTER it, the
    cadence should be unbroken -- one droplet per ~1.5 writes straight
    through, no flat spot.

    The startup banner also reports what was restored:

        Droplet state: restored from 0.3/uniform/dropletState
             dropOwed 3.21 mm^3 (0.568 droplet)   dropNum 9

    If it says 'fresh start' on a restart, the dictionary is not being found
    and nothing has been fixed -- check that <time>/uniform/dropletState exists
    in each processor* directory after a parallel write.""")


if __name__ == "__main__":
    revert() if "--revert" in sys.argv else apply()
