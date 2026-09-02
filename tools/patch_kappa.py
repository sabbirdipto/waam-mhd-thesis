#!/usr/bin/env python3
"""Two-branch metal thermal conductivity for laserbeamFoam.

WHY
    Aluminium's thermal conductivity halves on melting -- 162 -> 81 W/mK across
    a 91 K window. laserbeamFoam reads ONE Polynomial<8> per phase, and no
    smooth degree-7 polynomial can represent that step: the best full-range fit
    is 18% wrong in the solid and 84% wrong at the liquidus. That error sits at
    the solidification front, where Darcy damping kills velocity and conduction
    locally dominates -- the worst possible place for it.

HOW
    Two optional dictionary entries, poly_kappa_solid and poly_kappa_liquid,
    are blended by LIQUID FRACTION inside the metal term, leaving the existing
    metal/gas blend by alpha1 untouched.

    The liquid fraction is derived from temperature (Zhao Eq. 8), not read from
    the epsilon1 field, because updateKappaCp.H is first included at
    createFields.H:282 -- before epsilon1, TSolidus and TLiquidus are declared.
    It is the same linear relation the solver uses at createFields.H:349, and
    it keeps kappa a pure function of T, matching the polynomial fits.

    Measured result: worst error anywhere falls from 84% to 2.6%, and the
    remaining 2.6% is confined to the mushy window where it reflects a
    difference of mixing rule, not a fitting failure.

SAFETY
    Both entries are OPTIONAL. When absent the solver takes a branch that
    evaluates the ORIGINAL expression verbatim, so the legacy path is
    BIT-identical, not merely algebraically equal. That matters: on a chaotic
    case a 1-ULP difference grows to tens of percent within milliseconds, which
    would make every regression test unable to prove anything.

USAGE
    python3 patch_kappa.py            # apply
    python3 patch_kappa.py --revert   # restore the .orig backups
"""
import sys, shutil
from pathlib import Path

SOLVER = Path.home() / "thesis/src/lbf-com/applications/solvers/laserbeamFoam"
CF, UK = SOLVER / "createFields.H", SOLVER / "updateKappaCp.H"

# ---------------------------------------------------------------- createFields
CF_ANCHOR = 'const Polynomial<8> polykappa_g(transportPropertiesGas.lookup("poly_kappa"));'

CF_ADD = '''

// --- Two-branch metal conductivity (WAAM-MHD thesis addition) ---------------
// Aluminium's thermal conductivity halves on melting (162 -> 81 W/mK across a
// 91 K window). One Polynomial<8> cannot represent that step. These optional
// entries give the solid and liquid branches separately; updateKappaCp.H
// blends them by liquid fraction, derived from T via Zhao Eq. 8.
// Absent from the dictionary, both fall back to poly_kappa, so every existing
// case runs unchanged.
auto readMetalKappa = [&](const word& key) -> Polynomial<8>
{
    return transportPropertiesMetal.found(key)
         ? Polynomial<8>(transportPropertiesMetal.lookup(key))
         : polykappa_m;
};
const Polynomial<8> polykappa_mS(readMetalKappa("poly_kappa_solid"));
const Polynomial<8> polykappa_mL(readMetalKappa("poly_kappa_liquid"));

// The METAL's own solidus/liquidus, as scalars. NOT the TSolidus/TLiquidus
// volScalarFields declared further down -- those are blended between metal and
// gas by alpha1, and they do not exist yet at the point updateKappaCp.H is
// first included (line 282, before their declarations at 285 and 299).
const scalar TsolidusMetal(readScalar(transportPropertiesMetal.lookup("Tsolidus")));
const scalar TliquidusMetal(readScalar(transportPropertiesMetal.lookup("Tliquidus")));
const bool twoBranchKappa = transportPropertiesMetal.found("poly_kappa_solid");
Info<< "Metal kappa: "
    << (transportPropertiesMetal.found("poly_kappa_solid")
        ? "two-branch (solid/liquid blended by liquid fraction)"
        : "single poly_kappa (legacy)")
    << endl;
'''

# ---------------------------------------------------------------- updateKappaCp
UK_OLD = """        kappaI[celli] =
            (
                max(min(alpha1I[celli], 1.0), 0.0)
               *(polykappa_m.value(TI[celli]))
            )
          + (
                (1.0 - max(min(alpha1I[celli], 1.0), 0.0))
               *polykappa_g.value(TI[celli])
            );"""

UK_NEW = """        // alpha1 is the METAL fraction (VOF). The LIQUID fraction is taken
        // from temperature via Zhao Eq. 8 -- the same linear relation the
        // solver itself uses at createFields.H:349 and inverts in TEqn.H:89.
        // The epsilon1 FIELD is not used here: this file is first included at
        // createFields.H:282, before epsilon1 exists. Deriving it from T also
        // keeps kappa a pure function of temperature, which is what the
        // polynomial fits assume.
        const scalar aM = max(min(alpha1I[celli], 1.0), 0.0);

        scalar kappaMetal;
        if (twoBranchKappa)
        {
            const scalar eL =
                max(min((TI[celli] - TsolidusMetal)
                       /(TliquidusMetal - TsolidusMetal), 1.0), 0.0);

            kappaMetal =
                eL*polykappa_mL.value(TI[celli])
              + (1.0 - eL)*polykappa_mS.value(TI[celli]);
        }
        else
        {
            // EXACTLY the original expression. Writing eL*k + (1-eL)*k with
            // identical branches is algebraically equal but differs by ~1 ULP
            // about 5% of the time, and a chaotic case amplifies that to tens
            // of percent. Branching keeps the legacy path BIT-identical, so a
            // regression test can distinguish "unchanged" from "barely changed".
            kappaMetal = polykappa_m.value(TI[celli]);
        }

        kappaI[celli] =
            (aM*kappaMetal)
          + ((1.0 - aM)*polykappa_g.value(TI[celli]));"""

EDITS = [(CF, [(CF_ANCHOR, CF_ANCHOR + CF_ADD)]),
         (UK, [(UK_OLD, UK_NEW)])]


def revert():
    n = 0
    for path, _ in EDITS:
        bak = path.with_suffix(path.suffix + ".orig")
        if bak.exists():
            shutil.copy2(bak, path); bak.unlink(); n += 1
            print(f"  reverted {path.name}")
    print(f"  {n} file(s) restored" if n else "  nothing to revert")


def apply():
    for path, subs in EDITS:
        if not path.exists():
            sys.exit(f"ERROR: {path} not found. Is the clone where you expect?")
        text = path.read_text()
        if "WAAM-MHD thesis addition" in text or "epsilon1I" in text:
            sys.exit(f"ERROR: {path.name} looks already patched. "
                     f"Run --revert first if you want to reapply.")
        for old, _ in subs:
            if text.count(old) != 1:
                sys.exit(f"ERROR: in {path.name}, expected exactly one match for:\n"
                         f"-----\n{old}\n-----\nfound {text.count(old)}.\n"
                         f"The upstream file has changed. Nothing was modified.")

    cf = CF.read_text().splitlines()
    anchor = next(i for i, l in enumerate(cf) if CF_ANCHOR in l)
    include = next((i for i, l in enumerate(cf) if '#include "updateKappaCp.H"' in l), None)
    if include is None:
        sys.exit("ERROR: no #include \"updateKappaCp.H\" found in createFields.H")
    if anchor >= include:
        sys.exit(f"ERROR: the anchor (line {anchor+1}) is not before the "
                 f"updateKappaCp.H include (line {include+1}). The new symbols "
                 f"would not be in scope. Nothing was modified.")
    print(f"  scope check OK: anchor line {anchor+1} < include line {include+1}")

    for path, subs in EDITS:                       # only now, after ALL checks
        shutil.copy2(path, path.with_suffix(path.suffix + ".orig"))
        text = path.read_text()
        for old, new in subs:
            text = text.replace(old, new, 1)
        path.write_text(text)
        print(f"  patched {path.name}   (backup: {path.name}.orig)")

    print("\n  Next:")
    print("    of2506")
    print("    cd ~/thesis/src/lbf-com && ./Allwmake")
    print("  Then re-run Gate A. It uses no poly_kappa_solid, so it must give")
    print("  the SAME 1.39% as before -- that is the regression test.")


if __name__ == "__main__":
    revert() if "--revert" in sys.argv else apply()
