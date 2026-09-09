#!/usr/bin/env python3
"""Make hartmann_check.py read the case instead of assuming unity.

THE PROBLEM
    hartmann_check.py defaults to --sigma 1 --mu 1. The hartmann1d case has
    elec_conductivity_air_ = 1000, hidden behind

        #include "propertiesReference"

    in constant/transportProperties, which the script never reads. So the
    default invocation computes

        Ha = B*L*sqrt(sigma/mu) = 0.99993 * 0.5 * sqrt(1/1) = 0.500

    instead of the true 15.810, compares a genuine plug profile against a
    parabola, and reports

        max |u/u0 - theory| : 6.2903e-01
        VERDICT : *** FAIL - investigate before building on this module ***

    on a case that actually passes at 1.3869e-02. The correct invocation,
    --sigma 1000 --mu 1, is recorded nowhere in the repository: Gate A was
    only reproducible by remembering a command line.

    This is the same failure pattern as sigma_e_gas_ratio (Phase 2) and the
    kappa boundary defect (Edit 8): a parameter that appears configured,
    silently isn't, and produces a confident wrong answer.

THE FIX
    1. Read constant/propertiesReference for elec_conductivity_air_ and
       magnetic_permeability_air_, locating the case by walking up from the
       csv path until a directory containing constant/ is found.
    2. Read nu and rho from constant/transportProperties for mu = rho*nu.
       If the phases disagree, REFUSE rather than pick one -- in hartmann1d
       'metal' is commented "Unused Phase" and only 'gas' flows, and guessing
       which is live is exactly the kind of assumption that caused this.
    3. Defaults become None. A value that can be neither read nor given is a
       hard error, not an assumed 1.
    4. Print the provenance of every number used.

USAGE
    python3 patch_hartmann_check.py            # apply
    python3 patch_hartmann_check.py --revert
"""
import shutil, sys
from pathlib import Path

TARGET = Path.home() / "thesis/cases/00-tutorials/hartmann1d/hartmann_check.py"
BAK = TARGET.with_suffix(".py.hbak")

OLD_ARGS = '''    p.add_argument("--sigma", type=float, default=1.0)
    p.add_argument("--mu", type=float, default=1.0)'''

NEW_ARGS = '''    p.add_argument("--sigma", type=float, default=None,
                   help="electrical conductivity. Default: read from the case")
    p.add_argument("--mu", type=float, default=None,
                   help="dynamic viscosity rho*nu. Default: read from the case")'''

OLD_MUMAG = '''    p.add_argument("--mu-mag", type=float, default=1.0, dest="mumag",
                   help="magnetic permeability, to turn sampled H into B (default 1)")'''

NEW_MUMAG = '''    p.add_argument("--mu-mag", type=float, default=None, dest="mumag",
                   help="magnetic permeability, to turn sampled H into B. "
                        "Default: read from the case")'''

OLD_PARSE = "    args = p.parse_args()"

NEW_PARSE = '''    args = p.parse_args()

    # ------------------------------------------------------------------
    # Resolve material properties FROM THE CASE. Assuming unity turned a
    # passing benchmark (1.3869e-02) into a 63% FAIL, because hartmann1d
    # hides elec_conductivity_air_ = 1000 behind an #include.
    # ------------------------------------------------------------------
    def _scalars(path):
        """Every 'key value;' pair in a foam dict. Comments stripped."""
        out = {}
        if not path.exists():
            return out
        for raw in path.read_text().splitlines():
            line = raw.split("//")[0].strip()
            if not line.endswith(";"):
                continue
            bits = line[:-1].split()
            if len(bits) == 2:
                try:
                    out[bits[0]] = float(bits[1])
                except ValueError:
                    pass
        return out

    def _all_values(path, key):
        """Every value assigned to `key`, so disagreement can be detected."""
        vals = []
        if not path.exists():
            return vals
        for raw in path.read_text().splitlines():
            line = raw.split("//")[0].strip()
            if not line.endswith(";"):
                continue
            bits = line[:-1].split()
            if len(bits) == 2 and bits[0] == key:
                try:
                    vals.append(float(bits[1]))
                except ValueError:
                    pass
        return vals

    case = args.csv.resolve().parent
    for _ in range(6):
        if (case / "constant").is_dir():
            break
        case = case.parent
    else:
        case = None

    src = {}
    if case is not None:
        ref = _scalars(case / "constant" / "propertiesReference")
        tp = case / "constant" / "transportProperties"

        if args.sigma is None and "elec_conductivity_air_" in ref:
            args.sigma = ref["elec_conductivity_air_"]
            src["sigma"] = "constant/propertiesReference"
        if args.mumag is None and "magnetic_permeability_air_" in ref:
            args.mumag = ref["magnetic_permeability_air_"]
            src["mu_mag"] = "constant/propertiesReference"

        if args.mu is None:
            nus, rhos = _all_values(tp, "nu"), _all_values(tp, "rho")
            if nus and rhos and len(set(nus)) == 1 and len(set(rhos)) == 1:
                args.mu = rhos[0] * nus[0]
                src["mu"] = "constant/transportProperties (rho*nu)"
            elif nus and rhos:
                sys.exit(
                    "ERROR: the phases in constant/transportProperties disagree "
                    f"(nu={sorted(set(nus))}, rho={sorted(set(rhos))}).\\n"
                    "  Which phase is flowing is not something this script will "
                    "guess -- in hartmann1d 'metal' is marked \\"Unused Phase\\".\\n"
                    "  Pass --mu explicitly.")

    missing = [n for n, v in (("--sigma", args.sigma), ("--mu", args.mu)) if v is None]
    if missing:
        sys.exit(
            f"ERROR: could not determine {', '.join(missing)} from the case"
            + (f" at {case}" if case else " (no constant/ directory found)")
            + ".\\n  Refusing to assume 1: doing so scores a passing Hartmann "
              "benchmark\\n  as a 63% FAIL. Pass the value(s) explicitly.")
    if args.mumag is None:
        args.mumag = 1.0
        src["mu_mag"] = "assumed 1 (not found in the case)"

    for name, val in (("sigma", args.sigma), ("mu", args.mu),
                      ("mu_mag", args.mumag)):
        print(f"  {name:<16} : {val:<12g} [{src.get(name, 'given on command line')}]")
'''


def revert():
    if not BAK.exists():
        print("  nothing to revert"); return
    shutil.copy2(BAK, TARGET); BAK.unlink()
    print(f"  reverted {TARGET.name}")


def apply():
    if not TARGET.exists():
        sys.exit(f"ERROR: {TARGET} not found")
    text = TARGET.read_text()

    if "Resolve material properties FROM THE CASE" in text:
        sys.exit("ERROR: already patched. Run --revert first to reapply.")

    for label, old in (("--sigma/--mu block", OLD_ARGS),
                       ("--mu-mag block", OLD_MUMAG),
                       ("args = p.parse_args()", OLD_PARSE)):
        n = text.count(old)
        if n != 1:
            sys.exit(f"ERROR: expected exactly one {label}, found {n}. "
                     f"Nothing was modified.")

    shutil.copy2(TARGET, BAK)
    for old, new in ((OLD_ARGS, NEW_ARGS), (OLD_MUMAG, NEW_MUMAG),
                     (OLD_PARSE, NEW_PARSE)):
        text = text.replace(old, new, 1)
    TARGET.write_text(text)
    print(f"  patched {TARGET.name}   (backup: {BAK.name})")
    print("""
  Verify -- all three must hold:

    cd ~/thesis/cases/00-tutorials/hartmann1d

    # 1. reads the case, no flags, and still scores 1.3869e-02
    python3 hartmann_check.py postProcessing/sample/10/line_centreProfile_H_U.csv \\
        | grep -E 'sigma|mu |mu_mag|Hartmann|max \\||VERDICT'

    # 2. explicit flags still work and agree
    python3 hartmann_check.py postProcessing/sample/10/line_centreProfile_H_U.csv \\
        --sigma 1000 --mu 1 | grep -E 'Hartmann|max \\|'

    # 3. refuses rather than assuming, when the case cannot be found
    cp postProcessing/sample/10/line_centreProfile_H_U.csv /tmp/orphan.csv
    python3 hartmann_check.py /tmp/orphan.csv ; echo "exit $?" """)


if __name__ == "__main__":
    revert() if "--revert" in sys.argv else apply()
