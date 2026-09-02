#!/usr/bin/env python3
"""Wire the Phase 2 property set into generate_cases.py.

WHAT IT ADDS
    1. transportProperties is COPIED into every case's constant/ at generation
       time -- not symlinked. A symlink survives locally but arrives dangling
       when a case is tarred to the cluster, and worse, it means regenerating
       the properties silently rewrites the inputs of every case you already
       ran. A copy makes each finished case a self-contained record.

    2. case_meta.json records properties_sha, a 12-char fingerprint of the file
       that case actually used.

    3. The generator REFUSES to skip a case whose stored results were built
       with a different property set.

WHY (3) IS THE POINT
    The FROZEN columns are compared within one generation run. properties_sha
    cannot work that way -- every case built in a single run gets the same file
    by construction, so a within-run check would always pass and prove nothing.

    The real failure is across runs: regenerate the properties in month four,
    rebuild M04-M07, and compare them against M00 which still holds results
    from the old set. Nothing warns you, and Phase 2 measured the field effect
    at 0.14% against a 0.4-0.7% noise floor -- so a silent property change is
    more than large enough to manufacture a result.

    Overridable with the existing --ignore-audit flag.

USAGE
    python3 patch_generate.py            # apply
    python3 patch_generate.py --revert
"""
import sys, shutil
from pathlib import Path

TARGET = Path.home() / "thesis/cases/generate_cases.py"

IMPORT_ANCHOR = "import shutil"
IMPORT_ADD = "\nimport hashlib"

FN_ANCHOR = "def build_case(row, templates_dir, out_root, args, repo_root):"
FN_ADD = '''# --------------------------------------------------------------------------
# Phase 2 property set. One source of truth, copied into every case.

MATERIALS_DIR = Path.home() / "thesis" / "materials"
PROPERTIES_SRC = MATERIALS_DIR / "transportProperties"


def properties_sha():
    """12-char fingerprint of the property set currently on disk."""
    if not PROPERTIES_SRC.exists():
        die(f"missing {PROPERTIES_SRC}\\n"
            f"  build it with:  python3 {MATERIALS_DIR}/make_transport.py")
    return hashlib.sha256(PROPERTIES_SRC.read_bytes()).hexdigest()[:12]


def install_properties(case_dir):
    """Copy the material dictionary into a case. Returns its fingerprint.

    Copied, never symlinked: a symlink is stored as a link by tar and arrives
    dangling on the cluster, and it would make every past case silently point
    at whatever the properties are today rather than what they ran with.
    """
    dst = case_dir / "constant" / "transportProperties"
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(PROPERTIES_SRC, dst)     # copy2 preserves mtime
    return properties_sha()


def stored_properties_sha(case_dir):
    """The fingerprint recorded when this case was generated, or None."""
    meta = case_dir / "case_meta.json"
    if not meta.exists():
        return None
    try:
        return json.loads(meta.read_text()).get("properties_sha")
    except (ValueError, OSError):
        return None


'''

SKIP_OLD = '''    if dest.exists():
        if has_results(dest) and not args.force:
            print(f"      SKIPPED -- folder holds run output. Use --force to overwrite.")
            return "skipped"
        shutil.rmtree(dest)'''

SKIP_NEW = '''    if dest.exists():
        if has_results(dest) and not args.force:
            # This case keeps its old results while others are rebuilt. If it
            # ran on a different property set, comparing them is invalid --
            # and nothing else in the pipeline would notice.
            was = stored_properties_sha(dest)
            now = properties_sha()
            if was is not None and was != now:
                msg = (f"{case_id} holds results built with a DIFFERENT property "
                       f"set.\\n      recorded {was}, current {now}\\n"
                       f"      Comparing it against newly built cases is not valid. "
                       f"Re-run it with --force, or restore the old properties.")
                if getattr(args, "ignore_audit", False):
                    print(f"      *** {msg}\\n      (--ignore-audit: continuing)")
                else:
                    die(msg)
            print(f"      SKIPPED -- folder holds run output. Use --force to overwrite.")
            return "skipped"
        shutil.rmtree(dest)'''

CALL_OLD = "    unknown = substitute_tree(dest, values)"
CALL_NEW = ('''    props_sha = install_properties(dest)
    unknown = substitute_tree(dest, values)''')

META_OLD = '''        "csv_row": dict(row),
        "derived": {k: v for k, v in values.items() if k not in row},
    }'''
META_NEW = '''        "csv_row": dict(row),
        "derived": {k: v for k, v in values.items() if k not in row},
        "properties_sha": props_sha,
    }'''

EDITS = [(IMPORT_ANCHOR, IMPORT_ANCHOR + IMPORT_ADD),
         (FN_ANCHOR, FN_ADD + FN_ANCHOR),
         (SKIP_OLD, SKIP_NEW),
         (CALL_OLD, CALL_NEW),
         (META_OLD, META_NEW)]


def revert():
    bak = TARGET.with_suffix(".py.orig")
    if not bak.exists():
        print("  nothing to revert"); return
    shutil.copy2(bak, TARGET); bak.unlink()
    print(f"  reverted {TARGET.name}")


def apply():
    if not TARGET.exists():
        sys.exit(f"ERROR: {TARGET} not found")
    text = TARGET.read_text()
    if "properties_sha" in text:
        sys.exit("ERROR: already patched. Run --revert first to reapply.")
    if "import json" not in text:
        sys.exit("ERROR: expected 'import json' -- stored_properties_sha needs it. "
                 "Nothing modified.")
    for old, _ in EDITS:
        if text.count(old) != 1:
            sys.exit(f"ERROR: expected exactly one match for:\n-----\n{old}\n-----\n"
                     f"found {text.count(old)}. Nothing modified.")

    shutil.copy2(TARGET, TARGET.with_suffix(".py.orig"))
    for old, new in EDITS:
        text = text.replace(old, new, 1)
    TARGET.write_text(text)
    print(f"  patched {TARGET.name}   (backup: {TARGET.name}.orig)")
    print("\n  Check it with a dry run first:")
    print("    python3 ~/thesis/cases/generate_cases.py --dry-run --all")


if __name__ == "__main__":
    revert() if "--revert" in sys.argv else apply()
