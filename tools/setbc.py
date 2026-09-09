#!/usr/bin/env python3
"""Set the thermal condition on ONE patch of 0/T. Nothing else.

WHY THIS EXISTS AS A SCRIPT AND NOT A sed ONE-LINER
    The domain-independence test is only meaningful if exactly one thing
    changes between run A and run B. A hand edit cannot demonstrate that: it
    leaves no record of what the file looked like before, and if a stray
    keystroke also touches 'substrateBottom' the comparison silently stops
    measuring what it claims to measure. That is the same failure that made
    the first attempt at this test worthless -- 'walls' was the bottom and the
    far side at once, so varying it varied a physical modelling choice and a
    numerical truncation together, and both effects landed in one number.

    So this script does three things a manual edit does not:
      1. edits ONE named patch and refuses if that patch is not unique
      2. verifies the OTHER patch is still zeroGradient, and stops if it is
         not -- the control must actually be held
      3. prints the resulting boundaryField in full, so the run log contains
         proof of what was solved

WHAT THE TWO CONDITIONS MEAN, physically
    zeroGradient    dT/dn = 0 at the face. By Fourier's law q = -k dT/dn, the
                    heat flux through that face is exactly zero: adiabatic.
                    For a TRUNCATION patch this says "the plate continues and
                    carries no net heat across this imaginary line", which is
                    true only while no heat has arrived there yet.

    fixedValue 300  the face is pinned at ambient. The near-wall cell then
                    sees a gradient (T_cell - 300)/(dn/2) and loses
                    q = k (T_cell - 300)/(dn/2) through it: an infinite heat
                    sink. This says "beyond this line is a body so large it
                    cannot warm up".

    These are the two OPPOSITE extremes -- no heat leaves, versus as much heat
    leaves as the gradient allows. Any real condition lies between them. That
    is exactly why they make a good test: if the domain is large enough that
    the arc has not yet heated the far side, the near-wall cell is still AT
    300 K, the gradient is zero either way, and both conditions pass zero heat.
    The two runs must then agree to round-off. If they disagree, heat HAS
    reached the boundary and the answer depends on an arbitrary choice, which
    means the domain is too small -- not that one condition is wrong.

USAGE, from inside the case directory
    python3 setbc.py farField zeroGradient
    python3 setbc.py farField fixedValue 300
    python3 setbc.py --show

WHICH FILE IT EDITS, and why that is not obvious
    The builder writes the fields to initial/, and Allrun does

        rm -rf 0 ; cp -r initial 0

    every time it runs. So initial/T is the source of truth and 0/T is a
    copy that Allrun overwrites. Editing 0/T and then calling Allrun would
    throw the edit away WITHOUT SAYING SO -- the run would proceed with the
    original condition and the two "different" runs would be the same run.

    This script therefore writes every copy that exists, initial/T first,
    so the case is consistent whichever way it is started.
"""
import re, sys
from pathlib import Path

# initial/ first: it is the one Allrun copies FROM.
TARGETS = [p for p in (Path("initial/T"), Path("0/T")) if p.exists()]
CONTROL = {"substrateBottom": "zeroGradient"}   # patches that must NOT move


def block(text, patch):
    """The (start, end, body) of one patch entry in boundaryField.

    Brace-counted rather than regexed to a closing '}', because a patch body
    can itself contain braces and a lazy regex would stop at the first one.
    """
    m = re.search(r"^\s*" + re.escape(patch) + r"\s*$", text, re.M)
    if not m:
        return None
    i = text.index("{", m.end())
    depth, j = 1, i + 1
    while depth and j < len(text):
        if text[j] == "{":
            depth += 1
        elif text[j] == "}":
            depth -= 1
        j += 1
    return m.start(), j, text[i + 1:j - 1]


def show(text):
    k = text.index("boundaryField")
    print(text[k:].split("// *")[0].rstrip())


def main():
    if not TARGETS:
        sys.exit("ERROR: neither initial/T nor 0/T is here. "
                 "Run this from inside the case directory.")
    text = TARGETS[0].read_text()

    if "--show" in sys.argv:
        print(f"  {TARGETS[0]}")
        show(text)
        # If both exist and disagree, say so -- that is precisely the state
        # in which a run solves something other than what you just read.
        for other in TARGETS[1:]:
            if other.read_text() != text:
                print(f"\n  WARNING: {other} DIFFERS from {TARGETS[0]}.")
                print(f"  {other}")
                show(other.read_text())
                print("\n  Allrun overwrites 0/ from initial/, so the run will "
                      "use initial/T.\n  Re-run this script with an explicit "
                      "condition to bring them back into line.")
        return

    if len(sys.argv) < 3:
        sys.exit(__doc__.split("USAGE")[1])

    patch, kind = sys.argv[1], sys.argv[2]
    val = sys.argv[3] if len(sys.argv) > 3 else None

    if kind == "fixedValue" and val is None:
        sys.exit("ERROR: fixedValue needs a temperature, e.g. "
                 "'python3 setbc.py farField fixedValue 300'")
    if kind not in ("zeroGradient", "fixedValue"):
        sys.exit(f"ERROR: unknown type '{kind}'. This script writes only "
                 f"zeroGradient or fixedValue -- the two extremes the test "
                 f"needs. Anything else, edit 0/T by hand and say so.")

    # ---- the target patch must exist, exactly once ----
    if len(re.findall(r"^\s*" + re.escape(patch) + r"\s*$", text, re.M)) != 1:
        found = ", ".join(re.findall(r"^    (\w+)$", text, re.M))
        sys.exit(f"ERROR: '{patch}' does not appear exactly once in "
                 f"{TARGETS[0]}.\n  Patches found: {found}\n"
                 f"  If you expected substrateBottom/farField and see 'walls',"
                 f" this case predates the patch split -- rebuild it.")

    # ---- the control patches must genuinely be held ----
    for name, want in CONTROL.items():
        if name == patch:
            continue
        b = block(text, name)
        if b is None:
            sys.exit(f"ERROR: control patch '{name}' is missing from "
                     f"{TARGETS[0]}. Refusing to run a one-variable test on "
                     f"a case that does not have the variable held.")
        if want not in b[2]:
            sys.exit(f"ERROR: control patch '{name}' is not '{want}':\n"
                     f"{b[2].rstrip()}\n"
                     f"  Two things would change between runs and the "
                     f"comparison would mean nothing. Fix it, then re-run.")

    s, e, old = block(text, patch)
    body = ("        type            zeroGradient;\n" if kind == "zeroGradient"
            else f"        type            fixedValue;\n"
                 f"        value           uniform {val};\n")
    new = f"    {patch}\n    {{\n{body}    }}"   # no trailing \n: text[e:] has it

    was = old.split(';')[0].split()[-1]
    for path in TARGETS:
        t = path.read_text()
        bs, be, _ = block(t, patch)
        path.write_text(t[:bs] + new + t[be:])
    print(f"  {patch}: {was} -> {kind}{' ' + val if val else ''}"
          f"   written to {', '.join(str(p) for p in TARGETS)}")

    print()
    show(TARGETS[0].read_text())
    print("\n  substrateBottom is the control and is held at zeroGradient in "
          "both runs.\n  Only farField differs, so any difference in the "
          "result is attributable to it.")


if __name__ == "__main__":
    main()
