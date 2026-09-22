#!/usr/bin/env python3
"""patch_shieldgas.py -- Zhao et al. (2021) Eq. 28 shielding gas, as boundary
conditions only. No solver change, so Gates A and B are untouched.

WHAT THIS MODELS
    20 L/min of argon leaves the annulus between the wire (Rw = 0.6 mm) and the
    nozzle (Rn = 10 mm) as fully developed laminar flow. For that geometry the
    Navier-Stokes momentum equation collapses to

        (mu/r) d/dr( r du/dr ) = dp/dz,

    and integrating twice with no-slip at both walls (u(Rw) = u(Rn) = 0) gives

        u(r) = K [ Rn^2 - r^2 + (Rn^2 - Rw^2) ln(r/Rn) / ln(Rn/Rw) ].

    K is NOT free. Requiring the profile to carry the measured flow rate,
    integral of u.2.pi.r dr from Rw to Rn = Q, fixes

        K = (2Q/pi) / [ Rn^4 - Rw^4 - (Rn^2 - Rw^2)^2 / ln(Rn/Rw) ].

    With Q = 20 L/min: mean 1.0649 m/s, peak 1.6933 m/s at r = 4.208 mm, and
    Re = rho.U.Dh/mu = 1431 on the annulus hydraulic diameter, so the laminar
    assumption Zhao makes is self-consistent.

WHAT IT CHANGES
    0.orig/U       inlet (the top patch): pressureInletOutletVelocity
                   -> codedFixedValue carrying the profile above, centred on
                      the moving torch, Vz = 0 outside the nozzle and inside
                      the wire, exactly as the paper states.
    0.orig/p_rgh   inlet: totalPressure -> fixedFluxPressure. Mandatory: once
                   the velocity is prescribed on a patch, the pressure equation
                   needs a gradient consistent with that prescribed flux.
    constant/caseParameters   gains gas_Q_m3_s, gas_Rn_m, gas_Rw_m. The BC
                   reads arc_x0_m, arc_y0_m and travel_speed_m_s from the same
                   file, so the nozzle and the arc can never disagree.
    system/controlDict   gains a surfaceFieldValue function object summing phi
                   over the inlet, so every write records the volumetric flow
                   actually delivered (should be -Q while the nozzle is fully
                   inside the domain; see the clipping window printed below).

    0.orig/T and 0.orig/alpha.metal are NOT touched: inletOutlet with 300 K and
    0 already behaves as a Dirichlet condition wherever gas flows in and as
    zero-gradient where it does not.

USAGE
    python3 patch_shieldgas.py --clone <src-case> <dst-case>   clone, then patch
    python3 patch_shieldgas.py <case>                          patch in place
    python3 patch_shieldgas.py --check <case>                  report only
    python3 patch_shieldgas.py --revert <case>                 restore .gbak

Every edit writes a <file>.gbak first and refuses to run twice.
"""
import argparse
import math
import os
import re
import shutil
import sys

Q_DEFAULT = 20.0 / 60000.0      # 20 L/min in m^3/s
RN_DEFAULT = 0.010              # nozzle inner radius [m]   Zhao Table 3
RW_DEFAULT = 0.0006             # wire radius [m]           Zhao Table 2, d=1.2 mm
TAG = "shieldGasEq28"

# --------------------------------------------------------------- the BC entry
U_ENTRY = """    inlet
    {{
        // Zhao et al. (2021) Eq. 28 -- shielding gas.
        // Laminar annular Poiseuille profile between the wire (gas_Rw_m) and
        // the nozzle (gas_Rn_m), normalised to carry gas_Q_m3_s, centred on
        // the moving torch. Vz = 0 outside the nozzle and inside the wire.
        // Patched in by tools/patch_shieldgas.py -- {tag}
        type            codedFixedValue;
        value           uniform (0 0 0);
        name            {tag};

        codeInclude
        #{{
            #include "IOdictionary.H"
            #include "mathematicalConstants.H"
        #}};

        code
        #{{
            const fvPatch& pp = patch();
            const scalar t = this->db().time().value();

            // caseParameters never changes during a run, and updateCoeffs()
            // is called every timestep, so read it once per process.
            static bool init = false;
            static scalar x0 = 0, y0 = 0, V = 0, Rn = 0, Rw = 0, Q = 0;
            static scalar K = 0, Ln = 0;

            if (!init)
            {{
                IOdictionary cp
                (
                    IOobject
                    (
                        "caseParameters",
                        this->db().time().constant(),
                        this->db(),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE,
                        false
                    )
                );

                x0 = cp.get<scalar>("arc_x0_m");
                y0 = cp.getOrDefault<scalar>("arc_y0_m", 0.0);
                V  = cp.get<scalar>("travel_speed_m_s");
                Rn = cp.get<scalar>("gas_Rn_m");
                Rw = cp.get<scalar>("gas_Rw_m");
                Q  = cp.get<scalar>("gas_Q_m3_s");

                Ln = Foam::log(Rn/Rw);
                const scalar den =
                    Foam::pow4(Rn) - Foam::pow4(Rw)
                  - Foam::sqr(Foam::sqr(Rn) - Foam::sqr(Rw))/Ln;
                K = 2.0*Q/(Foam::constant::mathematical::pi*den);

                const scalar rpk =
                    Foam::sqrt((Foam::sqr(Rn) - Foam::sqr(Rw))/(2.0*Ln));
                const scalar upk =
                    K*( Foam::sqr(Rn) - Foam::sqr(rpk)
                      + (Foam::sqr(Rn) - Foam::sqr(Rw))
                       *Foam::log(rpk/Rn)/Ln );

                Info<< "{tag}: " << Q*6.0e4 << " L/min through Rw "
                    << Rw*1e3 << " mm to Rn " << Rn*1e3
                    << " mm, nozzle centre x = " << x0 << " + " << V
                    << "*t, peak " << upk << " m/s at r = " << rpk*1e3
                    << " mm, mean "
                    << Q/(Foam::constant::mathematical::pi
                         *(Foam::sqr(Rn) - Foam::sqr(Rw)))
                    << " m/s" << endl;

                init = true;
            }}

            const scalar xc = x0 + V*t;
            const vectorField& Cf = pp.Cf();
            vectorField& Up = *this;

            forAll(Up, i)
            {{
                const scalar dx = Cf[i].x() - xc;
                const scalar dy = Cf[i].y() - y0;
                const scalar r = Foam::sqrt(dx*dx + dy*dy);

                scalar u = 0.0;
                if (r > Rw && r < Rn)
                {{
                    u = K*( Foam::sqr(Rn) - Foam::sqr(r)
                          + (Foam::sqr(Rn) - Foam::sqr(Rw))
                           *Foam::log(r/Rn)/Ln );
                    u = Foam::max(u, scalar(0));
                }}

                Up[i] = vector(0, 0, -u);
            }}
        #}};
    }}
"""

PRGH_ENTRY = """    inlet
    {
        // Velocity is prescribed on this patch by the Eq. 28 jet, so the
        // pressure equation needs the flux-consistent gradient here, not a
        // pressure value. -- tools/patch_shieldgas.py
        type            fixedFluxPressure;
        value           uniform 0;
    }
"""

PARAMS_BLOCK = """
// --- Shielding gas, Zhao Eq. 28 (patch_shieldgas.py) ---
// Laminar annular Poiseuille flow between the wire and the nozzle, carrying
// gas_Q_m3_s and centred on the torch (arc_x0_m + travel_speed_m_s*t).
// The profile is imposed by 0.orig/U, patch 'inlet'. Delete these three
// entries and the BC will fail on a missing key rather than run silently
// with a default.
gas_Q_m3_s        {q:.6e};    // {qlpm:g} L/min, Zhao Table 2
gas_Rn_m          {rn:g};        // nozzle inner radius, Zhao Table 3
gas_Rw_m          {rw:g};       // wire radius (1.2 mm wire, Zhao Table 2)
"""

FUNC_ENTRY = """    gasFlux
    {
        // Volumetric flow through the top patch, every write. While the
        // nozzle is fully inside the domain this must equal -gas_Q_m3_s;
        // it drops when the nozzle overhangs an x boundary.
        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");
        fields          (phi);
        operation       sum;
        regionType      patch;
        name            inlet;
        writeControl    writeTime;
        writeFields     false;
        log             true;
    }
"""

SKIP_DIRS = ("processor", "postProcessing", "dynamicCode", "VTK")


# ------------------------------------------------------------------ utilities
def die(msg):
    sys.exit("ERROR: " + msg)


def read(p):
    with open(p) as fh:
        return fh.read()


def write(p, s):
    with open(p, "w") as fh:
        fh.write(s)


def backup(p):
    b = p + ".gbak"
    if os.path.exists(b):
        die(f"{b} already exists -- this case looks patched. Use --revert first.")
    shutil.copy2(p, b)


def patch_block(text, patch_name):
    """Return (start, end) of the '<patch_name> { ... }' entry, by brace count."""
    m = re.search(r"^\s*" + re.escape(patch_name) + r"\s*$", text, re.M)
    if not m:
        m = re.search(r"^\s*" + re.escape(patch_name) + r"\s*\{", text, re.M)
        if not m:
            return None
    i = text.index("{", m.start())
    depth, j = 0, i
    while j < len(text):
        if text[j] == "{":
            depth += 1
        elif text[j] == "}":
            depth -= 1
            if depth == 0:
                break
        j += 1
    if depth != 0:
        die(f"unbalanced braces around patch '{patch_name}'")
    return m.start(), j + 1


def dict_value(text, key):
    m = re.search(r"^\s*" + re.escape(key) + r"\s+([^;]+);", text, re.M)
    return m.group(1).strip() if m else None


def vertices(bm):
    m = re.search(r"vertices\s*\((.*?)\)\s*;", bm, re.S)
    if not m:
        die("no vertices block in system/blockMeshDict")
    pts = re.findall(r"\(\s*([-\d.eE+]+)\s+([-\d.eE+]+)\s+([-\d.eE+]+)\s*\)", m.group(1))
    return [tuple(float(v) for v in p) for p in pts]


# --------------------------------------------------------------------- checks
def inspect(case, q, rn, rw):
    """Verify the case is what this patch assumes, and print the geometry."""
    need = ["constant/caseParameters", "system/blockMeshDict", "system/controlDict",
            "0.orig/U", "0.orig/p_rgh", "0.orig/T", "0.orig/alpha.metal"]
    for f in need:
        if not os.path.isfile(os.path.join(case, f)):
            die(f"{case}: missing {f}")

    cp = read(os.path.join(case, "constant/caseParameters"))
    if dict_value(cp, "goldak_af_m") is None:
        die("constant/caseParameters has no goldak_af_m -- not a WAAM arc case")
    x0 = float(dict_value(cp, "arc_x0_m"))
    y0 = float(dict_value(cp, "arc_y0_m") or 0.0)
    V = float(dict_value(cp, "travel_speed_m_s"))

    bm = read(os.path.join(case, "system/blockMeshDict"))
    pts = vertices(bm)
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]; zs = [p[2] for p in pts]
    x0d, x1d, y0d, y1d, z0d, z1d = min(xs), max(xs), min(ys), max(ys), min(zs), max(zs)
    if "inlet" not in bm:
        die("system/blockMeshDict has no patch named 'inlet'")
    m = re.search(r"inlet\s*\{[^}]*faces\s*\(\s*\(([^)]*)\)", bm)
    if m:
        ids = [int(v) for v in m.group(1).split()]
        if not all(abs(pts[i][2] - z1d) < 1e-12 for i in ids):
            die("patch 'inlet' is not the top (z = zmax) face -- check blockMeshDict")

    u = read(os.path.join(case, "0.orig/U"))
    if TAG in u:
        die("0.orig/U already carries the shielding-gas BC")
    blk = patch_block(u, "inlet")
    if blk is None or "pressureInletOutletVelocity" not in u[blk[0]:blk[1]]:
        die("0.orig/U patch 'inlet' is not pressureInletOutletVelocity -- "
            "this case is not the expected starting point")

    p = read(os.path.join(case, "0.orig/p_rgh"))
    blk = patch_block(p, "inlet")
    if blk is None or "totalPressure" not in p[blk[0]:blk[1]]:
        die("0.orig/p_rgh patch 'inlet' is not totalPressure")

    # geometry consequences
    area = math.pi * (rn**2 - rw**2)
    umean = q / area
    Ln = math.log(rn / rw)
    den = rn**4 - rw**4 - (rn**2 - rw**2) ** 2 / Ln
    K = 2 * q / (math.pi * den)
    rpk = math.sqrt((rn**2 - rw**2) / (2 * Ln))
    upk = K * (rn**2 - rpk**2 + (rn**2 - rw**2) * math.log(rpk / rn) / Ln)

    t_in = (x0d + rn - x0) / V          # nozzle fully inside from this time
    t_out = (x1d - rn - x0) / V         # ... until this time
    gas_h = z1d - 0.0                   # gas height above the plate top (z = 0)
    a_out = 2 * (y1d - y0d) * gas_h     # two x-face outlets, gas part only

    print(f"case             {case}")
    print(f"  domain mm      x {x0d*1e3:.1f}..{x1d*1e3:.1f}  y {y0d*1e3:.1f}..{y1d*1e3:.1f}"
          f"  z {z0d*1e3:.1f}..{z1d*1e3:.1f}")
    print(f"  torch          x = {x0*1e3:.1f} mm + {V*1e3:.1f} mm/s * t,  y = {y0*1e3:.1f} mm")
    print(f"  jet            Q {q*6e4:.3f} L/min, mean {umean:.4f} m/s, "
          f"peak {upk:.4f} m/s at r = {rpk*1e3:.3f} mm")
    print(f"  nozzle fully inside the domain for t = {max(t_in,0):.3f} .. {t_out:.3f} s")
    print(f"  outlet area    {a_out*1e6:.0f} mm^2 of gas -> mean exit {q/a_out:.3f} m/s")

    tp = os.path.join(case, "constant/transportProperties")
    if os.path.isfile(tp):
        t = read(tp)
        g = patch_block(t, "gas") or patch_block(t, "air")
        rho = dict_value(t[g[0]:g[1]], "rho") if g else None
        if rho:
            try:
                rho = float(re.findall(r"[-\d.eE+]+", rho)[-1])
                print(f"  gas density    {rho:g} kg/m^3 -> stagnation 0.5*rho*u_peak^2 = "
                      f"{0.5*rho*upk**2:.3f} Pa   (arc pressure peak 228 Pa)")
            except (ValueError, IndexError):
                pass
    return dict(x0=x0, V=V, x0d=x0d, x1d=x1d, t_in=t_in, t_out=t_out)


# ---------------------------------------------------------------------- edits
def apply(case, q, rn, rw):
    info = inspect(case, q, rn, rw)

    fU = os.path.join(case, "0.orig/U")
    u = read(fU)
    s, e = patch_block(u, "inlet")
    backup(fU)
    write(fU, u[:s] + U_ENTRY.format(tag=TAG) + u[e + 1:].lstrip("\n"))

    fP = os.path.join(case, "0.orig/p_rgh")
    p = read(fP)
    s, e = patch_block(p, "inlet")
    backup(fP)
    write(fP, p[:s] + PRGH_ENTRY + p[e + 1:].lstrip("\n"))

    fC = os.path.join(case, "constant/caseParameters")
    c = read(fC)
    if "gas_Q_m3_s" in c:
        die("constant/caseParameters already has gas_Q_m3_s")
    backup(fC)
    write(fC, c.rstrip() + "\n" + PARAMS_BLOCK.format(q=q, qlpm=q * 6e4, rn=rn, rw=rw))

    fD = os.path.join(case, "system/controlDict")
    d = read(fD)
    backup(fD)
    if re.search(r"^functions\s*$|^functions\s*\{", d, re.M):
        s, e = patch_block(d, "functions")
        inner = d[s:e]
        close = inner.rindex("}")
        d = d[:s] + inner[:close] + FUNC_ENTRY + inner[close:] + d[e:]
    else:
        blk = "\nfunctions\n{\n" + FUNC_ENTRY + "}\n"
        m = re.search(r"\n//\s\*+\s//\s*$", d)          # trailing separator
        d = (d[:m.start()] + "\n" + blk + d[m.start():]) if m else (d.rstrip() + "\n" + blk)
    write(fD, d)

    print("\npatched:")
    for f in ("0.orig/U", "0.orig/p_rgh", "constant/caseParameters", "system/controlDict"):
        print(f"  {f}   (backup {f}.gbak)")
    print("  0.orig/T and 0.orig/alpha.metal unchanged (inletOutlet 300 K / 0 already correct)")
    return info


def revert(case):
    n = 0
    for f in ("0.orig/U", "0.orig/p_rgh", "constant/caseParameters", "system/controlDict"):
        b = os.path.join(case, f + ".gbak")
        if os.path.exists(b):
            shutil.move(b, os.path.join(case, f))
            print(f"  restored {f}")
            n += 1
    if not n:
        die("no .gbak files found")


def clone(src, dst):
    if os.path.exists(dst):
        die(f"{dst} already exists")
    src = src.rstrip("/")
    os.makedirs(dst)
    for name in sorted(os.listdir(src)):
        p = os.path.join(src, name)
        if name.startswith(SKIP_DIRS) or name.startswith("log") or name.endswith(".gbak"):
            continue
        if os.path.isdir(p):
            if name not in ("0.orig", "constant", "system"):
                continue                      # skip time directories
            shutil.copytree(p, os.path.join(dst, name),
                            ignore=shutil.ignore_patterns("polyMesh", "*.gbak"))
        else:
            shutil.copy2(p, os.path.join(dst, name))
    print(f"cloned {src} -> {dst} (0.orig, constant, system and scripts only)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", nargs="*")
    ap.add_argument("--clone", action="store_true", help="first argument is the source case")
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--revert", action="store_true")
    ap.add_argument("--q_lpm", type=float, default=20.0)
    ap.add_argument("--Rn_mm", type=float, default=RN_DEFAULT * 1e3)
    ap.add_argument("--Rw_mm", type=float, default=RW_DEFAULT * 1e3)
    a = ap.parse_args()

    q = a.q_lpm / 60000.0
    rn, rw = a.Rn_mm * 1e-3, a.Rw_mm * 1e-3

    if a.clone:
        if len(a.case) != 2:
            die("--clone needs a source and a destination case")
        clone(a.case[0], a.case[1])
        apply(a.case[1], q, rn, rw)
    elif a.revert:
        if len(a.case) != 1:
            die("--revert needs one case")
        revert(a.case[0])
    elif a.check:
        if len(a.case) != 1:
            die("--check needs one case")
        inspect(a.case[0], q, rn, rw)
    else:
        if len(a.case) != 1:
            die("give one case directory, or --clone SRC DST")
        apply(a.case[0], q, rn, rw)


if __name__ == "__main__":
    main()
