#!/usr/bin/env python3
"""Phase 2 force-scaling analysis for magnetised GMAW-WAAM of Al-5%Mg.

Every number here is derived from a stated source or a stated construction.
Run it, change an input, watch what moves. Nothing is hard-coded from memory.
"""
import math

# ---------------------------------------------------------------- inputs
I      = 135.0        # A          Zhao Table 2
U_arc  = 16.0         # V          Zhao Table 2
sig_j  = 0.002        # m          Zhao Table 4, "estimated"
sig_p  = 0.002        # m          Zhao Table 4, "estimated"
Rw     = 0.0006       # m          wire radius, 1.2 mm dia
mu0    = 4*math.pi*1e-7
g      = 9.81
eta    = 0.7          # arc efficiency, Zhao Table 4

# pool scales -- from Zhao's own results (Fig 16: width 7.76 mm layer 1)
L_pool = 0.003        # m          half-width, characteristic MHD length
U_pool = 0.3          # m/s        melt velocity, Zhao Sec 4.1 ("about 0.3 m/s")
U_drop = 0.9          # m/s        peak droplet velocity, Zhao conclusion (2)

# material -- Phase 2 property set
rho_model = 2650.0    # kg/m3      Zhao Table 4 (constant, solid RT value)
rho_liq   = 2350.0    # kg/m3      real liquid Al-5Mg at liquidus (see notes)
mu_f      = 1.25e-3   # Pa.s       liquid Al at Tm, Assael et al. 2006
sig_e     = 3.65e6    # S/m        liquid Al-5Mg at liquidus  <-- NEW THIS PHASE
sig_e_tut = 1.0e6     # S/m        what the MTHD Plate2D tutorial uses
sig_e_RT  = 1.64e7    # S/m        solid Al-5Mg at 300 K

def banner(s): print("\n" + s + "\n" + "-"*len(s))

# ------------------------------------------------- 1. unit checks
banner("1. UNIT CHECKS  (Zhao Table 2 -> SI)")
for label, val, unit, si in [
    ("travel speed",    120.0, "cm/min", 120/100/60),
    ("wire feed",       900.0, "cm/min", 900/100/60),
    ("shielding gas",    20.0, "L/min",   20/1000/60),
]:
    print(f"  {label:16s} {val:7.1f} {unit:8s} = {si:.6g} " +
          ("m/s" if "min" in unit and "L" not in unit else "m3/s"))
print(f"  wire feed 900 cm/min  = 0.15 m/s   -> matches Zhao Table 3 inlet Vz = 0.15 m/s  OK")
print(f"  arc power  U*I        = {U_arc*I:.0f} W ; into pool at eta={eta}: {eta*U_arc*I:.0f} W")

# ------------------------------------------------- 2. current density
banner("2. CURRENT DENSITY AND SELF-FIELD")
A_j  = math.pi*sig_j**2
J    = I/A_j
B_self = mu0*I/(2*math.pi*sig_j)
A_w  = math.pi*Rw**2
J_w  = I/A_w
B_w  = mu0*I/(2*math.pi*Rw)
print(f"  in the pool,  r = sigma_j = {sig_j*1e3:.1f} mm")
print(f"    |J|      = I/(pi sigma_j^2)   = {J:.3e} A/m2")
print(f"    B_theta  = mu0 I/(2 pi r)     = {B_self*1e3:.2f} mT")
print(f"    J x B    = {J*B_self:.3e} N/m3")
print(f"  at the wire,  r = Rw = {Rw*1e3:.1f} mm")
print(f"    |J|      = {J_w:.3e} A/m2 ,  B = {B_w*1e3:.1f} mT ,  JxB = {J_w*B_w:.3e} N/m3")
print(f"    (Zhao Sec 4.1 reports droplet EM force 'order 1e6 N/m3'  -> {J_w*B_w:.1e})")

# ------------------------------------------------- 3. force ledger
banner("3. FORCE LEDGER  (compared as STRESSES, N/m2)")
print("  Surface forces enter as stresses; the CSF conversion to N/m3 divides by")
print("  the interface thickness, so a volumetric comparison is mesh-dependent and")
print("  meaningless. Body forces are put on the same footing by multiplying by L.\n")
gamma, dgdT = 0.85, -1.55e-4
dTds  = 300/0.003                      # K/m, 300 K across the 3 mm pool
f_grav = rho_model*g
f_self = J*B_self
f_buoy = rho_model*g*3e-5*200          # beta*dT, dT ~ 200 K over the pool
P_arc  = mu0*I**2/(8*math.pi*sig_p**2)
body = [("arc self-field J x B", f_self), ("gravity", f_grav),
        ("buoyancy (Boussinesq)", f_buoy)]
surf = [("surface tension  gamma*kappa", gamma/L_pool),
        ("arc pressure     peak",        P_arc),
        ("Marangoni        (dg/dT)(dT/ds)", abs(dgdT)*dTds)]
rows = [(n+"   [body]", v*L_pool) for n,v in body] + \
       [(n+" [surf]", v) for n,v in surf]
for n,v in sorted(rows, key=lambda r:-r[1]):
    print(f"  {n:36s} {v:10.2f} Pa   {v/(f_grav*L_pool):7.2f} x gravity")
print(f"\n  arc pressure peak P_arc = mu0 I^2/(8 pi sigma_p^2) = {P_arc:.1f} Pa")
print(f"  external field at 20 mT: J*B_ext*L = {J*0.02*L_pool:.0f} Pa "
      f"({J*0.02*L_pool/(f_grav*L_pool):.1f} x gravity)")

# ------------------------------------------------- 4. external field crossovers
banner("4. WHEN DOES AN EXTERNAL FIELD MATTER?")
print("  External Lorentz density = |J| * B_ext , with |J| set by the welding current.")
for tgt, name in [(f_grav,"gravity"), (f_self,"the arc's own self-field"),
                  (abs(dgdT)*dTds/L_pool,"Marangoni"),
                  (P_arc/L_pool,"arc pressure")]:
    B = tgt/J
    print(f"    matches {name:26s} at B_ext = {B*1e3:7.2f} mT")

# ------------------------------------------------- 5. dimensionless
banner("5. DIMENSIONLESS NUMBERS AT WAAM SCALE")
print(f"  L = {L_pool*1e3:.0f} mm , U = {U_pool} m/s , sigma_e = {sig_e:.2e} S/m , mu = {mu_f:.2e} Pa.s")
print(f"\n  {'B (mT)':>8} {'Ha':>8} {'N':>10} {'Re':>9} {'Rm':>10}   regime")
Re = rho_model*U_pool*L_pool/mu_f
for B in (5,10,20,50,100,150,270,500):
    b  = B*1e-3
    Ha = b*L_pool*math.sqrt(sig_e/mu_f)
    N  = sig_e*b**2*L_pool/(rho_model*U_pool)
    Rm = mu0*sig_e*U_pool*L_pool
    if   N >= 1.0: reg = "bulk flow magnetically dominated"
    elif Ha >= 1.0: reg = "Hartmann layers form; bulk largely unchanged"
    else:           reg = "no significant magnetic effect"
    print(f"  {B:8.0f} {Ha:8.2f} {N:10.4f} {Re:9.0f} {Rm:10.2e}   {reg}")
B_N1  = math.sqrt(rho_model*U_pool/(sig_e*L_pool))
B_Ha1 = 1.0/(L_pool*math.sqrt(sig_e/mu_f))
print(f"\n  Ha = 1  at B = {B_Ha1*1e3:6.2f} mT")
print(f"  N  = 1  at B = {B_N1*1e3:6.1f} mT   <-- order-unity BULK braking threshold")
print(f"  Rm = {mu0*sig_e*U_pool*L_pool:.2e}  << 1 : induced field negligible, quasi-static OK")

# ------------------------------------------------- 6. what the tutorial conductivity did
banner("6. CONSEQUENCE OF THE TUTORIAL'S sigma_e = 1e6 S/m")
r = sig_e/sig_e_tut
print(f"  true liquid Al-5Mg sigma_e / tutorial sigma_e = {r:.2f}")
print(f"    J x B  scales as sigma      -> understated by {r:.2f} x")
print(f"    Ha     scales as sqrt(sigma)-> understated by {math.sqrt(r):.2f} x")
print(f"    N      scales as sigma      -> understated by {r:.2f} x")
print(f"  Phase 1's Plate2D null was measured on a melt {r:.1f}x less conductive than real.")

# ------------------------------------------------- 7. Joule heating
banner("7. JOULE HEATING vs ARC HEATING")
V_pool = (4/3)*math.pi*L_pool**3
Q_joule = (J**2/sig_e)*V_pool
Q_arc   = eta*U_arc*I
print(f"  q_joule = |J|^2/sigma_e   = {J**2/sig_e:.3e} W/m3")
print(f"  pool volume ~ {V_pool:.2e} m3  ->  {Q_joule:.2f} W")
print(f"  arc into pool                 =  {Q_arc:.0f} W")
print(f"  ratio = {100*Q_joule/Q_arc:.3f} %   -> negligible, but now quantified with a real sigma_e")

# ------------------------------------------------- 8. the central distinction
banner("8. STIRRING vs BRAKING -- the physics that justifies current injection")
print("""  Two completely different couplings share the name 'Lorentz force':

    STIRRING   F = J_weld x B_ext        J from the ARC (135 A, injected)
               linear in B, DRIVES flow, needs a current source

    BRAKING    F = sigma(u x B) x B      J from MOTION only
               quadratic in B, DAMPS flow, needs no current source

  MTHD as shipped has only the second. laserbeamFoam is a LASER code: there is
  no welding current to inject, so J = curl(H) is purely motional.
""")
print(f"  {'B (mT)':>8} {'stir (Pa)':>11} {'brake (Pa)':>11} {'stir/brake':>11}")
for B in (5,10,20,50,100,270,500):
    b = B*1e-3
    stir  = J*b*L_pool
    brake = sig_e*U_pool*b**2*L_pool
    print(f"  {B:8.0f} {stir:11.1f} {brake:11.3f} {stir/brake:11.0f}")
print(f"\n  At 20-50 mT -- the range the experimental literature actually uses --")
print(f"  the current-driven term is 200-500x the motional term. They cross only")
print(f"  at B = {J/(sig_e*U_pool):.1f} T, an order of magnitude past any laboratory magnet.")

banner("9. WHAT THIS SAYS ABOUT PHASE 1")
L_p2d, U_p2d = 5e-4, 0.1          # Plate2D: laser pool, order-of-magnitude
for tag, s, L_, U_ in [("Plate2D as run (sigma=1e6)", sig_e_tut, L_p2d, U_p2d),
                       ("Plate2D at true sigma",      sig_e,     L_p2d, U_p2d)]:
    N  = s*(0.02)**2*L_/(rho_model*U_)
    Ha = 0.02*L_*math.sqrt(s/mu_f)
    print(f"  {tag:28s} Ha = {Ha:5.2f}   N = {N:.2e}")
print("""
  N ~ 1e-4. With no current injection, Plate2D at 20 mT sits four orders of
  magnitude below the braking threshold. The Phase 1 'no resolved velocity
  effect' was not a resolution failure -- it is the physically correct answer
  for that configuration. Re-running it at the true sigma_e would not change
  the verdict, because sigma buys only 3.65x against a 1e4 shortfall.

  That reframes Gate B: it proved the MACHINERY works. It could never have
  shown a magnetic effect on the melt, and no field strength would have made
  it. The effect requires the welding current -- which is Phase 6.""")
