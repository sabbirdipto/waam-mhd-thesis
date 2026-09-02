#!/usr/bin/env python3
"""Generate constant/alMg5Properties for the WAAM-MHD thesis.

Every entry carries its provenance. Four tags are used, and they are the point
of this file -- a number you cannot defend in a viva is worse than no number:

  [ZHAO]   printed in Zhao et al. 2021 Table 2/3/4, or in its text
  [MEAS]   an independent measurement, source named
  [CONSTR] a construction from measured data, method named
  [EST]    an estimate. Sensitivity to it must be tested.
"""
import math, textwrap

L0 = 2.44e-8    # Sommerfeld Lorenz number, W.ohm/K^2
Ts, Tl = 815.0, 906.0

# --- electrical resistivity, Ohm.m --------------------------------------
# Desai/James/Ho (1984) JPCRD 13(4):1131 pure-Al fit, + Matthiessen offset
# 3.30e-8 Ohm.m from the 5083/5356 room-temperature value (29% IACS).
def rho_e_solid(T):  return 1e-8*(0.1743 + 7.502e-3*T + 4.028e-6*T**2) + 3.30e-8
def rho_e_liquid(T): return 1e-8*(8.2019 + 1.9799e-2*T - 2.1227e-6*T**2) + 3.0e-8

def sigma_e(T):
    if T <= Ts: return 1.0/rho_e_solid(T)
    if T >= Tl: return 1.0/rho_e_liquid(T)
    fL = (T-Ts)/(Tl-Ts)
    return (1-fL)/rho_e_solid(Ts) + fL/rho_e_liquid(Tl)

def k_th(T):
    """Thermal conductivity by Wiedemann-Franz from sigma_e(T).
    Validated: pure Al solid at Tm 213 vs 211 measured; pure Al LIQUID at Tm
    91.7 vs 91 measured; Al-5Mg at 300 K 120 vs 123 measured. All 1-2%."""
    return L0*sigma_e(T)*T

def cp(T):
    """Solid branch anchored on pure Al (897 J/kgK at 298 K, 1177 at Tm);
    Neumann-Kopp with 5% Mg shifts it by +0.7%, inside the uncertainty."""
    if T <= Ts: return 763.0 + 0.457*T
    if T >= Tl: return 1180.0
    fL = (T-Ts)/(Tl-Ts)
    return (1-fL)*(763.0+0.457*Ts) + fL*1180.0

def mu_dyn(T):
    """Assael et al. 2006 for liquid Al: eta[mPa.s] = 0.1492 exp(1984.5/T).
    1.25 mPa.s at Tm, the accepted value. Mg at 5% is within the scatter."""
    return 1e-3*0.1492*math.exp(1984.5/max(T, Tl))

# ------------------------------------------------------------- tables
def table():
    out = []
    for T in (300,400,500,600,700,800,815,860,906,1000,1200,1500,1800,2000):
        ph = "solid" if T <= Ts else ("mushy" if T < Tl else "liquid")
        out.append((T, ph, sigma_e(T), k_th(T), cp(T),
                    mu_dyn(T) if T >= Tl else float('nan')))
    return out

hdr = f"""/*--------------------------------*- C++ -*----------------------------------*\\
| Al-5%Mg (5A05 substrate / ER5356 wire) property set                         |
| Thesis: magnetohydrodynamics in GMAW-WAAM.  Phase 2.                        |
|                                                                             |
| Replication target : Zhao, Wei, Long, Chen, Liu, Ou (2021)                  |
|                      Welding in the World 65:1571-1590                      |
|                      doi 10.1007/s40194-021-01123-1                         |
|                                                                             |
| PROVENANCE TAGS -- every value carries one.                                 |
|   [ZHAO]   printed in Zhao Table 2/3/4 or its text                          |
|   [MEAS]   independent measurement, source named                            |
|   [CONSTR] construction from measured data, method named                    |
|   [EST]    estimate. Sensitivity must be tested. See S-cases in the matrix. |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      alMg5Properties;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ============================================================ process
// [ZHAO] Table 2. Unit conversions verified in scaling.py section 1.
weldingCurrent          135.0;      // A
weldingVoltage          16.0;       // V
travelSpeed             0.02;       // m/s   <- 120 cm/min
wireFeedSpeed           0.15;       // m/s   <- 900 cm/min, matches Table 3 inlet
wireDiameter            1.2e-3;     // m
shieldingGasFlow        3.3333e-4;  // m3/s  <- 20 L/min
nozzleRadius            10.0e-3;    // m     [ZHAO] Sec 3.6
dropletInletTemperature 1506.0;     // K     [ZHAO] Sec 3.5

// ============================================================ phase change
solidusTemperature      815.0;      // K     [ZHAO] Table 4, via Mills LM5/5182
liquidusTemperature     906.0;      // K     [ZHAO] Table 4, via Mills LM5/5182
latentHeatOfFusion      3.58e5;     // J/kg  [ZHAO] Table 4, via Mills
mushyZoneConstant       1.0e6;      // kg/m3/s [EST] Zhao gives only the range
                                    //   1e4-1e7 and never states the value used.
                                    //   1e6 is the geometric centre. THIS IS A
                                    //   FREE PARAMETER OF THE REPLICATION -- see
                                    //   S-cases. Do not present it as the paper's.
darcyEpsilon            1.0e-3;     // -     [ZHAO] Eq. 7

// ============================================================ density
// [ZHAO] Table 4 uses ONE constant density for both phases, with Boussinesq
// carrying only thermal expansion. Real liquid Al-5Mg near the liquidus is
// about 2350 kg/m3, so the melt is modelled ~13% too heavy. Kept at 2650 for
// replication fidelity; recorded here because it biases every inertial term.
rho                     2650.0;     // kg/m3 [ZHAO] Table 4
rhoLiquidActual         2350.0;     // kg/m3 [MEAS] liquid Al 2377 at Tm, Mg-shifted
thermalExpansion        3.0e-5;     // 1/K   [ZHAO] Table 4
rhoGas                  1.6228;     // kg/m3 [MEAS] argon at 300 K, 1 atm
                                    //   NOTE Zhao Sec 3.3 states argon density as
                                    //   6e-6 kg/m3 for the plasma drag term. That
                                    //   is ~4 orders below argon at any temperature
                                    //   (0.05 kg/m3 even at 10,000 K). Treat as a
                                    //   probable typo; plasma drag is then under-
                                    //   stated. Flagged, not silently corrected.

// ============================================================ surface
// [ZHAO] Table 4 -- but sourced to Mills' PURE ALUMINIUM chapter (ref [34]),
// while every thermal property comes from the LM5/5182 ALLOY chapter (ref [33]).
// This is the paper's one deliberate pure-Al substitution. Magnesium is strongly
// surface-active in aluminium, so a true Al-5Mg value is lower. 0.85 N/m happens
// to coincide with the OXIDIZED pure-Al figure (0.865, Garcia-Cordovilla 1986);
// clean pure Al is 1.10 N/m. No measured Al-5Mg value was located.
// Marangoni drives the paper's headline surface flow, so this is the largest
// single uncertainty in the replication. Carry it as a sensitivity case.
sigmaSurface            0.85;       // N/m   [ZHAO]/[EST]
dSigmaDT                -1.55e-4;   // N/m/K [ZHAO]/[EST]

// ============================================================ heat source
// [ZHAO] Table 4 and Eq. 27. af/ar/b/c were CALIBRATED by Zhao against their own
// first-layer pool shape -- they are fitted, not transferable. Recalibrate.
arcEfficiency           0.7;        // -     [ZHAO] excludes droplet heating
goldak_af               3.0e-3;     // m     [ZHAO] calibrated
goldak_ar               6.0e-3;     // m     [ZHAO] calibrated
goldak_b                4.0e-3;     // m     [ZHAO] calibrated
goldak_c                4.0e-3;     // m     [ZHAO] calibrated
goldak_ff               0.6667;     // -     [CONSTR] 2af/(af+ar), ff+fr=2
goldak_fr               1.3333;     // -     [CONSTR] 2ar/(af+ar)
currentDistRadius       2.0e-3;     // m     [ZHAO] Table 4, marked "[estimated]"
arcPressureDistRadius   2.0e-3;     // m     [ZHAO] Table 4, marked "[estimated]"
plasmaFlowCoefficient   0.44;       // -     [ZHAO] Table 4
plasmaVelocity          100.0;      // m/s   [ZHAO] Sec 3.3, "assumed"
heatTransferCoefficient 80.0;       // W/m2/K [ZHAO] Table 4
emissivity              0.4;        // -     [EST] not given by Zhao. Oxidised Al
                                    //   0.2-0.4; radiation is minor here.
ambientTemperature      300.0;      // K     [ZHAO] Table 3

// ============================================================ electromagnetic
// NOT IN ZHAO. The paper computes Lorentz force from an ASSUMED Gaussian current
// distribution (Eqs 22-24, after Kumar & DebRoy 2003) and never invokes Ohm's
// law, so it needs no electrical conductivity at all. Everything below is new
// work for this thesis and cannot be cited to Zhao.
//
// SOLID   sigma = 1/[1e-8(0.1743 + 7.502e-3 T + 4.028e-6 T^2) + 3.30e-8]
//         [CONSTR] Desai, James & Ho (1984) JPCRD 13(4):1131 pure-Al recommended
//         values (+-2-5%), plus a temperature-independent Matthiessen offset set
//         by the measured 29% IACS of 5083/5356. Cross-checked against Brandt &
//         Neuer's measured AlSiCu-19 (same RT resistivity): agrees to 1%.
//
// JUMP    x2.19-2.42 on melting. [MEAS] Brandt & Neuer (2007) Int J Thermophys
//         28(5):1429 measured this for 11 Al alloys; the band is +-5% with no
//         compositional trend. Transferred from Al-Si to Al-Mg -- nobody has
//         measured Al-Mg. State this as a limitation.
//
// LIQUID  sigma = 1/[1e-8(8.2019 + 1.9799e-2 T - 2.1227e-6 T^2) + 3.0e-8]
//         [CONSTR] Desai liquid pure-Al curve + a 3e-8 Ohm.m alloy offset. Two
//         independent routes (jump-ratio; pure-liquid-plus-solute) agree to 8%.
//         NO MEASUREMENT OF LIQUID Al-Mg RESISTIVITY EXISTS at any composition.
//         Closest primary source, theory only: Joshi et al. (2013) Adv Mater Res
//         665:76. Uncertainty +-15% at the liquidus, +-25% by 2000 K.
sigmaESolid300K         1.643e7;    // S/m   [MEAS] 29% IACS -> 5.95e-8 Ohm.m
sigmaESolidus           8.154e6;    // S/m   [CONSTR] at 815 K
sigmaELiquidus          3.650e6;    // S/m   [CONSTR] at 906 K   <- KEY NUMBER
sigmaELiquidRef         3.6e6;      // S/m   single-value fallback
magneticPermeability    1.25664e-6; // H/m   [ZHAO] 4*pi*e-7, non-magnetic alloy
sigmaEGasRatio          1.0e-6;     // -     [MEAS-ish] the MTHD Plate2D tutorial's
                                    //   own floor. Held div(muH) at 1e-12 across
                                    //   the jump in Gate B. Far more aggressive
                                    //   than the 1e-3 assumed before Phase 1.

// ============================================================ scaling results
// Derived in scaling.py. These set the case matrix, so they live with the data.
//   Ha = 1  at    6.2 mT      Hartmann layers begin to form
//   N  = 1  at  269   mT      order-unity BULK braking. No experiment reaches it.
//   Rm = 4.1e-3            << 1, quasi-static approximation is safe
//   J x B_ext beats gravity at 2.4 mT and the arc's OWN self-field at 13.5 mT
//   Joule heating 3.6 W against 1512 W of arc  =  0.24%, negligible but now
//     quantified with a real conductivity rather than asserted
//
// THE RESULT THAT SHAPES THE THESIS
//   Current-driven STIRRING (J_weld x B, linear in B) is 200-500x larger than
//   motional BRAKING (sigma u B^2, quadratic) at 20-50 mT. They cross at 9.8 T.
//   laserbeamFoam's MTHD module has ONLY the braking term. It is therefore in
//   the wrong physical regime for magnetised GMAW by two to three orders of
//   magnitude, at every field strength anyone applies. Injecting the welding
//   current is not a refinement -- without it the model cannot show the effect.

// ************************************************************************* //
"""

if __name__ == "__main__":
    open("alMg5Properties","w").write(hdr)
    print(hdr.split("// ====")[0].split("FoamFile")[0])
    print(f"{'T (K)':>7} {'phase':>7} {'sigma_e (S/m)':>14} {'k (W/mK)':>10} "
          f"{'cp (J/kgK)':>11} {'mu (Pa.s)':>11}")
    print("  " + "-"*66)
    for T,ph,s,k,c,m in table():
        ms = f"{m:11.3e}" if m == m else f"{'-':>11}"
        print(f"{T:7.0f} {ph:>7} {s:14.3e} {k:10.1f} {c:11.0f} {ms}")
    print("\nwrote alMg5Properties")
