# First liquid metal: the 0.1 s full-domain droplet run

Case `waam-arc-01`, full domain, 143,616 cells, 4 ranks, Edit 4 (droplet
injection) active, `electromagnetics false`, adiabatic walls. Three droplets
at 30 Hz. This is the run that closes out Edit 4 and supersedes
`Arc-Alone-Cannot-Melt.pdf` on the question of whether a pool is reachable.

## What was established

**Momentum transfer works.** `U metal` read `1.4164` m/s at the first
injection. That column had been exactly `0.0000` in every run this project
had ever produced — a true zero, not a rounded one, because the Darcy source

    S = -C (1-eps1)^2 / (eps1^3 + q) * U

is unbounded as the liquid fraction eps1 -> 0 and drives velocity to machine
zero in solid cells. It decays back to `0.0000` by t = 0.045 as the droplet
freezes, which is the enthalpy-porosity coupling behaving correctly, not a
failure. Mass, energy and momentum all now arrive with the droplet; the three
old-time fields (`alpha1`, `T`, `epsilon1`) plus `rho.oldTime()` and
`rhoCp.oldTime()` are what made it work.

**Mass is conserved to 3.7%.** The bead cross-section is fixed before the run
by mass conservation, A = Vdot_wire / v_travel = 169.6 / 20 = **8.482 mm^2**.
17.6 mm^3 deposited while the torch travelled 20 mm/s x 0.1 s = 2.00 mm gives
**8.80 mm^2**. Droplet frequency also checks: V_d = (4/3)pi(1.105 mm)^3 =
5.652 mm^3 at 169.6 mm^3/s is 30.02 Hz, period 33.3 ms, matching injections
observed at 0.035 / 0.070 / 0.100 within the 5 ms write resolution.

**Energy books balance.** `E kept` reached 138.2%, and exceeding 100% is
correct: the percentage is metal thermal energy over *arc* energy only, and
droplets carry their own. Arc 151.2 J in 0.1 s; droplet sensible enthalpy
52.9 J (a mean of **529 W**, a third again on top of the 1512 W arc); latent
heat up to 18.5 J released on freezing, which `pool_stats` counts as sensible
because it only knows rho*cp*dT. Sensible alone predicts 135.0%, with two of
three droplets solidified 143.1%. Observed 138.2 sits between.

**Domain independence still holds.** `far` = 300.0 at every written time.

## What failed: there is no melt pool

Six consecutive writes (0.045 through 0.065) show zero cells above solidus
anywhere in the domain. Each droplet arrives, freezes within ~10 ms, and the
next lands on cold solid metal 33 ms later. `depth` reads 0.25 mm, which on a
0.5 mm z-cell is **exactly one cell layer** — the number is the mesh, not the
physics.

The 126,000 K/s cooling rate is physical. For two semi-infinite bodies in
contact the interface jumps immediately to a temperature weighted by thermal
effusivity b = sqrt(k rho cp):

    T_c = (b1 T1 + b2 T2) / (b1 + b2)

Droplet and substrate are the same alloy, so b1 = b2 = 19,650 and this
collapses to the arithmetic mean, T_c = (1506 + 300)/2 = **903 K** — three
kelvin below the 906 K liquidus. A droplet of Al-5Mg landing on a 300 K plate
of Al-5Mg is mushy at the contact plane the instant it touches. And it is
fully penetrated at once: sqrt(4 a t) = 1.12 mm at 5 ms against a 1.105 mm
droplet radius, so no hot core survives.

## The prediction the 0.5 s run tests

0.1 s is pure startup. The arc-only portion of the table (369 K at 0.005
rising to 549 K at 0.030) tracks the arc-alone trajectory, which needed ~0.5 s
to reach its 832 K saturation — 74 K short of liquidus. But the droplets add
529 W to the arc's 1512 W. For a diffusion-limited quasi-steady source the
saturation rise scales with power:

    dT_sat = 532 K x 2041/1512 = 718 K   =>   T_sat ~ 1018 K

That is 112 K above liquidus. **The droplets are what make a pool possible** —
the answer to what has blocked this project since the arc-alone result. It is
a linear extrapolation of a nonlinear system, but it is falsifiable, and the
test is `endTime 0.5`.

The 0.5 s run deliberately keeps the **adiabatic** walls rather than adding
Edit 10 (Zhao Eq. 29, convective+radiative h = 80), so that only time changes.
The bias is favourable and bounded: even if the entire 4,320 mm^2 wetted
substrate boundary sat at 517 K, the missing loss is 75 W = 3.7% of input. If
no pool forms *with* that bias, the negative result is strong.

## Watch items

- `bottom` climbed to 364.8 K by 0.1 s. sqrt(4 a t) = 4.99 mm at 0.1 s against
  a 6 mm plate, and 11.16 mm at 0.5 s. Edit 10 moves from deferred to needed
  once this run is interpreted.
- `pool_stats` rev 3 and `check_droplet` rev 2 fix two measurement defects this
  run exposed: `max T` maximised over argon as well as metal (1277.8 K vs
  872.5 K at t = 0.040, a 400 K error on the runaway gate), and `metal all`
  counted whole cells above alpha > 0.5 rather than sum(alpha*V), so interface
  smearing looked like 1.5 mm^3 of mass loss that never happened.
