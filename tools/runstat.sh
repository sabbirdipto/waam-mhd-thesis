#!/bin/sh
# Status of a running (or finished) laserbeamFoam case, in one command.
#
# WHY THIS IS A FILE AND NOT A PASTED COMMAND
#     The awk one-liners this replaces are multi-line, and a multi-line paste
#     that loses its first line leaves bash halfway through parsing a quote,
#     which looks like a syntax error in the middle of a run and is alarming
#     for no reason. A file cannot half-arrive.
#
# USAGE
#     sh ~/thesis/tools/runstat.sh                 # waam-arc-01
#     sh ~/thesis/tools/runstat.sh ../some-case    # any case directory
#
# Needs no OpenFOAM environment -- it only reads the log.

CASE=${1:-$HOME/thesis/cases/00-tutorials/waam-arc-01}
LOG=$CASE/log.laserbeamFoam

if [ ! -f "$LOG" ]; then
    echo "  no log at $LOG"
    echo "  cases present:"
    ls -d "$(dirname "$CASE")"/waam-arc-01* 2>/dev/null | sed 's/^/    /'
    exit 1
fi

END=$(awk '/^endTime/{print $2}' "$CASE/system/controlDict" 2>/dev/null | tr -d ';')
DT=$(awk '/^maxDeltaT/{print $2}'  "$CASE/system/controlDict" 2>/dev/null | tr -d ';')
[ -z "$END" ] && END=0.5
[ -z "$DT"  ] && DT=1e-5

echo "  case      $CASE"
echo "  endTime   $END s      maxDeltaT $DT s"

# -x matches the process NAME exactly. -f would match any command line
# containing the string, which includes this script's own grep and awk
# children reading log.laserbeamFoam -- so it always reported RUNNING.
if pgrep -x laserbeamFoam >/dev/null 2>&1; then
    echo "  state     RUNNING"
else
    echo "  state     not running (finished, or killed)"
fi
echo

# ---- progress, and the rate measured as a SLOPE rather than an average -----
# Averaging from t=0 folds in solver startup -- field reading, mesh
# construction, the initial pcorr solve -- which makes an early reading
# pessimistic. The slope over recent steps drops all of that out.
awk -v end="$END" -v dt="$DT" '
/^Time = /       { t = $3 }
/^ExecutionTime/ { n++; T[n] = t; E[n] = $3 }
END {
  if (n == 0) { print "  no timesteps logged yet -- still meshing or decomposing"; exit }
  # int(x) truncates, and end/dt is 49999.999... in floating point, so the
  # step total printed as 49999. Round instead.
  steps  = int(t / dt + 0.5)
  total  = int(end / dt + 0.5)
  printf "  sim       %.4f / %s s   %.1f%%   (%d of %d steps)\n",
         t, end, 100*t/end, steps, total
  printf "  cpu       %.0f s elapsed\n", E[n]
  i = n - 500; if (i < 1) i = 1
  ds = (T[n] - T[i]) / dt
  if (ds >= 50) {
    r = (E[n] - E[i]) / ds
    printf "  rate      %.4f s/step   (slope over the last %d steps)\n", r, ds
    printf "  eta       %.0f min remaining   total %.0f min\n",
           r*(total - steps)/60, r*total/60
  } else {
    printf "  rate      %.4f s/step   (average from t=0; too few steps for a slope)\n",
           E[n]/steps
  }
}' "$LOG"

# ---- what the solver is actually spending its solves on --------------------
# Hx/Hy/Hz and pH mean the magnetic module is live. With B = 0 that is 15 of
# 19 solves per timestep doing nothing, and pH is elliptic, so it is two
# thirds of the EXPENSIVE solves. Measured cost on this case: 1.79x.
echo
echo "  linear solves per timestep:"
NSTEP=$(grep -c '^Time = ' "$LOG")
grep -oE 'Solving for [A-Za-z_.]+' "$LOG" | sort | uniq -c | sort -rn |
while read -r c f; do
    [ "$NSTEP" -gt 0 ] && per=$(awk -v c="$c" -v n="$NSTEP" 'BEGIN{printf "%.1f", c/n}') || per="?"
    printf "    %-22s %8d total   %s per step\n" "${f#Solving for }" "$c" "$per"
done

# ---- the startup gates, printed once at the top of every run ---------------
echo
echo "  startup gates:"
grep -iE 'Goldak|renormalis|two-branch|no metal' "$LOG" | head -5 | sed 's/^/    /'

# ---- has the timestep collapsed --------------------------------------------
echo
echo "  last deltaT:"
grep '^deltaT' "$LOG" | tail -2 | sed 's/^/    /'
