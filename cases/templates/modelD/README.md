# Model P template -- melt pool

Generated cases MUST NOT be edited by hand. Change caseMatrix.csv and rerun
generate_cases.py, or edit this template if the change applies to every case.

## Token substitution

Write placeholders as at-sign, name, at-sign. Valid names are any column in
caseMatrix.csv, plus these derived values:

    Rw_m                 wire radius
    Bx_T By_T Bz_T       applied field resolved onto axes (tesla)
    B_omega_rad_s        2*pi*f, zero for a static field
    B_is_alternating     0 or 1
    x_min_m .. z_max_m   domain box
    nx ny nz n_cells     uniform mesh counts
    end_time_s           filled in from length/speed if left blank
    arc_power_W          eta * U * I

## The other mechanism

Dictionaries can instead do

    #include "$FOAM_CASE/constant/caseParameters"

and reference $I_A, $dt_s, $Bext directly. Prefer this where OpenFOAM's own
expansion works -- it keeps values live and hand-tweakable. Fall back to token
substitution inside verbatim code blocks, where $-expansion does not reach.
