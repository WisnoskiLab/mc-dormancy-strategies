# The bet-hedging control for a cued seed bank, matched in realized emergence.
#
# Every seed bank has the same mean emergence fraction, but a cued species
# emerges most where conditions suit it, which is also where most of its
# propagules are, so it realizes a higher fraction than a bet-hedger at the
# same nominal value. Comparing the two at the same nominal fraction would mix
# what the cue does with what a smaller seed bank does. The main sweep
# therefore runs bet-hedgers at a ladder of fractions, and each cued run is
# compared with the bet-hedger interpolated to the fraction it realized: log
# values are interpolated linearly in realized emergence between the two
# bet-hedgers that bracket it, in the same replicate, regime, and dispersal
# rate. Interpolating on the log scale keeps sums of log terms (as in the
# pathway split) exact.
#
# Sourced by seedbank_stability.R, stability_pathways.R, and seedbank_checks.R.

# The interpolated value for one cued run. `g_controls` and `values` are the
# realized emergence and the values of the bet-hedgers in its cell. Values that
# are not positive and finite (a lost run) are left out. Outside the ladder the
# nearest end is used, and the caller flags it.
interpolate_control <- function(g_out, g_controls, values) {
  usable <- is.finite(log(values)) & is.finite(g_controls)
  if (sum(usable) < 2) {
    return(NA_real_)
  }
  exp(approx(g_controls[usable], log(values[usable]), xout = g_out, rule = 2,
             ties = mean)$y)
}

# Each cued run with the value of its matched bet-hedger.
#
# cued      one row per cued run and response, with the columns in `keys`,
#           `g_realized`, and `value`
# controls  one row per bet-hedging run and response, with the columns in
#           `keys`, `g_control` (its realized emergence), and `value_control`
# Returns `cued` with `value_reference`, the interpolated bet-hedger, and
# `outside_ladder`, TRUE where the run realized a fraction beyond the ladder.
match_bet_hedging <- function(cued, controls, keys) {
  cued |>
    mutate(row_id = row_number()) |>
    left_join(controls, by = keys, relationship = "many-to-many") |>
    group_by(row_id) |>
    summarise(value_reference = interpolate_control(first(g_realized), g_control,
                                                    value_control),
              outside_ladder  = first(g_realized) < min(g_control, na.rm = TRUE) |
                                first(g_realized) > max(g_control, na.rm = TRUE),
              .groups = "drop") |>
    right_join(mutate(cued, row_id = row_number()), by = "row_id") |>
    arrange(row_id) |>
    select(-row_id)
}
