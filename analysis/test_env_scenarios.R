# Check the environmental regimes.
#
# Generates the four regimes with the settings the simulation drivers use and
# measures what each one actually did. The point is to confirm
# that predictability and synchrony vary independently before any simulation
# is run on them.
#
# What to look for in the diagnostics table:
#   - mean_patch_cor is near 1 in the synchronous regimes and well below 1 in
#     the asynchronous ones
#   - predictability is high in the predictable regimes and low in the
#     unpredictable ones, and does not change with synchrony
#   - temporal_cv_mean is the same in all four regimes
#   - frac_clamped is small; a large value means the scaling is cutting off
#     the extremes
#
# Writes docs/env_regime_diagnostics.csv.

library(tidyr)
library(dplyr)
library(here)

source(here("simulations/metacom_functions.R"))
source(here("simulations/env_functions.R"))

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

patches   <- 100
x_dim     <- 100
y_dim     <- 100
timesteps <- 300
cycle_length <- 10
w <- (2 * pi) / cycle_length

# ---------------------------------------------------------------------------
# Generate each regime and measure it
# ---------------------------------------------------------------------------

set.seed(20260807)
landscape <- init_landscape(patches = patches, x_dim = x_dim, y_dim = y_dim)

# The regimes at their default settings, which are the ones the drivers use
regimes <- env_scenarios(extent = x_dim)

diagnostics  <- list()

for (i in 1:nrow(regimes)) {
  regime <- regimes[i, ]
  env_df <- env_generate(
    landscape          = landscape,
    timesteps          = timesteps,
    spat_heterogeneity = regime$spat_heterogeneity,
    amplitude          = regime$amplitude,
    w                  = w,
    phase_mode         = regime$phase_mode,
    phase_spread       = regime$phase_spread,
    temp_noise_sd      = regime$temp_noise_sd,
    temp_noise_color   = regime$temp_noise_color,
    spat_range         = regime$spat_range,
    torus              = TRUE,
    scale_half         = regime$scale_half,   # the same scaling for every regime
    seed               = 1000 + i
  )

  diagnostics[[regime$scenario]] <- cbind(
    scenario = regime$scenario,
    env_diagnostics(env_df, w = w),
    frac_clamped = attr(env_df, "regime")$frac_clamped
  )
}

diagnostics <- bind_rows(diagnostics)
print(diagnostics, digits = 3)

write.csv(diagnostics, here("docs/env_regime_diagnostics.csv"), row.names = FALSE)

if (any(diagnostics$frac_clamped >= 0.05)) {
  warning("Some regimes clamp 5% or more of their values: ",
          paste(diagnostics$scenario[diagnostics$frac_clamped >= 0.05], collapse = ", "))
}
