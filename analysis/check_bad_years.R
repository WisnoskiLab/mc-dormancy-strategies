# How often do conditions fall below replacement for the best-adapted species?
#
# Rebuilds the environments of a regimes sweep exactly, from the replicates'
# saved patches (patches.csv) and the seeds the driver gives env_generate()
# (4000 + 97 * regime + replicate), and measures, over the recorded steps, the
# share of time steps in which the species best adapted to a patch could not
# replace itself there (seeds per germinant below 1). "Best adapted" is the
# species with the highest mean log growth in the patch over the recorded
# steps, as in the driver's `suitability` matrix.
#
# This is the number behind the Methods sentence that justifies the size of
# the fluctuations (sigma_theta = 0.13).
#
# Writes docs/env_bad_years.csv: one row per regime, with the mean, minimum,
# and maximum over replicates, in percent of time steps.

library(dplyr)
library(here)

source(here("simulations/metacom_functions.R"))   # dist.torus
source(here("simulations/env_functions.R"))       # env_generate, env_scenarios

# ---------------------------------------------------------------------------
# Settings: those of simulations/metacom_regimes.R
# ---------------------------------------------------------------------------

# The sweep to rebuild. The most recent one is used unless a directory is named.
sweep_dir <- NULL
if (is.null(sweep_dir)) {
  sweep_dir <- sort(Sys.glob(here("sim_output", "regimes_*")), decreasing = TRUE)[1]
}
message("Rebuilding the environments of ", sweep_dir)

x_dim          <- 100
niche_breadth  <- 0.1
niche_optima   <- seq(0.25, 0.75, length.out = 10)
f_max          <- 5
burn_in        <- 300
timesteps      <- 300
w              <- (2 * pi) / 10
noise_colors   <- c(-0.6, 0, 0.6)

# ---------------------------------------------------------------------------
# Rebuild each environment and count the steps below replacement
# ---------------------------------------------------------------------------

patches_all <- read.csv(file.path(sweep_dir, "patches.csv"))
nreps <- max(patches_all$rep)

# The regimes in the driver's order, so that regime i gets the driver's seed
regimes <- env_scenarios(extent = x_dim, colors = noise_colors) |>
  filter(pred_frac < 0.5 | temp_noise_color == 0)

results <- list()
for (i in seq_len(nrow(regimes))) {
  for (rep_i in 1:nreps) {
    regime <- regimes[i, ]
    p <- patches_all[patches_all$rep == rep_i, ]
    p <- p[order(p$patch), ]

    env_df <- env_generate(
      landscape        = p[, c("x", "y")],
      timesteps        = burn_in + timesteps,
      patch_means      = p$patch_mean,
      amplitude        = regime$amplitude,
      w                = w,
      phase_mode       = regime$phase_mode,
      phase_spread     = regime$phase_spread,
      temp_noise_sd    = regime$temp_noise_sd,
      temp_noise_color = regime$temp_noise_color,
      spat_range       = regime$spat_range,
      torus            = TRUE,
      scale_half       = regime$scale_half,
      seed             = 4000 + 97 * i + rep_i
    )

    # The recorded steps, as in the driver
    recorded <- env_df$time > burn_in & env_df$time <= burn_in + timesteps
    env   <- env_df$env[recorded]
    patch <- as.integer(env_df$patch[recorded])

    # Log seeds per germinant for every species (rows: patch-steps)
    log_r <- sapply(niche_optima, function(z) {
      log(f_max) - ((z - env) / (2 * niche_breadth))^2
    })

    # The best-adapted species in each patch, and the steps it falls below 1
    mean_log_r <- apply(log_r, 2, function(v) tapply(v, patch, mean))
    best <- apply(mean_log_r, 1, which.max)
    below <- log_r[cbind(seq_along(env), best[patch])] < 0

    results[[length(results) + 1]] <- data.frame(
      scenario = regime$scenario,
      rep = rep_i,
      share_below = mean(below)
    )
  }
}

bad_years <- bind_rows(results) |>
  group_by(scenario) |>
  summarize(percent_mean = 100 * mean(share_below),
            percent_min  = 100 * min(share_below),
            percent_max  = 100 * max(share_below),
            replicates   = n(),
            .groups = "drop")

print(as.data.frame(bad_years), digits = 3)
write.csv(bad_years, here("docs", "env_bad_years.csv"), row.names = FALSE)
message("Wrote docs/env_bad_years.csv")
