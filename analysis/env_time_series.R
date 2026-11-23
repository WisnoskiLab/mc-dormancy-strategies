# Figure: the environment through time, the lower panel of Fig. 1.
#
# Draws the conditions in a subset of patches over time, as produced by the
# environments of simulations/metacom_regimes.R, with the species' niche
# optima marked on the same axis. It shows how far patches differ in their
# long-term means, how far each one fluctuates around its mean, and whether
# patches fluctuate together.
#
# Writes figures/main/panels/env_asynchronous.pdf: asynchronous patches,
# predictable and red noise side by side, which analysis/model_figure.tex
# places under the life cycle. Prints, for each regime, how the variance in
# conditions divides among patches and through time, and the share of
# conditions that lie beyond the range of the niche optima.

library(tidyr)
library(dplyr)
library(ggplot2)
library(here)

source(here("simulations/metacom_functions.R"))   # landscape
source(here("simulations/env_functions.R"))       # environments and regimes

# ---------------------------------------------------------------------------
# 1. Settings (the same values as simulations/metacom_regimes.R)
# ---------------------------------------------------------------------------

x_dim   <- 100
y_dim   <- 100
patches <- 100
species <- 10
optima_span <- c(0.25, 0.75)
cycle_length <- 10
timesteps    <- 600   # burn_in + timesteps in the driver
patch_mean_span  <- c(0.25, 0.75)

# The driver builds replicate `k` from the seed `replicate_seed + k`. This
# figure shows the first replicate.
replicate_seed <- 20260807

# Noise colors in the driver, and the one shown here for unpredictable regimes
noise_colors     <- c(-0.6, 0, 0.6)
main_noise_color <- 0.6
patch_mean_range <- 10

# What to draw: this many patches, evenly spread from the lowest long-term
# mean to the highest, over this many time steps. Few enough patches that
# asynchronous ones can be told apart.
patches_shown_main <- 5
steps_shown        <- 50

panel_dir <- here("figures", "main", "panels")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 2. Landscape and patch means, drawn as the driver draws them
# ---------------------------------------------------------------------------

set.seed(replicate_seed + 1)
landscape <- init_landscape(patches = patches, x_dim = x_dim, y_dim = y_dim)

# The driver draws a replicate's competition coefficients before its patch
# means, and the same draw is made here so that the patch means come out the
# same.
alpha <- competition_matrix(species = species, intra = 1, min_inter = 0,
                            max_inter = 1, scale = 0.05)
patch_means <- patch_mean_field(landscape, lower = patch_mean_span[1],
                                upper = patch_mean_span[2],
                                spat_range = patch_mean_range, torus = TRUE)

# ---------------------------------------------------------------------------
# 3. One environment per regime
# ---------------------------------------------------------------------------

# The regimes, built and numbered as the driver builds them, so that the seeds
# below give the same environments. The figure keeps the regimes of the main
# text: the predictable ones, and the unpredictable ones with red noise.
regimes <- env_scenarios(extent = x_dim, colors = noise_colors) |>
  filter(pred_frac < 0.5 | temp_noise_color == 0)
regimes$shown <- regimes$pred_frac > 0.5 | regimes$temp_noise_color == main_noise_color

conditions <- list()
for (i in which(regimes$shown)) {
  regime <- regimes[i, ]
  env_df <- env_generate(
    landscape          = landscape,
    timesteps          = timesteps,
    patch_means        = patch_means,
    amplitude          = regime$amplitude,
    w                  = (2 * pi) / cycle_length,
    phase_mode         = regime$phase_mode,
    phase_spread       = regime$phase_spread,
    temp_noise_sd      = regime$temp_noise_sd,
    temp_noise_color   = regime$temp_noise_color,
    spat_range         = regime$spat_range,
    torus              = TRUE,
    scale_half         = regime$scale_half,
    seed               = 4000 + 97 * i + 1
  )
  conditions[[i]] <- data.frame(scenario = regime$scenario, time = env_df$time,
                                patch = as.integer(env_df$patch), env = env_df$env)
}
conditions <- bind_rows(conditions) |>
  mutate(patch_mean = patch_means[patch],
         predictability = ifelse(grepl("^predictable", scenario),
                                 "Predictable", "Red noise"),
         synchrony = ifelse(grepl("_synchronous", scenario),
                            "Synchronous patches", "Asynchronous patches"),
         synchrony = factor(synchrony, levels = c("Synchronous patches", "Asynchronous patches")))

# ---------------------------------------------------------------------------
# 4. How the variance divides, and how much lies beyond the niche optima
# ---------------------------------------------------------------------------

# Among patches: the variance of the patches' long-term means. Through time:
# the mean variance within a patch. Their ratio says which kind of variation
# a species meets more of.
variance_split <- conditions |>
  group_by(scenario, patch) |>
  summarize(mean_in_patch = mean(env), var_in_patch = var(env), .groups = "drop") |>
  group_by(scenario) |>
  summarize(sd_among_patches = sd(mean_in_patch),
            sd_through_time  = sqrt(mean(var_in_patch)),
            variance_ratio   = var(mean_in_patch) / mean(var_in_patch))

beyond_optima <- conditions |>
  group_by(scenario) |>
  summarize(below_optima = mean(env < optima_span[1]),
            above_optima = mean(env > optima_span[2]),
            lowest = min(env), highest = max(env))

print(as.data.frame(left_join(variance_split, beyond_optima, by = "scenario")), digits = 3)

# ---------------------------------------------------------------------------
# 5. Figure
# ---------------------------------------------------------------------------

# Lines: conditions in one patch, colored by that patch's long-term mean on
# the viridis scale that colors niche optima in the niche-overlap figure, so a
# color stands for the same conditions in both. Marks on the right: the
# species' niche optima. Dashed lines: the lowest and highest optimum.
optima <- data.frame(env = seq(optima_span[1], optima_span[2], length.out = species))

# Patches to draw: evenly spread through the ranking of long-term means
shown_main <- order(patch_means)[round(seq(1, patches, length.out = patches_shown_main))]

# Text is sized near the life cycle above it (8-9 pt), and the width to its
# 16.75 cm, so the panel is placed at its natural size and Fig. 1 fits full
# page width (173 mm).
p <- conditions |>
  filter(synchrony == "Asynchronous patches", patch %in% shown_main, time <= steps_shown) |>
  ggplot(aes(x = time, y = env, color = patch_mean, group = patch)) +
  geom_hline(yintercept = optima_span, linetype = "dashed", linewidth = 0.3) +
  geom_line(linewidth = 0.45) +
  geom_rug(data = optima, aes(y = env), inherit.aes = FALSE, sides = "r",
           linewidth = 0.5, length = unit(0.03, "npc")) +
  facet_grid(. ~ predictability) +
  scale_color_viridis_c(end = 0.9) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  labs(x = "Time step", y = expression("Conditions " * theta[x](t)), color = "Patch mean") +
  theme_bw(base_size = 9) +
  theme(strip.background = element_rect(fill = "white", color = "white"),
        strip.text = element_text(size = 9),
        panel.grid.minor = element_blank())

ggsave(file.path(panel_dir, "env_asynchronous.pdf"), p,
       width = 16.75 / 2.54, height = 2.4, bg = "white")
message("Wrote ", file.path(panel_dir, "env_asynchronous.pdf"))
