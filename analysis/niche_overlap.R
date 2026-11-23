# Figure: the species' environmental niches against the conditions they meet.
#
# The number of seeds each plant sets is a Gaussian function of the match
# between the conditions in a patch and the species' niche optimum. This
# figure draws
# those curves along the environmental axis, with the line of replacement
# (one seed per plant, before competition), the long-term mean of
# every patch, and how far conditions in one patch move over time, as produced
# by the environments of simulations/metacom_regimes.R. It shows how much the
# species overlap, whether every species has patches that suit it, and whether
# fluctuations within a patch are large enough to carry a species below
# replacement.
#
# Writes figures/supplement/niche_overlap.pdf and prints, for each regime and
# each species, how much habitat it has and how often the best-adapted species
# in a patch falls below replacement.

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
niche_breadth <- 0.1
f_max   <- 5
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

fig_dir <- here("figures", "supplement")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Plotting theme
theme_set(theme_bw() +
            theme(strip.background = element_rect(fill = "white", color = "white"),
                  strip.text = element_text(size = 10),
                  panel.grid.minor = element_blank()))

# ---------------------------------------------------------------------------
# 2. Niches: seeds per plant of each species along the environmental axis
# ---------------------------------------------------------------------------

optima <- seq(optima_span[1], optima_span[2], length.out = species)

niches <- expand_grid(species = 1:species, env = seq(0, 1, length.out = 501)) |>
  mutate(optimum = optima[species],
         r = f_max * exp(-((optimum - env) / (2 * niche_breadth))^2))

# ---------------------------------------------------------------------------
# 3. Conditions: what the regimes of the main text produce
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
                                patch = env_df$patch, env = env_df$env)
}
conditions <- bind_rows(conditions)

# How far conditions move within a patch over time: the central 95% of its
# values, as a distance either side of its mean. The regimes share one
# temporal variance, so this is averaged over them.
within_patch <- conditions |>
  group_by(scenario, patch) |>
  summarize(half_width = (quantile(env, 0.975) - quantile(env, 0.025)) / 2,
            .groups = "drop")
half_width <- mean(within_patch$half_width)

# ---------------------------------------------------------------------------
# 4. How often is the best-adapted species below replacement?
# ---------------------------------------------------------------------------

# Seeds per plant of every species in every patch and time step
growth <- conditions |>
  expand_grid(species = 1:species) |>
  mutate(r = f_max * exp(-((optima[species] - env) / (2 * niche_breadth))^2))

# The best-adapted species in a patch has the highest geometric mean growth
by_patch <- growth |>
  group_by(scenario, patch, species) |>
  summarize(geom_mean_r = exp(mean(log(r))),
            sd_log_r    = sd(log(r)),
            below_replacement = mean(r < 1), .groups = "drop")

bad_years <- by_patch |>
  group_by(scenario, patch) |>
  slice_max(geom_mean_r, n = 1) |>
  group_by(scenario) |>
  summarize(geom_mean_r = mean(geom_mean_r),
            sd_log_r    = mean(sd_log_r),
            below_replacement = mean(below_replacement))
print(as.data.frame(bad_years), digits = 3)

# Habitat of each species: the number of patches where its geometric mean
# growth is above replacement, and the number where it is the best-adapted
habitat <- by_patch |>
  group_by(scenario, patch) |>
  mutate(best = geom_mean_r == max(geom_mean_r)) |>
  group_by(species, scenario) |>
  summarize(viable = sum(geom_mean_r > 1), best = sum(best), .groups = "drop") |>
  group_by(species) |>
  summarize(patches_viable = mean(viable), patches_best = mean(best))
print(as.data.frame(habitat), digits = 3)

# ---------------------------------------------------------------------------
# 5. Figure
# ---------------------------------------------------------------------------

# Marks along the bottom: the long-term mean of each patch. Grey band: the
# conditions that one patch, with a mean in the middle of the axis, moves
# through over time. Dashed line: replacement.
p <- ggplot(niches, aes(x = env, y = r)) +
  annotate("rect", xmin = 0.5 - half_width, xmax = 0.5 + half_width,
           ymin = -Inf, ymax = Inf, fill = "grey88") +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.3) +
  geom_line(aes(color = optimum, group = species), linewidth = 0.7) +
  geom_rug(data = data.frame(env = patch_means), aes(x = env), inherit.aes = FALSE,
           sides = "b", linewidth = 0.3, length = unit(0.03, "npc")) +
  scale_color_viridis_c(end = 0.9) +
  scale_x_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, f_max), expand = expansion(mult = c(0, 0.03))) +
  labs(x = "Environmental conditions", y = "Propagules per individual",
       color = "Niche optimum") +
  theme(panel.grid = element_blank())

ggsave(file.path(fig_dir, "niche_overlap.pdf"), p, width = 6, height = 4, bg = "white")
message("Wrote ", file.path(fig_dir, "niche_overlap.pdf"))
