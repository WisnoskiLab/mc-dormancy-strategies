# Figures: dormancy, dispersal, and stability across scales.
#
# Reads the sweep written by simulations/metacom_regimes.R. The figures ask two
# questions in order, along the gradient of dispersal rate (the probability
# that a new propagule leaves its patch):
#
#   What does a seed bank do?             every kind of dormancy, and none
#   Does the way emergence is controlled  the cued seed banks against the
#   matter?                               bet-hedger that emerges as much
#
# Dormancy (after Lennon et al. 2021):
#   none          every propagule emerges
#   bet-hedging   a fixed fraction emerges each step
#   anticipatory  emergence follows the predictable part of the coming step
#   responsive    emergence follows the conditions of the previous step
# all at the same mean emergence fraction, set below.
#
# Communities, censused once per step after emergence:
#   active    the active individuals; the main figures show them
#   dormant   the dormant propagules
#   pool      active and dormant combined (A + D, which at the census equals
#             the pool P); labeled "active + dormant" in figures and text
#   arriving  the propagules arriving at the end of a step (not drawn)
#
# Kinds of variation through time:
#   predictable   a seasonal cycle with a little white noise
#   red noise     unpredictable, but conditions persist from step to step
#   white noise   unpredictable, and each step is independent of the last
#   blue noise    unpredictable, and a bad step tends to follow a good one
#
# Main figures, in figures/main/, for asynchronous patches in a
# predictable environment and under red noise:
#   diversity          every kind of dormancy and none: alpha, spatial beta,
#                      gamma (active)
#   growth_variance    variance of log growth of the propagule pool, within
#                      patches and over all patches
#   stability          stability at four levels (active)
#   synchrony          the four synchrony terms (active)
#   strategy_layers    cued against the matched bet-hedger: stability of the
#                      active and of the dormant community
# (Fig. 1, the model, is analysis/model_figure.tex with a panel from
# analysis/env_time_series.R.)
#
# Supplementary figures, in figures/supplement/:
#   synchronous_diversity, synchronous_stability, synchronous_strategy_layers
#                      synchronous patches
#   noise_strategy_layers
#                      cued against the matched bet-hedger under red, white,
#                      and blue noise
#   combined_stability active and dormant combined
#   dormant_diversity  the dormant community

library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
library(here)

source(here("analysis/matched_control.R"))   # the matched bet-hedger

# ---------------------------------------------------------------------------
# 1. Settings
# ---------------------------------------------------------------------------

# The sweep to draw. The most recent one is used unless a directory is named.
sweep_dir <- NULL
if (is.null(sweep_dir)) {
  sweep_dir <- sort(Sys.glob(here("sim_output", "regimes_*")), decreasing = TRUE)[1]
}
message("Reading ", sweep_dir)

# The mean emergence fraction of the seed banks that are drawn
focal_g_mean <- 0.1

# The kinds of variation through time shown in the main figures
main_kinds <- c("Predictable", "Red noise")

# The community drawn in every figure but those of the dormant community and
# of active and dormant combined
focal_community <- "active"

# Main-text figures go in one directory, supporting figures in another
fig_dir  <- here("figures")
main_dir <- file.path(fig_dir, "main")
supp_dir <- file.path(fig_dir, "supplement")
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

# Plotting theme
theme_set(theme_bw() +
            theme(strip.background = element_rect(fill = "white", color = "white"),
                  strip.text = element_text(size = 10),
                  strip.placement = "outside",
                  panel.grid.minor = element_blank()))

# Colors for dormancy (Okabe-Ito, colorblind-safe)
dormancy_colors <- c("None" = "black", "Bet-hedging" = "#0072B2",
                     "Anticipatory" = "#D55E00", "Responsive" = "#009E73")

# ---------------------------------------------------------------------------
# 2. Read the sweep and label the treatments
# ---------------------------------------------------------------------------

runs    <- read_csv(file.path(sweep_dir, "summary.csv"), show_col_types = FALSE)
species <- read_csv(file.path(sweep_dir, "per_species.csv"), show_col_types = FALSE)

# Variance of log growth of the propagule pool, averaged over the species in a
# run, so that each run is one value
growth_variance <- species |>
  group_by(run) |>
  summarize(across(c(loglam_var, loglam_var_local), \(x) mean(x, na.rm = TRUE)))
runs <- left_join(runs, growth_variance, by = "run")

# A regime is two things, how conditions vary through time and how synchronous
# they are among patches, and each gets its own column so they can be the rows
# and columns of a figure. Stability is the inverse of variability (CV) at each
# level of the partition.
runs <- runs |>
  mutate(
    dormancy_type = factor(dormancy,
                           levels = c("none", "bet-hedging", "anticipatory", "responsive"),
                           labels = c("None", "Bet-hedging", "Anticipatory", "Responsive")),
    predictability = case_when(
      grepl("^predictable", scenario) ~ "Predictable",
      noise_color > 0                 ~ "Red noise",
      noise_color == 0                ~ "White noise",
      noise_color < 0                 ~ "Blue noise"
    ),
    predictability = factor(predictability,
                            levels = c("Predictable", "Red noise", "White noise", "Blue noise")),
    synchrony = ifelse(grepl("_synchronous", scenario),
                       "Synchronous patches", "Asynchronous patches"),
    synchrony = factor(synchrony, levels = c("Synchronous patches", "Asynchronous patches")),
    stab_population     = 1 / CV_S_L,
    stab_community      = 1 / CV_C_L,
    stab_metapopulation = 1 / CV_S_R,
    stab_metacommunity  = 1 / CV_C_R
  )

# Readable names for every response that gets drawn
response_labels <- c(
  alpha_div     = "Alpha diversity",
  beta_spatial  = "Spatial beta diversity",
  gamma_div     = "Gamma diversity",
  stab_population     = "Population stability",
  stab_community      = "Community stability",
  stab_metapopulation = "Metapopulation stability",
  stab_metacommunity  = "Metacommunity stability",
  phi_S2C_L = "Local species synchrony",
  phi_S2C_R = "Regional species synchrony",
  phi_S_L2R = "Spatial species synchrony",
  phi_C_L2R = "Spatial community synchrony"
)
diversity_vars <- c("alpha_div", "beta_spatial", "gamma_div")
stability_vars <- c("stab_population", "stab_community", "stab_metapopulation",
                    "stab_metacommunity")
synchrony_vars <- c("phi_S2C_L", "phi_S2C_R", "phi_S_L2R", "phi_C_L2R")

# One row per run, community, and response, which is the shape every figure
# below uses. `runs_long_every` holds every run; `runs_long` the runs drawn as
# treatments: no seed bank, and every kind at the focal mean emergence.
runs_long_every <- runs |>
  pivot_longer(all_of(names(response_labels)), names_to = "response", values_to = "value") |>
  mutate(response_label = factor(response_labels[response], levels = response_labels))
runs_long <- filter(runs_long_every, g_mean %in% c(1, focal_g_mean))

# ---------------------------------------------------------------------------
# 3. The cued seed banks against bet-hedging
# ---------------------------------------------------------------------------

# Each cued run is compared with bet-hedging in the same replicate, regime,
# dispersal rate, and community. Those runs saw the same realized environment.
pair_keys <- c("rep", "scenario", "dispersal_rate", "community", "response")

cued <- runs_long |>
  filter(dormancy_type %in% c("Anticipatory", "Responsive"))

# The bet-hedger interpolated to the cued run's realized emergence, along the
# ladder of bet-hedging fractions (analysis/matched_control.R)
ladder <- runs_long_every |>
  filter(dormancy_type == "Bet-hedging") |>
  select(all_of(pair_keys), g_control = g_realized, value_control = value)
ratio_to_matched_every <- cued |>
  match_bet_hedging(ladder, keys = pair_keys) |>
  mutate(value = value / value_reference)

# The ladder should bracket every cued run, so that the match interpolates
# rather than extrapolates
message(sprintf("%.1f%% of cued runs fall outside the bet-hedging ladder",
                100 * mean(ratio_to_matched_every$outside_ladder)))

# ---------------------------------------------------------------------------
# 4. The layouts
# ---------------------------------------------------------------------------

# Every figure is a response along the gradient of dispersal rate. Lines are
# means over replicates and ribbons are +/- 1 SE. Color is the kind of dormancy.
plot_gradient <- function(data, y_label = NULL) {
  ggplot(data, aes(x = dispersal_rate, y = value, color = dormancy_type,
                   fill = dormancy_type)) +
    stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.15, color = NA) +
    stat_summary(fun = mean, geom = "line", linewidth = 0.7) +
    scale_x_log10(labels = scales::label_log()) +
    scale_color_manual(values = dormancy_colors) +
    scale_fill_manual(values = dormancy_colors) +
    labs(x = "Dispersal rate", y = y_label, color = "Dormancy", fill = "Dormancy")
}

# Rows are responses, named on the left beside their values, where a y-axis
# title would be. Columns are the kinds of variation through time. The ratio
# figures share one y axis across rows, so that the same distance means the
# same change in every row.
wrap_rows <- labeller(response_label = label_wrap_gen(width = 18))
facet_time           <- facet_grid(response_label ~ predictability, scales = "free_y",
                                   switch = "y", labeller = wrap_rows)
facet_time_shared    <- facet_grid(response_label ~ predictability,
                                   switch = "y", labeller = wrap_rows)

# A reference line at 1 and a log axis for the ratio figures, so that halving
# and doubling are the same distance from the line
ratio_line <- geom_hline(yintercept = 1, linewidth = 0.3, color = "grey40")

# What following a cue does to the active and the dormant community: the cued
# seed banks as ratios to bet-hedging, for the active (solid) and the dormant
# community (dashed), each against the same community of its control.
plot_layers <- function(ratios, y_label) {
  ratios |>
    filter(community %in% c("active", "dormant"), response %in% stability_vars) |>
    mutate(layer = factor(community, levels = c("active", "dormant"),
                          labels = c("Active", "Dormant"))) |>
    ggplot(aes(x = dispersal_rate, y = value, color = dormancy_type,
               fill = dormancy_type, linetype = layer)) +
    stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.15, color = NA) +
    stat_summary(fun = mean, geom = "line", linewidth = 0.7) +
    scale_x_log10(labels = scales::label_log()) +
    scale_y_log10() +
    scale_color_manual(values = dormancy_colors) +
    scale_fill_manual(values = dormancy_colors) +
    scale_linetype_manual(values = c("Active" = "solid", "Dormant" = "22")) +
    # the community key above the dormancy key, in a fixed order
    guides(linetype = guide_legend(order = 1, override.aes = list(fill = NA)),
           color = guide_legend(order = 2), fill = guide_legend(order = 2)) +
    ratio_line + facet_time_shared +
    labs(x = "Dispersal rate", y = y_label, color = "Dormancy", fill = "Dormancy",
         linetype = NULL)
}

# Diversity and stability for one community and set of regimes, as absolute
# values for every kind of dormancy and none
plot_diversity <- function(data) {
  data |>
    filter(response %in% diversity_vars) |>
    plot_gradient() + facet_time
}
plot_stability <- function(data) {
  data |>
    filter(response %in% stability_vars) |>
    plot_gradient() + facet_time + scale_y_log10()
}

# The regimes of the main figures, and of the supplementary sets
in_main <- function(data) {
  filter(data, synchrony == "Asynchronous patches", predictability %in% main_kinds)
}
in_noise <- function(data) {
  filter(data, synchrony == "Asynchronous patches", predictability != "Predictable")
}
in_synchronous <- function(data) {
  filter(data, synchrony == "Synchronous patches", predictability %in% main_kinds)
}

active      <- filter(runs_long, community == focal_community)
matched     <- ratio_to_matched_every
matched_lab <- "Ratio to bet-hedging"

# ---------------------------------------------------------------------------
# 5. Main figures
# ---------------------------------------------------------------------------

# What does a seed bank do to diversity, and where does dispersal do the same?
ggsave(file.path(main_dir, "diversity.pdf"), plot_diversity(in_main(active)),
       width = 6.5, height = 6.5, bg = "white")

# How much of the environment's variation reaches growth, within patches and
# over the landscape? The variance through time of log growth of the propagule
# pool, in each patch (population) and over all patches (metapopulation), on
# one log axis so that the two levels can be compared.
growth_labels <- c(loglam_var_local = "Variance of log\npopulation growth",
                   loglam_var       = "Variance of log\nmetapopulation growth")
p <- runs |>
  filter(community == focal_community, g_mean %in% c(1, focal_g_mean)) |>
  in_main() |>
  pivot_longer(names(growth_labels), names_to = "response", values_to = "value") |>
  mutate(response_label = factor(growth_labels[response], levels = growth_labels)) |>
  plot_gradient() +
  facet_grid(response_label ~ predictability, switch = "y") +
  scale_y_log10(labels = scales::label_log())
ggsave(file.path(main_dir, "growth_variance.pdf"), p, width = 6.5, height = 4.75, bg = "white")

# What does it do to stability, and at which level?
ggsave(file.path(main_dir, "stability.pdf"), plot_stability(in_main(active)),
       width = 6.5, height = 8, bg = "white")

# Which synchrony does it change?
p <- in_main(active) |>
  filter(response %in% synchrony_vars) |>
  plot_gradient() + facet_time
ggsave(file.path(main_dir, "synchrony.pdf"), p, width = 6.5, height = 8, bg = "white")

# Does the way emergence is controlled matter, for the active and the dormant
# community? Each cued seed bank against the bet-hedger that emerges as much.
ggsave(file.path(main_dir, "strategy_layers.pdf"), plot_layers(in_main(matched), matched_lab),
       width = 6.5, height = 8, bg = "white")

# ---------------------------------------------------------------------------
# 6. Supplement: other regimes
# ---------------------------------------------------------------------------

# The color of the noise: the unpredictable regimes with red, white, and blue
# noise. A responsive seed bank reads the step just past, which tells it most
# under red noise, nothing under white, and the opposite of what follows under
# blue.
ggsave(file.path(supp_dir, "noise_strategy_layers.pdf"),
       plot_layers(in_noise(matched), matched_lab),
       width = 6.5, height = 8, bg = "white")

# Synchronous patches
ggsave(file.path(supp_dir, "synchronous_diversity.pdf"), plot_diversity(in_synchronous(active)),
       width = 6.5, height = 6.5, bg = "white")
ggsave(file.path(supp_dir, "synchronous_stability.pdf"), plot_stability(in_synchronous(active)),
       width = 6.5, height = 8, bg = "white")
ggsave(file.path(supp_dir, "synchronous_strategy_layers.pdf"),
       plot_layers(in_synchronous(matched), matched_lab),
       width = 6.5, height = 8, bg = "white")

# ---------------------------------------------------------------------------
# 7. Supplement: active and dormant combined, and the dormant community
# ---------------------------------------------------------------------------

# Runs without a seed bank have no dormant community, so only the seed banks
# are drawn for it. Without a seed bank, active and dormant combined is the
# active community.
dormant  <- filter(runs_long, community == "dormant", dormancy_type != "None")
combined <- filter(runs_long, community == "pool")

ggsave(file.path(supp_dir, "combined_stability.pdf"), plot_stability(in_main(combined)),
       width = 6.5, height = 8, bg = "white")
ggsave(file.path(supp_dir, "dormant_diversity.pdf"), plot_diversity(in_main(dormant)),
       width = 6.5, height = 6.5, bg = "white")

message("Wrote figures to ", fig_dir)
