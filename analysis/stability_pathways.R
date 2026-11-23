# Stability pathways: through which route does dormancy change metacommunity
# stability?
#
# Reads the sweep written by simulations/metacom_regimes.R. Wang et al.'s
# partition (Wang & Loreau 2016; Wang et al. 2019) divides every CV by the same
# mean total abundance, so in every run, exactly,
#
#   CV_C_R = CV_S_L * phi_S2C_L * phi_C_L2R      (species first, then patches)
#          = CV_S_L * phi_S_L2R * phi_S2C_R      (patches first, then species)
#
# Take the log ratio of a treatment to its reference. Then the log ratio of
# metacommunity stability (1 / CV_C_R) is the sum of three terms, with no
# remainder:
#
#   population     -dlog(CV_S_L)
#   among species  -(dlog(phi_S2C_L) + dlog(phi_S2C_R)) / 2
#   among patches  -(dlog(phi_C_L2R) + dlog(phi_S_L2R)) / 2
#
# Averaging the two paths removes the dependence on which aggregation comes
# first. A term above 1 (as a factor) means the treatment gains stability by
# that route relative to its reference.
#
# Contrasts:
#   each seed bank against none, for the active community. Runs without a seed
#     bank have no dormant community.
#   each cued seed bank against the bet-hedger matched in realized emergence
#     (analysis/matched_control.R), for the active and the dormant community.
#     Each log quantity is interpolated with the same weights, so the terms
#     still sum to the total.
#
# Outputs:
#   docs/seedbank_stability_pathways.csv                term means by contrast,
#                                                       community, regime, and
#                                                       dispersal rate
#   figures/supplement/pathways_seedbank.pdf  against none
#   figures/supplement/pathways_strategy.pdf  cued against the matched
#                                                       bet-hedger

library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
library(here)

source(here("analysis/matched_control.R"))   # the matched bet-hedger

# ---------------------------------------------------------------------------
# 1. Settings
# ---------------------------------------------------------------------------

# The sweep to read. The most recent one is used unless a directory is named.
sweep_dir <- NULL
if (is.null(sweep_dir)) {
  sweep_dir <- sort(Sys.glob(here("sim_output", "regimes_*")), decreasing = TRUE)[1]
}
message("Reading ", sweep_dir)

# The frame held fixed in the main figures
focal_g_mean   <- 0.1
main_scenarios <- c(Predictable = "predictable_asynchronous_color0",
                    `Red noise` = "unpredictable_asynchronous_color0.6")

fig_dir <- here("figures", "supplement")
doc_dir <- here("docs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Plotting theme, the same as in analysis/seedbank_stability.R
theme_set(theme_bw() +
            theme(strip.background = element_rect(fill = "white", color = "white"),
                  strip.text = element_text(size = 10),
                  strip.placement = "outside",
                  panel.grid.minor = element_blank()))

# Figure sizes, as in analysis/seedbank_stability.R: 6.5 in wide, and a
# height set by the number of rows
fig_width   <- 6.5
fig_heights <- c(`3` = 6.5, `4` = 8)

# Names and colors for the terms (Okabe-Ito, colorblind-safe). Each term is
# read as stability gained, so asynchrony rather than synchrony.
term_labels <- c(total         = "Metacommunity stability",
                 population    = "Population stability",
                 among_species = "Asynchrony among species",
                 among_patches = "Asynchrony among patches")
term_colors <- c("black", "#E69F00", "#56B4E9", "#CC79A7")
names(term_colors) <- term_labels

dormancy_names  <- c(none = "None", `bet-hedging` = "Bet-hedging",
                     anticipatory = "Anticipatory", responsive = "Responsive")
community_names <- c(active = "active", dormant = "dormant")

# ---------------------------------------------------------------------------
# 2. Read the sweep
# ---------------------------------------------------------------------------

runs <- read_csv(file.path(sweep_dir, "summary.csv"), show_col_types = FALSE)

# The quantities whose log ratios make up the partition, one row per run,
# community, and quantity
log_vars  <- c("CV_S_L", "CV_C_R", "phi_S2C_L", "phi_S2C_R", "phi_C_L2R", "phi_S_L2R")
pair_keys <- c("rep", "scenario", "dispersal_rate", "community", "quantity")

runs_long <- runs |>
  pivot_longer(all_of(log_vars), names_to = "quantity", values_to = "value")
treatments <- filter(runs_long, g_mean %in% c(1, focal_g_mean))

# ---------------------------------------------------------------------------
# 3. Pair each run with its reference and take log ratios
# ---------------------------------------------------------------------------

# Against no seed bank: the run with none in the same replicate, regime, and
# dispersal rate, which saw the same realized environment
no_bank <- treatments |>
  filter(dormancy == "none") |>
  select(all_of(pair_keys), value_reference = value)
to_none <- treatments |>
  filter(dormancy != "none", community == "active") |>
  left_join(no_bank, by = pair_keys) |>
  mutate(reference = "none")

# Against the matched bet-hedger, interpolated along the ladder of bet-hedging
# fractions to the cued run's realized emergence
ladder <- runs_long |>
  filter(dormancy == "bet-hedging") |>
  select(all_of(pair_keys), g_control = g_realized, value_control = value)
to_matched <- treatments |>
  filter(dormancy %in% c("anticipatory", "responsive"),
         community %in% c("active", "dormant")) |>
  match_bet_hedging(ladder, keys = pair_keys) |>
  select(-outside_ladder) |>
  mutate(reference = "bet-hedging")

# The log ratio of each quantity, one row per pair, and the terms of the
# partition. Each term is a log factor of stability, so a positive value is a
# gain.
paths <- bind_rows(to_none, to_matched) |>
  mutate(d = log(value / value_reference)) |>
  select(reference, dormancy, rep, scenario, dispersal_rate, community, quantity, d) |>
  pivot_wider(names_from = quantity, values_from = d, names_prefix = "d_") |>
  mutate(
    total         = -d_CV_C_R,
    population    = -d_CV_S_L,
    among_species = -(d_phi_S2C_L + d_phi_S2C_R) / 2,
    among_patches = -(d_phi_C_L2R + d_phi_S_L2R) / 2,
    residual = total - (population + among_species + among_patches),
    defined  = is.finite(total) & is.finite(among_species) & is.finite(among_patches)
  )

# ---------------------------------------------------------------------------
# 4. Check: the terms sum to the total in every pair
# ---------------------------------------------------------------------------

max_residual <- max(abs(paths$residual[paths$defined]))
message("Largest residual, total - sum of terms: ", signif(max_residual, 3))
if (max_residual > 1e-8) {
  stop("The terms do not sum to the total.")
}

# Pairs in which a log ratio is undefined (a run in which a community was lost)
lost <- paths |>
  filter(!defined) |>
  count(reference, dormancy, community, scenario, name = "pairs_undefined")
print(lost)

# ---------------------------------------------------------------------------
# 5. Means over replicates, written to docs/
# ---------------------------------------------------------------------------

# Log ratios are averaged over replicates (a geometric mean of the ratios), so
# the term means still sum to the total mean. Every path-specific log ratio is
# kept too, so the two orders of aggregation can be compared.
term_vars <- c("total", "population", "among_species", "among_patches",
               "d_phi_S2C_L", "d_phi_S2C_R", "d_phi_C_L2R", "d_phi_S_L2R")

pathway_means <- paths |>
  group_by(reference, dormancy, community, scenario, dispersal_rate) |>
  summarize(pairs = n(),
            pairs_defined = sum(defined),
            max_abs_residual = max(abs(residual[defined]), -Inf),
            across(all_of(term_vars), \(x) mean(x[defined])),
            .groups = "drop")
write_csv(pathway_means, file.path(doc_dir, "seedbank_stability_pathways.csv"))

# ---------------------------------------------------------------------------
# 6. Figures
# ---------------------------------------------------------------------------

# Long table of terms in the main frame: mean and SE of each log factor over
# replicates, drawn back on the factor scale against a log axis
main_terms <- paths |>
  filter(defined, scenario %in% main_scenarios) |>
  mutate(predictability = factor(names(main_scenarios)[match(scenario, main_scenarios)],
                                 levels = names(main_scenarios))) |>
  pivot_longer(names(term_labels), names_to = "term", values_to = "log_factor") |>
  group_by(reference, dormancy, community, predictability, dispersal_rate, term) |>
  summarize(mean = mean(log_factor),
            se   = sd(log_factor) / sqrt(n()),
            .groups = "drop") |>
  mutate(term = factor(term_labels[term], levels = term_labels),
         factor_mean = exp(mean),
         factor_low  = exp(mean - se),
         factor_high = exp(mean + se))

# Rows are contrasts, named on the left beside their values; the y-axis title
# names the reference. One y axis across rows, so that a given distance is the
# same change in every row, as in the other ratio figures.
plot_pathways <- function(data, y_label) {
  ggplot(data, aes(x = dispersal_rate, y = factor_mean, color = term, fill = term)) +
    geom_hline(yintercept = 1, linewidth = 0.3, color = "grey40") +
    geom_ribbon(aes(ymin = factor_low, ymax = factor_high), alpha = 0.15, color = NA,
                show.legend = FALSE) +
    geom_line(linewidth = 0.7) +
    facet_grid(row_label ~ predictability, switch = "y",
               labeller = labeller(row_label = label_wrap_gen(width = 18))) +
    scale_x_log10(labels = scales::label_log()) +
    scale_y_log10() +
    scale_color_manual(values = term_colors) +
    scale_fill_manual(values = term_colors) +
    labs(x = "Dispersal rate", y = y_label, color = NULL)
}

# What does each seed bank change, and by which route? Active community,
# against no seed bank.
p <- main_terms |>
  filter(reference == "none") |>
  mutate(row_label = factor(dormancy_names[dormancy], levels = dormancy_names[-1])) |>
  plot_pathways("Ratio to no dormancy")
ggsave(file.path(fig_dir, "pathways_seedbank.pdf"), p,
       width = fig_width, height = fig_heights["3"], bg = "white")

# Does following a cue change it, and by which route? Active and dormant
# communities, against the matched bet-hedger.
row_levels <- c("Anticipatory, active", "Anticipatory, dormant",
                "Responsive, active", "Responsive, dormant")
p <- main_terms |>
  filter(reference == "bet-hedging") |>
  mutate(row_label = factor(paste0(dormancy_names[dormancy], ", ",
                                   community_names[community]),
                            levels = row_levels)) |>
  plot_pathways("Ratio to bet-hedging")
ggsave(file.path(fig_dir, "pathways_strategy.pdf"), p,
       width = fig_width, height = fig_heights["4"], bg = "white")
