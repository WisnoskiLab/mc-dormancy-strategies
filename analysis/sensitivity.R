# Sensitivity of metacommunity stability to the parameters held fixed in the
# main sweep.
#
# Reads the sweep written by simulations/metacom_sensitivity.R, in which one
# parameter at a time is set below and above its value in the main sweep, in
# the setting of the main figures. For each parameter set it asks whether the
# two contrasts of the main text hold: each seed bank against no dormancy, and
# each cued seed bank against bet-hedging.
#
# Output:
#   figures/supplement/sensitivity.pdf
#       metacommunity stability (A) against no dormancy and (B) cued against
#       bet-hedging, averaged over the dispersal gradient, one column per
#       parameter; the manuscript's supplementary figure S07
#   docs/seedbank_sensitivity_summary.csv
#       the points of that figure

library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)
library(here)

# ---------------------------------------------------------------------------
# 1. Settings
# ---------------------------------------------------------------------------

# The sweep to read. The most recent is used unless a directory is named.
sensitivity_dir <- NULL
if (is.null(sensitivity_dir)) {
  sensitivity_dir <- sort(Sys.glob(here("sim_output", "sensitivity_*")), decreasing = TRUE)[1]
}
message("Reading ", sensitivity_dir)

# The community the figure is computed on
focal_community <- "active"

# Names for the parameters varied in the sensitivity sweep
parameter_labels <- c(g_mean        = "Emergence",
                      beta          = "Responsiveness",
                      survival      = "Survival in the bank",
                      niche_breadth = "Niche breadth")

supp_dir <- here("figures", "supplement")
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
# 2. Metacommunity stability of every run
# ---------------------------------------------------------------------------

# The parameter sets, with their values named `set_<parameter>`
parameter_sets <- read_csv(file.path(sensitivity_dir, "parameter_sets.csv"),
                           show_col_types = FALSE) |>
  rename_with(\(x) paste0("set_", x), all_of(names(parameter_labels)))

# One row per run. Metacommunity stability is the inverse of the
# metacommunity CV. The sets that vary emergence at low survival in the bank
# are left out; analysis/seedbank_checks.R examines them.
sens_meta <- read_csv(file.path(sensitivity_dir, "summary.csv"), show_col_types = FALSE) |>
  filter(community == focal_community, varied != "g_mean_low_survival") |>
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
    value = 1 / CV_C_R
  ) |>
  left_join(select(parameter_sets, set, starts_with("set_")), by = "set")

# ---------------------------------------------------------------------------
# 3. Paired ratios, averaged over the dispersal gradient
# ---------------------------------------------------------------------------

# Metacommunity stability as a ratio, paired within parameter set, replicate,
# regime, and dispersal rate: each seed bank against none, and each cued seed
# bank against the bet-hedging run of the same set. The sensitivity sweep has
# one bet-hedger per set, at the set's mean emergence fraction, so the second
# contrast is nominal; for the active community it differs little from the
# matched one, because the stability of a bet-hedger's active community
# changes little across the ladder of emergence fractions.
summary_keys <- c("set", "rep", "scenario", "dispersal_rate")

meta_none <- sens_meta |>
  filter(dormancy == "none") |>
  select(all_of(summary_keys), value_reference = value)
meta_bet_hedging <- sens_meta |>
  filter(dormancy == "bet-hedging") |>
  select(all_of(summary_keys), value_reference = value)

sens_pairs <- bind_rows(
  sens_meta |>
    filter(dormancy != "none") |>
    left_join(meta_none, by = summary_keys, relationship = "many-to-one") |>
    mutate(contrast = "to none"),
  sens_meta |>
    filter(dormancy %in% c("anticipatory", "responsive")) |>
    left_join(meta_bet_hedging, by = summary_keys, relationship = "many-to-one") |>
    mutate(contrast = "to bet-hedging")
) |>
  mutate(log_ratio = log(value / value_reference))

# Pairs in which either run lost every species
sens_pairs |>
  filter(!is.finite(log_ratio)) |>
  count(contrast, dormancy, predictability, set, name = "pairs_undefined") |>
  print()

# Every value a parameter takes, in increasing order, so that each column
# reads from low to high
value_levels <- sort(unique(unlist(select(parameter_sets, paste0("set_", names(parameter_labels))))))

# One row per pair and parameter, for the baseline and the sets that vary that
# parameter, with the value it took. The log ratio is averaged over the
# dispersal gradient within each replicate, then over replicates, so each
# point is a geometric mean with its standard error.
sens_summary <- sens_pairs |>
  filter(is.finite(log_ratio)) |>
  pivot_longer(paste0("set_", names(parameter_labels)), names_to = "parameter",
               names_prefix = "set_", values_to = "parameter_value") |>
  filter(varied == "baseline" | varied == parameter) |>
  group_by(contrast, dormancy_type, predictability, parameter, parameter_value, rep) |>
  summarize(log_ratio = mean(log_ratio), .groups = "drop") |>
  group_by(contrast, dormancy_type, predictability, parameter, parameter_value) |>
  summarize(replicates = n(),
            ratio      = exp(mean(log_ratio)),
            ratio_low  = exp(mean(log_ratio) - sd(log_ratio) / sqrt(n())),
            ratio_high = exp(mean(log_ratio) + sd(log_ratio) / sqrt(n())),
            .groups = "drop") |>
  mutate(parameter_label = factor(parameter_labels[parameter], levels = parameter_labels),
         parameter_value = factor(parameter_value, levels = value_levels))
write_csv(sens_summary, here("docs", "seedbank_sensitivity_summary.csv"))

# ---------------------------------------------------------------------------
# 4. Figure
# ---------------------------------------------------------------------------

# The value of each parameter in the main sweep, marked by a dashed line
baseline_values <- parameter_sets |>
  filter(varied == "baseline") |>
  pivot_longer(starts_with("set_"), names_to = "parameter", names_prefix = "set_",
               values_to = "parameter_value") |>
  filter(parameter %in% names(parameter_labels)) |>
  mutate(parameter_label = factor(parameter_labels[parameter], levels = parameter_labels),
         parameter_value = factor(parameter_value, levels = value_levels))

# One panel per contrast: columns are the parameters, each with the baseline
# and the two values around it
plot_sensitivity <- function(data, y_label) {
  dodge <- position_dodge(width = 0.4)
  ggplot(data, aes(x = parameter_value, y = ratio, color = dormancy_type,
                   group = dormancy_type)) +
    geom_hline(yintercept = 1, linewidth = 0.3, color = "grey40") +
    geom_line(position = dodge, linewidth = 0.5) +
    geom_pointrange(aes(ymin = ratio_low, ymax = ratio_high), position = dodge,
                    size = 0.3) +
    # the baseline line comes after the data, which set the order of the values
    geom_vline(data = baseline_values, aes(xintercept = parameter_value),
               linewidth = 0.3, linetype = "22", color = "grey60") +
    facet_grid(predictability ~ parameter_label, scales = "free_x", switch = "both") +
    scale_y_log10() +
    scale_color_manual(values = dormancy_colors) +
    labs(x = NULL, y = y_label, color = "Dormancy") +
    theme(strip.placement = "outside")
}

p_seedbank <- sens_summary |>
  filter(contrast == "to none") |>
  plot_sensitivity("Metacommunity stability,\nratio to no dormancy")
p_strategy <- sens_summary |>
  filter(contrast == "to bet-hedging") |>
  plot_sensitivity("Metacommunity stability,\nratio to bet-hedging")

# Both contrasts in one figure: (A) each seed bank against none, (B) each cued
# bank against bet-hedging. The colors are shared, so panel A's legend, which
# lists every kind of dormancy, serves both. It goes below, so the four columns
# have the full width.
p <- (p_seedbank + theme(legend.direction = "horizontal")) /
  (p_strategy + theme(legend.position = "none")) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A", theme = theme(legend.position = "bottom"))
ggsave(file.path(supp_dir, "sensitivity.pdf"), p, width = 6.5, height = 8, bg = "white")

message("Wrote ", file.path(supp_dir, "sensitivity.pdf"))
