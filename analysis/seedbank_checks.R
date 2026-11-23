# Checks on the results for dormancy, dispersal, and stability.
#
# Three questions about the figures drawn by analysis/seedbank_stability.R and
# analysis/sensitivity.R:
#
#   Do lost runs bias the means?  A run in which every species is lost has no
#     variability to measure, so it drops out of the stability and synchrony
#     means. Section 2 counts those runs for every treatment.
#   Do the contrasts hold within replicates?  Section 3 gives each contrast as
#     a ratio, paired within replicate, with its standard error and t value,
#     in three bands of the dispersal gradient.
#   Why do cued seed banks do better when survival in the bank is low?
#     Section 4 compares them, at low survival, with the bet-hedger matched in
#     realized emergence, from the sweep written by
#     simulations/metacom_sensitivity.R.
#
# Output:
#   docs/seedbank_runs_lost.csv         runs with no survivors, by treatment
#   docs/seedbank_paired_contrasts.csv  the paired contrasts (community "pool" is
#                                       active and dormant combined)
#   docs/seedbank_yield_per_active.csv  propagules per active individual, cued
#                                       against the matched bet-hedger
#   figures/supplement/sensitivity_low_survival.pdf
#                                       every kind of dormancy at each emergence
#                                       fraction, at low survival
#   docs/seedbank_low_survival_matched.csv
#                                       cued against the matched bet-hedger, at
#                                       low survival

library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
library(here)

source(here("analysis/matched_control.R"))   # the matched bet-hedger

# ---------------------------------------------------------------------------
# 1. Settings
# ---------------------------------------------------------------------------

# The sweeps to read. The most recent of each is used unless a directory is
# named.
sweep_dir <- NULL
if (is.null(sweep_dir)) {
  sweep_dir <- sort(Sys.glob(here("sim_output", "regimes_*")), decreasing = TRUE)[1]
}
sensitivity_dir <- NULL
if (is.null(sensitivity_dir)) {
  sensitivity_dir <- sort(Sys.glob(here("sim_output", "sensitivity_*")), decreasing = TRUE)[1]
}
message("Reading ", sweep_dir)
message("Reading ", sensitivity_dir)

# The setting of the main figures: asynchronous patches, one mean emergence
# fraction, and two kinds of variation through time
focal_g_mean   <- 0.1
main_scenarios <- c("predictable_asynchronous_color0",
                    "unpredictable_asynchronous_color0.6")

# The community the runs-lost count and the low-survival figure are computed on
focal_community <- "active"

# The communities whose paired contrasts are tabulated in section 3
contrast_communities <- c("active", "dormant", "pool")

# Bands of the dispersal gradient over which the contrasts are averaged
band_breaks <- c(0, 1e-3, 5e-2, 1)
band_names  <- c("low", "intermediate", "high")

# The parameters varied in the sensitivity sweep
parameter_labels <- c(g_mean        = "Emergence",
                      beta          = "Responsiveness",
                      survival      = "Survival in the bank",
                      niche_breadth = "Niche breadth")

# The survival at which the emergence fraction was also varied
low_survival <- 0.8

x_label <- "Dispersal rate"

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

# Readable names for the responses
response_labels <- c(
  alpha_div           = "Alpha diversity",
  gamma_div           = "Gamma diversity",
  stab_population     = "Population stability",
  stab_community      = "Community stability",
  stab_metapopulation = "Metapopulation stability",
  stab_metacommunity  = "Metacommunity stability"
)

# Labels shared by the two sweeps: the kind of dormancy, the kind of variation
# through time, the band of the dispersal gradient, and stability as the
# inverse of variability (CV) at each level
label_runs <- function(runs) {
  runs |>
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
      band = cut(dispersal_rate, breaks = band_breaks, labels = band_names),
      stab_population     = 1 / CV_S_L,
      stab_community      = 1 / CV_C_L,
      stab_metapopulation = 1 / CV_S_R,
      stab_metacommunity  = 1 / CV_C_R
    )
}

# A table of paired runs, with `value` and `value_reference`, as one ratio for
# each group: the mean log ratio over the dispersal rates in a band for each
# replicate, then the mean over replicates, its standard error, and t. A ratio
# is undefined when either run has no survivors, and `replicates` says how
# many replicates each row rests on.
summarize_paired <- function(pairs) {
  pairs |>
    mutate(log_ratio = log(value / value_reference)) |>
    filter(is.finite(log_ratio)) |>
    group_by(contrast, community, predictability, dormancy, response, band, rep) |>
    summarize(log_ratio = mean(log_ratio), .groups = "drop") |>
    group_by(contrast, community, predictability, dormancy, response, band) |>
    summarize(replicates = n(),
              ratio  = exp(mean(log_ratio)),
              se_log = sd(log_ratio) / sqrt(n()),
              t      = mean(log_ratio) / se_log,
              .groups = "drop") |>
    arrange(contrast, community, predictability, dormancy, response, band)
}

# ---------------------------------------------------------------------------
# 2. Runs with no survivors
# ---------------------------------------------------------------------------

runs_every <- read_csv(file.path(sweep_dir, "summary.csv"), show_col_types = FALSE) |>
  label_runs()
runs <- filter(runs_every, community == focal_community)

# A run is lost when no species is present in the recorded steps. Its
# diversity is zero and is part of every mean, but its stability and synchrony
# are undefined, so those means describe the surviving runs only.
runs_lost <- runs |>
  group_by(predictability, synchrony, dormancy, g_mean, band) |>
  summarize(runs = n(),
            lost = sum(gamma_div == 0),
            share_lost = lost / runs,
            .groups = "drop")
write_csv(runs_lost, here("docs", "seedbank_runs_lost.csv"))

message("Treatments that lost runs:")
runs_lost |>
  filter(lost > 0) |>
  as.data.frame() |>
  print(digits = 2)

# ---------------------------------------------------------------------------
# 3. Paired contrasts in the setting of the main figures
# ---------------------------------------------------------------------------

main_long <- runs_every |>
  filter(community %in% contrast_communities, scenario %in% main_scenarios) |>
  pivot_longer(all_of(names(response_labels)), names_to = "response", values_to = "value")

# Runs are paired when they share a replicate, a regime, and a dispersal rate,
# and so saw the same realized environment, and are compared within a community
pair_keys <- c("rep", "scenario", "dispersal_rate", "community", "response")

# The seed banks at the focal mean emergence fraction
focal_banks <- filter(main_long, dormancy != "none", g_mean == focal_g_mean)
cued_banks  <- filter(focal_banks, dormancy %in% c("anticipatory", "responsive"))

# Contrast 1: each seed bank against no seed bank
no_bank <- main_long |>
  filter(dormancy == "none") |>
  select(all_of(pair_keys), value_reference = value)
to_none <- focal_banks |>
  left_join(no_bank, by = pair_keys) |>
  mutate(contrast = "to no dormancy")

# Contrast 2: each cued seed bank against the bet-hedger matched in realized
# emergence, interpolated along the ladder of bet-hedging fractions. This is
# the comparison in the main text.
ladder <- main_long |>
  filter(dormancy == "bet-hedging") |>
  select(all_of(pair_keys), g_control = g_realized, value_control = value)
to_matched <- cued_banks |>
  match_bet_hedging(ladder, keys = pair_keys) |>
  select(-outside_ladder) |>
  mutate(contrast = "to bet-hedging, matched emergence")

# Contrast 3: each cued seed bank against bet-hedging at the same mean
# emergence fraction
nominal <- main_long |>
  filter(dormancy == "bet-hedging", g_mean == focal_g_mean) |>
  select(all_of(pair_keys), value_reference = value)
to_nominal <- cued_banks |>
  left_join(nominal, by = pair_keys) |>
  mutate(contrast = "to bet-hedging, same mean emergence")

# Contrast 4: active and dormant combined against the active community of the
# same run: how much more stable the seed bank and the active individuals are
# together than the active individuals alone
active_same_run <- main_long |>
  filter(community == "active") |>
  select(run, response, value_reference = value)
to_active <- main_long |>
  filter(community == "pool", dormancy != "none", g_mean == focal_g_mean) |>
  left_join(active_same_run, by = c("run", "response")) |>
  mutate(contrast = "active + dormant to active")

paired_contrasts <- summarize_paired(bind_rows(to_none, to_matched, to_nominal, to_active))
write_csv(paired_contrasts, here("docs", "seedbank_paired_contrasts.csv"))

# Propagules each active individual contributes to the next step's pool, per
# species and averaged over species, cued against the matched bet-hedger: do
# cued individuals, which emerge in better conditions, end up producing more?
yield <- read_csv(file.path(sweep_dir, "per_species.csv"), show_col_types = FALSE,
                  col_select = c(run, yield_per_active)) |>
  group_by(run) |>
  summarize(value = mean(yield_per_active, na.rm = TRUE))
yield_runs <- runs |>
  filter(scenario %in% main_scenarios) |>
  select(run, rep, scenario, dispersal_rate, dormancy, g_mean, g_realized,
         predictability, band) |>
  left_join(yield, by = "run") |>
  mutate(community = "active", response = "yield_per_active")
yield_pairs <- yield_runs |>
  filter(dormancy %in% c("anticipatory", "responsive"), g_mean == focal_g_mean) |>
  match_bet_hedging(yield_runs |>
                      filter(dormancy == "bet-hedging") |>
                      select(all_of(pair_keys), g_control = g_realized, value_control = value),
                    keys = pair_keys) |>
  mutate(contrast = "to bet-hedging, matched emergence")
write_csv(summarize_paired(yield_pairs), here("docs", "seedbank_yield_per_active.csv"))

# ---------------------------------------------------------------------------
# 4. Low survival: following a cue, or emerging more?
# ---------------------------------------------------------------------------

# The parameter sets, with their values named `set_<parameter>`. In the table
# of runs, `g_mean` and `beta` are the values a run used: 1 with no seed bank,
# and 0 for a seed bank that reads no cue.
parameter_sets <- read_csv(file.path(sensitivity_dir, "parameter_sets.csv"),
                           show_col_types = FALSE) |>
  rename_with(\(x) paste0("set_", x), all_of(names(parameter_labels)))

# One row per run and response. Realized emergence is kept beside the
# responses, because it is what the bet-hedger is matched on.
sens_long <- read_csv(file.path(sensitivity_dir, "summary.csv"), show_col_types = FALSE) |>
  filter(community == focal_community) |>
  label_runs() |>
  left_join(select(parameter_sets, set, starts_with("set_")), by = "set") |>
  pivot_longer(all_of(names(response_labels)), names_to = "response", values_to = "value") |>
  mutate(response_label = factor(response_labels[response], levels = response_labels))

# When few propagules survive in the bank, a bet-hedging seed bank that
# emerges a small fraction loses most of its propagules before they emerge,
# and the cued seed banks do better. A cued species also emerges more than its
# mean fraction, so it could do better for either reason. The sets at low
# survival separate the two.
low_long <- filter(sens_long, set_survival == low_survival)

# Every kind of dormancy at each emergence fraction, with realized emergence
# as the last row (without the runs that have no seed bank, where it is 1)
low_emergence <- low_long |>
  filter(response == "alpha_div", dormancy != "none") |>
  mutate(response = "g_realized", value = g_realized,
         response_label = "Realized emergence fraction")
low_values <- low_long |>
  filter(response %in% c("alpha_div", "gamma_div",
                         "stab_population", "stab_metacommunity")) |>
  mutate(response_label = as.character(response_label)) |>
  bind_rows(low_emergence) |>
  mutate(response_label = factor(response_label, levels = unique(response_label)),
         level = paste("Emergence", set_g_mean))

# A response along the dispersal gradient, as in the main figures. Lines are
# means over replicates and ribbons are +/- 1 SE.
plot_gradient <- function(data, y_label = NULL) {
  ggplot(data, aes(x = dispersal_rate, y = value, color = dormancy_type,
                   fill = dormancy_type)) +
    stat_summary(fun.data = mean_se, geom = "ribbon", alpha = 0.15, color = NA) +
    stat_summary(fun = mean, geom = "line", linewidth = 0.7) +
    scale_x_log10(labels = scales::label_log()) +
    scale_color_manual(values = dormancy_colors) +
    scale_fill_manual(values = dormancy_colors) +
    labs(x = x_label, y = y_label, color = "Dormancy", fill = "Dormancy")
}
wrap_rows <- labeller(response_label = label_wrap_gen(width = 18))

p <- plot_gradient(low_values) +
  facet_grid(response_label ~ predictability + level, scales = "free_y",
             switch = "y", labeller = wrap_rows)
ggsave(file.path(supp_dir, "sensitivity_low_survival.pdf"), p,
       width = 14, height = 9.5, bg = "white")

# Each cued seed bank at the focal mean emergence fraction against the
# bet-hedger interpolated to its realized emergence, along the bet-hedgers at
# every emergence fraction of the low-survival sets
low_keys <- c("rep", "scenario", "dispersal_rate", "response")
low_ladder <- low_long |>
  filter(dormancy == "bet-hedging") |>
  select(all_of(low_keys), g_control = g_realized, value_control = value)
low_matched <- low_long |>
  filter(dormancy %in% c("anticipatory", "responsive"), set_g_mean == focal_g_mean) |>
  match_bet_hedging(low_ladder, keys = low_keys) |>
  mutate(contrast = "to bet-hedging, matched emergence, low survival")

# What was matched: the cued seed bank's realized emergence, and the share of
# its runs that fall outside the ladder
low_matching <- low_matched |>
  filter(response == "alpha_div") |>
  group_by(predictability, dormancy, band) |>
  summarize(g_realized     = mean(g_realized),
            outside_ladder = mean(outside_ladder),
            .groups = "drop")

low_contrasts <- summarize_paired(select(low_matched, -outside_ladder)) |>
  left_join(low_matching, by = c("predictability", "dormancy", "band"))
write_csv(low_contrasts, here("docs", "seedbank_low_survival_matched.csv"))

message("Wrote checks to ", here("docs"), " and figures to ", supp_dir)
