# Dormancy, dispersal, and stability: sensitivity to parameters held fixed.
#
# The main sweep (simulations/metacom_regimes.R) holds four parameters at one
# value each. Do its results depend on those values? Here each one is varied in
# turn while the others stay at the value of the main sweep:
#   - mean emergence fraction of the seed banks (how large the seed bank is)
#   - responsiveness of the cued seed banks to their cue
#   - survival in the bank
#   - niche breadth (how often conditions fall below replacement)
#
# The bet-hedging seed bank fails when survival in the bank is low. To ask
# whether the cued seed banks do better there because they follow a cue, or
# only because they emerge more, the emergence fraction is also varied at the
# lowest survival. The bet-hedging runs of those sets are the controls: a cued
# seed bank can be compared with the one that emerges as much as it does.
#
# Everything else is the main sweep at a lower resolution, in the setting of
# its main figures: asynchronous patches and two kinds of variation through
# time (predictable, and red noise). The replicates are the first replicates of
# the main sweep, with the same landscapes, patch means, competition
# coefficients, and realized environments, and each run takes the seed of the
# main-sweep runs in the same cell, so the baseline repeats the main sweep run
# for run and every parameter set is run on the same draws. The life cycle in
# section 5 is that of the main sweep, and the same four communities are
# recorded: active, dormant, the pool, and arriving propagules.
#
# Output, in sim_output/sensitivity_<timestamp>/:
#   summary.csv         one row per run and community: diversity and
#                       variability partitions, and the run's realized
#                       emergence
#   parameter_sets.csv  the parameter sets; in summary.csv, `g_mean` and `beta`
#                       are the values a run used (1 with no seed bank, 0 for
#                       a seed bank that reads no cue)
#   run_settings.csv    the settings below

library(tidyr)
library(dplyr)
library(tibble)
library(readr)
library(foreach)
library(doParallel)
library(here)

source(here("simulations/metacom_functions.R"))    # landscape, species, reproduction
source(here("simulations/env_functions.R"))        # environments and regimes
source(here("simulations/emergence_functions.R"))  # fixed and cued emergence
source(here("simulations/dispersal_functions.R"))  # dispersal kernels
source(here("simulations/seedbank_functions.R"))   # survival in the bank
source(here("analysis/diversity_partitioning.R"))  # alpha, beta, gamma
source(here("analysis/metacommunity_variability_partitioning.R"))  # CV and synchrony

# ---------------------------------------------------------------------------
# 1. Settings
# ---------------------------------------------------------------------------

# TRUE runs a small version that finishes in under a minute, to check that the
# script works. FALSE runs the full sweep.
test_run <- FALSE

# The parameters that are varied. Each row is one set of values. The first is
# the baseline (the main sweep), and each of the next eight differs from it in
# one parameter, named in `varied`. The last three vary the emergence fraction
# at the lowest survival.
#   g_mean         mean emergence fraction of the seed banks
#   beta           how strongly a cued seed bank responds to its cue
#   survival       survival in the bank, s
#   niche_breadth  width of each species' niche, sigma
parameter_sets <- tribble(
  ~varied,                ~g_mean, ~beta, ~survival, ~niche_breadth,
  "baseline",              0.1,     1.5,   0.95,      0.1,
  "g_mean",                0.3,     1.5,   0.95,      0.1,
  "g_mean",                0.5,     1.5,   0.95,      0.1,
  "beta",                  0.1,     0.5,   0.95,      0.1,
  "beta",                  0.1,     3.0,   0.95,      0.1,
  "survival",              0.1,     1.5,   0.80,      0.1,
  "survival",              0.1,     1.5,   0.99,      0.1,
  "niche_breadth",         0.1,     1.5,   0.95,      0.075,
  "niche_breadth",         0.1,     1.5,   0.95,      0.15,
  "g_mean_low_survival",   0.2,     1.5,   0.80,      0.1,
  "g_mean_low_survival",   0.3,     1.5,   0.80,      0.1,
  "g_mean_low_survival",   0.5,     1.5,   0.80,      0.1
)
parameter_sets$set <- 1:nrow(parameter_sets)

# Dormancy (after Lennon et al. 2021), and what each type reads:
#   "signal"    the predictable part of the coming step's conditions
#   "previous"  the realized conditions of the previous step
# Every type is run under every parameter set, so each set has its own run
# with no seed bank and its own bet-hedging run to compare against. With no
# seed bank every propagule emerges, and only the cued types use `beta`.
dormancy_types <- tribble(
  ~dormancy,       ~cue,
  "none",          "none",
  "bet-hedging",   "none",
  "anticipatory",  "signal",
  "responsive",    "previous"
)

# The regimes of the main figures
main_scenarios <- c("predictable_asynchronous_color0",
                    "unpredictable_asynchronous_color0.6")

# Landscape and community, as in the main sweep
x_dim   <- 100
y_dim   <- 100
patches <- 100
species <- 10
optima_span  <- c(0.25, 0.75)
niche_optima <- seq(optima_span[1], optima_span[2], length.out = species)
f_max        <- 5   # propagules per active individual at the niche optimum

# Length of a run, in time steps, as in the main sweep. The community first
# assembles under constant conditions (initialization), then experiences the
# fluctuating environment. Only the last `timesteps` steps are recorded.
initialization <- 100
burn_in        <- 300
timesteps      <- 300

# Replicate `k` is built from the seed `replicate_seed + k`, as in the main
# sweep, so these are its first `nreps` replicates
nreps          <- 5
replicate_seed <- 20260807

# Environment, as in the main sweep
cycle_length <- 10
w <- (2 * pi) / cycle_length
patch_mean_span  <- c(0.25, 0.75)
patch_mean_range <- 10
noise_colors <- c(-0.6, 0, 0.6)

# Dispersal: every second dispersal rate of the main sweep, and its kernel
dispersal_grid  <- 10^seq(-5, 0, length.out = 20)   # the main sweep's rates
dispersal_rates <- dispersal_grid[seq(2, 20, by = 2)]
kernel_decay    <- 0.1
exclude_self    <- TRUE   # emigrants cannot land in their home patch

if (test_run) {
  nreps           <- 2
  dispersal_rates <- dispersal_rates[c(2, 6, 10)]
  parameter_sets  <- parameter_sets[c(1, 3, 6, 11), ]
}

# Number of parallel workers: all but one local core, or what SLURM allocated
ncores <- parallel::detectCores() - 1
if (Sys.getenv("SLURM_CPUS_PER_TASK") != "") {
  ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
}

# ---------------------------------------------------------------------------
# 2. Replicates: landscape, patch means, and competition
# ---------------------------------------------------------------------------

# A replicate is one draw of everything random about the setting. Every
# parameter set and dormancy type is run on every replicate, so they are
# compared against the same background.
landscapes       <- list()
alphas           <- list()
patch_means_list <- list()

for (rep_i in 1:nreps) {
  set.seed(replicate_seed + rep_i)

  landscapes[[rep_i]] <- init_landscape(patches = patches, x_dim = x_dim, y_dim = y_dim)

  # Competition coefficients: intraspecific 1, interspecific drawn between 0
  # and 1, all scaled by 0.05, so species limit themselves more than they limit
  # each other
  alphas[[rep_i]] <- competition_matrix(species = species, intra = 1, min_inter = 0,
                                        max_inter = 1, scale = 0.05)

  # Long-term mean conditions of each patch
  patch_means_list[[rep_i]] <- patch_mean_field(landscapes[[rep_i]],
                                                lower = patch_mean_span[1],
                                                upper = patch_mean_span[2],
                                                spat_range = patch_mean_range,
                                                torus = TRUE)
}

# ---------------------------------------------------------------------------
# 3. Environments
# ---------------------------------------------------------------------------

# The table of regimes is built as in the main sweep and numbered before the
# main-figure regimes are picked out, because the number sets the seed of each
# environment and of each run. The environments here are therefore those of
# the main sweep.
regimes <- env_scenarios(extent = x_dim, colors = noise_colors) |>
  filter(pred_frac < 0.5 | temp_noise_color == 0)
regimes$regime_id <- 1:nrow(regimes)
regimes <- filter(regimes, scenario %in% main_scenarios)

env_list <- list()

for (i in 1:nrow(regimes)) {
  for (rep_i in 1:nreps) {
    regime <- regimes[i, ]
    env_list[[paste(regime$scenario, rep_i)]] <- env_generate(
      landscape          = landscapes[[rep_i]],
      timesteps          = burn_in + timesteps,
      patch_means        = patch_means_list[[rep_i]],
      amplitude          = regime$amplitude,
      w                  = w,
      phase_mode         = regime$phase_mode,
      phase_spread       = regime$phase_spread,
      temp_noise_sd      = regime$temp_noise_sd,
      temp_noise_color   = regime$temp_noise_color,
      spat_range         = regime$spat_range,
      torus              = TRUE,
      scale_half         = regime$scale_half,
      seed               = 4000 + 97 * regime$regime_id + rep_i
    )
  }
}

# ---------------------------------------------------------------------------
# 4. Design: one row per run
# ---------------------------------------------------------------------------

design <- expand_grid(
  rep            = 1:nreps,
  scenario       = regimes$scenario,
  parameter_sets,
  dormancy_types,
  dispersal_rate = dispersal_rates
)

# With no seed bank every propagule emerges, and a seed bank that reads no cue
# does not respond to one
design$g_mean[design$dormancy == "none"] <- 1
design$beta[design$cue == "none"]        <- 0

design$kernel_decay <- kernel_decay
design$noise_color  <- regimes$temp_noise_color[match(design$scenario, regimes$scenario)]
design$run <- 1:nrow(design)

# The seed of the main-sweep runs in the same cell (replicate, regime, and
# dispersal rate), so the baseline repeats them and every parameter set is
# run on the same draws
design$regime_id       <- regimes$regime_id[match(design$scenario, regimes$scenario)]
design$dispersal_index <- match(design$dispersal_rate, dispersal_grid)
design$seed <- 100000 * design$rep + 100 * design$regime_id + design$dispersal_index
design$replicate_seed <- replicate_seed + design$rep   # rebuilds the replicate's setting

# How long a propagule stays in the bank
design$bank_residence <- bank_residence(design$g_mean, design$survival)

message("Running ", nrow(design), " simulations on ", ncores, " cores")

# ---------------------------------------------------------------------------
# 5. Simulations
# ---------------------------------------------------------------------------

registerDoParallel(cores = ncores)

results <- foreach(run = 1:nrow(design), .errorhandling = "pass") %dopar% {

  set.seed(design$seed[run])

  dispersal_rate <- design$dispersal_rate[run]
  g_mean         <- design$g_mean[run]
  beta           <- design$beta[run]
  survival       <- design$survival[run]
  niche_breadth  <- design$niche_breadth[run]
  cue            <- design$cue[run]
  env_df         <- env_list[[paste(design$scenario[run], design$rep[run])]]

  # The replicate's landscape, patch means, and competition coefficients
  landscape   <- landscapes[[design$rep[run]]]
  patch_means <- patch_means_list[[design$rep[run]]]
  alpha       <- alphas[[design$rep[run]]]

  # --- species: identical dormancy and dispersal, niche optima spread evenly --
  species_traits <- init_species(niche_optima   = niche_optima,
                                 niche_breadth  = niche_breadth,
                                 f_max          = f_max,
                                 dispersal_rate = dispersal_rate,
                                 kernel_decay   = kernel_decay,
                                 g_mean         = g_mean,
                                 survival       = survival)

  disp_array <- generate_dispersal_matrices(landscape, species, patches,
                                            species_traits, torus = TRUE,
                                            exclude_self = exclude_self)

  # --- emergence cue ------------------------------------------------------------
  # As in the main sweep: the cue is the species' fecundity under the
  # predictable part of the coming step's conditions (anticipatory) or under
  # the previous step's conditions (responsive), standardized over every patch
  # and step of the fluctuating environment, with mean emergence fixed at
  # `g_mean` for any `beta`.
  cue_values <- env_df$env_signal
  if (cue == "previous") {
    cue_values <- env_df$env
  }
  cue_all <- list()
  for (i in 1:species) {
    mismatch <- (species_traits$niche_optimum[i] - cue_values) /
      (2 * species_traits$niche_breadth[i])
    cue_all[[i]] <- species_traits$f_max[i] * exp(-mismatch^2)
  }
  cue_setup <- emergence_setup(cue_all, g_mean = species_traits$g_mean, beta = beta)

  # --- starting state: arriving propagules and an empty seed bank ---------------
  R      <- init_community(species, patches) + 1
  D_star <- R * 0

  # --- storage for the recorded steps --------------------------------------------
  # One census per step, after emergence: active individuals and dormant
  # propagules divide the pool without counting any propagule twice. The
  # arriving propagules are recorded as a fourth community.
  active_record   <- array(0, dim = c(species, timesteps, patches))
  dormant_record  <- array(0, dim = c(species, timesteps, patches))
  arriving_record <- array(0, dim = c(species, timesteps, patches))
  pool_record     <- array(0, dim = c(species, timesteps, patches))

  # Conditions of the step just past, which a responsive species reads
  env_previous <- patch_means

  for (step in 1:(initialization + burn_in + timesteps)) {

    # While the community assembles, conditions in each patch are held at its
    # long-term mean. After that they fluctuate.
    t_env <- max(1, step - initialization)
    env_now <- env_df$env[env_df$time == t_env]
    cue_now <- env_df$env_signal[env_df$time == t_env]
    if (cue == "previous") {
      cue_now <- env_previous
    }
    if (step <= initialization) {
      env_now <- patch_means
      cue_now <- patch_means
    }

    # Propagules per active individual without competition, under the
    # realized conditions (f) and under the cue (f_cue)
    f     <- fecundity(species_traits, env_now)
    f_cue <- fecundity(species_traits, cue_now)

    # 1. Emergence: a fraction of the pool becomes active; the rest stay dormant
    P <- R + D_star
    A <- emergence(P, f_cue, cue_setup)
    D <- P - A

    # 2. Reproduction, limited by the environment and then by competition
    B_expected <- expected_births(A, f, alpha)

    # 3. Dormant propagules survive in the bank
    D_star <- bank_survival(D, survival)

    # 4. Demographic stochasticity in reproduction
    B <- matrix(rpois(species * patches, lambda = pmax(B_expected, 0)),
                nrow = patches, ncol = species)

    # 5. Dispersal of the new propagules. Dormant propagules do not disperse.
    E <- matrix(0, patches, species)
    for (i in 1:species) {
      E[, i] <- rbinom(patches, size = B[, i], prob = dispersal_rate)
    }
    I <- matrix(0, patches, species)
    for (i in 1:species) {
      if (sum(E[, i]) > 0) {
        expected_arrivals <- disp_array[, , i] %*% E[, i]
        I[, i] <- rmultinom(1, size = sum(E[, i]), prob = as.vector(expected_arrivals))
      }
    }
    R <- B - E + I
    env_previous <- env_now

    # Record once the burn-in is over
    t_rec <- step - initialization - burn_in
    if (t_rec > 0) {
      active_record[, t_rec, ]   <- t(A)
      dormant_record[, t_rec, ]  <- t(D)
      arriving_record[, t_rec, ] <- t(R)
      pool_record[, t_rec, ]     <- t(A + D)
    }
  }

  # --- community summaries: diversity, variability, synchrony --------------------
  communities <- list(active   = active_record,
                      dormant  = dormant_record,
                      arriving = arriving_record,
                      pool     = pool_record)

  # Realized emergence of the run, as in the main sweep
  g_realized_run <- sum(active_record) / sum(pool_record)

  run_summary <- list()
  for (community in names(communities)) {
    N_record    <- communities[[community]]
    diversity   <- diversity.partition(N_record)   # alpha, beta, gamma
    beta_parts  <- beta.div.calc(N_record)         # spatial and temporal beta
    variability <- var.partition(N_record)         # CVs and synchrony (Wang et al. 2019)

    run_summary[[community]] <- tibble(
      run           = run,
      community     = community,
      alpha_div     = diversity[1],
      beta_div      = diversity[2],
      gamma_div     = diversity[3],
      beta_spatial  = beta_parts$mean_beta_spatial,
      beta_temporal = beta_parts$mean_beta_temporal,
      CV_S_L    = variability["CV_S_L"],      # variability: species, local
      CV_C_L    = variability["CV_C_L"],      #              community, local
      CV_S_R    = variability["CV_S_R"],      #              species, regional
      CV_C_R    = variability["CV_C_R"],      #              community, regional
      phi_S_L2R = variability["phi_S_L2R"],   # synchrony: spatial, of species
      phi_C_L2R = variability["phi_C_L2R"],   #            spatial, of communities
      phi_S2C_L = variability["phi_S2C_L"],   #            among species, local
      phi_S2C_R = variability["phi_S2C_R"],   #            among species, regional
      g_realized = g_realized_run
    )
  }
  bind_rows(run_summary)
}

stopImplicitCluster()

# ---------------------------------------------------------------------------
# 6. Collect and save
# ---------------------------------------------------------------------------

# A run that failed returns its error and not a table of results. Keep its row
# in the output, flagged, so a missing result is never mistaken for a run that
# was not in the design.
ran <- !sapply(results, inherits, what = "error")
message("Completed ", sum(ran), " of ", nrow(design), " simulations")
if (any(!ran)) {
  print(results[!ran][[1]])
}

# One row per run and community, so a failed run is a flagged row in each
design$ran <- ran
summaries  <- expand_grid(design, community = c("active", "dormant", "arriving", "pool")) |>
  left_join(bind_rows(results[ran]), by = c("run", "community"))

outdir <- here("sim_output", paste0("sensitivity_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(outdir, recursive = TRUE)

write_csv(summaries,      file.path(outdir, "summary.csv"))
write_csv(parameter_sets, file.path(outdir, "parameter_sets.csv"))

settings <- tibble(test_run, patches, species,
                   optima_lower = optima_span[1], optima_upper = optima_span[2],
                   initialization, burn_in, timesteps, nreps, cycle_length,
                   exclude_self, kernel_decay, replicate_seed, f_max,
                   patch_mean_lower = patch_mean_span[1],
                   patch_mean_upper = patch_mean_span[2], patch_mean_range)
write_csv(settings, file.path(outdir, "run_settings.csv"))

message("Wrote output to ", outdir)
