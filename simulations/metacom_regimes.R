# Dormancy, dispersal, and stability across scales: the main sweep.
#
# How do seed banks of different types change diversity, stability, and
# synchrony along a gradient of dispersal rate, and does the answer depend on
# how predictable and how spatially synchronous the environment is?
#
# Within a run, all species share one dormancy strategy and one dispersal
# kernel. They differ only in their niche optimum. Runs cross:
#   - dispersal rate (the probability that a new propagule emigrates)
#   - dormancy (after Lennon et al. 2021): none; bet-hedging (a fixed emergence
#     fraction); anticipatory (emergence follows the predictable part of the
#     coming step's conditions); responsive (emergence follows the conditions
#     of the previous step). All at a mean emergence fraction of 0.1, with
#     bet-hedging also at a ladder of fractions, the controls that emerge as
#     much as a cued seed bank does.
#   - environmental regime: predictable or unpredictable in time, crossed with
#     synchronous or asynchronous among patches; the unpredictable regimes
#     with blue, white, or red noise
#
# Notation follows the manuscript. In each step, the propagule pool P of a
# species in a patch divides into active individuals A and dormant propagules
# D; the active individuals produce propagules B, of which E emigrate and I
# immigrate, giving the arriving propagules R; dormant propagules survive in
# the bank (D_star); and R + D_star is the next step's pool.
#
# Every run is censused once per step, after emergence. Diversity, variability,
# and synchrony are computed for four communities: the active individuals, the
# dormant propagules, the two together (the pool at the census), and the
# arriving propagules.
#
# Output, in sim_output/regimes_<timestamp>/:
#   summary.csv          one row per run and community: diversity and
#                        variability partitions, and the run's realized
#                        emergence
#   per_species.csv      one row per species per run: Cohen's mixture, the
#                        variance of log growth of the propagule pool in each
#                        patch and over all patches, and abundance and
#                        occupancy in each community
#   env_diagnostics.csv  realized properties of each environment
#   patches.csv          each replicate's patches: position and long-term mean
#   competition.csv      each replicate's competition coefficients
#   run_settings.csv     the settings below

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
source(here("simulations/seedbank_functions.R"))   # survival in the bank, Cohen's mixture
source(here("analysis/diversity_partitioning.R"))  # alpha, beta, gamma
source(here("analysis/metacommunity_variability_partitioning.R"))  # CV and synchrony

# ---------------------------------------------------------------------------
# 1. Settings
# ---------------------------------------------------------------------------

# TRUE runs a small version that finishes in under a minute, to check that the
# script works. FALSE runs the full sweep.
test_run <- FALSE

# Landscape and community
x_dim   <- 100
y_dim   <- 100
patches <- 100
species <- 10

# Niches. Optima are spread evenly over the same span as the patches'
# long-term means (below), so every species is the best-adapted somewhere. At
# a breadth of 0.1, an active individual produces fewer than one propagule when
# conditions are more than 0.25 from its optimum, which is about two standard
# deviations of the fluctuations within a patch.
niche_breadth <- 0.1
optima_span   <- c(0.25, 0.75)
niche_optima  <- seq(optima_span[1], optima_span[2], length.out = species)
f_max         <- 5   # propagules per active individual at the niche optimum

# Length of a run, in time steps. The community first assembles under constant
# conditions (initialization), then experiences the fluctuating environment.
# Only the last `timesteps` steps are recorded.
initialization <- 100
burn_in        <- 300
timesteps      <- 300
nreps          <- 10

# Replicate `k` is built from the seed `replicate_seed + k`
replicate_seed <- 20260807

# Environment: the periodic layer has a cycle of 10 time steps
cycle_length <- 10
w <- (2 * pi) / cycle_length

# Each patch has its own long-term mean conditions, and fluctuates around it.
# The means are spread evenly over `patch_mean_span`, chosen so that conditions
# vary about as much among patches as they do through time within a patch.
# They are arranged in space as a random field in which patches closer than
# about `patch_mean_range` are alike.
patch_mean_span  <- c(0.25, 0.75)
patch_mean_range <- 10

# Color of the noise layer: the correlation between its values in successive
# time steps. 0 is white noise, in which a bad year says nothing about the
# next; 0.6 is red noise, in which bad years come in runs; -0.6 is blue noise,
# in which a bad year tends to follow a good one. The color is varied in the
# unpredictable regimes, where the noise is most of the variation. The
# predictable regimes keep white noise.
noise_colors <- c(-0.6, 0, 0.6)

# Dispersal: dispersal rates, and the rate at which the kernel decays with
# distance. Seeds are indexed on `dispersal_grid`, so a test run, which uses a
# subset of it, draws the same numbers as the full sweep.
dispersal_grid  <- 10^seq(-5, 0, length.out = 20)
dispersal_rates <- dispersal_grid
kernel_decay    <- 0.1
exclude_self    <- TRUE   # emigrants cannot land in their home patch

# Dormancy. `g_mean` is the mean emergence fraction; with g_mean = 1 every
# propagule emerges, so there is no seed bank. `beta` is how strongly emergence
# responds to a cue (0 = not at all), and `cue` is what the species reads:
#   "signal"    the predictable part of the coming step's conditions
#   "previous"  the realized conditions of the previous step
# A cued species emerges most where conditions suit it, so it realizes a
# different fraction than its mean of 0.1: more at low and intermediate
# dispersal rates, less at the highest. The bet-hedgers from 0.05 to 0.25
# bracket the realized fractions, so that each cued run can be compared with
# the bet-hedger interpolated to the fraction it realized
# (analysis/matched_control.R).
dormancy_types <- tribble(
  ~dormancy,       ~g_mean, ~beta, ~cue,
  "none",           1.0,     0,    "none",
  "bet-hedging",    0.05,    0,    "none",
  "bet-hedging",    0.075,   0,    "none",
  "bet-hedging",    0.1,     0,    "none",
  "bet-hedging",    0.125,   0,    "none",
  "bet-hedging",    0.15,    0,    "none",
  "bet-hedging",    0.175,   0,    "none",
  "bet-hedging",    0.2,     0,    "none",
  "bet-hedging",    0.25,    0,    "none",
  "anticipatory",   0.1,     1.5,  "signal",
  "responsive",     0.1,     1.5,  "previous"
)
survival <- 0.95   # survival in the bank, s

# Variance of log growth in each patch. Growth is undefined in a step where a
# patch holds none of a species' propagules, so a patch counts toward a
# species' local growth variance only if its growth is defined in at least
# this share of the recorded steps.
local_growth_min_share <- 0.5

if (test_run) {
  patches <- 25
  nreps   <- 2
  initialization  <- 20
  burn_in         <- 50
  timesteps       <- 30
  dispersal_rates <- dispersal_grid[c(4, 9, 14, 20)]
}

# Number of parallel workers: all but one local core, or what SLURM allocated
ncores <- parallel::detectCores() - 1
if (Sys.getenv("SLURM_CPUS_PER_TASK") != "") {
  ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
}

# ---------------------------------------------------------------------------
# 2. Replicates: landscape, patch means, and competition
# ---------------------------------------------------------------------------

# A replicate is one draw of everything random about the setting: where the
# patches are, their long-term mean conditions, how strongly the species
# compete, and (in section 3) how conditions fluctuate. Every treatment is run
# on every replicate, so dispersal and dormancy are compared against the same
# background, and replicates show whether a pattern depends on one draw.
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

  # Long-term mean conditions of each patch, the same in every regime, so that
  # regimes differ only in how conditions fluctuate
  patch_means_list[[rep_i]] <- patch_mean_field(landscapes[[rep_i]],
                                                lower = patch_mean_span[1],
                                                upper = patch_mean_span[2],
                                                spat_range = patch_mean_range,
                                                torus = TRUE)
}

# ---------------------------------------------------------------------------
# 3. Environments
# ---------------------------------------------------------------------------

# The regimes: predictable or unpredictable, synchronous or asynchronous. The
# unpredictable regimes are built at each noise color, the predictable ones
# with white noise only.
regimes <- env_scenarios(extent = x_dim, colors = noise_colors) |>
  filter(pred_frac < 0.5 | temp_noise_color == 0)

# One environment per regime and replicate, built on that replicate's
# landscape and patch means. Every run in a replicate sees the same realized
# environment, so treatments are compared on identical conditions. Conditions
# are held at the patch means while the community assembles, so the
# environment covers the burn-in and the recorded steps.
env_list <- list()
env_diag <- list()

for (i in 1:nrow(regimes)) {
  for (rep_i in 1:nreps) {
    regime <- regimes[i, ]
    env_df <- env_generate(
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
      seed               = 4000 + 97 * i + rep_i
    )

    env_name <- paste(regime$scenario, rep_i)
    env_list[[env_name]] <- env_df
    env_diag[[env_name]] <- cbind(scenario = regime$scenario, rep = rep_i,
                                  env_diagnostics(env_df, w = w))
  }
}
env_diag <- bind_rows(env_diag)

# Check that the regimes are what they claim to be
env_diag |>
  group_by(scenario) |>
  summarize(across(c(predictability, mean_patch_cor, rank_persistence,
                     decorr_time, temporal_cv_mean, spatial_cv_mean),
                   \(x) mean(x, na.rm = TRUE))) |>
  as.data.frame() |>
  print(digits = 3)

# ---------------------------------------------------------------------------
# 4. Design: one row per run
# ---------------------------------------------------------------------------

design <- expand_grid(
  rep            = 1:nreps,
  scenario       = regimes$scenario,
  dormancy_types,
  dispersal_rate = dispersal_rates
)
design$kernel_decay <- kernel_decay
design$noise_color  <- regimes$temp_noise_color[match(design$scenario, regimes$scenario)]
design$run <- 1:nrow(design)

# Each run's seed is set by its replicate, regime, and dispersal rate, so every
# dormancy type in a cell starts from the same draws, and the sensitivity
# sweep (metacom_sensitivity.R) reuses them: its baseline repeats these runs.
design$regime_id       <- match(design$scenario, regimes$scenario)
design$dispersal_index <- match(design$dispersal_rate, dispersal_grid)
design$seed <- 100000 * design$rep + 100 * design$regime_id + design$dispersal_index
design$replicate_seed <- replicate_seed + design$rep   # rebuilds the replicate's setting

# How long a propagule stays in the bank, and how far emigrants actually travel
# on the replicate's set of patches. Both describe the treatments better than
# the raw parameters do.
design$bank_residence <- bank_residence(design$g_mean, survival)
design$distance_realized <- NA_real_
for (rep_i in 1:nreps) {
  in_rep <- design$rep == rep_i
  design$distance_realized[in_rep] <-
    realized_dispersal_distance(landscapes[[rep_i]], design$kernel_decay[in_rep],
                                torus = TRUE, exclude_self = exclude_self)
}

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
  # A seed bank that responds to a cue emerges more when the cue says its
  # species would do well: the cue is the species' fecundity under some reading
  # of the conditions. An anticipatory species reads the predictable part of
  # the coming step's conditions (the seasonal signal, without the noise). A
  # responsive species reads the conditions of the step just past, noise
  # included, which say something about the coming step only if conditions
  # persist. The setup standardizes the cue over every patch and step of the
  # fluctuating environment and fixes mean emergence at `g_mean` for any
  # `beta`, so every kind of seed bank emerges the same fraction on average
  # over the landscape.
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
  # Each step is censused once, after emergence. Every propagule in a patch is
  # then either active or dormant, so the two communities divide the pool
  # without counting any propagule twice. The arriving propagules are recorded
  # as a fourth community; they are counted again, as active or dormant, at the
  # next step's census.
  active_record   <- array(0, dim = c(species, timesteps, patches))
  dormant_record  <- array(0, dim = c(species, timesteps, patches))
  arriving_record <- array(0, dim = c(species, timesteps, patches))
  pool_record     <- array(0, dim = c(species, timesteps, patches))
  cohen_record <- list()
  # Log growth of each patch's propagule pool (this step's pool to the next)
  # and the pool itself, for the variance of local growth
  loglam_local_record <- array(NA_real_, dim = c(timesteps, patches, species))
  P_record            <- array(0, dim = c(timesteps, patches, species))

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

    # 1. Emergence. The pool is the propagules that arrived at the end of the
    #    last step plus those that survived in the bank. A fraction becomes
    #    active; the rest stay dormant.
    P <- R + D_star
    A <- emergence(P, f_cue, cue_setup)
    D <- P - A

    # 2. Reproduction, limited by the environment and then by competition
    #    among the active individuals of every species in the patch
    B_expected <- expected_births(A, f, alpha)

    # 3. Dormant propagules survive in the bank
    D_star <- bank_survival(D, survival)

    # 4. Demographic stochasticity in reproduction
    B <- matrix(rpois(species * patches, lambda = pmax(B_expected, 0)),
                nrow = patches, ncol = species)

    # 5. Dispersal of the new propagules. Emigrants (E) leave each patch, and
    #    all of a species' emigrants are distributed among patches (I)
    #    according to its kernel. Dormant propagules do not disperse.
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
      loglam_local_record[t_rec, , ] <- safe_logratio(R + D_star, P)
      P_record[t_rec, , ]            <- P
      cohen_record[[t_rec]] <- cbind(species = 1:species,
                                     cohen_partition(P, A, D_star, R))
    }
  }

  # --- community summaries: diversity, variability, synchrony --------------------
  communities <- list(active   = active_record,
                      dormant  = dormant_record,
                      arriving = arriving_record,
                      pool     = pool_record)

  # Realized emergence of the run: all active individuals over the whole pool,
  # summed over species, patches, and recorded steps
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
  run_summary <- bind_rows(run_summary)

  # --- per-species summaries: growth of the propagule pool over time --------------
  # Cohen's mixture is computed over all patches together, so it describes the
  # metapopulation. The mean of log pool growth is the log geometric mean
  # growth rate, which is what bet-hedging acts on; its variance is what a seed
  # bank reduces.
  run_species <- bind_rows(cohen_record) |>
    group_by(species) |>
    summarize(
      g_realized       = mean(g_realized, na.rm = TRUE),
      yield_per_active = mean(yield_per_active, na.rm = TRUE),
      s_realized       = mean(s_realized, na.rm = TRUE),
      loglam_mean      = mean(log_lambda_pool, na.rm = TRUE),
      loglam_var       = var(log_lambda_pool, na.rm = TRUE),
      lambda_arith     = mean(lambda_pool, na.rm = TRUE)
    )

  # The same variance of log pool growth within each patch, then averaged over
  # the patches that count, weighted by the patch's mean pool, so that it is
  # the local variance met by an average propagule. Local growth includes
  # arrivals from other patches. The coverage columns say how much of the
  # species the average describes.
  run_species$loglam_var_local            <- NA_real_   # pool-weighted
  run_species$loglam_var_local_unweighted <- NA_real_
  run_species$local_patches_counted       <- NA_real_
  run_species$local_pool_share_counted    <- NA_real_   # share of the pool in those patches
  run_species$local_steps_defined         <- NA_real_   # share of patch-steps with growth defined
  for (i in 1:species) {
    var_in_patch <- rep(NA_real_, patches)
    steps_defined <- rep(0, patches)
    for (x in 1:patches) {
      loglam_x <- loglam_local_record[, x, i]
      steps_defined[x] <- sum(is.finite(loglam_x))
      if (steps_defined[x] >= local_growth_min_share * timesteps) {
        var_in_patch[x] <- var(loglam_x, na.rm = TRUE)
      }
    }
    mean_pool <- colMeans(P_record[, , i])
    counted <- is.finite(var_in_patch)
    run_species$loglam_var_local[i] <-
      sum(mean_pool[counted] * var_in_patch[counted]) / sum(mean_pool[counted])
    run_species$loglam_var_local_unweighted[i] <- mean(var_in_patch[counted])
    run_species$local_patches_counted[i]       <- sum(counted)
    run_species$local_pool_share_counted[i]    <- sum(mean_pool[counted]) / sum(mean_pool)
    run_species$local_steps_defined[i]         <- sum(steps_defined) / (timesteps * patches)
  }

  # Mean abundance, and the share of patch-steps occupied, in each community
  for (community in names(communities)) {
    N_record <- communities[[community]]
    run_species[[paste0("abundance_mean_", community)]] <- apply(N_record, 1, mean)
    run_species[[paste0("occupancy_", community)]]      <- apply(N_record > 0, 1, mean)
  }
  run_species$run <- run

  list(summary = run_summary, per_species = run_species)
}

stopImplicitCluster()

# ---------------------------------------------------------------------------
# 6. Collect and save
# ---------------------------------------------------------------------------

# A run that failed returns its error and not a list of results. Keep its row
# in the output, flagged, so a missing result is never mistaken for a run that
# was not in the design.
ran <- !sapply(results, inherits, what = "error")
message("Completed ", sum(ran), " of ", nrow(design), " simulations")
if (any(!ran)) {
  print(results[!ran][[1]])
}

summaries   <- bind_rows(lapply(results[ran], function(x) x$summary))
per_species <- bind_rows(lapply(results[ran], function(x) x$per_species))

# One row per run and community, so a failed run is a flagged row in each
design$ran  <- ran
summaries   <- expand_grid(design, community = c("active", "dormant", "arriving", "pool")) |>
  left_join(summaries, by = c("run", "community"))
per_species <- left_join(design, per_species, by = "run")

outdir <- here("sim_output", paste0("regimes_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(outdir, recursive = TRUE)

write_csv(summaries,   file.path(outdir, "summary.csv"))
write_csv(per_species, file.path(outdir, "per_species.csv"))
write_csv(env_diag,    file.path(outdir, "env_diagnostics.csv"))

# The setting of each replicate: where its patches are and their long-term
# means, and its competition coefficients
patch_table <- list()
competition_table <- list()
for (rep_i in 1:nreps) {
  patch_table[[rep_i]] <- tibble(rep = rep_i, patch = 1:patches,
                                 x = landscapes[[rep_i]]$x,
                                 y = landscapes[[rep_i]]$y,
                                 patch_mean = patch_means_list[[rep_i]])
  competition_table[[rep_i]] <- tibble(rep = rep_i,
                                       species_acting   = rep(1:species, times = species),
                                       species_affected = rep(1:species, each = species),
                                       coefficient      = as.vector(alphas[[rep_i]]))
}
write_csv(bind_rows(patch_table),       file.path(outdir, "patches.csv"))
write_csv(bind_rows(competition_table), file.path(outdir, "competition.csv"))

settings <- tibble(test_run, patches, species, niche_breadth,
                   optima_lower = optima_span[1], optima_upper = optima_span[2],
                   initialization, burn_in, timesteps, nreps, cycle_length,
                   kernel_decay, exclude_self, survival, replicate_seed, f_max,
                   local_growth_min_share,
                   noise_colors = paste(noise_colors, collapse = " "),
                   patch_mean_lower = patch_mean_span[1],
                   patch_mean_upper = patch_mean_span[2], patch_mean_range)
write_csv(settings, file.path(outdir, "run_settings.csv"))

message("Wrote output to ", outdir)
