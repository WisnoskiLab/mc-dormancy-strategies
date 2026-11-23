# Environments that vary in how predictable they are in time and how
# synchronous they are among patches.
#
# The condition in patch i at time t is the sum of three layers:
#
#     env_i(t) = mu_i  +  amplitude * sin(w * t + phi_i)  +  sigma * z_i(t)
#
#   mu_i              a static spatial gradient
#   amplitude, phi_i  a periodic, predictable layer. Giving patches different
#                     phases (phi_i) makes them asynchronous without making
#                     them unpredictable.
#   z_i(t)            a noise layer: autocorrelated in time (AR(1), coefficient
#                     `temp_noise_color`) and correlated among patches over a
#                     distance `spat_range`
#
# Predictability is the share of temporal variance in the periodic layer.
# Synchrony is set by the spread of phases and by `spat_range`. The two are
# set independently, which gives the four regimes in env_scenarios().
#
# env_generate()   build one environment
# env_diagnostics()   measure what an environment actually did
# env_scenarios()     the settings for the four regimes
#
# The manuscript's Methods describe the regimes and their settings.

# ---------------------------------------------------------------------------
# Spatial correlation
# ---------------------------------------------------------------------------

# Correlation among patches that declines exponentially with distance.
# spat_range = Inf makes every patch perfectly correlated.
patch_cor_matrix <- function(landscape, spat_range, torus = FALSE) {
  if (torus) {
    d <- as.matrix(dist.torus(landscape))
  } else {
    d <- as.matrix(dist(landscape))
  }

  if (is.infinite(spat_range)) {
    matrix(1, nrow = nrow(d), ncol = ncol(d))
  } else {
    exp(-d / spat_range)
  }
}

# Square root of a correlation matrix. Multiplying independent normal draws by
# it gives draws with that correlation among patches.
cor_sqrt <- function(C) {
  ev <- eigen(C, symmetric = TRUE)
  vals <- pmax(ev$values, 0)
  ev$vectors %*% diag(sqrt(vals), nrow = length(vals)) %*% t(ev$vectors)
}

# Long-term mean conditions for each patch: evenly spaced values between
# `lower` and `upper`, arranged in space as a random field that is correlated
# over a distance `spat_range`, so that nearby patches are alike.
patch_mean_field <- function(landscape, lower = 0, upper = 1, spat_range = 10,
                             torus = TRUE) {
  n_patches <- nrow(landscape)
  Cs <- cor_sqrt(patch_cor_matrix(landscape, spat_range, torus))
  field <- as.vector(Cs %*% rnorm(n_patches))
  lower + (upper - lower) * (rank(field) - 1) / (n_patches - 1)
}

# ---------------------------------------------------------------------------
# Generate one environment
# ---------------------------------------------------------------------------

#' @param landscape          data.frame with x, y (from init_landscape)
#' @param timesteps          length of the series
#' @param spat_heterogeneity slope of the static spatial gradient
#' @param patch_means        optional: the long-term mean of each patch, on the
#'                           scaled axis (0 to 1 is the range of niche optima).
#'                           When given, it replaces the gradient, and the
#'                           fluctuations are added on top of it without
#'                           clamping, so a patch at the end of the axis
#'                           fluctuates as far as one in the middle.
#' @param amplitude          amplitude of the periodic layer
#' @param w                  angular frequency (2 * pi / cycle length)
#' @param phi                baseline phase
#' @param phase_mode         "uniform": all patches in phase.
#'                           "wave": phase increases across the landscape.
#'                           "random": phases are a random field, correlated
#'                           among patches over `spat_range`.
#' @param phase_spread       0 to 1; the fraction of a full cycle the phases span
#' @param temp_noise_sd      standard deviation of the noise layer
#' @param temp_noise_color   AR(1) coefficient of the noise layer: negative is
#'                           blue noise, 0 white, positive red
#' @param spat_range         distance over which the noise layer (and "random"
#'                           phases) are correlated among patches
#' @param torus              wrap distances around the landscape edges
#' @param range_sd           used only when `scale_half` is not given
#' @param scale_half         half-width of the range mapped onto [0, 1]. Use
#'                           the same value for every regime being compared
#'                           (env_scenarios() supplies it); otherwise regimes
#'                           with more noise are compressed more.
#' @param seed               optional seed
#'
#' @return long data.frame with one row per time and patch: `env`, the
#'   realized conditions, and `env_signal`, the gradient plus the periodic
#'   layer only. Both are scaled to [0, 1] with the same map. `env_signal` is
#'   what a cued species can read.
env_generate <- function(landscape,
                            timesteps = 1000,
                            spat_heterogeneity = 0.1,
                            patch_means = NULL,
                            amplitude = 5,
                            w = (2 * pi) / 10,
                            phi = 0,
                            phase_mode = c("uniform", "wave", "random"),
                            phase_spread = 0,
                            temp_noise_sd = 0,
                            temp_noise_color = 0,
                            spat_range = Inf,
                            torus = FALSE,
                            range_sd = 3,
                            scale_half = NULL,
                            seed = NULL) {

  phase_mode <- match.arg(phase_mode)
  if (!is.null(seed)) set.seed(seed)

  n_patches <- nrow(landscape)

  # --- static spatial gradient ------------------------------------------------
  mu <- spat_heterogeneity * landscape$x + spat_heterogeneity * landscape$y
  if (!is.null(patch_means)) {
    # patch means are added after scaling, below
    mu <- rep(0, n_patches)
  }

  # --- phase of the periodic layer in each patch ------------------------------
  if (phase_mode == "uniform") {
    phi_i <- rep(phi, n_patches)
  }
  if (phase_mode == "wave") {
    position <- landscape$x + landscape$y
    position_range <- diff(range(position))
    if (position_range == 0) {
      phi_i <- rep(phi, n_patches)
    } else {
      phi_i <- phi + phase_spread * 2 * pi * (position - min(position)) / position_range
    }
  }
  if (phase_mode == "random") {
    # a random field correlated over the same distance as the noise layer, so
    # asynchrony has one spatial scale
    Cs <- cor_sqrt(patch_cor_matrix(landscape, spat_range, torus))
    field <- as.vector(Cs %*% rnorm(n_patches))
    phi_i <- phi + phase_spread * 2 * pi * (rank(field) - 1) / max(1, n_patches - 1)
  }

  # --- noise layer: AR(1) in time, correlated among patches -------------------
  z <- matrix(0, nrow = timesteps, ncol = n_patches)
  if (temp_noise_sd > 0) {
    Cs <- cor_sqrt(patch_cor_matrix(landscape, spat_range, torus))
    a <- temp_noise_color
    b <- sqrt(1 - a^2)   # keeps the variance of z at 1 for any a
    z[1, ] <- as.vector(Cs %*% rnorm(n_patches))
    if (timesteps > 1) {
      for (t in 2:timesteps) {
        z[t, ] <- a * z[t - 1, ] + b * as.vector(Cs %*% rnorm(n_patches))
      }
    }
  }

  # --- add the layers ----------------------------------------------------------
  times <- 1:timesteps
  gradient <- matrix(mu, nrow = timesteps, ncol = n_patches, byrow = TRUE)
  periodic <- outer(times, phi_i, function(t, p) amplitude * sin(w * t + p))

  env_mat    <- gradient + periodic + temp_noise_sd * z
  signal_mat <- gradient + periodic

  # --- scale to [0, 1] ----------------------------------------------------------
  # One linear map, centered on the gradient, for both matrices. Values beyond
  # the range are clamped, and the clamped fraction is reported.
  center <- mean(range(mu))
  if (!is.null(scale_half)) {
    half <- scale_half
  } else {
    half <- diff(range(mu)) / 2 + amplitude + range_sd * temp_noise_sd
  }
  if (half <= 0) half <- 1

  env_scaled <- (env_mat - (center - half)) / (2 * half)
  frac_clamped <- mean(env_scaled < 0 | env_scaled > 1)
  env_scaled[env_scaled < 0] <- 0
  env_scaled[env_scaled > 1] <- 1

  sig_scaled <- (signal_mat - (center - half)) / (2 * half)
  sig_scaled[sig_scaled < 0] <- 0
  sig_scaled[sig_scaled > 1] <- 1

  # --- or: patch means, with the fluctuations on top ------------------------------
  # The fluctuations are scaled by the same factor as above, so their size on
  # the scaled axis does not depend on how the patch means are spread.
  if (!is.null(patch_means)) {
    mean_mat <- matrix(patch_means, nrow = timesteps, ncol = n_patches, byrow = TRUE)
    env_scaled <- mean_mat + env_mat / (2 * half)
    sig_scaled <- mean_mat + signal_mat / (2 * half)
    frac_clamped <- 0
  }

  # --- long format: one row per time and patch ---------------------------------
  colnames(env_scaled) <- 1:n_patches
  colnames(sig_scaled) <- 1:n_patches
  env_df <- tidyr::pivot_longer(
    cbind.data.frame(time = times, env_scaled),
    cols = (1:n_patches) + 1, names_to = "patch", values_to = "env"
  )
  sig_df <- tidyr::pivot_longer(
    cbind.data.frame(time = times, sig_scaled),
    cols = (1:n_patches) + 1, names_to = "patch", values_to = "env_signal"
  )
  env_df$env_signal <- sig_df$env_signal

  attr(env_df, "regime") <- list(
    spat_heterogeneity = spat_heterogeneity, amplitude = amplitude, w = w,
    phase_mode = phase_mode, phase_spread = phase_spread,
    temp_noise_sd = temp_noise_sd, temp_noise_color = temp_noise_color,
    spat_range = spat_range, torus = torus,
    frac_clamped = frac_clamped
  )
  env_df
}

# ---------------------------------------------------------------------------
# Measure what an environment actually did
# ---------------------------------------------------------------------------

#' Realized properties of a generated environment.
#'
#' Report these, and not the settings that produced the environment.
#'
#' @param env_df    output of env_generate
#' @param w         angular frequency of the periodic layer
#' @param value_col "env" (realized conditions) or "env_signal" (the cue)
#' @return one row:
#'   mean_patch_cor    mean correlation between patches; 1 = synchronous
#'   synchrony_ldm     Loreau & de Mazancourt synchrony among patches
#'   predictability    share of temporal variance explained by the cycle
#'   rank_persistence  correlation of the patch ranking from one step to the
#'                     next; 1 = the ranking never changes
#'   top_switch_rate   how often the best patch changes
#'   acf_lag1, acf_lag2  temporal autocorrelation at lags 1 and 2
#'   decorr_time       lag at which autocorrelation falls below 1/e
#'   temporal_cv_mean, spatial_cv_mean, env_range
env_diagnostics <- function(env_df, w = (2 * pi) / 10, value_col = "env") {

  # time x patch matrix of the chosen column
  d <- as.data.frame(env_df)[, c("time", "patch", value_col)]
  names(d)[3] <- ".value"
  wide <- tidyr::pivot_wider(d, names_from = "patch", values_from = ".value")
  times <- wide$time
  m <- as.matrix(wide[, -1, drop = FALSE])

  patch_sd <- apply(m, 2, sd)
  varies <- patch_sd > .Machine$double.eps

  # --- synchrony among patches --------------------------------------------------
  mean_cor <- NA_real_
  if (sum(varies) > 1) {
    cors <- suppressWarnings(cor(m[, varies, drop = FALSE]))
    mean_cor <- mean(cors[upper.tri(cors)])
  }

  synchrony_ldm <- NA_real_
  if (sum(patch_sd) > 0) {
    synchrony_ldm <- sd(rowSums(m))^2 / sum(patch_sd)^2
  }

  # --- predictability: R^2 of a regression on the known cycle -------------------
  X <- cbind(1, sin(w * times), cos(w * times))
  r2 <- rep(NA_real_, ncol(m))
  for (j in 1:ncol(m)) {
    y <- m[, j]
    if (sd(y) < .Machine$double.eps) next
    fit <- .lm.fit(X, y)
    r2[j] <- 1 - sum(fit$residuals^2) / sum((y - mean(y))^2)
  }

  # --- rank turnover: does the ordering of patches change? ----------------------
  # This decides whether dispersal can spread risk. If the best patch is always
  # the same patch, moving is relocation and not insurance.
  ranks <- t(apply(m, 1, rank))
  lag_cor <- rep(NA_real_, nrow(ranks) - 1)
  for (t in 1:(nrow(ranks) - 1)) {
    now  <- ranks[t, ]
    nxt  <- ranks[t + 1, ]
    if (sd(now) < .Machine$double.eps || sd(nxt) < .Machine$double.eps) next
    lag_cor[t] <- cor(now, nxt)
  }
  rank_persistence <- mean(lag_cor, na.rm = TRUE)

  best_patch <- apply(m, 1, which.max)
  top_switch_rate <- NA_real_
  if (length(best_patch) > 1) {
    top_switch_rate <- mean(best_patch[-1] != best_patch[-length(best_patch)])
  }

  # --- timescale: temporal autocorrelation ---------------------------------------
  # Predictability and timescale are different properties. Populations average
  # over fast fluctuations, so they have to be reported separately.
  max_lag <- max(2, min(50, floor(nrow(m) / 4)))
  acf_by_patch <- matrix(NA_real_, nrow = max_lag, ncol = ncol(m))
  for (j in 1:ncol(m)) {
    y <- m[, j]
    if (sd(y) < .Machine$double.eps) next
    acf_by_patch[, j] <- acf(y, lag.max = max_lag, plot = FALSE, demean = TRUE)$acf[-1]
  }
  acf_mean <- rowMeans(acf_by_patch, na.rm = TRUE)
  below <- which(acf_mean < exp(-1))
  decorr_time <- NA_real_
  if (length(below) > 0) decorr_time <- below[1]

  # --- variability in time and space ----------------------------------------------
  cv <- function(v) {
    if (abs(mean(v)) < .Machine$double.eps) return(NA_real_)
    sd(v) / mean(v)
  }

  data.frame(
    mean_patch_cor   = mean_cor,
    synchrony_ldm    = synchrony_ldm,
    predictability   = mean(r2, na.rm = TRUE),
    rank_persistence = rank_persistence,
    top_switch_rate  = top_switch_rate,
    acf_lag1         = acf_mean[1],
    acf_lag2         = acf_mean[2],
    decorr_time      = decorr_time,
    temporal_cv_mean = mean(apply(m, 2, cv), na.rm = TRUE),
    spatial_cv_mean  = mean(apply(m, 1, cv), na.rm = TRUE),
    env_range        = diff(range(m))
  )
}

# ---------------------------------------------------------------------------
# The four regimes
# ---------------------------------------------------------------------------

#' Settings for the four environmental regimes.
#'
#' Predictable or unpredictable in time, crossed with synchronous or
#' asynchronous among patches.
#'
#' Total temporal variance is the same in every regime. Predictability is the
#' fraction of it in the periodic layer (`pred_frac`):
#'
#'     amplitude = sqrt(2 * pred_frac * total_var)   amplitude of the cycle
#'     sd = sqrt((1 - pred_frac) * total_var)    sd of the noise
#'
#' Synchronous regimes correlate the noise over a distance far longer than the
#' landscape (`sync_range`); asynchronous regimes over a short one
#' (`async_range`). In the predictable and asynchronous regime, patches differ
#' in the phase of the cycle, with phases correlated over that same short
#' distance.
#'
#' Every regime keeps a shallow spatial gradient, so that there is always
#' somewhere better to disperse to, but shallow enough that the fluctuations
#' can reorder the patches.
#'
#' @param extent      width of the landscape
#' @param total_var   temporal variance, the same in every regime
#' @param gradient    slope of the spatial gradient
#' @param pred_hi,pred_lo predictable fraction in the predictable and
#'   unpredictable regimes
#' @param colors      AR(1) coefficients of the noise. More than one value
#'   repeats the four regimes for each color.
#' @return data.frame with one row per regime, holding the arguments for
#'   env_generate(), including the shared `scale_half`
env_scenarios <- function(extent = 100, total_var = 12.5,
                          gradient = 0.02,
                          pred_hi = 0.9, pred_lo = 0.1,
                          colors = 0, range_sd = 3,
                          sync_range = 50 * extent, async_range = 0.1 * extent) {

  four_regimes <- data.frame(
    scenario     = c("predictable_synchronous", "predictable_asynchronous",
                     "unpredictable_synchronous", "unpredictable_asynchronous"),
    pred_frac    = c(pred_hi, pred_hi, pred_lo, pred_lo),
    phase_mode   = c("uniform", "random", "uniform", "uniform"),
    phase_spread = c(0, 1, 0, 0),
    spat_range   = c(sync_range, async_range, sync_range, async_range),
    stringsAsFactors = FALSE
  )

  # the four regimes at each noise color
  regimes <- data.frame()
  for (color in colors) {
    this_color <- four_regimes
    this_color$temp_noise_color <- color
    if (length(colors) > 1) {
      this_color$scenario <- paste0(this_color$scenario, "_color", color)
    }
    regimes <- rbind(regimes, this_color)
  }

  regimes$spat_heterogeneity <- gradient
  regimes$amplitude     <- sqrt(2 * regimes$pred_frac * total_var)
  regimes$temp_noise_sd <- sqrt((1 - regimes$pred_frac) * total_var)

  out <- regimes

  # one scaling for all regimes: wide enough for the most variable of them
  gradient_half_range <- diff(range(c(0, gradient * 2 * extent))) / 2
  out$scale_half <- max(gradient_half_range + out$amplitude + range_sd * out$temp_noise_sd)
  rownames(out) <- NULL
  out
}
