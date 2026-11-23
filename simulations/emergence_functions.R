# Emergence from dormancy: a fixed fraction (bet-hedging) or a response to a cue.
#
# Both are one model. Each step, an individual of species i in the propagule
# pool of patch x becomes active with probability g, where
#
#     logit(g) = b_i + beta * c_tilde
#
#   c_tilde  the species' cue, standardized over every patch and step of the
#            fluctuating environment
#   beta     responsiveness to the cue
#   b_i      the intercept, solved so that mean emergence over every patch and
#            step equals g_mean, whatever beta is
#
# With beta = 0, g is the fixed fraction of a bet-hedger. With beta > 0 the
# species emerges more when the cue is favorable. Because mean emergence is
# held at g_mean for every beta, a cued and a bet-hedging seed bank differ in
# when they emerge, not in how much on average over the landscape.
#
# The cue is the species' fecundity f_i under some reading of the conditions:
# the predictable part of the coming step's conditions (anticipatory) or the
# realized conditions of the previous step (responsive). The drivers build it.

#' The intercept b that makes mean emergence equal `g_mean`.
#'
#' Mean emergence rises steadily with the intercept, so there is one solution.
#' With beta = 0 it is simply logit(g_mean).
#'
#' @param c_tilde  the standardized cue over every patch and step
#' @param beta     responsiveness to the cue
#' @param g_mean   the mean emergence fraction to match
solve_emergence_intercept <- function(c_tilde, beta, g_mean, lim = 40) {
  if (isTRUE(all.equal(beta, 0))) {
    return(qlogis(g_mean))
  }

  # mean emergence at a trial intercept, minus the target
  mismatch <- function(b) {
    mean(plogis(b + beta * c_tilde)) - g_mean
  }

  # if the target cannot be reached within the search limits, use the
  # intercept without a cue
  if (mismatch(-lim) > 0 || mismatch(lim) < 0) {
    return(qlogis(g_mean))
  }
  uniroot(mismatch, c(-lim, lim), tol = 1e-9)$root
}

#' Standardize each species' cue and solve its intercept. Call once per run,
#' before the time loop.
#'
#' @param cue_all  list with one element per species: its cue over every patch
#'                 and step of the fluctuating environment
#' @param g_mean   mean emergence fraction, per species or shared
#' @param beta     responsiveness, per species or shared; 0 = bet-hedger
#' @return list with the per-species cue mean, cue sd, intercept, and beta
emergence_setup <- function(cue_all, g_mean, beta) {
  n_species <- length(cue_all)
  g_mean <- rep_len(g_mean, n_species)
  beta   <- rep_len(beta, n_species)

  cue_mean  <- numeric(n_species)
  cue_sd    <- numeric(n_species)
  intercept <- numeric(n_species)

  for (i in 1:n_species) {
    cue <- cue_all[[i]]
    cue_mean[i] <- mean(cue, na.rm = TRUE)
    cue_sd[i]   <- sd(cue, na.rm = TRUE)

    # a cue that never varies carries no information; avoid dividing by zero
    if (!is.finite(cue_sd[i]) || cue_sd[i] < .Machine$double.eps) {
      cue_sd[i] <- 1
    }

    c_tilde <- (cue - cue_mean[i]) / cue_sd[i]
    intercept[i] <- solve_emergence_intercept(c_tilde, beta[i], g_mean[i])
  }

  list(mean = cue_mean, sd = cue_sd, intercept = intercept, beta = beta)
}

#' Emergence for one time step: A ~ Binomial(P, g).
#'
#' @param P      patches x species propagule pool
#' @param cue    patches x species cue this step
#' @param setup  output of emergence_setup()
#' @return patches x species active individuals
emergence <- function(P, cue, setup) {
  A <- P * 0
  for (i in 1:ncol(P)) {
    c_tilde <- (cue[, i] - setup$mean[i]) / setup$sd[i]
    g <- plogis(setup$intercept[i] + setup$beta[i] * c_tilde)
    A[, i] <- rbinom(nrow(P), P[, i], g)
  }
  A
}
