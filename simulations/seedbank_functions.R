# Survival in the seed bank, and growth of the propagule pool.
#
#   bank_survival()    survival of propagules that stay dormant
#   cohen_partition()  growth of the propagule pool as a mixture of two fates
#                      (Cohen 1966), per species at the metacommunity scale
#                      (summed over patches), which avoids dividing by zero
#                      whenever a patch is empty
#   bank_residence()   expected time a propagule spends in the bank

# Divide, returning NA where the denominator is zero or the result is not finite.
safe_div <- function(num, den) {
  out <- num / den
  out[den == 0 | !is.finite(out)] <- NA_real_
  out
}

# Log of a ratio, returning NA where either side is zero.
safe_logratio <- function(num, den) {
  out <- log(num / den)
  out[!is.finite(out)] <- NA_real_
  out
}

#' Survival in the bank: each dormant propagule survives with probability s,
#' D_star ~ Binomial(D, s).
#'
#' @param D  patches x species dormant propagules
#' @param s  survival, per species or shared
#' @return patches x species surviving propagules
bank_survival <- function(D, s) {
  s <- rep_len(s, ncol(D))
  D_star <- D * 0
  for (i in 1:ncol(D)) {
    D_star[, i] <- rbinom(n = nrow(D), size = D[, i], prob = s[i])
  }
  D_star
}

#' Growth of the propagule pool as a mixture of emerging and staying dormant.
#'
#'     lambda_pool = g * y + (1 - g) * s
#'
#' g is the realized emergence fraction, y the propagules each active
#' individual contributes to the next step's pool (after competition,
#' demographic stochasticity, and dispersal), and s the realized survival of
#' propagules that stayed dormant. The identity is exact.
#'
#' Dormancy sets the weights of the mixture. Its benefit is a lower temporal
#' variance of log growth, so summarize `log_lambda_pool` over time: its mean
#' is the log geometric mean growth rate, and its variance is what a seed bank
#' reduces.
#'
#' @param P       patches x species propagule pool at the start of the step
#' @param A       patches x species active individuals
#' @param D_star  patches x species dormant propagules that survived
#' @param R       patches x species arriving propagules at the end of the step
#' @return one row per species
cohen_partition <- function(P, A, D_star, R) {
  pool      <- colSums(P, na.rm = TRUE)
  active    <- colSums(A, na.rm = TRUE)
  survivors <- colSums(D_star, na.rm = TRUE)
  arriving  <- colSums(R, na.rm = TRUE)

  data.frame(
    g_realized       = safe_div(active, pool),
    yield_per_active = safe_div(arriving, active),
    s_realized       = safe_div(survivors, pool - active),
    lambda_pool      = safe_div(arriving + survivors, pool),
    log_lambda_pool  = safe_logratio(arriving + survivors, pool)
  )
}

#' Expected time a propagule spends in the bank, in time steps.
#'
#' A propagule stays dormant and survives with probability (1 - g_mean) * s
#' each step, so its expected residence is 1 / (1 - that probability). With
#' g_mean = 1 there is no bank and the residence time is 1.
bank_residence <- function(g_mean, s) {
  retention <- (1 - g_mean) * s
  ifelse(retention >= 1, Inf, 1 / (1 - retention))
}
