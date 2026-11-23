# Landscape, species, competition, and reproduction.
#
# The other model functions are in env_functions.R (environments and regimes),
# emergence_functions.R (leaving dormancy), dispersal_functions.R (kernels),
# and seedbank_functions.R (survival in the bank and Cohen's mixture).
#
# Notation follows the manuscript: A active individuals, f_i(theta) the
# propagules an active individual of species i produces under conditions
# theta without competition, alpha the competition coefficients.

# Distances between patches on a torus, so the landscape has no edges
# (dist.torus from the package som.nn)
dist.torus <- function(coors) {
  x <- coors[, 1]
  dx <- stats::dist(x, diag = TRUE, upper = TRUE)
  y <- coors[, 2]
  dy <- stats::dist(y, diag = TRUE, upper = TRUE)
  max.x <- max(dx) + 1
  max.y <- max(dy) + 1
  mdx <- max.x - dx
  mdy <- max.y - dy
  dx <- pmin(dx, mdx)
  dy <- pmin(dy, mdy)
  sqrt(dx^2 + dy^2)
}

# Patches at distinct integer coordinates, drawn uniformly on the grid
init_landscape <- function(patches, x_dim = 100, y_dim = 100) {
  repeat {
    landscape <- data.frame(x = sample(1:x_dim, size = patches, replace = TRUE),
                            y = sample(1:y_dim, size = patches, replace = TRUE))
    if (nrow(unique(landscape)) == patches) {
      break
    }
  }
  landscape
}

# One row per species with its traits. Every argument but `niche_optima` is
# one value shared by all species, or one value per species.
init_species <- function(niche_optima, niche_breadth, f_max, dispersal_rate,
                         kernel_decay, g_mean, survival) {
  data.frame(species        = seq_along(niche_optima),
             niche_optimum  = niche_optima,
             niche_breadth  = niche_breadth,
             f_max          = f_max,
             dispersal_rate = dispersal_rate,
             kernel_decay   = kernel_decay,
             g_mean         = g_mean,
             survival       = survival)
}

# The starting propagules in each patch (patches x species). The drivers add
# one, so every patch starts with at least one propagule of every species.
init_community <- function(species, patches) {
  matrix(rpois(n = species * patches, lambda = 0.5), nrow = patches, ncol = species)
}

# Competition coefficients (species x species): `intra` on the diagonal,
# interspecific coefficients drawn uniformly between `min_inter` and
# `max_inter`, all multiplied by `scale`. Element [j, i] is the per capita
# effect of species j on species i (alpha_ij in the manuscript), so that
# A %*% alpha sums the competition each species meets in each patch.
competition_matrix <- function(species, intra = 1, min_inter = 0, max_inter = 1,
                               scale = 0.05) {
  alpha <- matrix(runif(n = species * species, min = min_inter, max = max_inter),
                  nrow = species, ncol = species)
  diag(alpha) <- intra
  alpha * scale
}

# Propagules per active individual without competition, f_i(theta_x), for every
# patch and species (patches x species): f_max reduced by a Gaussian function
# of the mismatch between the conditions and the species' niche optimum.
#
# env  conditions in each patch, one value per patch
fecundity <- function(species_traits, env) {
  n_patches <- length(env)
  n_species <- nrow(species_traits)
  env_mat <- matrix(rep(env, each = n_species), nrow = n_species, ncol = n_patches)
  suitability <- exp(-((species_traits$niche_optimum - env_mat) /
                         (2 * species_traits$niche_breadth))^2)
  species_traits$f_max * t(suitability)
}

# Expected propagules produced (patches x species): each species' active
# individuals times their fecundity, reduced by Beverton-Holt competition among
# the active individuals of every species in the patch,
#   A_ix f_i / (1 + sum_j alpha_ij A_jx)
# A and f are patches x species; alpha is from competition_matrix().
expected_births <- function(A, f, alpha) {
  (A * f) / (1 + A %*% alpha)
}
