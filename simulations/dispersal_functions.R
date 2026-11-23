# Dispersal kernels in which emigrants leave their home patch.
#
# An emigrant lands in another patch with a probability that declines
# exponentially with distance, exp(-kernel_decay * distance). The home patch
# gets zero weight, so `dispersal_rate` is the fraction of individuals that
# leaves the home patch for every species, and the kernel only decides where
# those individuals land.
#
# If the home patch kept its weight of exp(0) = 1, the largest in the kernel,
# a share of each species' emigrants would land back where it started, and
# that share would depend on the kernel. `exclude_self = FALSE` gives that
# kernel.

#' Dispersal kernels for every species, as a patches x patches x species array.
#'
#' Each slice has columns that sum to 1: disp_array[j, i, k] is the probability
#' that an emigrant of species k leaving patch i arrives in patch j, so
#' disp_array[, , k] %*% E gives expected arrivals from emigrants E.
#'
#' @param species_traits  data frame from init_species(); uses `kernel_decay`
#' @param exclude_self    TRUE = emigrants cannot land in their home patch
generate_dispersal_matrices <- function(landscape, species, patches,
                                        species_traits, torus = TRUE,
                                        exclude_self = TRUE) {

  if (torus) {
    dist_mat <- as.matrix(dist.torus(coors = landscape))
  } else {
    dist_mat <- as.matrix(dist(landscape))
  }

  disp_array <- array(dim = c(patches, patches, species))
  for (k in 1:species) {
    weights <- exp(-species_traits[k, "kernel_decay"] * dist_mat)
    if (exclude_self) diag(weights) <- 0
    # apply(., 1, ...) returns each normalized row as a column of the result
    disp_array[, , k] <- apply(weights, 1, function(x) x / sum(x))
  }
  disp_array
}

#' Mean distance an emigrant travels, averaged over home patches.
#'
#' The mean of the kernel over a continuous field of patches is
#' 2 / kernel_decay. On a finite set of patches the distance
#' actually traveled cannot be shorter than the spacing to the nearest patch or
#' longer than the landscape allows. Report this value alongside the kernel.
#'
#' @param kernel_decay  one or more rates at which the kernel decays with distance
#' @return realized mean distance, one value per element of `kernel_decay`
realized_dispersal_distance <- function(landscape, kernel_decay, torus = TRUE,
                                        exclude_self = TRUE) {

  if (torus) {
    dist_mat <- as.matrix(dist.torus(coors = landscape))
  } else {
    dist_mat <- as.matrix(dist(landscape))
  }

  distance_realized <- numeric(length(kernel_decay))
  for (k in seq_along(kernel_decay)) {
    weights <- exp(-kernel_decay[k] * dist_mat)
    if (exclude_self) diag(weights) <- 0
    distance_realized[k] <- mean(rowSums(weights * dist_mat) / rowSums(weights))
  }
  distance_realized
}
