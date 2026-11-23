# partition alpha, beta, gamma div
diversity.partition <- function(metacomm_tsdata){
  ts_dims <- dim(metacomm_tsdata)
  names(ts_dims) <- c("species", "time", "sites")
  
  div_dynamics <- matrix(ncol = 3, nrow = ts_dims["time"])
  
  for(time in 1:ts_dims["time"]){
    snapshot <- metacomm_tsdata[,time,]
    
    occupancy <- rowSums((snapshot > 0) * 1)
    gamma_div <- sum(occupancy > 0)
    
    alpha_div <- mean(colSums(snapshot > 0))
    
    beta_div <- gamma_div / alpha_div
    
    div_dynamics[time,] <- c(alpha_div, beta_div, gamma_div)
  }
  
  return(colMeans(div_dynamics))
}



# compute spatial and temporal beta
#
# Beta diversity is the total variance of the Hellinger-transformed community
# table (Legendre & De Caceres 2013). The transformation divides by a
# community's total abundance, so it is undefined for an empty community.
# Empty patches and empty time steps are therefore left out: spatial beta is
# measured among the patches occupied at a time step, and temporal beta over
# the time steps in which a patch was occupied. At least two communities are
# needed to measure either.
beta.div.calc <- function(metacomm_tsdata){
  ts_dims <- dim(metacomm_tsdata)
  names(ts_dims) <- c("species", "time", "sites")

  # how does spatial beta diversity change over time
  spat_beta_dynamics <- rep(NA_real_, ts_dims["time"])

  for(time in 1:ts_dims["time"]){
    snapshot <- metacomm_tsdata[,time,]
    sbs <- t(snapshot)
    sbs <- sbs[rowSums(sbs) > 0, , drop = FALSE]   # occupied patches only
    if(nrow(sbs) < 2) next
    sbs <- vegan::decostand(sbs, "hellinger")
    tot_ss <- sum((scale(sbs, center = TRUE, scale = FALSE))^2)
    spat_beta_dynamics[time] <- tot_ss / (nrow(sbs)-1)
  }

  # how does temporal beta diversity vary across space
  temp_beta_distribution <- rep(NA_real_, ts_dims["sites"])
  for(site in 1:ts_dims["sites"]){
    site_ts <- metacomm_tsdata[,,site]
    tbs <- t(site_ts)
    tbs <- tbs[rowSums(tbs) > 0, , drop = FALSE]   # occupied time steps only
    if(nrow(tbs) < 2) next
    tbs <- vegan::decostand(tbs, "hellinger")
    tot_ss <- sum((scale(tbs, center = TRUE, scale = FALSE))^2)
    temp_beta_distribution[site] <- tot_ss / (nrow(tbs)-1)
  }

  beta_div <- data.frame(mean_beta_spatial = mean(spat_beta_dynamics, na.rm = TRUE),
                         mean_beta_temporal = mean(temp_beta_distribution, na.rm = TRUE))
  return(beta_div)
}
