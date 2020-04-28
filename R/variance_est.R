# Variance estimators for systematic samples with unequal probability sampling
library(spsurvey)
library(ggplot2)

#' Estimates the sampling variance of the mean under SRSWoR
var_srs <- function(samp, N) {
  n <- nrow(samp)
  return ((n - 1)^(-1) * n^(-1) * (1 - (n / N)) * sum((samp$z - mean(samp$z))^2))
}

#' Estimates the sampling variance of the mean under STSRSWoR
var_strs <- function(samp, N) {
  strs_var <- 0
  strata_ids <- unique(samp$strat_id)
  
  for(stratum_id in strata_ids) {
    samp_i <- samp[samp$strat_id == stratum_id,]
    n_i <- samp_i[[1, 'n_i']]
    N_i <- samp_i[[1, 'N_i']]
    var_mean <- var_srs(samp_i, N_i)
    strs_var <- strs_var + (N_i / N)^2 * var_mean
  }
  
  return(strs_var)
  
}

#' Estimates the sampling varaince of the mean using 
#' the Stevens and Olsen (2003) estimator.
var_so <- function(samp, N) {
  wt <- localmean.weight(samp$x, samp$y, samp$pi_i)
  var_total <- localmean.var(samp$z, wt)
  return((1/N) * var_total) # TODO Why is it 1/N and not 1/N^2?
}

#' Estimates the sampling variance of the mean using 
#' a Monte Carlo simulation. Generally treated as the 
#' "true" sampling variance.
var_sys_monte_carlo <- function(pop_mat, pop_df, strat_spec, M) {
  N <- nrow(pop_df)
  mu <- sum(pop_mat) / N
  
  var_hat <- 0
  vars <- c()
  
  for(i in 1:M) {
    print(i)
    samp <- samp_tri_grid_pi(pop_mat, pop_df, strat_spec)
    mu_hat <- (1/N) * sum(samp$z / samp$pi_i)
    vars <- c(vars, ((mu - mu_hat)^2))
  }
  
  var_sys <- mean(vars)
  return(var_sys)
}
