# Variance estimators for systematic samples with unequal probability sampling
library(spsurvey)
library(ggplot2)
library(FNN)
library(tidyr)

##' Estimates the sampling variance of the mean using
##' finite population block kriging
#var_fpbk <- function(samp, N) {
#  centroids <- gCentroid(samp, byid=TRUE)@coords
#  slm_df <- cbind(samp@data, centroids)
#  
#  # TODO sptotal does a lot of work computing the distance matrix each time,
#  # it could be possible to compute a distance matrix in the outermost simulation
#  # procedure and then simply index it based on what subsample is selected.
#  # sptotal is not well modularized so it will take some "ripping out" of functions
#  # to do this, but it will save a lot of simulation time
#  
#  # maybe move on to other estimators before getting too involved here, but I like 
#  # the idea of using FPBK as a very basic model-based variance estimator
#  slm <- slmfit(VOL ~ 1, data=slm_df, xcoordcol='x', ycoordcol='y')
#  slm.pred <- predict(slm)
#  
#  var_total <- slm.pred$PredVar[1,1]
#  n <- sum(!is.na(slm_df$VOL))
#  fpc <- 1 - n/N
#  
#  return(var_total * fpc * 1/N^2)
#}
#
#
##' Calculate the variance using a denominator of n
##' instead of n-1
#pop_var <- function(z) {
#  n <- length(z)
#  ssq <- sum((z - mean(z))^2)
#  ssq/n
#}
#
##' Non-overlapping neighborhood variance
##' estimator for hexagonal nearest neighbors.
#var_nnbh_hex <- function(samp, a) {
# N <- nrow(samp)
# neighborhoods <- get_hex_neighborhoods(samp[,c('c', 'r')], a)
# neighborhoods <- merge(neighborhoods, samp, by.x=c('c_n', 'r_n'), by.y=c('c', 'r'), all.x=TRUE)
# 
# N_neighbs <- neighborhoods %>%
#   group_by(r, c) %>%
#   summarize(n=n())
# N_neighbs <- nrow(N_neighbs)
# 
# neighborhoods <- neighborhoods %>%
#   filter(!is.na(z_1)) %>%
#   group_by(r, c) %>%
#   summarize(pop_var = pop_var(z_1), q_j = n()) %>%
#   mutate(N_j = q_j * a^2) %>%
#   mutate(w_j_sq = (N_j / N)^2, fpc = ((N_j - q_j) / N_j)) %>%
#   mutate(nbh_var = (1/N_neighbs)^2 *  (pop_var / q_j) * fpc) # TODO N_neighbs needs to be fixed
# 
# sum(neighborhoods$nbh_var)
#}
#
##' The Dorazio (2003) variance estimator.
#var_dorazio <- function(samp, z, a) {
#  N <- length(z)
#  C <- gearys_C(samp[,c('c', 'r')], z, a)
#  var_srs(samp$z_1, N) * C
#}
#
##' Estimates the sampling variance of the mean using 
##' a Monte Carlo simulation. Generally treated as the 
##' "true" sampling variance.
#var_sys_monte_carlo <- function(pop_mat, pop_df, strat_spec, M) {
#  N <- nrow(pop_df)
#  mu <- sum(pop_mat) / N
#  
#  var_hat <- 0
#  vars <- c()
#  
#  for(i in 1:M) {
#    print(i)
#    samp <- samp_tri_grid_pi(pop_mat, pop_df, strat_spec)
#    mu_hat <- (1/N) * sum(samp$z / samp$pi_i)
#    vars <- c(vars, ((mu - mu_hat)^2))
#  }
#  
#  var_sys <- mean(vars)
#  return(var_sys)
#}
#