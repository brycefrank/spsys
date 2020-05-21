# Variance estimators for systematic samples with unequal probability sampling
library(spsurvey)
library(ggplot2)
library(FNN)
library(tidyr)

#' Estimates the sampling variance of the mean under SRSWoR
var_srs <- function(samp, N) {
  n <- nrow(samp)
  fpc <- 1 - n/N
  var_mu <- (n - 1)^(-1) * n^(-1) *  sum((samp$z - mean(samp$z))^2)
  return(var_mu * fpc)
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

#' Estimates the sampling variance of the mean using 
#' the Stevens and Olsen (2003) estimator.
var_so <- function(samp, N) {
  wt <- localmean.weight(samp$x, samp$y, samp$pi_i)
  var_total <- localmean.var(samp$z/samp$pi_i, wt)
  n <- nrow(samp)
  fpc <- 1 - n/N
  return((1/N^2) * (1-(n/N))  * var_total) 
}

#' Estimates the sampling variance of the mean using
#' finite population block kriging
var_fpbk <- function(samp, N) {
  centroids <- gCentroid(samp, byid=TRUE)@coords
  slm_df <- cbind(samp@data, centroids)
  
  # TODO sptotal does a lot of work computing the distance matrix each time,
  # it could be possible to compute a distance matrix in the outermost simulation
  # procedure and then simply index it based on what subsample is selected.
  # sptotal is not well modularized so it will take some "ripping out" of functions
  # to do this, but it will save a lot of simulation time
  
  # maybe move on to other estimators before getting too involved here, but I like 
  # the idea of using FPBK as a very basic model-based variance estimator
  slm <- slmfit(VOL ~ 1, data=slm_df, xcoordcol='x', ycoordcol='y')
  slm.pred <- predict(slm)
  
  var_total <- slm.pred$PredVar[1,1]
  n <- sum(!is.na(slm_df$VOL))
  fpc <- 1 - n/N
  
  return(var_total * fpc * 1/N^2)
}

#' For a given set of (subsampled) hexagonal indices, translate
#' to the base coordinate system.
translate_hex_ix <- function(hex_ix, a) {
  r <- hex_ix$r
  c <- hex_ix$c
  
  min_r <- min(r)
  min_c <- min(c)
  
  r_t <- (r - min_r) / a + 1
  c_t <- (c - min_c) / a + 1
  
  return(data.frame(r_t , c_t))
}

#' Matern's variance estimator modified for hexagonal grids
var_mat_hex <- function(samp, a, attributes) {
  # Translate to the base coordinate system
  t_ix <-  translate_hex_ix(samp[,c('r', 'c')], a)
  
  z_bar <- mean(samp$z_1)
 
  # Subsample these translated indices at a=3 at 1,1
  # TODO may make sense to do this at all possible indices
  # and take the mean??
  centers <- subsample_hex(t_ix, c(1,1), 3)
  
  # TODO umm probably slow but works for now.
  # TODO something is causing zero estimates in some configurations
  # TODO used in other functions, take out and move to neighborhood.r (will also be the home of the rectangular versions)
  neighborhoods <- list()
  for(i in 1:nrow(centers)) {
    row <- centers[i, 1]
    col <- centers[i, 2]
    neighborhood <- matrix(NA, nrow=6, ncol=3)
    
    neighborhood[,1] <- c(row-1, row-1, row, row, row+1, row+1)
    neighborhood[,2] <- c(col-1, col+1, col-2, col+2, col-1, col+1)
    neighborhood[,3] <- c(1, 1, 1, -1, -1, -1)
    
    neighborhood <- data.frame(neighborhood)
    neighborhood$r_t <- row
    neighborhood$c_t <- col
    neighborhoods[[i]] <- neighborhood
  }
  
 neighborhoods <- bind_rows(neighborhoods)
 colnames(neighborhoods)[1:3]  <- c('r_n', 'c_n', 'grp')
 samp[,c('r_t', 'c_t')] <- t_ix
 
 # A dataframe such that each neighbor point is paired with its center
 neighborhoods <- merge(neighborhoods, samp, by.x=c('r_n', 'c_n'), by.y=c('r_t', 'c_t'), all.x=TRUE)
 
 # Get the number of neighborhoods with at least one point in Q
 in_Q <- neighborhoods %>%
   mutate(z_ct = z_1 - z_bar) %>%
   replace_na(list(z_ct=0))
 
 # What is the size of each group?
 n_groups <- neighborhoods %>%
   filter(is.na(r)) %>%
   group_by(r_t, c_t) %>%
   summarize(n = n())
  
 q <- 9
 n <- nrow(samp)
 
 Ti <- in_Q %>%
   group_by(r_t, c_t) %>%
   summarize(test = sum(grp * z_ct)^2 / 4)
 
 var <- (q * sum(Ti$test)) / n^2
 
 return(var)
}

#' Calculate the variance using a denominator of n
#' instead of n-1
pop_var <- function(z) {
  n <- length(z)
  ssq <- sum((z - mean(z))^2)
  ssq/n
}

#' Non-overlapping neighborhood variance
#' estimator for hexagonal nearest neighbors.
var_nnbh_hex <- function(samp, a) {
  # Translate to the base coordinate system
  t_ix <-  translate_hex_ix(samp[,c('r', 'c')], a)
  centers <- subsample_hex(t_ix, c(1,1), 3)
  samp[,c('r_t', 'c_t')] <- t_ix
  
  neighborhoods <- list()
  for(i in 1:nrow(centers)) {
    row <- centers[i, 1]
    col <- centers[i, 2]
    neighborhood <- matrix(NA, nrow=7, ncol=2)
    
    neighborhood[,1] <- c(row-1, row-1, row, row, row, row+1, row+1)
    neighborhood[,2] <- c(col-1, col+1, col-2, col, col+2, col-1, col+1)
    
    neighborhood <- data.frame(neighborhood)
    neighborhood$r_t <- row
    neighborhood$c_t <- col
    neighborhoods[[i]] <- neighborhood
  }
  
 neighborhoods <- bind_rows(neighborhoods)
 colnames(neighborhoods)[1:2]  <- c('r_n', 'c_n')
 
 N <- nrow(samp)
 N_neighbs <- neighborhoods %>%
   group_by(r_t, c_t) %>%
   summarize(n=n())
 
 N_neighbs <- nrow(N_neighbs)
 
 neighborhoods <- merge(neighborhoods, samp, by.x=c('r_n', 'c_n'), by.y=c('r_t', 'c_t'), all.x=TRUE)
 neighborhoods <- neighborhoods %>%
   filter(!is.na(z_1)) %>%
   group_by(r_t, c_t) %>%
   summarize(pop_var = pop_var(z_1), q_j = n()) %>%
   mutate(N_j = q_j * a^2) %>%
   mutate(w_j_sq = (N_j / N)^2, fpc = ((N_j - q_j) / N_j)) %>%
   mutate(nbh_var = (1/N_neighbs)^2 *  (pop_var / q_j) * fpc) # TODO N_neighbs needs to be fixed
 
 sum(neighborhoods$nbh_var)
 
 #neighborhoods %>%
#   group_by()
 
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
