# Functions for dealing with spatial components of some
# variance estimators
library(Rfast)

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

#' Gets the hexagonal neighborhoods
#' 
#' @param hex_ix a dataframe with columns r and c corresponding to the row and column
#' of the hexagonal indices
#' @param a the sampling interval used to produce hex_ix
#' @param grouping a vector of 7 elements assigning contrast coefficients to each neighborhood point
#' @return a dataframe with columns r and c (the original sample positions) and
#' columns r_n and c_n that are the coordinates of the neighbors.
get_hex_neighborhoods <- function(hex_ix, a, contrasts=NA) {
  t_ix <-  translate_hex_ix(hex_ix, a)
  centers <- subsample_hex(t_ix, c(1,1), 3)
  hex_ix[,c('r_t', 'c_t')] <- t_ix
  
  neighborhoods <- list()
  for(i in 1:nrow(centers)) {
    row <- centers[i, 1]
    col <- centers[i, 2]
    neighborhood <- matrix(NA, nrow=7, ncol=2)
    
    neighborhood[,1] <- c(row-1, row-1, row, row, row, row+1, row+1)
    neighborhood[,2] <- c(col-1, col+1, col-2, col, col+2, col-1, col+1)
    
    if(!is.na(grouping)) {
      cbind(neighborhood, contrasts)
    }
    
    neighborhood <- data.frame(neighborhood)
    neighborhood$r <- row
    neighborhood$c <- col
    neighborhoods[[i]] <- neighborhood
  }
  
 neighborhoods <- bind_rows(neighborhoods)
 colnames(neighborhoods)[1:2]  <- c('r_n', 'c_n')
 
 # "Untranslate" the neighborhoods
 neighborhoods[,c('r', 'r_n')] <- a*(neighborhoods[,c('r', 'r_n')] - 1) + min(hex_ix[,'r'])
 neighborhoods[,c('c', 'c_n')] <- a*(neighborhoods[,c('c', 'c_n')] - 1) + min(hex_ix[,'c'])
 return(neighborhoods)
}

gearys_C <- function(hex_ix, z, a) {
  N <- nrow(hex_ix)
  W <- Dist(hex_ix) <= 2*a
  
  # TODO a bit of a slow solution here...
  # probably a way to do this with a matrix computation
  # e.g. outer(z, z, "-")
  numerator <- 0 
  for(i in 1:N) {
    for(j in 1:N) {
      e_ij <- z[[i]] - z[[j]]
      w_ij <- W[[i, j]]
      w_ij * e_ij^2
      numerator <- numerator + w_ij
    }
  }
  
  denominator <- 2 * sum(W) * sum((z - mean(z))^2)
  C <- ((N-1) * numerator) / denominator
  return(C)
}