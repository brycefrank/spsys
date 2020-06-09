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

setGeneric('neighborhood_matrix', function(sys_frame, ...) {
  standardGeneric('neighborhood_matrix')
})

setMethod('neighborhood_matrix', signature(sys_frame='HexFrame'), 
  function(sys_frame, order=1) {
    W <- Dist(sys_frame@data[,c('r', 'c')]) <= 2*order
    diag(W) <- 0
    return(W)
  }
)

setMethod('neighborhood_matrix', signature(sys_frame='RectFrame'), 
  function(sys_frame, order=1) {
    # TODO Check if this is ok, it is just grabbing the NSEW neighbors, not
    # the diagonals
    W <- Dist(sys_frame@data[,c('r', 'c')]) <= 2*order
    diag(W) <- 0
    return(W)
  }
)

setGeneric('gearys_c', function(sys_frame, ...) {
  standardGeneric('gearys_c')
})

setMethod('gearys_c', signature(sys_frame='SysFrame'),
  function(sys_frame, order=1) {
    n <- nrow(sys_frame@data)
    W <- neighborhood_matrix(sys_frame, order=order)
    
    att_df <- as.matrix(sys_frame@data[, sys_frame@attributes, drop=FALSE])
    p <- length(names(sys_frame[,sys_frame@attributes]))
    
    # TODO could be optimized further
    neighbor_inds <- which(W == 1, arr.ind=TRUE)
    numerator <- rep(0, p)
    for(k in 1:nrow(neighbor_inds)) {
      i <- neighbor_inds[[k, 1]]
      j <- neighbor_inds[[k, 2]]
      
      e_ij <- (att_df[i, , drop=FALSE] - att_df[j, , drop=FALSE])^2
      numerator <- numerator + e_ij
    }
    
    denominator <- 2 * sum(W) * colSums(((att_df - colMeans(att_df))^2))
    C <- ((n-1) * numerator) / denominator
    
    return(C)
  }
)

setGeneric('morans_i', function(sys_frame, ...) {
  standardGeneric('morans_i')
})

setMethod('morans_i', signature(sys_frame='SysFrame'),
  function(sys_frame, order=1) {
    n <- nrow(sys_frame@data)
    W <- neighborhood_matrix(sys_frame, order=order)
    
    att_df <- as.matrix(sys_frame@data[, sys_frame@attributes, drop=FALSE])
    p <- length(names(sys_frame[,sys_frame@attributes]))
    
    # TODO could be optimized further
    neighbor_inds <- which(W == 1, arr.ind=TRUE)
    numerator <- rep(0, p)
    att_means <- colMeans(att_df)
    for(k in 1:nrow(neighbor_inds)) {
      i <- neighbor_inds[[k, 1]]
      j <- neighbor_inds[[k, 2]]
      
      e_ij <- (att_df[i, , drop=FALSE] - att_means) * (att_df[j, , drop=FALSE] - att_means)
      numerator <- numerator + e_ij
    }
    
    denominator <- colSums(((att_df - colMeans(att_df))^2))
    morans_I <- (n / sum(W)) * (numerator / denominator)
    
    return(morans_I)
  }
)