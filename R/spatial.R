# Functions for dealing with spatial components of some
# variance estimators
library(Rfast)

#' Create the binary neigborhood matrix for D'Orazio's estimators
#' 
#' @param sys_frame A `SysFrame` object
#' @param order An integer representing the neighborhood order. For `HexFrame` 1 represents 
#' the immediate neighbors, 2 represents the next 'layer' of neighbors and so on. For `RectFrame` 1 represents
#' the NESW neighbors, 2 represents the 8 nearest neighbors plus four next-nearest NESW neighbors and so on.
#' @return
#' @keywords internal
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

#' Retrieve the Geary's C for a given `SysFrame`
#' 
#' @param sys_frame A `SysFrame` object
#' @param order An integer representing the neighborhood order. For `HexFrame` 1 represents 
#' the immediate neighbors, 2 represents the next 'layer' of neighbors and so on. For `RectFrame` 1 represents
#' the NESW neighbors, 2 represents the 8 nearest neighbors plus four next-nearest NESW neighbors and so on.
#' @return A named vector of Geary's C values for each attribute
#' @keywords internal
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
    
    # If all elements are 0 then C will be equal to 1 for all attributes
    if(sum(numerator == 0) == p) {
      return(rep(1, p))
    }
    
    denominator <- 2 * sum(W) * colSums((sweep(att_df, 2, colMeans(att_df))^2))
    C <- ((n-1) * numerator) / denominator
    C <- as.vector(C)
    names(C) <- sys_frame@attributes
    
    return(C)
  }
)


#' Retrieve the Moran's I for a given `SysFrame`
#' 
#' @param sys_frame A `SysFrame` object
#' @param order An integer representing the neighborhood order. For `HexFrame` 1 represents 
#' the immediate neighbors, 2 represents the next 'layer' of neighbors and so on. For `RectFrame` 1 represents
#' the NESW neighbors, 2 represents the 8 nearest neighbors plus four next-nearest NESW neighbors and so on.
#' @return A named vector of Moran's C values for each attribute
#' @keywords internal
setGeneric('morans_i', function(sys_frame, ...) {
  standardGeneric('morans_i')
})

setMethod('morans_i', signature(sys_frame='SysFrame'),
  function(sys_frame, order=1) {
    n <- nrow(sys_frame@data)
    W <- neighborhood_matrix(sys_frame, order=order)
    
    att_df <- as.matrix(sys_frame@data[, sys_frame@attributes, drop=FALSE])
    p <- length(names(sys_frame[,sys_frame@attributes]))
    
    neighbor_inds <- which(W == 1, arr.ind=TRUE)
    numerator <- rep(0, p)
    att_means <- colMeans(att_df)
    for(k in 1:nrow(neighbor_inds)) {
      i <- neighbor_inds[[k, 1]]
      j <- neighbor_inds[[k, 2]]
      
      e_ij <- (att_df[i, , drop=FALSE] - att_means) * (att_df[j, , drop=FALSE] - att_means)
      numerator <- numerator + e_ij
    }
    
    denominator <- colSums(((sweep(att_df, 2, colMeans(att_df)))^2))
    morans_I <- (n / sum(W)) * (numerator / denominator)
    morans_I <- as.vector(morans_I)
    names(morans_I) <- sys_frame@attributes
    
    return(morans_I)
  }
)