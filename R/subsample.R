# Various functions for creating subsets of sampling frames used for assessment
# and to create neighborhood centers and anchors.

#' Retrieve the set of all possible starting positions
#' 
#' In two-dimensional systematic samples, all possible samples can be obtained
#' merely by moving the starting position within a "starting region" contained
#' by the top-left-most sample group. For a sampling interveal `a` the hexagonal
#' and rectangular configurations will contain `a^2` possible samples. This function
#' retrieves those starting positions and is used internally as part of `compare_estimators()`.
#' 
#' @param sys_frame A `SysFrame` object
#' @param a The sampling interval
#' @return  A dataframe of starting positions for all possible systematic samples 
#' for the interval a.
setGeneric('subsample_starts', function(sys_frame, a) {
  standardGeneric('subsample_starts')
})

setMethod('subsample_starts', signature = list(sys_frame='RectFrame', a='numeric'), 
  function(sys_frame, a) {
    r_starts <- seq(1, a)
    c_starts <- seq(1, a)
    
    starts <- expand.grid(r_starts, c_starts)
    starts <- data.frame(r = starts[,1], c = starts[,2])
    return(starts)
})

setMethod('subsample_starts', signature = list(sys_frame='HexFrame', a='numeric'), 
  function(sys_frame, a) {
    r_starts <- c()
    c_starts <- c()
    
    for(r in 1:a) {
     if(r %% 2 == 0)  {
       c_starts <- c(c_starts, seq(1, 2*a, 2) + 1)
     } else {
       c_starts <- c(c_starts, seq(1, 2*a, 2))
     }
       r_starts <- c(r_starts, rep(r, a))
    }
    return(data.frame(r=r_starts, c=c_starts))
})


#' Transforms a set of population unit coordinates to an integer index.
#' 
#' Note that input coordinates must be spatially aligned vertically and horizontally. This 
#' may not always be the case, e.g. for some spatial projections that curve geographic space.
#' 
#' @param coords A dataframe of x,y coordinates of population unit locations.
#' @param d_x The distance between points in the x dimension
#' @param d_y The distance between points in the y dimension
#' @return A dataframe with two columns, r and c, that correspond to row and column indices
#' @keywords internal
transform_coords <- function(coords, d_x=NA, d_y=NA) {
  x_shift <- coords[,1] - min(coords[,1])
  y_shift <- coords[,2] - max(coords[,2])
  
  if(is.na(d_x)) {
    u_x <- unique(x_shift)
    u_x <- u_x[order(u_x)]
    d_x <- u_x[[2]] - u_x[[1]]
  }
  
  if(is.na(d_y)) {
    u_y <- unique(y_shift)
    u_y <- u_y[order(u_y)]
    d_y <- u_y[[2]] - u_y[[1]]
  }
  
  
  # Create a dataframe of integers representing a grid of hexagons
  c <-    x_shift / d_x + 1
  r <-    abs(y_shift / d_y) + 1
  
  # There may be some floating point error
  c <- round(c)
  r <- round(r)
  
  return(data.frame(r=r, c=c))
  
}

#' Retrieves the subsample for a given set of hexagonal indices
#' 
#' @param hex_ix A dataframe of hexagonal indices
#' @param start_pos A starting position
#' @param a A sampling interval
#' @return A dataframe of sample indices
#' @keywords internal
subsample_hex_ix <- function(hex_ix, start_pos, a) {
  max_r <- max(hex_ix$r)
  max_c <- max(hex_ix$c)
  
  r_seq <- seq(0, max_r-1, a)
  r_samp <- c()
  c_samp <- c()
  
  j <- 0
  for(r in r_seq) {
    if (j %% 2 == 0) {
      add <- seq(-a*2, max_c-1, a*2)
    } else {
      add <- seq(-a, max_c-1, a*2)
    }
    
    c_samp <- c(c_samp, add)
    r_samp <- c(r_samp, rep(r, length(add)))
    
    j <- j + 1
  }
  
  r_samp <- r_samp + start_pos[[1]]
  c_samp <- c_samp + start_pos[[2]]
  
  samp_ix   <- data.frame(r = r_samp, c = c_samp)
  samp_ix
}


#' Retrieves a subset of indices that correspond to the center of a compact set 
#' of hexagonal neighborhoods
#' 
#' @param ix A dataframe of hexagonal indices
#' @return A dataframe of hexagonal neigborhood centers
#' @keywords internal
subsample_hex_ix_compact <- function(ix) {
  max_r <- max(max(ix$r), 14)
  max_c <- max(max(ix$c), 14)
  
  samp_ix <- list()
  
  # Make a grid for each row
  for(j in 0:14) {
    col_shift <- -j * 5
    
    r_seq <- seq(j, max_r, 14)
    c_seq <- seq(col_shift, max_c, 14)
    c_seq <- c_seq[c_seq>=0]
    
    samp_ix[[j+1]] <- expand.grid(r_seq, c_seq)
  }
  
  samp_ix <- bind_rows(samp_ix)
  
  # Bump over to (1,1 origin)
  samp_ix <- samp_ix+1
  colnames(samp_ix) <- c('r', 'c')
  
  samp_ix
}

#' Retrieves the subsample for a given set of rectangular indices
#' 
#' @param ix A dataframe of hexagonal indices
#' @param start_pos A starting position
#' @param a A sampling interval
#' @return A dataframe of sample indices
#' @keywords internal
subsample_rect_ix <- function(ix, start_pos, a) {
  max_r <- max(ix$r)
  max_c <- max(ix$c)
  
  r_start <- start_pos[[1]]
  c_start <- start_pos[[2]]
  
  r_seq <- seq(r_start, max_r, a)
  c_seq <- seq(c_start, max_c, a)
  
  a_grid <- expand.grid(r_seq, c_seq)
  colnames(a_grid) <- c('r', 'c')
  return(a_grid)
}