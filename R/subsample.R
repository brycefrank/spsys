#' For a given sampling interval, return the set of all possible starting positions
#' 
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


#' Transforms a set of population unit coordinates to an integer index. Used in the
#' `index_hex` and `index_sq` functions.
#' 
#' Note that input coordinates must be spatially aligned vertically and horizontally. This 
#' may not always be the case, e.g. for some spatial projections that curve geographic space.
#' 
#' @param coords A dataframe of x,y coordinates of population unit locations.
#' @return A dataframe with two columns, r and c, that correspond to row and column indices
transform_coords <- function(coords) {
  x_shift <- coords[,1] - min(coords[,1])
  y_shift <- coords[,2] - max(coords[,2])
  
  u_x <- unique(x_shift)
  u_x <- u_x[order(u_x)]
  h_dist <- u_x[[2]] - u_x[[1]]
  
  u_y <- unique(y_shift)
  u_y <- u_y[order(u_y)]
  v_dist <- u_y[[2]] - u_y[[1]]
  
  # Create a dataframe of integers representing a grid of hexagons
  c <-    x_shift / h_dist + 1
  r <-    abs(y_shift / v_dist) + 1
  
  # There may be some floating point error
  c <- round(c)
  r <- round(r)
  
  return(data.frame(r=r, c=c))
  
}

#' Creates a dataframe of hexagon indices for use in `subsample_sq`.
#' 
#' @param square_polys A SpatialPolygonsDataFrame of a square grid.
#' @return A dataframe with two columns, r and c, that correspond to row and column indices
#' of the hexagon grid.
index_sq <- function(square_polys) {
  square_polys@polygons[[1]]@Polygons
  
  top_points <- matrix(0, nrow=length(square_polys@polygons), ncol=2)
  i <- 1
  for(polygon in square_polys@polygons) {
    poly_coords <- polygon@Polygons[[1]]@coords
    max_y_ix <- which(poly_coords[,2] == max(poly_coords[,2]))[[1]]
    top_points[i,] <- poly_coords[max_y_ix,]
    i <- i + 1
  }
  
  sq_ix <- transform_coords(top_points)
  return(sq_ix)
  
}

# TODO just make this a method
subsample_hex_ix <- function(hex_ix, start_pos, a) {
  max_r <- max(hex_ix$r)
  max_c <- max(hex_ix$c)
  
  r_seq <- seq(0, max_r-1, a)
  r_samp <- c()
  c_samp <- c()
  
  j <- 0
  for(r in r_seq) {
    if (j %% 2 == 0) {
      add <- seq(0, max_c-1, a*2)
    } else {
      add <- seq(-a, max_c-1, a*2)
      add <- add[add>0]
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


#' There is a way to subsample a set of hexagonal
#' indices such that we obtain a compact set of 
#' neighborhoods. Implemented here.
subsample_hex_ix_compact <- function(ix) {
  max_r <- max(max(ix$r), 14)
  max_c <- max(max(ix$c), 14)
  
  samp_ix <- list()
  
  browser()
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