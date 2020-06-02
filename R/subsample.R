library(sp)
library(rgeos)


#' For a given sampling interval, return the set of all possible starting positions
#' 
#' @param a The sampling interval
#' @return  A dataframe of starting positions for all possible systematic samples 
#' for the interval a.
subsample_starts <- function(a) {
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
}

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
  
  return(data.frame(c=c, r=r))
  
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

#' Subsamples a dataframe of square indices produced by the `index_sq` function.
#' For some specified starting index and sampling interval, subsamples a square index 
#' set such that the spatial structure is preserved, i.e. subsamples also retain a
#' square structure.
#' 
#' @param sq_ix A dataframe of square indices from the `index_sq` function.
#' @param start_pos the starting position
#' @param a The order of the subset: 2 represents sampling every other
#' @return A dataframe of row and column indices that have been sampled.
subsample_sq <- function(sq_ix, start_pos, a) {
  max_row <- max(sq_ix$r)
  max_col <- max(sq_ix$c)
  
  r_start <- start_pos[[1]]
  c_start <- start_pos[[2]]
  
  r_seq <- seq(r_start, max_row, a)
  c_seq <- seq(c_start, max_col, a)
  
  a_grid <- expand.grid(r_seq, c_seq)
  colnames(a_grid) <- c('r', 'c')
  a_grid <- a_grid[a_grid$r <= max_row || a_grid$c <= max_col,]
  return(a_grid)
}