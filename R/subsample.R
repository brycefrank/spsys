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
#' @param coords A dataframe of x,y coordinates of population unit locations.
#' @return A dataframe with two columns, r and c, that correspond to row and column indices
transform_coords <- function(coords) {
  # The top left point (min x max y) needs to be transformed linearly such that is (1,1)
  top_left_ix <- which(coords[,1] == min(coords[,1]) & coords[,2] == max(coords[,2]))
  top_left <- coords[top_left_ix,]
  
  x_shift <- (coords[,1] - min(coords[,1]))
  h_dist <- min(x_shift[x_shift!=0])
  
  y_shift <- (coords[,2] - min(coords[,2]))
  v_dist <- min(y_shift[y_shift!=0])
  
  # Create a dataframe of integers representing a grid of hexagons
  c <- ((coords[,1] - min(coords[,1]))  / h_dist) + 1
  r <- abs((coords[,2] - max(coords[,2])) / v_dist) + 1
  
  c <- round(c, 3)
  r <- round(r, 3)
  
  return(data.frame(c=c, r=r))
  
}

#' Creates a dataframe of hexagon indices for use in `subsample_hex`
#' 
#' @param hex_polys A SpatialPolygonsDataFrame of a hexagon grid.
#' @return A dataframe with two columns, r and c, that correspond to row and column indices
#' of the hexagon grid.
index_hex <- function(hex_polys) {
  top_points <- matrix(0, nrow=length(hex_polys@polygons), ncol=2)
  i <- 1
  for(polygon in hex_polys@polygons) {
    poly_coords <- polygon@Polygons[[1]]@coords
    max_y_ix <- which(poly_coords[,2] == max(poly_coords[,2]))
    top_points[i,] <- poly_coords[max_y_ix,]
    i <- i + 1
  }
  
  hex_ix <- transform_coords(top_points)
  return(hex_ix)
  
}

#' Subsamples a dataframe of hexagonal indices produced by the `index_hex` function.
#' For some specified starting index and sampling interval, subsamples a hexagonal index 
#' set such that the spatial structure is preserved, i.e. subsamples also retain a
#' hexagonal structure.
#' 
#' @param hex_ix A dataframe of hex indices from the `index_hex` function.
#' @param start_pos the starting position
#' @param a The order of the subset: 2 represents sampling every other
#' @return A dataframe of row and column indices that have been sampled.
subsample_hex <- function(hex_ix, start_pos, a) {
  max_row <- max(hex_ix$r)
  max_col <- max(hex_ix$c)
  
  r_start <- start_pos[[1]]
  c_start <- start_pos[[2]]
  
  r_seq <- seq(r_start, max_row, a)
  r_samp <- c()
  c_samp <- c()
  
  j <- 1
  for(r in r_seq) {
    if (j %% 2 == 0) {
      add <- seq(c_start-a, max_col, a*2)
    } else {
      add <- seq(c_start-2*a, max_col, a*2)
    }
    
    add <- add[add>0]
    c_samp <- c(c_samp, add)
    r_samp <- c(r_samp, rep(r, length(add)))
    
    
    j <- j + 1
  }
  
  return(data.frame(r = r_samp, c = c_samp))
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