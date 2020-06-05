#' @include SysFrame.R


setClass('HexFrame', contains='SysFrame',
         slots = list(
           a='numeric'
         ))

# TODO not sure why this needs to be repeated from sysframe
# is it possible to do something like callNextMethod but for
# instantiation?
HexFrame <- function(splydf, attributes=character(), index=NA, standardize=TRUE) {
  sys_frame <- SysFrame(splydf, attributes)
  
  hex_frame <- new('HexFrame')
  hex_frame@data <- sys_frame@data
  hex_frame@bbox <- sys_frame@bbox
  hex_frame@proj4string <- sys_frame@proj4string
  hex_frame@coords <- sys_frame@coords
  hex_frame@attributes <- sys_frame@attributes
  
  # TODO make this condition a bit cleaner
  if(length(index) > 1) {
    hex_frame@data[,c('r', 'c')] <- index
  } else {
    hex_frame@data[,c('r', 'c')] <- transform_coords(hex_frame@coords)
  }
  
  hex_frame <- hex_frame[complete.cases(hex_frame@data[,c('r', 'c')]),]
  
  if(standardize) {
    
    hex_frame@data$c <- ceiling( (hex_frame@data$c - min(hex_frame@data$c))) + 1
    hex_frame@data$r <- ((hex_frame@data$r - min(hex_frame@data$r))) + 1
    
    # If odd rows correspond to even columns, add 1 to column indices
    first_r <- hex_frame@data$r[[1]]
    first_c <- hex_frame@data$c[hex_frame@data$r == first_r][[1]]
    
    print(first_r)
    print(first_c)
    
    odd_row_even_col <- (first_r %% 2 != 0) & (first_c %% 2 == 0)
    even_row_odd_col <- (first_r %% 2 == 0) & (first_c %% 2 != 0)
    
    if(odd_row_even_col | even_row_odd_col) {
      hex_frame@data$c <- hex_frame@data$c + 1
    }
    
  }
  
  hex_frame
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
setGeneric('subsample', function(object, start_pos, a, standardize=TRUE){
  standardGeneric('subsample')
})

setMethod('merge', signature(x='HexFrame', y='data.frame'),
  function(x, y, ...) {
    merged <- callNextMethod(x, y, ...)
    HexFrame(merged, attributes=x@attributes, index=merged@data[,c('r', 'c')])
  }
)

setMethod('subsample', 'HexFrame', function(object, start_pos, a, standardize) {
  if(nrow(object) <= a) {
    stop('Attempting to subsample SysFrame that is the same size
          or smaller than a.')
  }
  
  max_r <- max(object@data$r)
  max_c <- max(object@data$c)
  
  r_samp <- seq(1, max_r, a)
  c_samp <- seq(1, max_c, a)
  
  samp_ix <- expand.grid(r_samp, c_samp)
  colnames(samp_ix) <- c('r', 'c')
  
  # Even-row odd-columns and odd-row even-columns should not exist
  # Get only even-row even-column indices and odd-row odd-column indices
  even_row_even_column <- (samp_ix$r %% 2 == 0) & (samp_ix$c %% 2 == 0)
  odd_row_odd_column   <- (samp_ix$r %% 2 != 0) & (samp_ix$c %% 2 != 0)
  
  samp_ix <- rbind(samp_ix[even_row_even_column,], samp_ix[odd_row_odd_column,])
  
  samp_ix[,1] <- samp_ix[,1] + start_pos[[1]]
  samp_ix[,2] <- samp_ix[,2] + start_pos[[2]]
  
  samp <- merge(object, samp_ix, by=c('r', 'c'), all.x=FALSE)
  
  if(standardize) {
    samp@data$c <- ceiling( (samp@data$c - min(samp@data$c)) / a) + 1
    samp@data$r <- ((samp@data$r - min(samp@data$r)) / a) + 1
  }
  
  samp@a <- a
  return(samp)
})

#' Creates a dataframe of hexagon indices for use in `subsample_hex`
#' 
#' @param hex_polys A SpatialPolygonsDataFrame of a hexagon grid.
#' @return A dataframe with two columns, r and c, that correspond to row and column indices
#' of the hexagon grid.
index_hex_polys <- function(hex_polys) {
  top_points <- matrix(0, nrow=length(hex_polys), ncol=2)
  i <- 1
  for(polygon in hex_polys) {
    poly_coords <- polygon@Polygons[[1]]@coords
    max_y_ix <- which(poly_coords[,2] == max(poly_coords[,2]))
    
    if(length(max_y_ix) > 1) {
      max_y_ix <- max_y_ix[[1]]
    }
    
    top_points[i,] <- poly_coords[max_y_ix,]
    i <- i + 1
  }
  
  hex_ix <- transform_coords(top_points)
  return(hex_ix)
}