#' @include SysFrame.R


setClass('HexFrame', contains='SysFrame',
         slots = list(
           a='numeric'
         ))

# TODO not sure why this needs to be repeated from sysframe
# is it possible to do something like callNextMethod but for
# instantiation?
HexFrame <- function(splydf, attributes=character(), index=NA) {
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
setGeneric('subsample', function(object, start_pos, a){
  standardGeneric('subsample')
})

setMethod('merge', signature(x='HexFrame', y='data.frame'),
  function(x, y, ...) {
    HexFrame(callNextMethod(x,y, ...), attributes=x@attributes)
  }
)

setMethod('subsample', 'HexFrame', function(object, start_pos, a) {
  max_row <- max(object@data$r[!is.na(object@data$r)])
  max_col <- max(object@data$c[!is.na(object@data$c)])
  
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
 
  # I am wondering if this is correct? 
  # There was a problem with RectFrame where it was returning
  # the same sample...do these need to be standardized?
  samp_ix   <- data.frame(r = r_samp, c = c_samp)
  samp <- merge(object, samp_ix, by=c('r', 'c'), all.x=FALSE)
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