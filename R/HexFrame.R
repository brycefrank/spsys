#' @include SysFrame.R
#' @include subsample.R

#validate_hex_frame <- function(object) {
#  ix <- object@data[,c('r', 'c')]
#  er_oc <- sum(ix$c[ix$r %% 2 == 0] %% 2 != 0)
#  or_ec <- sum(ix$c[ix$r %% 2 != 0] %% 2 == 0)
#  
#  if(er_oc + or_ec > 0) {
#    error_msg <- sprintf('This index had %i even-row odd-row matchings and %i odd-row even-row matchings,
#              there should be 0 of each', er_oc, or_ec)
#    return(error_msg)
#  } else {
#    return(TRUE)
#  }
#}

setClass('HexFrame', contains='SysFrame',
         slots = list(
           a='numeric',
           N='numeric'
         ))

# TODO not sure why this needs to be repeated from sysframe
# is it possible to do something like callNextMethod but for
# instantiation?
HexFrame <- function(splydf, attributes=character(), index=NA, standardize=TRUE, a=1, N=Inf) {
  sys_frame <- SysFrame(splydf, attributes)
  
  hex_frame <- new('HexFrame')
  hex_frame@data <- sys_frame@data
  hex_frame@bbox <- sys_frame@bbox
  hex_frame@proj4string <- sys_frame@proj4string
  hex_frame@coords <- sys_frame@coords
  hex_frame@attributes <- sys_frame@attributes
  hex_frame@N <- N
  
  
  # TODO make this condition a bit cleaner
  if(length(index) > 1) {
    hex_frame@data[,c('r', 'c')] <- index
  } else {
    hex_frame@data[,c('r', 'c')] <- transform_coords(hex_frame@coords)
  }
  
  hex_frame <- hex_frame[complete.cases(hex_frame@data[,c('r', 'c')]),]
  
  # Ensure the rows start from 1 and the columns start from 1
  min_c <- min(hex_frame@data$c)
  min_r <- min(hex_frame@data$r)
  
  
  hex_frame@data$r <- hex_frame@data$r - min_r + 1
  min_r_new <- min(hex_frame@data$r)
  # For columns, if the min_r and min_c "disagree" we have a 'beta'
  # configured index set, and need to bump the columns over by one index
  if (min_r %% 2 == 0 & min_c %% 2!=0) {
    hex_frame@data$c <- hex_frame@data$c - min_c + 2
  } else if (min_r %% 2 !=0 & min_c %% 2 == 0) {
    hex_frame@data$c <- hex_frame@data$c - min_c + 2
  } else {
    # In some cases the grid is still misaligned from the standard set
    # and that is checked for here and corrected.
    min_c_at_min_r <- min(hex_frame@data$c[hex_frame@data$r == min_r_new])
    if(min_c %%2 == 0 & min_c_at_min_r %% 2 != 0) {
      hex_frame@data$c <- hex_frame@data$c - min_c + 2
    } else if (min_c %% 2 !=0 & min_c_at_min_r %% 2 == 0) {
      hex_frame@data$c <- hex_frame@data$c - min_c + 2
    } else {
      hex_frame@data$c <- hex_frame@data$c - min_c + 1
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

setMethod('subsample', 'HexFrame', function(object, start_pos, a) {
  if(nrow(object) <= a) {
    stop('Attempting to subsample SysFrame that is the same size
          or smaller than a.')
  }
  
  samp_ix <- subsample_hex_ix(object@data[,c('r', 'c')], start_pos, a)
  object@data$TEMP <- 1:nrow(object@data)
  keep_ix <- merge(object@data, samp_ix, by=c('r', 'c'), all.x=FALSE)[,c('r', 'c', 'TEMP')]
  
  # Now "collapse" keep_ix to the standard index set
  keep_ix$r <- (keep_ix$r - start_pos[[1]]) / a + start_pos[[1]]
  keep_ix$c <- (keep_ix$c - start_pos[[2]]) / a + start_pos[[2]]
  
  new_spdf <- SpatialPointsDataFrame(coords = object@coords[keep_ix$TEMP,], data=object@data[keep_ix$TEMP,])
  crs(new_spdf) <- crs(object)
  
  samp <- HexFrame(new_spdf, attributes=object@attributes, index=keep_ix[,c('r', 'c')])
  samp@data$TEMP <- NA
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