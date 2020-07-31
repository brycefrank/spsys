#' @include SysFrame.R
#' @include subsample.R

setClass('HexFrame', contains='SysFrame',
         slots = list(
           a='numeric',
           N='numeric'
         ))


#' Hexagonal sampling frame
#' 
#' `HexFrame` represents sampling positions oriented in a hexagonal (or triangular) pattern. 
#' All `HexFrame`s must be constructed from an existing `sp::SpatialPointsDataFrame` that represents
#' the sample positions and contains their observation values.
#' 
#' @param attributes A character vector indicating the columns that will be treated as attributes, i.e.
#' the response variables at each sample position
#' @param index An optional dataframe that manually specifies the row and column indices
#' @param a The sampling interval at which the systematic sample was collected originally. Typically set
#' internally by `subsample`
#' @param N The population size
#' @return An object of class `HexFrame`
#' @example 
#' hex_frame <- HexFrame(hex_pts, attributes=c('z_1'))
#' @export
HexFrame <- function(splydf, attributes=character(), index=NA, a=1, N=Inf) {
  sys_frame <- SysFrame(splydf, attributes)
  
  hex_frame <- new('HexFrame')
  hex_frame@data <- sys_frame@data
  hex_frame@bbox <- sys_frame@bbox
  hex_frame@proj4string <- sys_frame@proj4string
  hex_frame@coords <- sys_frame@coords
  hex_frame@attributes <- sys_frame@attributes
  hex_frame@N <- N
  
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

#' Create a systematic sample from an existing `SysFrame`
#' 
#' Creates a systematic sample from a `SysFrame` for a specified starting position
#' and sampling interval. Typically used to conduct variance assessments, e.g. as a part of
#' `compare_estimators()`.
#' 
#' @param object A `SysFrame` from which to sample
#' @param start_pos A numeric vector of two elements indicating the starting position
#' @param a The sampling interval
#' @return A new `SysFrame` object that contains only the sampled elements
#' @example 
#' samp <- subsample(hex_frame, c(1,1), 3)
#' @export
setGeneric('subsample', function(object, start_pos, a){
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