#' @include SysFrame.R

setClass('RectFrame', contains='SysFrame',
         slots = list(
           a='numeric',
           N='numeric'
         ))

#' Rectangular sampling frame
#' 
#' `RectFrame` represents sampling positions oriented in a rectangular pattern. 
#' All `RectFrame`s must be constructed from an existing `sp::SpatialPointsDataFrame` that represents
#' the sample positions and contains their observation values.
#' 
#' @param attributes A character vector indicating the columns that will be treated as attributes, i.e.
#' the response variables at each sample position
#' @param index An optional dataframe that manually specifies the row and column indices
#' @param a The sampling interval at which the systematic sample was collected originally. Typically set
#' internally by `subsample`
#' @param N The population size
#' @param d_x The distance between points in the x dimension
#' @param d_y The distance between points in the y dimension
#' @return An object of class `RectFrame`
#' @example 
#' hex_frame <- RectFrame(rect_pts, attributes=c('z_1'))
#' @export
RectFrame <- function(splydf, attributes=character(), index=NA, N=Inf, a=1, d_x=NA, d_y=NA) {
  sys_frame <- SysFrame(splydf, attributes)
  
  rect_frame <- new('RectFrame')
  rect_frame@data <- sys_frame@data
  rect_frame@bbox <- sys_frame@bbox
  rect_frame@proj4string <- sys_frame@proj4string
  rect_frame@coords <- sys_frame@coords
  rect_frame@attributes <- sys_frame@attributes
  
  if(length(index) > 1) {
    rect_frame@data[,c('r', 'c')] <- index
  } else {
    rect_frame@data[,c('r', 'c')] <- transform_coords(rect_frame@coords, d_x, d_y)
  }
  
  # Ensure the rows start from 1 and the columns start from 1
  min_r <- min(rect_frame@data$r)
  min_c <- min(rect_frame@data$c)
  
  rect_frame@data$r <- rect_frame@data$r - min_r + 1
  rect_frame@data$c <- rect_frame@data$c - min_c + 1
  
  rect_frame
}

setMethod('subsample', 'RectFrame', function(object, start_pos, a) {
  if(nrow(object) <= a) {
    stop('Attempting to subsample SysFrame that is the same size
          or smaller than a.')
  }
  
  samp_ix <- subsample_rect_ix(object@data[,c('r', 'c')], start_pos, a)
  object@data$TEMP <- 1:nrow(object@data)
  keep_ix <- merge(object@data, samp_ix, by=c('r', 'c'), all.x=FALSE)[,c('r', 'c', 'TEMP')]
  
  # Now "collapse" keep_ix to the standard index set
  keep_ix$r <- (keep_ix$r - start_pos[[1]]) / a + start_pos[[1]]
  keep_ix$c <- (keep_ix$c - start_pos[[2]]) / a + start_pos[[2]]
  
  new_spdf <- SpatialPointsDataFrame(coords = object@coords[keep_ix$TEMP,], data=object@data[keep_ix$TEMP,])
  crs(new_spdf) <- crs(object)
  
  samp <- RectFrame(new_spdf, attributes=object@attributes, index=keep_ix[,c('r', 'c')])
  samp@data$TEMP <- NA
  samp@a <- a
  return(samp)
})

