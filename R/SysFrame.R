library(sp)

setClass('SysFrame', contains='SpatialPointsDataFrame',
         slots=list(attributes="character"))

#' Base class for systematic sampling frames
#' 
#' `SysFrame` represents a generic systematic sampling frame. Nearly all users will want to use the sub-classes `HexFrame()`
#'  or `RectFrame()` corresponding to hexagonal and rectangular sampling frames, respectively.
#' 
#' @param attributes A character vector indicating the columns that will be treated as attributes, i.e.
#' the response variables at each sample position
#' @return An object of class `SysFrame`
#' @export
SysFrame <- function(splydf, attributes=character()) {
  sys_frame <- new('SysFrame')
  sys_frame@data <- splydf@data
  sys_frame@bbox <- splydf@bbox
  sys_frame@proj4string <- splydf@proj4string
  sys_frame@coords <- splydf@coords
  sys_frame@attributes <- attributes
  sys_frame
}

setGeneric('merge', function(x, y, ...) {
  standardGeneric('merge')
})

setMethod('merge', signature(x='SysFrame', y='data.frame'), 
  function(x, y, by=intersect(names(x), names(y)), by.x=by, 
           by.y=by, all.x=TRUE, suffixes = c(".x",".y"), 
           incomparables = NULL, duplicateGeoms=FALSE, ...){
    # Originally from sp - merge.R
    # Author: Robert J. Hijmans
    # Date : November 2011 / October 2015
    # Version 2
    # Licence GPL v3
    if (!('data' %in% slotNames(x)))
      stop('x has no attributes')
    
    
    x$temp_ID <- 1:nrow(x)
    d <- merge(x@data, y, by=by, by.x=by.x, by.y=by.y, suffixes=suffixes, 
               incomparables=incomparables, all.x=all.x, all.y=FALSE)
    
    if(nrow(d) == 0) {
      stop('Empty merge')
    }
    
    if (!all.x) {
      # Spatial* objects cannot have NULL geometries
      if (nrow(d) == 0) {
        warning('no matching records')
        return(NULL)
      }
    }
    
    # sort the merged table
    d <- d[order(d$temp_ID), ]
    
    # Normally we want one-to-one joins with spatial data
    if (!duplicateGeoms) {
      if (any(table(d$DoNotUse_temp_sequential_ID_963) > 1)) {
        stop('non-unique matches detected')
      }
    } 
    
    # duplicate (duplicateGeoms = TRUE) or remove (all.x=FALSE) records if needed
    x <- x[d$temp_ID, ]
    
    d$temp_ID <- NULL
    x@data <- d
    x
}
)