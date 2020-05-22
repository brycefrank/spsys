setClass('SysFrame', contains='SpatialPolygonsDataFrame')

SysFrame <- function(splydf) {
  sys_frame <- new('SysFrame')
  sys_frame@data <- splydf@data
  sys_frame@bbox <- splydf@bbox
  sys_frame@proj4string <- splydf@proj4string
  sys_frame@plotOrder <- splydf@plotOrder
  sys_frame@polygons <- splydf@polygons
  sys_frame
}

setClass('HexFrame', contains='SysFrame',
         slots = list(
           attributes='character',
           a='integer'
         ))


HexFrame <- function(splydf, attributes) {
  sys_frame <- SysFrame(splydf)
  
  hex_frame <- new('HexFrame')
  hex_frame@attributes <- attributes
  hex_frame@data <- sys_frame@data
  hex_frame@bbox <- sys_frame@bbox
  hex_frame@proj4string <- sys_frame@proj4string
  hex_frame@plotOrder <- sys_frame@plotOrder
  hex_frame@polygons <- sys_frame@polygons
  hex_frame@data[,c('r', 'c')] <- index_hex(hex_frame@polygons)
  
  hex_frame
}

setGeneric('subsample', function(object, start_pos, a){
  standardGeneric('subsample')
})

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
  
  samp_ix   <- data.frame(r = r_samp, c = c_samp)
  #samp <- HexFrame(sp::merge(object, samp_ix, by=c('r', 'c'), all.x=FALSE))
  #fmerge(object, samp_ix, by=c('r', 'c'), all.x=FALSE)
  foo(object, 2)
})

#' Creates a dataframe of hexagon indices for use in `subsample_hex`
#' 
#' @param hex_polys A SpatialPolygonsDataFrame of a hexagon grid.
#' @return A dataframe with two columns, r and c, that correspond to row and column indices
#' of the hexagon grid.
index_hex <- function(hex_polys) {
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