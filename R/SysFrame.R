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