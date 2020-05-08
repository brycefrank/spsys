library(sp)
#library(rgdal)
library(raster)
library(dplyr)
library(maptools)

cast_hex_grid <- function(plots) {
  bbox <- bbox(plots)
  min_x <- bbox[1,1]
  max_x <- bbox[1,2]
  min_y <- bbox[2,1]
  max_y <- bbox[2,2]
  bbox <- data.frame(x=c(min_x, min_x, max_x, max_x), y=c(min_y, max_y, max_y, min_y))
  bbox <- Polygons(list(Polygon(bbox)), 1)
  bbox <- SpatialPolygons(list(bbox))
  
  side_to_side <- 5187 #meters
  HexPts <- spsample(bbox, type="hexagonal", cellsize=side_to_side)
  HexPols <- HexPoints2SpatialPolygons(HexPts)
  row.names(HexPols) <- as.character(seq(1, length(HexPols)))
  
  hex_df <- data.frame(id=row.names(HexPols))
  row.names(hex_df) <- row.names(HexPols)
  HexPols <- SpatialPolygonsDataFrame(HexPols, data=hex_df)
  crs(HexPols) <- newcrs
  
  ov <- over(plots, HexPols)
  colnames(ov) <- 'hex_ix'
  ov$VOL <- plots@data$VOL
  ov <- data.frame(plt_ix = as.numeric(row.names(ov)), hex_ix = ov[,'hex_ix'], VOL=ov[,'VOL'])
  
  first <- ov %>% group_by(hex_ix) %>% slice(1)
  HexPols <- sp::merge(HexPols, first, by.y='hex_ix', by.x='id')
  
  return(HexPols)
}