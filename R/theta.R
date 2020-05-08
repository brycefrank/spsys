#library(rgdal)
#library(raster)
#library(sp)
#
#northsouth <- readOGR('D:\\test_grids\\northsouth.shp')
#utm <- CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs')
#
#northsouth <- spTransform(northsouth, utm)
#
#coords <- data.frame(northsouth@coords)
#
#mod <- lm(coords.x2 ~ coords.x1, data=coords)
#
#x_int <- -coefficients(mod)[[1]] / coefficients(mod)[[2]]
#x_int
#y_int <- coefficients(mod)[[1]]
#
#theta <- (atan(y_int/x_int) / (2*pi)) * 360
#