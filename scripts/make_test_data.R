library(sp)
library(raster)
library(devtools)
library(rgdal)
library(sp)
devtools::load_all()
set.seed(12)

# A demonstration of variance estimation for a hexagonal grid

## Make a synthetic dataset ##
# TODO replace with real data?
plots_o <- read.csv('data/old/plots_safe.csv')
plots <- SpatialPoints(plots_o[,c('LON', 'LAT')])
crs(plots) <- CRS("+proj=longlat +datum=WGS84 +no_defs")

newcrs <- CRS("+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
plots <- SpatialPointsDataFrame(spTransform(plots, newcrs), data = plots_o)

# Fake data just to get the next part running
plots@data$z_1 <- rnorm(length(plots))
block <- cast_hex_grid(plots)

index_hex(block)
block@data <- cbind(block@data, index_hex(block))

n_missing <- sum(is.na(block@data$z_1))
block@data[is.na(block@data$z_1), 'z_1'] <- rnorm(n_missing)
saveRDS(block, file="data/block_hex.RDS")

or_hex <- block[complete.cases(block@data),]
saveRDS(or_hex, file="data/or_hex.RDS")

hex_pts <- readRDS('data/hex_pts.RDS')

d <- 7000

min_x <- bbox(hex_pts)[1,1]
min_y <- bbox(hex_pts)[2,1]
max_x <- bbox(hex_pts)[1,2]
max_y <- bbox(hex_pts)[2,2]

x <- seq(min_x, max_x, by=d)
y <- seq(min_y, max_y, by=d)

coords <- expand.grid(x,y)

coords$z_1 <- rnorm(nrow(coords))

rect_pts <- SpatialPointsDataFrame(coords[,c('Var1', 'Var2')], data=coords[,c('z_1'),drop=FALSE])
saveRDS(rect_pts, 'data/rect_pts.RDS')
