plots <- read.csv('D:\\plot040820.csv')
plots <- plots[complete.cases(plots[,c('LON_ACTUAL_NAD83', 'LAT_ACTUAL_NAD83')]),]

plots <- SpatialPointsDataFrame(plots[,c('LON_ACTUAL_NAD83', 'LAT_ACTUAL_NAD83')], data=plots)
crs(plots) <- CRS('+proj=longlat +datum=WGS84 +no_defs')
utm <- CRS('+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs')
plots <- spTransform(plots, utm)

find_hex <- function(M) {
  prev_loss <- 10000
  for (i in 1:M) {
    bbox <- bbox(plots)
    min_x <- bbox[1,1] - 20000 - runif(1, 0, 6000)
    max_x <- bbox[1,2] + 20000
    width <- max_x - min_x
    min_y <- bbox[2,1] - 20000 - runif(1, 0, 6000)
    max_y <- min_y + width + 20000
    #max_x <- min_x + 100000
    #max_y <- min_y + 100000
    
    bbox <- data.frame(x=c(min_x, min_x, max_x, max_x), y=c(min_y, max_y, max_y, min_y))
    bbox <- Polygons(list(Polygon(bbox)), 1)
    bbox <- SpatialPolygons(list(bbox))
    
    rand_angle <- runif(1, 0, 90)
    side_to_side <- 5187 # Length of hexagon from flat-flat side (avg distance between plots)
    HexPts <- spsample(bbox, type="hexagonal", cellsize=side_to_side)
    
    # Transform coords such that the origin is the center
    mid <- c(min_x + (width/2), min_y + (width/2))
    HexPols <- HexPoints2SpatialPolygons(HexPts)
    HexPols <- elide(HexPols, rotate=90-theta, center=mid)
    
    HexPols_df <- data.frame(x=rep(0, length(HexPols)))
    crs(HexPols) <- newcrs
    plot(HexPols)
    
    overlay <- data.frame(hex_id = over(plots, HexPols))
    
    ov <- overlay %>%
      group_by(hex_id) %>%
      summarize(n=n()) %>%
      group_by(n) %>%
      summarize(n1=n())
    
    loss <- ov[[2,2]]
    if(loss < prev_loss) {
      print(i)
      print(loss)
      prev_loss <- loss
      best <- HexPols
      plot(best)
      plot(plots, add=TRUE)
    }
    
  }
  return(best)
}

grd <- find_hex(1000)

HexPols_df <- data.frame(id=seq(1, length(HexPols)))
row.names(HexPols_df) <- row.names(HexPols)
HexPols_spdf <- SpatialPolygonsDataFrame(HexPols, data=HexPols_df)
writeOGR(HexPols_spdf, dsn='D:\\test_grids', layer='test_layer.shp', driver='ESRI Shapefile')
