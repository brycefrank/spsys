# Functions that make the simulated data sets
library(sp)
library(raster)
library(dplyr)
library(gstat)
library(ggplot2)
library(tidyr)
library(colorspace)
library(latex2exp)
set.seed(12)

make_hex_pops <- function() {
  min_x <- 0
  max_x <- 150
  min_y <- 0
  max_y <- 150
  bbox <- data.frame(x=c(min_x, min_x, max_x, max_x), y=c(min_y, max_y, max_y, min_y))
  bbox <- Polygons(list(Polygon(bbox)), 1)
  bbox <- SpatialPolygons(list(bbox))
  
  HexPts <- spsample(bbox, type="hexagonal", cellsize=5)
  
  # Now we want to simulate spatially correlated data
  coords <- HexPts@coords
  colnames(coords) <- c('x', 'y')
  coords <- data.frame(coords)
  
  slm1 <-  gstat(formula=z~1, locations=~x+y, model=vgm(psill=2, range=1, model='Exp', nugget = 0.5), dummy=T, beta=1, nmax=1)
  slm50 <-  gstat(formula=z~1, locations=~x+y, model=vgm(psill=2, range=50, model='Exp', nugget = 0.5), dummy=T, beta=1, nmax=1)
  slm100 <- gstat(formula=z~1, locations=~x+y, model=vgm(psill=2, range=100, model='Exp', nugget = 0.5), dummy=T, beta=1, nmax=1)
  
  pop1 <-  data.frame(predict(slm1, newdata=coords, nsim=1))
  pop50 <-  data.frame(predict(slm50, newdata=coords, nsim=1))
  pop100 <- data.frame(predict(slm100, newdata=coords, nsim=1))
  
  pop_df <- cbind(pop1, pop50[,c('sim1')], pop100[,c('sim1')])
  colnames(pop_df) <- c('s_1', 's_2', 'z_1', 'z_50', 'z_100')
  
  pop_df
}

hex_pts <- make_hex_pops()
usethis::use_data(hex_pts, overwrite = TRUE)

#' Makes a small dataset for testing the HexFrame functionality
make_hex_small <- function() {
  hex_pts_small <- hex_pts[hex_pts$s_1 < 20 & hex_pts$s_2 > 130,]
  return(hex_pts_small)
}

hex_pts_small <- make_hex_small()
usethis::use_data(hex_pts_small, overwrite = TRUE)

make_rect_pops <- function() {
  min_x <- 0
  max_x <- 150
  min_y <- 0
  max_y <- 150
  bbox <- data.frame(x=c(min_x, min_x, max_x, max_x), y=c(min_y, max_y, max_y, min_y))
  bbox <- Polygons(list(Polygon(bbox)), 1)
  bbox <- SpatialPolygons(list(bbox))
  
  RectPts <- spsample(bbox, type="regular", cellsize=5)
  
  # Now we want to simulate spatially correlated data
  coords <- RectPts@coords
  colnames(coords) <- c('x', 'y')
  coords <- data.frame(coords)
  
  slm1   <- gstat(formula=z~1, locations=~x+y, model=vgm(psill=2, range=1, model='Exp', nugget = 0.5), dummy=T, beta=1, nmax=1)
  slm50  <- gstat(formula=z~1, locations=~x+y, model=vgm(psill=2, range=50, model='Exp', nugget = 0.5), dummy=T, beta=1, nmax=1)
  slm100 <- gstat(formula=z~1, locations=~x+y, model=vgm(psill=2, range=100, model='Exp', nugget = 0.5), dummy=T, beta=1, nmax=1)
  
  pop1   <- data.frame(predict(slm1, newdata=coords, nsim=1))
  pop50  <- data.frame(predict(slm50, newdata=coords, nsim=1))
  pop100 <- data.frame(predict(slm100, newdata=coords, nsim=1))
  
  pop_df <- cbind(pop1, pop50[,c('sim1')], pop100[,c('sim1')])
  colnames(pop_df) <- c('s_1', 's_2', 'z_1', 'z_50', 'z_100')
  
  pop_df
}

rect_pts <- make_rect_pops()
usethis::use_data(rect_pts, overwrite = TRUE)

make_hex_fig <- function(pops) {
  plot_df <- pivot_longer(pops, c('z_1', 'z_50', 'z_100'))
  plot_df$name <- factor(plot_df$name, levels=c('z_1', 'z_50', 'z_100'))
  levels(plot_df$name) <- c(TeX('$\\phi = 1$'), TeX('$\\phi = 50$'), TeX('$\\phi = 100$'))
  
  ggplot(plot_df) +
    geom_hex(aes(x=s_1, y=s_2, color=value, fill=value), stat='identity') + 
    scale_fill_continuous_sequential(palette = "BluYl", name=TeX('z(s)')) +
    scale_color_continuous_sequential(palette = "BluYl", name=TeX('z(s)')) +
    xlab(unname(TeX('$s_1$'))) +
    ylab(unname(TeX('$s_2$'))) +
    facet_wrap(~name, labeller=label_parsed) +
    theme(strip.text.x = element_text(size = 16), legend.position='bottom')
}

#' Makes a small dataset for testing the RectFrame functionality
make_rect_small <- function() {
  rect_pts_small <- rect_pts[rect_pts$s_1 < 20 & rect_pts$s_2 > 130,]
  return(rect_pts_small)
}

rect_pts_small <- make_rect_small()
usethis::use_data(rect_pts_small, overwrite = TRUE)