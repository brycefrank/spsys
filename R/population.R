library(gstat)
library(sp)
library(reshape2)
library(data.table)

#' @param w width in pixels
#' @param h height in pixels
make_pop <- function(w, h) {
  xy <- expand.grid(1:w, 1:h)
  names(xy) <- c('x', 'y')
  g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, model=vgm(psill=0.025, range=5, model='Exp'), nmax=20)
  yy <- predict(g.dummy, newdata=xy, nsim=1)
  gridded(yy) = ~x+y
  
  pop_mat <- matrix(yy@data[,1], nrow=h, ncol=w)
  return(pop_mat)
}

#' @param a The length of a side of the equilateral triangle. Must be an
#' integer prepresenting the number of pixels (does not have to be)
samp_tri_grid <- function(pop_mat, pop_df, a) {
  colnames(pop_df) <- c('x', 'y', 'z')
  
  sx <- runif(1, 1, 1 + (sqrt(3)/2) * a)
  sy <- runif(1, 1, 1 + 2*a)
  
  samp_x <- seq(sx, 300, (sqrt(3)/2) * a)
  samp_y <- seq(sy, 300, 2*a)
  
  grid1 <- expand.grid(samp_x, samp_y)
  grid2 <- grid1
  
  grid2[,2] <- grid2[,2] + (sqrt(3)/2)*a
  grid2[,1] <- grid2[,1] + a/2 
  
  grid <- rbind(grid1, grid2)
  grid <- round(grid[(grid[,1] <= 300) & (grid[,2] <= 300),])
  names(grid) <- c('x', 'y')
  grid <- data.table(grid)
  
  samp <- merge(grid, pop_df, all.x=TRUE, all.y=FALSE, by=c('x', 'y'))
  names(samp) <- c('x', 'y', 'z')
  
  return(samp)
}


mean_pop <- mean(pop_mat)
pop_df <- melt(pop_mat)
pop_df <- data.table(pop_df)
colnames(pop_df) <- c('x', 'y', 'z')

sys_vars <- c()
for(i in 1:300) {
  samp <- samp_tri_grid(pop_mat, pop_df, 70)
  plot(samp$x, samp$y)
  z_bar <- mean(samp$z)
  print(z_bar - mean_pop)
}
