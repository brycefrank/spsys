library(gstat)
library(sp)
library(reshape2)
library(data.table)
library(dplyr)

#' @param w width in pixels
#' @param h height in pixels
make_pop <- function(w, h) {
  xy <- expand.grid(1:w, 1:h)
  names(xy) <- c('x', 'y')
  g.dummy <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=0, model=vgm(psill=0.025, range=5, model='Exp'), nmax=20)
  yy <- predict(g.dummy, newdata=xy, nsim=1)
  gridded(yy) = ~x+y
  
  pop_mat <- matrix(yy@data[,1], nrow=h, ncol=w)
  return(pop_mat)
}

#' @param a The length of a side of the equilateral triangle. Must be an
#' integer prepresenting the number of pixels (does not have to be)
samp_tri_grid <- function(pop_mat, pop_df, a) {
  colnames(pop_df) <- c('x', 'y', 'z')
  h <- dim(pop_mat)[[1]]
  w <- dim(pop_mat)[[2]]
  
  # Select a random starting position
  # TODO this induces a random nature to n...how big of a deal
  # is this?
  sx <- runif(1, 1, 1 + (sqrt(3)/2) * a)
  sy <- runif(1, 1, 1 + 2*a)
  
  samp_x <- seq(sx, 300, (sqrt(3)/2) * a)
  samp_y <- seq(sy, 300, 2*a)
  
  grid1 <- expand.grid(samp_x, samp_y)
  grid2 <- grid1
  
  grid2[,2] <- grid2[,2] + (sqrt(3)/2)*a
  grid2[,1] <- grid2[,1] + a/2 
  
  grid <- rbind(grid1, grid2)
  grid <- round(grid[(grid[,1] <= w) & (grid[,2] <= h),])
  names(grid) <- c('x', 'y')
  grid <- data.table(grid)
  
  samp <- merge(grid, pop_df, all.x=TRUE, all.y=FALSE, by=c('x', 'y'))
  names(samp) <- c('x', 'y', 'z')
  
  return(samp)
}

samp_tri_grid_pi <- function(pop_mat, pop_df, strat_spec) {
  strata_grids <- list()
  i <- 1
  
  strat_mat  <- strat_spec[[1]]
  strata_ids <- strat_spec[[2]]
  
  for(id in names(strata_ids)) {
    a <- strata_ids[[id]]
    a_ix <- data.table(which(strat_mat == as.numeric(id), arr.ind=TRUE))
    colnames(a_ix) <- c('y', 'x')
    
    full_grid <- samp_tri_grid(pop_mat, pop_df, a)
    sub_grid  <- merge(full_grid, a_ix, all.x=FALSE, all.y=FALSE)
    
    # TODO the pi_a are conditioned on the sample, should be fixed b/c n is random and a function
    # of the starting position. One option is to do a large number of simulations, since it may be dependent
    # on very complex geomeries.
    sub_grid$n_i      <- nrow(sub_grid)
    sub_grid$N_i      <- nrow(a_ix)
    sub_grid$pi_i     <- sub_grid$n_i / sub_grid$N_i
    sub_grid$strat_id <- id
    strata_grids[[i]] <- sub_grid
    
    i <- i + 1
  }
  
  return(bind_rows(strata_grids))
  
}
