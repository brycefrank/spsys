devtools::load_all()
library(ggplot2)

pop_df <- readRDS('data/pop_df.RDS')
pop_mat <- readRDS('data/pop_mat.RDS')
strat_mat <- matrix(0, 100, 100)

# Specifying ID of strata in space...
strat_mat[1:30, 1:50] <- 1
strat_mat[strat_mat == 0] <- 2

# Linking strata IDs to distances btwn samples
a <- list('1'=10, '2'=10)

# A strat spec is 2-element list containing the 
# strat_mat and the a list
strat_spec <- list(strat_mat, a)

var_sys <- var_sys_monte_carlo(pop_mat, pop_df, strat_spec, 1000)
simulate <- function(n_sim, pop_mat, pop_df, strat_spec) {
  N <- nrow(pop_df)
  
  # Specify the variance functions and their arguments
  var_funcs <- list(
    'var_strs' = list(func = var_strs, data=c()),
    'var_so'   = list(func = var_so, data=c())
  )
  
  for(i in 1:n_sim) {
    samp <- samp_tri_grid_pi(pop_mat, pop_df, strat_spec)
    samp <- data.frame(samp) # TODO would be nice to not have to do  this, some problem with data.table subsetting...
    
    for (key in names(var_funcs)) {
      func <- var_funcs[[key]][['func']]
      var_ <- func(samp, N)
      var_funcs[[key]][['data']] <- c(var_funcs[[key]][['data']], var_ / var_sys)
    }
  }
  
  sims <- data.frame(matrix(NA, nrow=n_sim, ncol=length(var_funcs)))
  colnames(sims) <- names(var_funcs)
  
  for(key in names(var_funcs)) {
    sims[,key] <- var_funcs[[key]]['data']
  }
  
  return(sims)
  
}

sim  <- simulate(1000, pop_mat, pop_df, strat_spec)
ggplot(sim) +
  geom_point(aes(x=var_strs, y=var_so))
