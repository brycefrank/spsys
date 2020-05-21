devtools::load_all()
block_hex <- readRDS('data/block_hex.RDS')

plot(subsample_hex(block_hex, c(1,1), 10))
mu <- mean(block_hex@data$z_1)

estimators <- list(
  'var_srs' = var_srs,
  'var_dorazio' = var_dorazio
)

compare_estimators <- function(a_vec, pop, estimators) {
  sys_vars <- data.frame(matrix(0, nrow=length(a_vec), ncol=2))
  colnames(sys_vars) <- c('a', 'v_sys')
  N <- nrow(pop)
  est_vars <- list()
  
  for(i in 1:length(a_vec)) {
    a <- a_vec[[i]]
    means <- c()
    
    starts <- subsample_starts(a)
    est_vars_i <- matrix(0, nrow=nrow(starts), ncol=length(estimators))
    
    for(j in 1:nrow(starts)) {
      print(j)
      start <- starts[j,]
      subsamp_ix <- subsample_hex(pop, start, a)
      subsamp <- sp::merge(pop, subsamp_ix, by=c('r', 'c'), all.x=FALSE, all.y=TRUE)
      n <- nrow(subsamp)
      means <- c(means, mean(subsamp$z_1))
      
      #Convert subsamp to a dataframe compatible with the variance estimators
      subsamp_df <- subsamp@data
      subsamp_centroids <- gCentroid(subsamp, byid=TRUE)@coords
      
      # TODO this is just temporary, we need actual plot coordinates here
      subsamp_df$x <- subsamp_centroids[,1]
      subsamp_df$y <- subsamp_centroids[,2]
      subsamp_df$z <- subsamp_df$VOL
      subsamp_df$pi_i <- n/N
      
      
      for(k in 1:length(estimators)) {
        v <- estimators[[k]]
        
        if(names(estimators)[[k]] == 'var_fpbk')  {
          subsamp_df <- merge(pop[,c('plt_ix'), drop=FALSE], subsamp)
        } else if (names(estimators)[[k]] %in% c('var_mat_hex', 'var_nnbh_hex')) {
          est_vars_i[j,k] <- v(subsamp@data, a)
        } else{
          est_vars_i[j,k] <- v(subsamp_df, N)
        }
      }
    }
    
    est_vars_i <- data.frame(est_vars_i)
    colnames(est_vars_i) <- names(estimators)
    est_vars_i$a <- a
    est_vars[[i]] <- est_vars_i
    sys_vars[i, 'a'] <- a
    sys_vars[i, 'v_sys'] <- mean((means - mu)^2)
    
  }
  
  est_vars <- bind_rows(est_vars)
  return(list(sys_vars, est_vars))
  
}

a_vec   <- c(4,5,6,7,8,9,10)
results <- compare_estimators(a_vec, block_hex, estimators)

ggplot() +
  geom_point(data=results[[2]], aes(x=a, y=var_dorazio), color='red') +
  geom_point(data=results[[2]], aes(x=a, y=var_srs), color='blue') +
  geom_line(data=results[[1]], aes(x=a, y=v_sys), color='black')
