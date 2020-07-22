#' Using subsampling, assess a series of variance estimators.
#'
#'@param sys_frame
#'@param a_vec
#'@param estimators
compare_estimators <- function(sys_frame, a_vec, estimators) {
  N <- nrow(sys_frame)
  mu <- colMeans(sys_frame@data[,sys_frame@attributes, drop=FALSE])
  sigma2 <- (1/N) * colSums(sweep(sys_frame@data[,sys_frame@attributes], 2, mu)^2)
  gearys_C <- gearys_c(sys_frame)
  morans_I <- morans_i(sys_frame)
  
  p <- length(sys_frame@attributes)
  sys_vars <- data.frame(matrix(0, nrow=length(a_vec), ncol=p+1))
  colnames(sys_vars) <- c('a', 'v_sys')
  h <- 1
  var_sys_df <- list()
  est_vars <- list()
  
  for(i in 1:length(a_vec)) {
    a <- a_vec[[i]]
    print(paste('Processing sampling interval a =', a))
    
    starts <- subsample_starts(sys_frame, a)
    means <- data.frame(matrix(0, nrow=nrow(starts), ncol=p))
    colnames(means) <- sys_frame@attributes
    
    for(j in 1:nrow(starts)) {
      print(paste('---- Processing sampling position =', j))
      start <- starts[j,]
      subsamp <- subsample(sys_frame, start, a)
      subsamp@data$pi_i <- 1/a^2
      subsamp@N <- N
      n <- nrow(subsamp)
      means[j,] <- colMeans(subsamp@data[,subsamp@attributes, drop=FALSE])
      
      for(k in 1:length(estimators)) {
        est_vars_k <- data.frame(matrix(0, nrow=p, ncol=7))
        est_vars_k[,1] <- names(estimators)[[k]]
        est_vars_k[,2] <- a
        est_vars_k[,3] <- sys_frame@attributes
        est_vars_k[,4] <- j
        est_vars_k[,6] <- n
        est_vars_k[,7] <- t(means[j,])
        
        v <- estimators[[k]]
        v_name <- names(estimators)[[k]]
        var_est <- v(subsamp)
        
        est_vars_k[,5] <- unlist(var_est)
        est_vars[[h]] <- data.frame(est_vars_k)
        h <- h + 1
      }
    }
    # Compute and store the systematic variance
    var_sys <- colMeans(sweep(means, 2, mu)^2)
    var_sys_i <- data.frame(matrix(0, nrow=p, ncol=4))
    var_sys_i[,1] <- 'var_sys'
    var_sys_i[,2] <- a
    var_sys_i[,3] <- sys_frame@attributes
    var_sys_i[,4] <- var_sys
    var_sys_df[[i]] <- var_sys_i
  }
  
  est_vars <- bind_rows(est_vars)
  sys_vars <- bind_rows(var_sys_df)
  colnames(sys_vars) <- c('estimator', 'a', 'attribute', 'variance')
  colnames(est_vars) <- c('estimator', 'a', 'attribute', 's', 'variance', 'n', 'mu_hat')
  pop_frame <- data.frame(N=N, mu=mu, sigma2=sigma2, attribute=sys_frame@attributes, gearys_C = gearys_C, morans_I = morans_I)
  return(list('est_frame' = est_vars, 'pop_frame' = pop_frame, 'sys_frame' = sys_vars))
}