library(ggplot2)
devtools::load_all()
block_hex <- readRDS('data/hex_pts.RDS')

hex_frame <- HexFrame(block_hex, attributes=c('z_1'))
p <- length(hex_frame@attributes)
mu <- colMeans(hex_frame@data[,hex_frame@attributes, drop=FALSE])

estimators <- list(
  'var_srs' = var_srs,
  'var_dorazio_i' = var_dorazio_i,
  'var_dorazio_c' = var_dorazio_c,
  'var_non_overlap'= var_non_overlap
)

compare_estimators <- function(a_vec, pop, estimators) {
  sys_vars <- data.frame(matrix(0, nrow=length(a_vec), ncol=p+1))
  colnames(sys_vars) <- c('a', 'v_sys')
  N <- nrow(pop)
  h <- 1
  var_sys_df <- list()
  est_vars <- list()
  
  for(i in 1:length(a_vec)) {
    a <- a_vec[[i]]
    print(paste('Processing sampling interval a =', a))
    
    starts <- subsample_starts(a)
    means <- data.frame(matrix(0, nrow=nrow(starts), ncol=p))
    colnames(means) <- pop@attributes
    
    
    for(j in 1:nrow(starts)) {
      start <- starts[j,]
      subsamp <- subsample(pop, start, a)
      n <- nrow(subsamp)
      means[j,] <- colMeans(subsamp@data[,subsamp@attributes, drop=FALSE])
      
      for(k in 1:length(estimators)) {
        est_vars_k <- data.frame(matrix(0, nrow=p, ncol=5))
        est_vars_k[,1] <- names(estimators)[[k]]
        est_vars_k[,2] <- a
        est_vars_k[,3] <- pop@attributes
        est_vars_k[,4] <- j
        
        v <- estimators[[k]]
        v_name <- names(estimators)[[k]]
        
        # Some variance estimators need particular arguments
        est_vars_k[,5] <- v(subsamp, fpc=TRUE, N=N)
        est_vars[[h]] <- data.frame(est_vars_k)
        h <- h + 1
      }
    }
    # Compute and store the systematic variance
    var_sys <- colMeans(sweep(means, 2, mu)^2)
    var_sys_i <- data.frame(matrix(0), nrow=p, ncol=5)
    var_sys_i[,1] <- 'var_sys'
    var_sys_i[,2] <- a
    var_sys_i[,3] <- pop@attributes
    var_sys_i[,4] <- 1
    var_sys_i[,5] <- var_sys
    var_sys_df[[i]] <- var_sys_i
  }
  
  est_vars <- bind_rows(est_vars)
  sys_vars <- bind_rows(var_sys_df)
  colnames(sys_vars) <- c('estimator', 'a', 'attribute', 's', 'variance')
  colnames(est_vars) <- c('estimator', 'a', 'attribute', 's', 'variance')
  rbind(sys_vars, est_vars)
}

a_vec   <- c(4:10)
results <- compare_estimators(a_vec, hex_frame, estimators)

ggplot(results) +
  geom_jitter(aes(x=as.numeric(a), variance, color=estimator)) +
  geom_line(aes(x=as.numeric(a),  variance), data=filter(results, estimator=='var_sys'))
  ylab("")

