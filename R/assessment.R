
wt_mean <- function(x, W_g) {return(x*W_g)}

#' Compare variance estimators via simulation
#' 
#' This function accepts a `SysFrame` that is treated as a population. All possible
#' systematic samples for a given sampling interval are made and all estimators specified
#' in `estimators` are used to estimate the variance for that sample. The variance estimates
#' for each sample are returned along with other information, e.g. the population means, variances
#' etc.
#' 
#' The outputs of this function can be used to create mean squared errors, biases and other
#' assessments of the estimators.
#' 
#' @param sys_frame A population from which to sample
#' @param a_vec A vector of integers representing the sampling intervals desired for assessment
#' @param estimators A named list of constructed estimator functions
#' @param covt_cols A character vector indicating columns of post-strata variables in the `@data` slot of `sys_frame`.
#' @return A named list containing three dataframes
#' \itemize{
#' \item `'est_frame'` - A dataframe of the variance estimates
#' \item `'pop_frame'` - A dataframe of population parameters
#' \item `'sys_frame'` - A dataframe of true systematic variances
#' }
#' @example 
#' estimators <- list(
#'    'var_srs' = VarSRS(),
#'    'var_so' = VarS()
#' )
#' compare_estimators(hex_frame, c(4), estimators)
#' @export
compare_estimators <- function(sys_frame, a_vec, estimators, covt_cols=NA) {
  # TODO this could be moved to SysFrame declaration (covt_cols argument, for example)
  
  if(is.na(covt_cols)) {
    covt_cols <- c('intercept')
    sys_frame@data$intercept <- 1
  }
  
  mu <- sys_frame@data %>%
    group_by_at(covt_cols) %>%
    summarize_at(.vars=sys_frame@attributes, .funs=~mean(.))
  
  N_g <- sys_frame@data %>%
    group_by_at(covt_cols) %>%
    summarize(N_g = n())
  
  N <- nrow(sys_frame)
  
  p <- length(sys_frame@attributes)
  sys_vars <- data.frame(matrix(0, nrow=length(a_vec), ncol=p+1))
  colnames(sys_vars) <- c('a', 'v_sys')
  h <- 1
  var_sys_df <- list()
  est_vars <- list()
  pop_resids <- sys_frame@data[,c(covt_cols, sys_frame@attributes, 'r', 'c')] %>%
    merge(mu, by=covt_cols, suffixes=c('_s', '_mu_hat'))
  pop_r <- pop_resids %>% 
    dplyr::select(dplyr::ends_with('_s'))
  pop_group_means <- pop_resids %>%
    dplyr::select(dplyr::ends_with('_mu_hat'))
  
  pop_resids <- cbind(pop_resids[,c(covt_cols, 'r', 'c')], pop_r - pop_group_means)
  colnames(pop_resids) <- c(covt_cols, 'r', 'c', sys_frame@attributes)
  resid_frame <- sys_frame
  resid_frame@data <- merge(resid_frame@data[,!(colnames(resid_frame@data) %in% c(covt_cols, resid_frame@attributes))], pop_resids, by=c('r', 'c'))
  gc <- gearys_c(resid_frame)
  mi <- morans_i(resid_frame)
  
  weights <- pop_resids %>%
    group_by_at(covt_cols) %>%
    summarize(w_h =  n() / nrow(pop_resids))
  
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
      
      # Calculate group-specific means
      mu_hat <- subsamp@data[,c(covt_cols, sys_frame@attributes)] %>%
        group_by_at(covt_cols) %>%
        summarize_at(.vars=sys_frame@attributes, .funs=~mean(.)) %>%
        merge(mu, by=covt_cols, all.y=TRUE) %>%
        dplyr::select(covt_cols, dplyr::ends_with('.x'))
        
      colnames(mu_hat) <- c(covt_cols, sys_frame@attributes)
      
      # Replace NAns of non-present groups with zero
      zero_replace <- as.list(rep(0, p))
      names(zero_replace) <- sys_frame@attributes
      mu_hat <- mu_hat %>%
        replace_na(zero_replace)
      
      # Compute the group-level residuals
      samp_resids <- subsamp@data[,c(covt_cols, sys_frame@attributes, 'r', 'c')] %>%
        merge(mu_hat, by=covt_cols, suffixes=c('_s', '_mu_hat'))
      
      r <- samp_resids %>% 
        dplyr::select(dplyr::ends_with('_s'))
      group_means <- samp_resids %>%
        dplyr::select(dplyr::ends_with('_mu_hat'))
      
      samp_resids <- cbind(samp_resids[,c(covt_cols, 'r', 'c')], r - group_means)
      colnames(samp_resids) <- c(covt_cols, 'r', 'c', sys_frame@attributes)
      
      # TODO Calculate the extra variance term
      extra_var <- samp_resids %>%
        group_by_at(covt_cols) %>%
        summarize_at(.vars=sys_frame@attributes, .funs=~var(.)) %>%
        merge(weights, by=covt_cols) %>%
        replace_na(zero_replace) %>%
        mutate_at(.vars=sys_frame@attributes, .funs=~wt_mean(., 1 - w_h)) %>%
        summarize_at(.vars=sys_frame@attributes, .funs=~sum(. / n^2))
      
      # Now over-write with the residuals instead of the original data
      #subsamp@data <- cbind(samp_resids, subsamp@data[,c('r', 'c', 'pi_i')])
      subsamp@data <- merge(subsamp@data[,!(colnames(subsamp@data) %in% c(covt_cols, sys_frame@attributes))], samp_resids, by=c('r', 'c'))
      # Calculate the mean estimates for each attribute using the strata
      mu_hat_ps <- merge(mu_hat, N_g, by=covt_cols) %>%
        mutate(W_g = N_g/N) %>%
        mutate_at(.vars = sys_frame@attributes, .funs = ~wt_mean(., W_g)) %>%
        summarize_at(.vars = sys_frame@attributes, .funs = ~sum(.))
      
      means[j,] <- mu_hat_ps
        
      for(k in 1:length(estimators)) {
        # TODO there are two ways to do this:
        # 1. This is the current implementation: just send the residuals as the data
        # 2. Separately estimate the variances in each strata, V. discouraged this because of the way the neighborhoods are made
        #    we need access to the residual even if they are in separate strata.
        est_vars_k <- data.frame(matrix(0, nrow=p, ncol=7))
        est_vars_k[,1] <- names(estimators)[[k]]
        est_vars_k[,2] <- a
        est_vars_k[,3] <- sys_frame@attributes
        est_vars_k[,4] <- j
        est_vars_k[,6] <- n
        est_vars_k[,7] <- t(means[j,])
        
        v <- estimators[[k]]
        v_name <- names(estimators)[[k]]
        var_est <- v(subsamp) + extra_var
        
        
        est_vars_k[,5] <- unlist(var_est)
        est_vars[[h]] <- data.frame(est_vars_k)
        h <- h + 1
      }
    }
    # Compute and store the systematic variance
    mu  <- merge(mu, N_g, by=covt_cols) %>%
      mutate(W_g = N_g/N) %>%
      mutate_at(.vars = sys_frame@attributes, .funs = ~wt_mean(., W_g)) %>%
      summarize_at(.vars = sys_frame@attributes, .funs = ~sum(.))
      
    # TODO replace with the actual post-stratified variance here!
    var_sys <- colMeans(sweep(as.matrix(means), 2, as.numeric(mu))^2)
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
  pop_frame <- data.frame(N=N, mu=t(mu), attribute=sys_frame@attributes, gearys_C = gc, morans_I = mi)
  return(list('est_frame' = est_vars, 'pop_frame' = pop_frame, 'sys_frame' = sys_vars))
}