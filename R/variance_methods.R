# Variance estimators in method form
library(spsurvey)

setGeneric('var_srs', function(sys_frame, ...) {
  standardGeneric('var_srs')
})


# TODO could be possible to attach N in the case of subsamples
setMethod('var_srs', signature(sys_frame='SysFrame'),
  function(sys_frame, fpc=FALSE, N=NA_real_) {
    if(fpc == TRUE & is.na(N)) {
      stop('If fpc is set to true you must provide a population size N.')
    }
    
    att_df <- sys_frame@data[, sys_frame@attributes, drop=FALSE]
    n <- length(sys_frame)
    var_mu <- (n - 1)^(-1) * n^(-1) *  colSums((att_df - colMeans(att_df))^2)
    
    if(fpc) {
      var_mu <- var_mu * (1-n/N)
    }
    return(var_mu)
  }
)


#' Computes the variance estimate of a systematic 
#' sample using the Stevens and Olsen (2003) local mean estimator.
#' 
#' This is a wrapper for functions that exist in `spsurvey`
setGeneric('var_so', function(sys_frame, ...){
  standardGeneric('var_so')
})

setMethod('var_so', signature(sys_frame='SysFrame'),
  function(sys_frame, fpc=FALSE, N=NA_real_) {
    # TODO this should not strictly be needed
    if(is.na(N)) {
      stop('var_so requires a population size N')
    }
    
    n <- nrow(sys_frame)
    
    if(!'pi_i' %in% colnames(sys_frame@data)) {
      warning('Inclusion probability column - "pi_i" not found in @data, assuming
               inclusion probabilities are n/N.')
      pi_i <- rep(n/N, n)
    } else {
      pi_i <- sys_frame@data$pi_i
    }
    
    wt <- localmean.weight(sys_frame@coords[,1], sys_frame@coords[,2], pi_i)
    att_df <- sys_frame@data[, sys_frame@attributes, drop=FALSE]
    var_total <- sapply(colnames(att_df), so_att, att_df, pi_i, wt)
    
    if(fpc) {
      var_mu <- (1/N^2) * var_total * (1 - n/N)
    } else {
      var_mu <- (1/N^2) * var_total
    }
    return(var_mu)
  }
)

#' Computes the Stevens and Olsen estimator for a specific attribute. 
#' Used internally with var_so function.
so_att <- function(att, att_df, pi_i, wt) {
  z <- att_df[,att]
  localmean.var(z/pi_i, wt)
}