# Variance estimators in method form


setGeneric('var_srs', function(sys_frame, ...) {
  standardGeneric('var_srs')
})


# TODO could be possible to attach N in the case of subsamples
setMethod('var_srs', signature(sys_frame='SysFrame'),
  function(sys_frame, fpc=TRUE, N=NA_real_) {
    att_df <- sys_frame@data[, sys_frame@attributes, drop=FALSE]
    n <- length(sys_frame)
    var_mu <- (n - 1)^(-1) * n^(-1) *  colSums((att_df - colMeans(att_df))^2)
    
    if(fpc) {
      var_mu <- var_mu * (1-n/N)
    }
    return(var_mu)
  }
)