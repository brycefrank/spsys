#' Base class for all variance outputs
setClass('VarOut', slots = list(
  estimate = 'numeric',
  n = 'numeric',
  N = 'numeric',
  mu = 'numeric'
))

VarOut <- function(estimate, n, N, mu, diagnostic) {
  if(diagnostic) {
    var_out <- new('VarOut')
    var_out@estimate <- estimate
    var_out@n <- n
    var_out@N <- N
    var_out@mu <- mu
    return(var_out)
  } else {
    return(estimate)
  }
}

setGeneric('summary', function(var_out){
  standardGeneric('summary')
})

setMethod('summary', signature(var_out = 'VarOut'), 
  function(var_out) {
    cat(var_out@name, 'summary:\n\n')
    
    # TODO what to do for large number of variables?
    display_est <- t(var_out@estimate)
    colnames(display_est) <- c('Var. Estimate')
    print(display_est)
    
    if(!is.na(var_out@N)) {
      cat('\n\nn:', var_out@n, '  N:', var_out@N)
    } else {
      cat('\n\nn:', var_out@n)
    }
  }
)