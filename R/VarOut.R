#' Base class for all variance outputs
setClass('VarOut', slots = list(
  estimate = 'data.frame',
  n = 'numeric',
  N = 'numeric'
))

VarOut <- function(estimate, n, N, diagnostic) {
  if(diagnostic) {
    var_out <- new('VarOut')
    var_out@estimate <- estimate
    var_out@n <- n
    var_out@N <- N
    return(var_out)
  } else {
    return(estimate)
  }
}

