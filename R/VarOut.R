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
