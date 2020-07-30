#' Base class for all variance outputs
setClass('VarOut', slots = list(
  estimate = 'numeric',
  n = 'numeric',
  N = 'numeric',
  mu = 'numeric'
))

#' The base class for variance outputs. All estimators will contain
#' at least this information
#' 
#' @param estimate A named vector of variance estimates
#' @param n The sample size
#' @param N The population size
#' @param mu A named vector of mean estimates
#' @param diagnostic If TRUE then return the diagnostic information, if FALSE
#' return only the variance estimates.
#' @keywords internal
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
