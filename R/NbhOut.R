#' @include VarOut.R

setClass('NbhOut', contains='VarOut', slots=list(
  neighborhoods='data.frame',
  n='numeric',
  N='numeric',
  mu='numeric',
  n_groups='numeric',
  name='character'
))

#' Contains the diagnostic output of neighborhood-based estimators
#' 
#' @param estimate A named vector of variance estimates
#' @param neighborhoods A dataframe of neighborhood groupings and values, estimator dependent
#' @param n The sample size
#' @param N The population size
#' @param mu A named vector of mean estimates
#' @param diagnostic If TRUE then return the diagnostic information, if FALSE
#' return only the variance estimates.
#' @keywords internal
NbhOut <- function(estimate, neighborhoods, n, N, mu, name, diagnostic) {
  if(diagnostic) {
    nbh_out <- new('NbhOut')
    nbh_out@estimate <- estimate
    nbh_out@neighborhoods <- neighborhoods
    nbh_out@n_groups <- nrow(neighborhoods)
    nbh_out@n <- n
    nbh_out@N <- N
    nbh_out@mu <- mu
    nbh_out@name <- name
    return(nbh_out)
  } else{
    return(estimate)
  }
}
