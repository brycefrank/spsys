#' @include VarOut.R

#' Base class for all neighborhood variance
#' estimators
setClass('NbhOut', contains='VarOut', slots=list(
  neighborhoods='data.frame',
  n='numeric',
  N='numeric',
  mu='numeric',
  n_groups='numeric',
  name='character'
))

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
