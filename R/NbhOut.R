#' @include VarOut.R

#' Base class for all neighborhood variance
#' estimators
setClass('NbhOut', contains='VarOut', slots=list(
  neighborhoods='data.frame',
  n='numeric',
  N='numeric'
))

NbhOut <- function(estimate, neighborhoods, n, N, diagnostic) {
  if(diagnostic) {
    nbh_out <- new('NbhOut')
    nbh_out@estimate <- estimate
    nbh_out@neighborhoods <- neighborhoods
    nbh_out@n <- n
    nbh_out@N <- N
    return(nbh_out)
  } else{
    return(estimate)
  }
}