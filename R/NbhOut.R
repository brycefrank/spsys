#' @include VarOut.R

#' Base class for all neighborhood variance
#' estimators
setClass('NbhOut', contains='VarOut', slots=list(
  neighborhoods='data.frame',
  n='numeric',
  N='numeric',
  n_groups='numeric',
  name='character'
))

NbhOut <- function(estimate, neighborhoods, n, N, name, diagnostic) {
  if(diagnostic) {
    nbh_out <- new('NbhOut')
    nbh_out@estimate <- estimate
    nbh_out@neighborhoods <- neighborhoods
    nbh_out@n <- n
    nbh_out@N <- N
    nbh_out@name <- name
    return(nbh_out)
  } else{
    return(estimate)
  }
}

setMethod('summary', signature(var_out = 'NbhOut'), 
  function(var_out) {
    callNextMethod(var_out)
    cat('\nNumber of neighborhoods: ', nrow(var_out@neighborhoods), '\n')
    print(head(var_out@neighborhoods))
  }
)