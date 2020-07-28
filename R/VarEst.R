# Closures used to construct variance estimators, these are the entry points into
# the package for most users.

#' Variance estimator assuming simple random sampling
#' 
#' Constructor for the variance estimator assuming simple random sampling.
#' 
#' @param fpc If TRUE, finite population correction factor will be applied.
#' @param diagnostic If TRUE, diagnostic information, i.e. a VarOut, object will be returned.
#' @return Returns a function with the specified arguments that can be called on a 
#' SysFrame.
VarSRS <- function(fpc=FALSE, diagnostic=FALSE) {
  function(sys_frame) {
    return(var_srs(sys_frame, fpc=fpc, diagnostic=diagnostic))
  }
}

#' Stevens and Olsen (2003) local variance estimator
#' 
#' Constructor for a local variance estimator introduced in Stevens and Olsen (2003). This is a wrapper
#' for an implementation in the `spsurvey` package. A column named `pi_i` must be present in the dataframe 
#' contained in the `@data` slot that specifies the first order inclusion probabilities of each element in the sample.
#' 
#' @param fpc If TRUE, finite population correction factor will be applied.
#' @param diagnostic If TRUE, diagnostic information, i.e. a VarOut, object will be returned.
#' @param coord_cols A character vector of the names of the coordinates in the `sys_frame` that 
#' contain the x and y coordinates. The x column name should be listed first.
#' @param nbh The number of neighbors to use in the construction of local neighborhoods.
#' @return Returns a function with the specified arguments that can be called on a 
#' SysFrame.
VarSO <- function(fpc=FALSE, diagnostic=FALSE, coord_cols=NA, nbh=4) {
  function(sys_frame) {
    return(var_so(sys_frame, fpc=fpc, diagnostic=diagnostic, coord_cols=coord_cols, nbh=nbh))
  }
}

#' Matern (1984) local variance estimator
#' 
#' Constructor for a local variance estimator introduced by Matern in [] and [].
#' 
#' @param fpc If TRUE, finite population correction factor will be applied.
#' @param diagnostic If TRUE, diagnostic information, i.e. a VarOut, object will be returned.
#' @param nbh The neighborhood structure. For `HexFrame` one of either `par` or `hex`, for `RectFrame` only
#' `par` is supported and reflects the original structure implemented by Matern.
#' @return Returns a function with the specified arguments that can be called on a 
#' SysFrame.
VarMAT <- function(fpc=FALSE, diagnostic=FALSE, nbh='par') {
  function(sys_frame) {
    return(var_mat(sys_frame, fpc=fpc, diagnostic=diagnostic, nbh=nbh))
  }
}

#' Non-overlapping stratified random sampling estimator
#' 
#' Constructor for a local variance estimator that aggregates variance estimates 
#' from non-overlapping strata.
#' 
#' @param fpc If TRUE, finite population correction factor will be applied.
#' @param diagnostic If TRUE, diagnostic information, i.e. a VarOut, object will be returned.
#' @param nbh The neighborhood structure. For `HexFrame` one of either `tri` `par` or `hex`.
#' For `RectFrame` only `par` is supported.
VarNON <- function(fpc=FALSE, diagnostic=FALSE, nbh='par') {
  function(sys_frame) {
    return(var_non_overlap(sys_frame, fpc=fpc, diagnostic=diagnostic, nbh=nbh))
  }
}

VarDI <- function(fpc=FALSE, diagnostic=FALSE, order=1) {
  function(sys_frame) {
    return(var_dorazio_i(sys_frame, fpc=fpc, diagnostic=diagnostic, order=order))
  }
}

VarDC <- function(fpc=FALSE, diagnostic=FALSE, order=1) {
  function(sys_frame) {
    return(var_dorazio_c(sys_frame, fpc=fpc, diagnostic=diagnostic, order=order))
  }
}

VarSYS <- function(sys_frame, a) {
  function(sys_frame) {
    return(var_sys(sys_frame, a))
  }
}