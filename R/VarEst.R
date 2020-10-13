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
#' @examples
#' v_srs <- VarSRS(fpc=TRUE)
#' v_srs(hex_frame)
#' @export
VarSRS <- function(fpc=FALSE, diagnostic=FALSE) {
  function(sys_frame) {
    return(var_srs(sys_frame, fpc=fpc, diagnostic=diagnostic))
  }
}

#' Neighborhood-based estimator from \insertCite{stevens_variance_2003;textual}{spsys}
#' 
#' Constructor for a local variance estimator introduced in \insertCite{stevens_variance_2003;textual}{spsys}. 
#' This is a wrapper
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
#' @example 
#' v_so <- VarSO(fpc=FALSE, coord_cols=c('s_1', 's_2'))
#' v_so(hex_frame)
#' @importFrom Rdpack reprompt
#' @references {
#'     \insertAllCited{}
#' }
#' @export
VarSO <- function(fpc=FALSE, diagnostic=FALSE, coord_cols=NA, nbh=4, wt_fun=localmean.weight) {
  function(sys_frame) {
    return(var_so(sys_frame, fpc=fpc, diagnostic=diagnostic, coord_cols=coord_cols, nbh=nbh, wt_fun=wt_fun))
  }
}

#' \insertCite{matern_spatial_1986;textual}{spsys} local variance estimator
#' 
#' Constructor for a local variance estimator introduced in \insertCite{matern_spatial_1986;textual}{spsys}.
#' 
#' @param fpc If TRUE, finite population correction factor will be applied.
#' @param diagnostic If TRUE, diagnostic information, i.e. a VarOut, object will be returned.
#' @param nbh The neighborhood structure. For `HexFrame` one of either `'par'` or `'hex'`, for `RectFrame` only
#' `par` is supported and reflects the original structure implemented by Matèrn
#' @return Returns a function with the specified arguments that can be called on a `SysFrame``
#' @example 
#' v_mat <- VarMAT(fpc=FALSE, nbh='par')
#' v_mat(hex_frame)
#' @references {
#'     \insertAllCited{}
#' }
#' @export
VarMAT <- function(fpc=FALSE, diagnostic=FALSE, nbh='par') {
  function(sys_frame) {
    return(var_mat(sys_frame, fpc=fpc, diagnostic=diagnostic, nbh=nbh))
  }
}

#' Non-overlapping stratified random sampling estimator
#' 
#' Constructor for a local variance estimator that aggregates variance estimates 
#' from non-overlapping strata. An example of this estimator is given in \insertCite{aune-lundberg_comparison_2014;textual}{spsys}.
#' 
#' @param fpc If TRUE, finite population correction factor will be applied
#' @param diagnostic If TRUE, diagnostic information, i.e. a VarOut, object will be returned
#' @param nbh The neighborhood structure. For `HexFrame` one of either `tri` `par` or `hex`
#' For `RectFrame` only `par` is supported.
#' @return Returns a function with the specified arguments that can be called on a `SysFrame`
#' @example 
#' v_non <- VarNON(fpc=FALSE, nbh='hex')
#' v_non(hex_frame)
#' @references {
#'     \insertAllCited{}
#' }
#' @export
VarNON <- function(fpc=FALSE, diagnostic=FALSE, nbh='par') {
  function(sys_frame) {
    return(var_non_overlap(sys_frame, fpc=fpc, diagnostic=diagnostic, nbh=nbh))
  }
}

#' D'Orazio's estimator based on Moran's I
#' 
#' Constructor for an estimator that uses Moran's I to adjust an initial estimate via the simple random
#' sampling variance estimator, introduced in \insertCite{d'orazio_estimating_2003}.
#' 
#' @param fpc If TRUE, finite population correction factor will be applied
#' @param diagnostic If TRUE, diagnostic information, i.e. a VarOut, object will be returned
#' @param nbh The neighborhood structure. For `HexFrame` one of either `tri` `par` or `hex`
#' For `RectFrame` only `par` is supported.
#' @return Returns a function with the specified arguments that can be called on a `SysFrame`
#' @example 
#' v_di <- VarDI(fpc=FALSE, order=1)
#' v_di(hex_frame)
#' @references {
#'     \insertAllCited{}
#' }
#' @export
VarDI <- function(fpc=FALSE, diagnostic=FALSE, order=1) {
  function(sys_frame) {
    return(var_dorazio_i(sys_frame, fpc=fpc, diagnostic=diagnostic, order=order))
  }
}

#' D'Orazio's estimator based on Geary's C
#' 
#' Constructor for an estimator that uses Geary's C to adjust an initial estimate via the simple random
#' sampling variance estimator, introduced in \insertCite{d'orazio_estimating_2003}.
#' 
#' @param fpc If TRUE, finite population correction factor will be applied
#' @param diagnostic If TRUE, diagnostic information, i.e. a VarOut, object will be returned
#' @param nbh The neighborhood structure. For `HexFrame` one of either `tri` `par` or `hex`
#' For `RectFrame` only `par` is supported.
#' @return Returns a function with the specified arguments that can be called on a `SysFrame`
#' @example 
#' v_dc <- VarDc(fpc=FALSE, order=1)
#' v_dc(hex_frame)
#' @references {
#'     \insertAllCited{}
#' }
#' @export
VarDC <- function(fpc=FALSE, diagnostic=FALSE, order=1) {
  function(sys_frame) {
    return(var_dorazio_c(sys_frame, fpc=fpc, diagnostic=diagnostic, order=order))
  }
}

#' Systematic variance
#' 
#' Constructor for a function that retrieves the true systematic variance.
#' 
#' @param a The sampling interval
#' @return Returns a function with the specified arguments that can be called on a `SysFrame`
#' @example 
#' v_sys_3 <- VarSYS(a=3)
#' v_sys_3(hex_frame)
#' @export
VarSYS <- function(sys_frame, a) {
  function(sys_frame) {
    return(var_sys(sys_frame, a))
  }
}