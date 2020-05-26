# Variance estimators in method form
library(spsurvey)

setGeneric('var_srs', function(sys_frame, ...) {
  standardGeneric('var_srs')
})

# TODO could be possible to attach N in the case of subsamples
setMethod('var_srs', signature(sys_frame='SysFrame'),
  function(sys_frame, fpc=FALSE, N=NA_real_) {
    if(fpc == TRUE & is.na(N)) {
      stop('If fpc is set to true you must provide a population size N.')
    }
    
    att_df <- sys_frame@data[, sys_frame@attributes, drop=FALSE]
    n <- length(sys_frame)
    var_mu <- (n - 1)^(-1) * n^(-1) *  colSums((att_df - colMeans(att_df))^2)
    
    if(fpc) {
      var_mu <- var_mu * (1-n/N)
    }
    return(var_mu)
  }
)

#' Computes the variance estimate of a systematic 
#' sample using the Stevens and Olsen (2003) local mean estimator.
#' 
#' This is a wrapper for functions that exist in `spsurvey`
setGeneric('var_so', function(sys_frame, ...){
  standardGeneric('var_so')
})

setMethod('var_so', signature(sys_frame='SysFrame'),
  function(sys_frame, fpc=FALSE, N=NA_real_) {
    # TODO this should not strictly be needed
    if(is.na(N)) {
      stop('var_so requires a population size N')
    }
    
    n <- nrow(sys_frame)
    
    if(!'pi_i' %in% colnames(sys_frame@data)) {
      warning('Inclusion probability column - "pi_i" not found in @data, assuming
               inclusion probabilities are n/N.')
      pi_i <- rep(n/N, n)
    } else {
      pi_i <- sys_frame@data$pi_i
    }
    
    wt <- localmean.weight(sys_frame@coords[,1], sys_frame@coords[,2], pi_i)
    att_df <- sys_frame@data[, sys_frame@attributes, drop=FALSE]
    var_total <- sapply(colnames(att_df), so_att, att_df, pi_i, wt)
    
    if(fpc) {
      var_mu <- (1/N^2) * var_total * (1 - n/N)
    } else {
      var_mu <- (1/N^2) * var_total
    }
    return(var_mu)
  }
)

#' Computes the Stevens and Olsen estimator for a specific attribute. 
#' Used internally with var_so function.
so_att <- function(att, att_df, pi_i, wt) {
  z <- att_df[,att]
  localmean.var(z/pi_i, wt)
}

setGeneric('var_mat', function(sys_frame, ...) {
  standardGeneric('var_mat')
})

# TODO make sure in documentation that we explain the assumption
# that points provided to this function are ONLY those points that 
# are inside Q (which is the study area in Matern's notatio)
setMethod('var_mat', signature(sys_frame='HexFrame'),
  # TODO fix contrast defaults to something that makes sense
  function(sys_frame, fpc=FALSE, N=NA_real_, contrasts = c(-1,-1,-1, 0, 1, 1, 1)) {
    neighborhoods <- get_hex_neighborhoods(sys_frame@data[,c('r','c')], contrasts=contrasts)
    
    # Mean-center the attributes
    att_df <- sys_frame@data[, sys_frame@attributes, drop=FALSE]
    att_ct <- att_df  - rep(colMeans(att_df), rep.int(nrow(att_df), ncol(att_df)))
    d_ct <- cbind(sys_frame@data[,c('r', 'c')], att_ct)
    neighborhoods <- merge(neighborhoods, d_ct, by.x=c('r_n', 'c_n'), by.y=c('r', 'c'), all.x=TRUE)
    
    # Get the number of neighborhoods with at least one point in Q
    p <- length(colnames(att_df))
    
    # Set the NA values of neighborhood attributes to zero
    # this implies that they are outside of the study area
    zero_list <- as.list(rep(0,p))
    names(zero_list) <- colnames(att_df)
    in_Q <- neighborhoods %>%
      replace_na(zero_list)
    
    # What is the size of each group?
    # TODO how to handle semi-complete cases? e.g. att1=30 and att2=NA, etc? Maybe force the issue
    # when SysFrame is created
    n_groups <- neighborhoods %>%
      na.omit() %>%
      group_by(r, c) %>%
      summarize(n = n())
    
    q <- 9
    n <- nrow(sys_frame@data)
    summ_df <- function(x) {return(sum(x * in_Q$contrasts)^2  / 4)}
    
    Ti <- in_Q %>%
      group_by(r, c) %>%
      summarize_at(.vars=c('z_1'), .funs='summ_df')
    
    var <- (q * sum(Ti[,sys_frame@attributes])) / n^2
    return(var)
  }
)