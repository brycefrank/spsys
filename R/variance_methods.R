# Variance estimators in method form
library(spsurvey)
library(tidyr)

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
# TODO check appropriateness of FPC
setMethod('var_mat', signature(sys_frame='HexFrame'),
  # TODO fix contrast defaults to something that makes sense
  function(sys_frame, fpc=FALSE, N=NA_real_, contrasts = c(1,-1, 0, 0, 0, 1, -1)) {
    neighborhoods <- get_hex_neighborhoods(sys_frame@data[,c('r','c')], contrasts=contrasts)
    
    # Mean-center the attributes
    att_df <- sys_frame@data[, sys_frame@attributes, drop=FALSE]
    att_ct <- att_df  - rep(colMeans(att_df), rep.int(nrow(att_df), ncol(att_df)))
    d_ct <- cbind(sys_frame@data[,c('r', 'c')], att_ct)
    neighborhoods <- merge(neighborhoods, d_ct, by.x=c('r_n', 'c_n'), by.y=c('r', 'c'), all.x=TRUE)
    
    # Get the number of neighborhoods with at least one point in Q
    p <- length(colnames(att_df))
    
    # Set the NA values of neighborhood attributes (i.e. they are out of the area) 
    # to zero if they are
    zero_list <- as.list(rep(0,p))
    names(zero_list) <- colnames(att_df)
    in_Q <- neighborhoods %>%
      replace_na(zero_list)
    
    q <- 9
    n <- nrow(sys_frame@data)
    
    # A function that generates the T value for each neighborhood
    summ_df <- function(x) {return(sum(x * in_Q$contrasts)^2  / 4)}
    
    Ti <- in_Q %>%
      group_by(r, c) %>%
      mutate_at(.vars=colnames(att_df), .funs=function(x){return(x*contrasts)}) %>%
      summarize_at(.vars=colnames(att_df), .funs=function(x){return(sum(x)^2/4)})
    
    var <- (q * sum(Ti[,sys_frame@attributes,drop=FALSE])) / n^2
    return(var)
  }
)

#' Calculate the variance using a denominator of n
#' instead of n-1
pop_var <- function(z) {
  n <- length(z)
  ssq <- sum((z - mean(z))^2)
  ssq/n
}


# TODO give this a more descriptive name
# TODO implement fpc in a better way upstream
weight_var <- function(var, q_j, fpc, N_neighbs) {
  (1/N_neighbs)^2 * (var / q_j) * fpc
}

setGeneric('var_non_overlap', function(sys_frame, ...) {
  standardGeneric('var_non_overlap')
})

setMethod('var_non_overlap', signature(sys_frame = 'HexFrame'), 
  function(sys_frame, fpc=FALSE, N=NA_real_) {
    neighborhoods <- get_hex_neighborhoods(sys_frame@data[,c('r','c')])
    neighborhoods <- merge(neighborhoods, sys_frame@data, by.x=c('r_n', 'c_n'), by.y=c('r', 'c'), all.x=TRUE)
    att_df <- sys_frame@data[, sys_frame@attributes, drop=FALSE]
    n <- nrow(sys_frame@data)
    
    N_neighbs <- neighborhoods %>%
      group_by(r, c) %>%
      summarize(n=n()) %>%
      nrow()
    
    # Prepare a vector specifying the pop variance function for each attribute
    p <- length(colnames(att_df))
    funs <- rep('pop_var', p)
    names(funs) <- paste(colnames(att_df), 'var', sep='_')
    
    q <- neighborhoods %>%
      na.omit() %>%
      group_by(r, c) %>%
      summarize(q_j=n())
    
    # TODO some neighborhoods return a 0 variance
    neighborhoods <- neighborhoods %>%
      na.omit() %>%
      group_by(r, c) %>%
      summarize_at(.vars = colnames(att_df), .funs = funs) %>%
      merge(q) %>%
      mutate(N_j = q_j * sys_frame@a^2) %>%
      mutate(w_j_sq = (N_j / n)^2, fpc = ((N_j - q_j) / N_j)) %>%
      mutate_at(.vars = names(funs), .funs=~weight_var(., q_j, fpc, N_neighbs)) %>%
      summarize_at(.vars = names(funs), .funs=~sum(.))
    
    neighborhoods
  }
)

setGeneric('var_dorazio_c', function(sys_frame, ...) {
  standardGeneric('var_dorazio_c')
})


setMethod('var_dorazio_c', signature(sys_frame = 'HexFrame'), 
  function(sys_frame, fpc=FALSE, N=NA_real_, order=1) {
    v_srs <- var_srs(sys_frame, fpc=fpc, N=N)
    C <- gearys_c(sys_frame, order=order)
    print(paste('C:', C))
    return(v_srs * C)
  }
)

setGeneric('var_dorazio_i', function(sys_frame, ...) {
  standardGeneric('var_dorazio_i')
})


setMethod('var_dorazio_i', signature(sys_frame = 'HexFrame'), 
  function(sys_frame, fpc=FALSE, N=NA_real_, order=1) {
    v_srs <- var_srs(sys_frame, fpc=fpc, N=N)
    morans_I <- morans_i(sys_frame, order=order)
    
    if(morans_I > 0) {
      w <- 1 + 2/log(morans_I) + 2/(1/morans_I - 1)
    } else {
      w <- 1
    }
    print(paste('W:', w))
    return(v_srs * w)
  }
)