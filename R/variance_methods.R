library(spsurvey)
library(tidyr)

#' Internal function used to compute the simple random sampling variance estimate.
#' 
#' @param sys_frame A `SysFrame` object
#' @param fpc If TRUE, use the finite population correction factor
#' @param diagnostic If TRUE return diagnostic information, else return only the variance estimate
#' @keywords internal
setGeneric('var_srs', function(sys_frame, ...) {
  standardGeneric('var_srs')
})

setMethod('var_srs', signature(sys_frame='SysFrame'),
  function(sys_frame, fpc=FALSE, diagnostic=FALSE) {
    N <- sys_frame@N
    if(fpc == TRUE & N==Inf) {
      stop('If fpc is set to true you must provide a finite population size N.')
    }
    
    att_df <- sys_frame@data[, sys_frame@attributes, drop=FALSE]
    
    
    n <- length(sys_frame)
    
    att_means <- colMeans(att_df)
    ssq <- colSums(sweep(att_df, 2, att_means)^2)
    
    var_mu <- (n - 1)^(-1) * n^(-1) * ssq 
    
    if(fpc) {
      var_mu <- var_mu * (1-n/N)
    }
    
    mu <- colMeans(att_df)
    var_mu <- VarOut(var_mu, n, N, mu, diagnostic)
    return(var_mu)
  }
)

#' Internal function used to compute the systematic variance.
#' 
#' @param sys_frame A `SysFrame` object
#' @param a The sampling interval
#' @keywords internal
setGeneric('var_sys', function(sys_frame, ...){
  standardGeneric('var_sys')
})

setMethod('var_sys', signature(sys_frame='SysFrame'),
  function(sys_frame, a) {
    starts <- subsample_starts(sys_frame, a)
    K <- nrow(starts)
    mu <- colMeans(sys_frame@data[,sys_frame@attributes])
    sse <- 0
    
    # TODO diagnostic should be a dataframe with the sample indices and the
    # means of each sample
    # can't think of anything else to add?
    
    for(k in 1:K) {
      sub <- subsample(sys_frame, starts[k,], a)
      mu_hat <- colMeans(sub@data[,sys_frame@attributes])
      sq_err <- (mu - mu_hat)^2
      sse <- sse + sq_err
    }
    return(1/K * sse)
  }
)

#' Internal function used to compute the Stevens and Olsen (2003) variance estimate.
#' 
#' @param sys_frame A `SysFrame` object
#' @param fpc If TRUE, use the finite population correction factor
#' @param diagnostic If TRUE return diagnostic information, else return only the variance estimate
#' @param coord_cols A named vector of coordinate columns
#' @param nbh An integer representing the number of neighbors to use in the variance estimate
#' @keywords internal
setGeneric('var_so', function(sys_frame, ...){
  standardGeneric('var_so')
})

setMethod('var_so', signature(sys_frame='SysFrame'),
  function(sys_frame, fpc=FALSE, diagnostic=FALSE, coord_cols=NA, nbh=4) {
    N <- sys_frame@N
    if(N==Inf) {
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
    
    if(!is.na(coord_cols)) {
      coords <- sys_frame@data[,coord_cols]
    } else {
      coords <- sys_frame@coords
    }
    
    wt <- localmean.weight(coords[,1], coords[,2], pi_i, nbh=nbh)
    att_df <- sys_frame@data[, sys_frame@attributes, drop=FALSE]
    var_total <- sapply(colnames(att_df), so_att, att_df, pi_i, wt)
    mu <- colMeans(att_df)
    
    if(fpc) {
      var_mu <- (1/N^2) * var_total * (1 - n/N)
    } else {
      var_mu <- (1/N^2) * var_total
    }
    
    var_mu <- VarOut(var_mu, n, N, mu, diagnostic)
    return(var_mu)
  }
)

#' Computes the Stevens and Olsen estimator for a specific attribute. 
#' Used internally with var_so function.
so_att <- function(att, att_df, pi_i, wt) {
  z <- att_df[,att]
  localmean.var(z/pi_i, wt)
}

#' Internal function used to compute the Matèrn (1986) variance estimate.
#' 
#' @param sys_frame A `SysFrame` object
#' @param fpc If TRUE, use the finite population correction factor
#' @param diagnostic If TRUE return diagnostic information, else return only the variance estimate
#' @param nbh For `HexFrames` either 'par' or 'hex', for `RectFrames` only `par` is supported
#' @keywords internal
setGeneric('var_mat', function(sys_frame, ...) {
  standardGeneric('var_mat')
})

setMethod('var_mat', signature(sys_frame='HexFrame'),
  function(sys_frame, fpc=FALSE, diagnostic=FALSE, nbh='par') {
    if(nbh=='par') {
      neighborhoods <- neighborhoods_par(sys_frame)
      h <- 4
      q <- 4 # Number of elements per neighborhood
    } else if(nbh=='hex') {
      neighborhoods <- neighborhoods_non(sys_frame)
      h <- 9
      q <- 7
    } else {
      stop('Please specify either a "par" or "hex" neighborhood structure.')
    }
    
    N <- sys_frame@N
    
    # Mean-center the attributes
    att_df <- sys_frame@data[, sys_frame@attributes, drop=FALSE]
    att_ct <- att_df  - rep(colMeans(att_df), rep.int(nrow(att_df), ncol(att_df)))
    d_ct <- cbind(sys_frame@data[,c('r', 'c')], att_ct)
    neighborhoods <- merge(neighborhoods, d_ct, by.x=c('r_n', 'c_n'), by.y=c('r', 'c'), all.x=TRUE)
    atts <- colnames(att_df)
    
    # Get the number of neighborhoods with at least one point in Q
    p <- length(atts)
    
    # Set the NA values of neighborhood attributes (i.e. they are out of the area) 
    # to zero if they are
    zero_list <- as.list(rep(0,p))
    names(zero_list) <- atts
    
    in_Q <- neighborhoods %>%
      replace_na(zero_list)
    
    n <- nrow(sys_frame@data)
    
    # A function that generates the T value for each neighborhood
    get_Ti <- function(x, contr) {return(sum(x * contr)^2  / h)}
    
    Ti <- in_Q %>%
      group_by(r, c) %>%
      summarize_at(.vars=atts, .funs=~get_Ti(., contr))
    
    var_mat <- (q * colSums(Ti[,sys_frame@attributes,drop=FALSE])) / n^2
    
    if(fpc) {
      var_mat <- (1-n/N) * var_mat
    }
    
    mu <- colMeans(att_df)
    var_mat <- NbhOut(var_mat, Ti, n, N, mu, 'var_mat', diagnostic)
    
    return(var_mat)
  }
)

setMethod('var_mat', signature(sys_frame='RectFrame'),
  function(sys_frame, fpc=FALSE, diagnostic=FALSE, nbh='par') {
    if(nbh!='par') {
      stop('Only the "par" neigbhorhood structure is defined for RectFrame')
    }
    
    neighborhoods <- neighborhoods_par(sys_frame)
    h <- 4
    q <- 4 # Number of elements per neighborhood
    
    N <- sys_frame@N
    
    # Mean-center the attributes
    att_df <- sys_frame@data[, sys_frame@attributes, drop=FALSE]
    att_ct <- att_df  - rep(colMeans(att_df), rep.int(nrow(att_df), ncol(att_df)))
    d_ct <- cbind(sys_frame@data[,c('r', 'c')], att_ct)
    neighborhoods <- merge(neighborhoods, d_ct, by.x=c('r_n', 'c_n'), by.y=c('r', 'c'), all.x=TRUE)
    atts <- colnames(att_df)
    
    # Get the number of neighborhoods with at least one point in Q
    p <- length(atts)
    
    # Set the NA values of neighborhood attributes (i.e. they are out of the area) 
    # to zero if they are
    zero_list <- as.list(rep(0,p))
    names(zero_list) <- atts
    
    in_Q <- neighborhoods %>%
      replace_na(zero_list)
    
    n <- nrow(sys_frame@data)
    
    # A function that generates the T value for each neighborhood
    get_Ti <- function(x, contr) {return(sum(x * contr)^2  / h)}
    
    Ti <- in_Q %>%
      group_by(r, c) %>%
      summarize_at(.vars=atts, .funs=~get_Ti(., contr))
    
    var_mat <- (q * colSums(Ti[,sys_frame@attributes,drop=FALSE])) / n^2
    
    if(fpc) {
      var_mat <- (1-n/N) * var_mat
    }
    
    mu <- colMeans(att_df)
    var_mat <- NbhOut(var_mat, Ti, n, N, mu, 'var_mat', diagnostic)
    
    return(var_mat)
  }
)

#' Calculate the variance using a denominator of n
#' instead of n-1
pop_var <- function(z) {
  n <- length(z)
  ssq <- sum((z - mean(z))^2)
  ssq/n
}

weight_var <- function(var, q_j, fpc, N_neighbs) {
  (1/N_neighbs)^2 * (var / q_j) * fpc
}


#' Internal function used to compute the non-overlapping variance estimate.
#' 
#' @param sys_frame A `SysFrame` object
#' @param fpc If TRUE, use the finite population correction factor
#' @param diagnostic If TRUE return diagnostic information, else return only the variance estimate
#' @param nbh For `HexFrames` either 'par' or 'hex' or 'tri', for `RectFrames` only `par` is supported
#' @keywords internal
setGeneric('var_non_overlap', function(sys_frame, ...) {
  standardGeneric('var_non_overlap')
})

setMethod('var_non_overlap', signature(sys_frame = 'RectFrame'),
  function(sys_frame, fpc=FALSE, diagnostic=FALSE, nbh='par') {
    if(nbh!='par') {
      stop('Only the "par" neigbhorhood structure is defined for RectFrame')
    }
    
    nbh <- neighborhoods_non(sys_frame)
    nbh <- merge(nbh, sys_frame@data, by.x=c('r_n', 'c_n'), by.y=c('r', 'c'), all.x=TRUE)
    
    atts <- sys_frame@attributes
    att_df <- sys_frame@data[, atts, drop=FALSE]
    n <- nrow(sys_frame@data)
    N <- sys_frame@N
    
    neighbor_groups <- nbh %>%
      drop_na(c('r', 'c', atts)) %>%
      group_by(r,c)
    
    N_neighbs <- neighbor_groups %>%
      summarize(n=n()) %>%
      nrow()
    
    # Prepare a vector specifying the pop variance function for each attribute
    p <- length(colnames(att_df))
    
    q <- neighbor_groups %>%
      summarize(q_j=n())
    
    # TODO some neighborhoods return a 0 variance
    nbh_var <- neighbor_groups %>%
      summarize_at(.vars = atts, pop_var) %>%
      merge(q) %>%
      mutate(N_j = q_j * sys_frame@a^2) %>% # TODO is there someway to specify this without a?
      mutate(w_j_sq = (N_j / n)^2, fpc = ((N_j - q_j) / N_j)) %>%
      mutate_at(.vars = colnames(att_df), .funs=~weight_var(., q_j, fpc, N_neighbs))
    
    var_non <- nbh_var %>%
      summarize_at(.vars = colnames(att_df), .funs=~sum(.))
    
    if(fpc) {
      var_non <- (1-n/N) * var_non
    }
    
    mu <- colMeans(att_df)
    var_non <- as.numeric(var_non)
    names(var_non) <- sys_frame@attributes
    var_non <- NbhOut(var_non, nbh_var, n, N, mu, 'var_non_overlap', diagnostic)
    return(var_non)
  }
)

setMethod('var_non_overlap', signature(sys_frame = 'HexFrame'),
  function(sys_frame, fpc=FALSE, nbh='tri', diagnostic=FALSE) {
    if(nbh=='tri') {
      nbh <- neighborhoods_tri(sys_frame)
    } else if(nbh=='hex') {
      nbh <- neighborhoods_non(sys_frame)
    } else if(nbh=='par') {
      nbh <- neighborhoods_par(sys_frame)
    } else {
      stop('Please specify a neighborhood structure - either "tri", "hex" or "par"')
    }
    
    if(identical(sys_frame@a, numeric(0))) {
      stop('This estimator requires specification of the sampling interval. Please set one using the
           @a slot.')
    }
    
    N <- sys_frame@N
    nbh <- merge(nbh, sys_frame@data, by.x=c('r_n', 'c_n'), by.y=c('r', 'c'), all.x=TRUE)
    atts <- sys_frame@attributes
    att_df <- sys_frame@data[, atts, drop=FALSE]
    n <- nrow(sys_frame@data)
    
    neighbor_groups <- nbh %>%
      drop_na(c('r', 'c', atts)) %>%
      group_by(r,c)
    
    N_neighbs <- neighbor_groups %>%
      summarize(n=n()) %>%
      nrow()
    
    # Prepare a vector specifying the pop variance function for each attribute
    p <- length(colnames(att_df))
    
    q <- neighbor_groups %>%
      summarize(q_j=n())
    
    nbh_var <- neighbor_groups %>%
      summarize_at(.vars = atts, pop_var) %>%
      merge(q) %>%
      mutate(N_j = q_j * sys_frame@a^2) %>% # TODO is there someway to specify this without a?
      mutate(w_j_sq = (N_j / n)^2, fpc = ((N_j - q_j) / N_j)) %>%
      mutate_at(.vars = colnames(att_df), .funs=~weight_var(., q_j, fpc, N_neighbs))
    
    var_non <- nbh_var %>%
      summarize_at(.vars = colnames(att_df), .funs=~sum(.))
    
    if(fpc) {
      var_non <- (1-n/N) * var_non
    }
    
    mu <- colMeans(att_df)
    
    var_non <- as.numeric(var_non)
    names(var_non) <- sys_frame@attributes
    var_non <- NbhOut(var_non, nbh_var, n, N, mu, 'var_non_overlap', diagnostic)
    return(var_non)
  }
)

#' Internal function used to compute the D'Orazio's - C variance estimate.
#' 
#' @param sys_frame A `SysFrame` object
#' @param fpc If TRUE, use the finite population correction factor
#' @param diagnostic If TRUE return diagnostic information, else return only the variance estimate
#' @param order The neighborhood order. 1 represents the immediate neighbors
#' @keywords internal
setGeneric('var_dorazio_c', function(sys_frame, ...) {
  standardGeneric('var_dorazio_c')
})


setMethod('var_dorazio_c', signature(sys_frame = 'SysFrame'), 
  function(sys_frame, fpc=FALSE, order=1, diagnostic=FALSE) {
    att_df <- sys_frame@data[,sys_frame@attributes]
    v_srs <- var_srs(sys_frame, fpc=fpc)
    C <- gearys_c(sys_frame, order=order)
    var_c <- v_srs * C
    n <- nrow(sys_frame@data)
    mu <- colMeans(att_df)
    
    var_c <- AdjOut(var_c, n, sys_frame@N, mu, C, diagnostic)
  }
)


#' Internal function used to compute the D'Orazio's - I variance estimate.
#' 
#' @param sys_frame A `SysFrame` object
#' @param fpc If TRUE, use the finite population correction factor
#' @param diagnostic If TRUE return diagnostic information, else return only the variance estimate
#' @param order The neighborhood order. 1 represents the immediate neighbors
#' @keywords internal
setGeneric('var_dorazio_i', function(sys_frame, ...) {
  standardGeneric('var_dorazio_i')
})

setMethod('var_dorazio_i', signature(sys_frame = 'SysFrame'), 
  function(sys_frame, fpc=FALSE, order=1, diagnostic=FALSE) {
    v_srs <- var_srs(sys_frame, fpc=fpc)
    morans_I <- morans_i(sys_frame, order=order)
    n <- nrow(sys_frame@data)
    p <- length(sys_frame@attributes)
    att_df <- sys_frame@data[,sys_frame@attributes, drop=FALSE]
    mu <- colMeans(att_df)
    
    w <- rep(1, p)
    
    gt0 <- morans_I > 0
    gt0[is.na(gt0)] <- FALSE
    
    w[gt0]  <- 1 + (2/log(morans_I[gt0])) + (2/(1/morans_I[gt0] - 1))
    var_i <- v_srs * w
    names(w) <- sys_frame@attributes
    var_i <- AdjOut(var_i, n, sys_frame@N, mu, w, diagnostic)
    return(var_i)
  }
)

