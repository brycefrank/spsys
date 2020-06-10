# Variance estimators in method form
library(spsurvey)
library(tidyr)

setGeneric('var_srs', function(sys_frame, ...) {
  standardGeneric('var_srs')
})

setMethod('var_srs', signature(sys_frame='SysFrame'),
  function(sys_frame, fpc=FALSE, N=NA_real_) {
    if(fpc == TRUE & is.na(N)) {
      stop('If fpc is set to true you must provide a population size N.')
    }
    
    att_df <- sys_frame@data[, sys_frame@attributes, drop=FALSE]
    n <- length(sys_frame)
    
    att_means <- colMeans(att_df)
    ssq <- colSums(sweep(att_df, 2, att_means)^2)
    
    var_mu <- (n - 1)^(-1) * n^(-1) * ssq 
    
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

# FIXME is FPC appropriate for these?
setMethod('var_mat', signature(sys_frame='SysFrame'),
  function(sys_frame, fpc=FALSE, N=NA_real_) {
    neighborhoods <- neighborhoods_mat(sys_frame)
    
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
    
    q <- 4
    n <- nrow(sys_frame@data)
    
    # A function that generates the T value for each neighborhood
    get_Ti <- function(x, contr) {return(sum(x * contr)^2  / 4)}
    
    Ti <- in_Q %>%
      group_by(r, c) %>%
      summarize_at(.vars=atts, .funs=~get_Ti(., contr))
    
    var_mat <- (q * colSums(Ti[,sys_frame@attributes,drop=FALSE])) / n^2
    
    if(fpc) {
      var_mat <- (1-n/N) * var_mat
    }
    
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

# TODO give this a more descriptive name
# TODO implement fpc in a better way upstream
weight_var <- function(var, q_j, fpc, N_neighbs) {
  (1/N_neighbs)^2 * (var / q_j) * fpc
}

setGeneric('var_non_overlap', function(sys_frame, ...) {
  standardGeneric('var_non_overlap')
})

# FIXME is FPC appropriate for these?
setMethod('var_non_overlap', signature(sys_frame = 'SysFrame'),
  function(sys_frame, fpc=FALSE, N=NA_real_) {
    nbh <- neighborhoods_non(sys_frame)
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
    
    # TODO some neighborhoods return a 0 variance
    var_non <- neighbor_groups %>%
      summarize_at(.vars = atts, pop_var) %>%
      merge(q) %>%
      mutate(N_j = q_j * sys_frame@a^2) %>% # TODO is there someway to specify this without a?
      mutate(w_j_sq = (N_j / n)^2, fpc = ((N_j - q_j) / N_j)) %>%
      mutate_at(.vars = colnames(att_df), .funs=~weight_var(., q_j, fpc, N_neighbs)) %>%
      summarize_at(.vars = colnames(att_df), .funs=~sum(.))
    
    if(fpc) {
      var_non <- (1-n/N) * var_non
    }
    
    return(var_non)
  }
)

setGeneric('var_overlap', function(sys_frame, ...) {
  standardGeneric('var_overlap')
})

mse <- function(x) {
  m <- length(x)
  x_bar <- mean(x)
  (1/m) * sum((x - x_bar)^2)
}

# FIXME I have checked this estimator up and down and it seems
# to match the Aune-Lundberg description. Wait to see how it does for
# other populations / RectFrames before going further.
setMethod('var_overlap', signature(sys_frame = 'SysFrame'),
  function(sys_frame, fpc=FALSE, N=NA_real_) {
    browser()
    nbh <- neighborhoods_ov(sys_frame)
    nbh <- merge(nbh, sys_frame@data, by.x=c('r_n', 'c_n'), by.y=c('r', 'c'), all.x=TRUE)
    
    atts <- sys_frame@attributes
    n <- nrow(sys_frame@data)
    
    neighbor_groups <- nbh %>%
      drop_na(c('r', 'c', atts)) %>%
      group_by(r,c)
    
    grp_vars <- neighbor_groups %>%
      summarize_at(.vars = atts, .funs=~mse(.))
    
    var_ov <- (1/n^2) * colSums(grp_vars[,atts])
    
    if(fpc) {
      var_ov <- (1-n/N) * var_ov
    }
    
    return(var_ov)
    
  }
)


setGeneric('var_dorazio_c', function(sys_frame, ...) {
  standardGeneric('var_dorazio_c')
})


# FIXME does not work for multiple variables
setMethod('var_dorazio_c', signature(sys_frame = 'SysFrame'), 
  function(sys_frame, fpc=FALSE, N=NA_real_, order=1) {
    v_srs <- var_srs(sys_frame, fpc=fpc, N=N)
    C <- gearys_c(sys_frame, order=order)
    var_c <- v_srs * C
    n <- nrow(sys_frame@data)
    
    return(v_srs * C)
  }
)

setGeneric('var_dorazio_i', function(sys_frame, ...) {
  standardGeneric('var_dorazio_i')
})


# TODO this seems to consistently underestimate most of the variables
# but seems to be fine for uncorrelated. Could just be a poor estimator *shrug*
setMethod('var_dorazio_i', signature(sys_frame = 'SysFrame'), 
  function(sys_frame, fpc=FALSE, N=NA_real_, order=1) {
    v_srs <- var_srs(sys_frame, fpc=fpc, N=N)
    morans_I <- morans_i(sys_frame, order=order)
    n <- nrow(sys_frame@data)
    p <- length(sys_frame@attributes)
    
    w <- rep(1, p)
    
    gt0 <- morans_I[morans_I > 0]
    w[morans_I > 0]  <- 1 + (2/log(gt0)) + (2/(1/gt0 - 1))
    var_i <- v_srs * w
    
    return(var_i)
  }
)