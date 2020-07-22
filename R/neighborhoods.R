# Methods for constructing neighborhoods


#' For some input set of neighborhood center points
#' return all immediate hexagonal neighbors.
#' TODO could make this Hex/Rect method
#' TODO the default contrasts is a bit messy...
get_hex_neighborhoods <- function(centers, contrasts=c(1,1,1,0,-1,-1,-1)) {
  neighborhoods <- list()
  for(i in 1:nrow(centers)) {
    row <- centers[i, 1]
    col <- centers[i, 2]
    neighborhood <- matrix(NA, nrow=7, ncol=3)
    
    neighborhood[,1] <- c(row-1, row-1, row, row, row, row+1, row+1)
    neighborhood[,2] <- c(col-1, col+1, col-2, col, col+2, col-1, col+1)
    neighborhood[,3] <- contrasts
    
    neighborhood <- data.frame(neighborhood)
    neighborhood$r <- row
    neighborhood$c <- col
    neighborhoods[[i]] <- neighborhood
  }
  neighborhoods <- bind_rows(neighborhoods)
  colnames(neighborhoods)[1:3]  <- c('r_n', 'c_n', 'contr')
  
  return(neighborhoods)
}

#' The Matern estimator does not use 'centers' per se
#' but rather 'anchors'. Also, a set of contrasts is needed
#' for each neighborhood.
get_mat_hex_neighborhoods <- function(anchors, contrasts=c(1,1,-1,-1)) {
  neighborhoods <- list()
  for(i in 1:nrow(anchors)) {
    row <- anchors[i, 1]
    col <- anchors[i, 2]
    neighborhood <- matrix(NA, nrow=4, ncol=3)
    
    neighborhood[,1] <- c(row, row, row + 1, row+1)
    neighborhood[,2] <- c(col, col+2, col-1, col+1)
    neighborhood[,3] <- contrasts
    
    neighborhood <- data.frame(neighborhood)
    neighborhood$r <- row
    neighborhood$c <- col
    neighborhoods[[i]] <- neighborhood
  }
  neighborhoods <- bind_rows(neighborhoods)
  colnames(neighborhoods)[1:3]  <- c('r_n', 'c_n', 'contr')
  
  return(neighborhoods)
}

get_rect_neighborhoods <- function(centers) {
  neighborhoods <- list()
  for(i in 1:nrow(centers)) {
    row <- centers[i, 1]
    col <- centers[i, 2]
    neighborhood <- matrix(NA, nrow=9, ncol=2)
    
    neighborhood[,1] <- c(row-1, row-1, row-1, row, row, row, row+1, row+1, row+1)
    neighborhood[,2] <- c(col-1, col, col+1, col-1, col, col+1, col-1, col, col+1)
    
    neighborhood <- data.frame(neighborhood)
    neighborhood$r <- row
    neighborhood$c <- col
    neighborhoods[[i]] <- neighborhood
  }
  neighborhoods <- bind_rows(neighborhoods)
  colnames(neighborhoods)[1:2]  <- c('r_n', 'c_n')
  
  return(neighborhoods)
  
}

get_mat_rect_neighborhoods <- function(anchors, contrasts) {
  neighborhoods <- list()
  for(i in 1:nrow(anchors)) {
    row <- anchors[i, 1]
    col <- anchors[i, 2]
    neighborhood <- matrix(NA, nrow=4, ncol=3)
    
    neighborhood[,1] <- c(row, row, row+1, row+1)
    neighborhood[,2] <- c(col, col+1, col, col+1)
    neighborhood[,3] <- contrasts
    
    neighborhood <- data.frame(neighborhood)
    neighborhood$r <- row
    neighborhood$c <- col
    neighborhoods[[i]] <- neighborhood
  }
  neighborhoods <- bind_rows(neighborhoods)
  colnames(neighborhoods)[1:3]  <- c('r_n', 'c_n', 'contr')
  
  return(neighborhoods)
}

get_tri_neighborhoods <- function(tops, bottoms) {
  neighborhoods <- matrix(NA, nrow=nrow(tops)*6, ncol=4)
  for(i in 1:nrow(tops)) {
    t_row <- tops[i, 1]
    t_col <- tops[i, 2]
    
    b_row <- bottoms[i, 1]
    b_col <- bottoms[i, 2]
    
    lb <- (i - 1) * 6 + 1
    ub <- lb + 5
    neighborhoods[lb:ub, 1] <- c(t_row, t_row+1, t_row+1, b_row, b_row - 1, b_row - 1)
    neighborhoods[lb:ub, 2] <- c(t_col, t_col-1, t_col+1, b_col, b_col - 1, b_col +1)
    neighborhoods[lb:ub, 3] <- c(rep(t_row, 3), rep(b_row, 3))
    neighborhoods[lb:ub, 4] <- c(rep(t_col, 3), rep(b_col, 3))
  }
  
  colnames(neighborhoods) <- c('r_n', 'c_n', 'r', 'c')
  neighborhoods <- data.frame(neighborhoods)
  return(neighborhoods)
}

setGeneric('neighborhoods_tri', function(hex_frame) {
  standardGeneric('neighborhoods_tri')
})

setMethod('neighborhoods_tri', signature(hex_frame='HexFrame'),
  function(hex_frame) {
    # First we need the top points
    ix <- hex_frame@data[,c('r', 'c')]
    max_r <- max(ix$r)
    max_c <- max(ix$c)
    
    rows <- seq(1, max_r, 2)
    r_samp <- c()
    c_samp <- c()
    
    k <- 1
    for(row in rows) {
      if(k %% 2 != 0) {
        add <- seq(1, max_c, 6)
        c_samp <- c(c_samp, add)
        r_samp <- c(r_samp, rep(row, length(add)))
      } else {
        add <- seq(3, max_c, 6)
        c_samp <- c(c_samp, add)
        r_samp <- c(r_samp, rep(row, length(add)))
      }
      k <- k + 1
    }
    
    tops <- data.frame(r=r_samp, c=c_samp)
    
    bottoms <- tops
    bottoms[,'r'] <- bottoms[,'r'] + 1
    bottoms[,'c'] <- bottoms[,'c'] + 3
    
    neighbors <- get_tri_neighborhoods(tops, bottoms)
    return(neighbors)
  }
)


setGeneric('neighborhoods_ov', function(sys_frame, ...) {
  standardGeneric('neighborhoods_ov')
})


setMethod('neighborhoods_ov', signature(sys_frame='HexFrame'),
  function(sys_frame) {
    ix <- sys_frame@data[,c('r', 'c')]
    centers <- subsample_hex_ix_ov(ix)
    neighbors <- get_hex_neighborhoods(centers)
    return(neighbors)
  }
)


setGeneric('neighborhoods_non', function(sys_frame, ...) {
  standardGeneric('neighborhoods_non')
})


setMethod('neighborhoods_non', signature(sys_frame='HexFrame'),
  function(sys_frame, contrasts=c(1,1,1,0,-1,-1,-1)) {
    ix <- sys_frame@data[,c('r', 'c')]
    centers <- subsample_hex_ix_compact(ix)
    neighbors <- get_hex_neighborhoods(centers, contrasts)
    return(neighbors)
  }
)

setMethod('neighborhoods_non', signature(sys_frame='RectFrame'),
  function(sys_frame) {
    ix <- sys_frame@data[,c('r', 'c')]
    centers <- subsample_rect_ix(ix, c(1,1), 3)
    neighbors <- get_rect_neighborhoods(centers)
    return(neighbors)
  }
)

setGeneric('neighborhoods_par', function(sys_frame, ...) {
  standardGeneric('neighborhoods_par')
})


setMethod('neighborhoods_par', signature(sys_frame='HexFrame'),
  function(sys_frame, contrasts=c(1,1,-1,-1)) {
    ix <- sys_frame@data[,c('r', 'c')]
    anchors <- subsample_hex_ix(ix, c(1,1), 2)
    neighbors <- get_mat_hex_neighborhoods(anchors, contrasts)
    return(neighbors)
  }
)

setMethod('neighborhoods_par', signature(sys_frame='RectFrame'),
  function(sys_frame) {
    ix <- sys_frame@data[,c('r', 'c')]
    anchors <- subsample_rect_ix(ix, c(1,1), 2)
    neighbors <- get_mat_rect_neighborhoods(anchors)
    return(neighbors)
  }
)