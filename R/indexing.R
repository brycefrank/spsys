library(rgdal)
library(rgeos)

#' A function to index arbitrarily "warped" sets of hexagon points
#' this arose because the hexagonal tesselation used by the FIA
#' is warped over geographic space, which can break "naive"
#' indexing operations that use linear transformations.
#' 
#' The general approach is to iteratively map neighboring hexagons
#' to the index set, which allows for some "shifting" between points
#' without breaking across large distances.
#' 
#' The requirements are that: 1) the hexagons are all approximately the same
#' size (i.e. the hexagon points are approximately evenly spaced) and 
#' 2) that every hexagon that needs to be indexed has at least one neighbor
#' (i.e. there are no "holes"). The standard FIA hexagon system meets these
#' requirements and 3) the hexagons are all roughly oriented in the same direction.
#' 
#' (3) may not immediately be the case, but some projections may alleviate extreme
#' rotations.
neighbor_indexing <- function(tesselation) {
  # Step 0 - Build touching index and centroids
  centroids  <- gCentroid(tesselation, byid=TRUE)
  centroid_c     <- data.frame(centroids@coords)
  centroid_c$rix <- 1: nrow(centroid_c)
  touche_rix <- gTouches(tesselation, byid=TRUE, returnDense=FALSE)
  
  # Step 1 - Select the first hexagon that has 6 neighbors
  start_rix <- which(lengths(touche_rix) == 6)[[1]]
  start_c   <- centroids@coords[start_rix,]
  
  # Step 2  - Build a key that assigns angles to North/South/Northwest indices etc
  nbh_c  <- centroids@coords[touche_rix[[start_rix]],]
  dx <- nbh_c[,1] - start_c[[1]]
  dy <- nbh_c[,2] - start_c[[2]]
  
  # TODO used twice make a function
  theta <- atan2(dy, dx) * (180 / pi)
  theta[theta < 0] <- round(theta[theta < 0] + 360)
  theta <- round(theta/30) * 30
  
  theta_key <- data.frame(theta)
  theta_key <- theta_key[order(theta), ,drop=FALSE]
  theta_key$r <- c(0, -1, -1,  0,  1, 1)
  theta_key$c <- c(2,  1,  1, -2, -1, 1)
  
  to_visit <- c()
  visited  <- data.frame(matrix(NA, nrow=nrow(centroid_c), ncol=3))
  colnames(visited) <- c('rix', 'r', 'c')
  
  this_hex_rix <- start_rix
  n_visited <- 0
  
  while(n_visited < nrow(centroid_c)-1) {
    n_visited <- sum(!is.na(visited$rix))
    to_visit <- to_visit[2:length(to_visit)]
    nbh_c  <- centroid_c[touche_rix[[this_hex_rix]],]
    this_c <- centroid_c[this_hex_rix, , drop=FALSE]
    
    # Compute the angles of this neighborhood
    theta <- atan2(nbh_c[,2] - this_c[[2]], nbh_c[,1] - this_c[[1]]) * (180 / pi)
    theta[theta < 0] <- round(theta[theta < 0] + 360)
    theta <- round(theta/30) * 30
    nbh_c <- cbind(nbh_c, theta)
    
    nbh_v <- merge(nbh_c, visited, by='rix')
    visited[this_c$rix, 1] <- this_hex_rix
    
    if(n_visited == 0) {
      # This is the first point
      visited[this_c$rix, c('r', 'c')] <- c(0, 0)
    } else if(nrow(nbh_v) > 0) { # Then we have neighbors with thetas
      # Get the first non-NA neighbor
      non_na <- nbh_v[complete.cases(nbh_v[,c('r', 'c')]),]
      non_na <- non_na[1, ,drop=FALSE]
      
      # Now attach the theta key
      non_na <- merge(non_na, theta_key, by='theta', all.x=TRUE)
      visited[this_c$rix, c('r', 'c')] <- c(non_na$r.x - non_na$r.y, non_na$c.x - non_na$c.y)
    }
    
    # Prepare the next step, what indices do we need to visit?
    to_visit <- c(to_visit, nbh_c$rix[!(nbh_c$rix %in% to_visit | nbh_c$rix %in% visited[,1])])
    this_hex_rix <- to_visit[[1]]
  }
  
  # Standardize the index positions
  visited[,2] <- visited[,2] - min(visited[,2]) + 1
  visited[,3] <- visited[,3] - min(visited[,3]) + 1
  
  return(visited)
}