# A script to index arbitrarily "warped" sets of hexagon points
# this arose because the hexagonal tesselation used by the FIA
# is warped over geographic space, which can break "naive"
# indexing operations that use linear transformations.

# The general approach is to recursively map neighboring hexagons
# to the index set, which allows for some "shifting" between points
# without breaking across large distances.

# The requirements are that: 1) the hexagons are all approximately the same
# size (i.e. the hexagon points are approximately evenly spaced) and 
# 2) that every hexagon that needs to be indexed has at least one neighbor
# (i.e. there are no "holes"). The standard FIA hexagon system meets these
# requirements.
library(rgdal)
library(rgeos)

# Some formalisms. Anything ending in _c is a set of coordinates of hexagonal centroids
# anything ending in _ix is the desired 2dimensional index
# anything ending in _rix is the 'raw' one dimensional index of the original list of hexagons

block_hex <- readRDS('data/block_hex.RDS')

# Step 0 - Build touching index and centroids
touche_rix <- gTouches(block_hex, byid=TRUE, returnDense=FALSE)
centroids <- gCentroid(block_hex, byid=TRUE)


# Step 1 - Select the first hexagon that has 6 neighbors
start_rix <- which(lengths(touche_rix) == 6)[[1]]
start_c  <- centroids@coords[start_rix,]

# Step 2  - Determine the transformation angle. This will be used to rotate the set of hexagon centroids to a standard orientation
# Step 2a - Consider only those neighbors whose x coordinates are leq the starting x

# Step 2b - Of these coordinates, consider the one with the maximum y

# Step 2c - compute the angle between "north" and this point
dx <- top_c[[1]] - start_c[[1]]
dy <- top_c[[2]] - start_c[[2]]
theta_radians <- atan2(dy, dx) - (pi/2)

# Step 3 - Make a matrix to store the hexagonal indices, set the starting hexagon to 0,0
hex_ix <- matrix(NA, nrow=nrow(centroids@coords), ncol=2)
hex_ix[start_rix,] <- c(0,0)

# Step 4 - The Big One - Recursively Index the neighbors
recurse <- function(centroids, touche_ix, hex_ix, this_hex_rix, this_hex_ix) {
  print(this_hex_ix)
  # Find the neighbors of start_ix
  nbh_c <- centroids@coords[touche_ix[[this_hex_rix]],]
  print(nbh_c)
  this_c <- centroids@coords[this_hex_rix, , drop=FALSE]
  
  # Get neighbors to left of this point, order descending by y
  left_c   <- nbh_c[(nbh_c[,1] <= this_c[1]), ,drop=FALSE]
  left_c   <- left_c[order(left_c[,2], decreasing=TRUE),]
  
  # The first element is the "north" coordinate
  n_c   <- left_c[1, ,drop=FALSE]
  n_rix <- as.numeric(row.names(n_c))
  n_ix  <- c(this_hex_ix[[1]] - 2, this_hex_ix[[2]])
  hex_ix[n_rix, ] <- n_ix
  
  # The second element is the "northwest" coordinate
  nw_c  <- left_c[2, ,drop=FALSE]
  nw_rix <- as.numeric(row.names(nw_c))
  nw_ix <- c(this_hex_ix[[1]] - 1, this_hex_ix[[2]] - 1)
  hex_ix[nw_rix, ] <- nw_ix
  
  # The third element is the "southwest" coordinate
  sw_c  <- left_c[3, ,drop=FALSE]
  sw_rix <- as.numeric(row.names(nw_c))
  sw_ix <- c(this_hex_ix[[1]] + 1, this_hex_ix[[2]] - 1)
  hex_ix[sw_rix, ] <- sw_ix
  
  # Get neighbors to right of this point, order descending by y
  right_c   <- nbh_c[(nbh_c[,1] > this_c[1]), ,drop=FALSE]
  right_c   <- right_c[order(right_c[,2], decreasing=TRUE),]
  
  # The first element is the "northeast" coordinate
  ne_c   <- right_c[1, ,drop=FALSE]
  ne_rix <- as.numeric(row.names(ne_c))
  ne_ix  <- c(this_hex_ix[[1]] - 1, this_hex_ix[[2]] + 1)
  hex_ix[ne_rix, ] <- ne_ix
  
  # The second element is the "southeast" coordinate
  se_c  <- right_c[2, ,drop=FALSE]
  se_rix <- as.numeric(row.names(se_c))
  se_ix <- c(this_hex_ix[[1]] + 1, this_hex_ix[[2]] + 1)
  hex_ix[se_rix, ] <- se_ix
  
  # The third element is the "south" coordinate
  s_c  <- left_c[3, ,drop=FALSE]
  s_rix <- as.numeric(row.names(nw_c))
  s_ix <- c(this_hex_ix[[1]] + 2, this_hex_ix[[2]])
  hex_ix[s_rix, ] <- s_ix
  
  print(hex_ix)
  recurse(centroids, touche_ix, hex_ix, n_rix, n_ix)
  
}

recurse(centroids, touche_rix, hex_ix, start_rix, c(0,0))
