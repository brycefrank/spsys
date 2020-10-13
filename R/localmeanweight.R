#' This is a simpler version of the localmean.weight function used
#' by spsurvey to determine the neighborhood weights for the variance
#' computation. This version does not use a Moore-Penrose inversion to satisfy
#' a condition that the sum of the marginals of the weights matrix equals the
#' sample size.
localmean.weight.no_inv <- function(x, y, prb, nbh=4) {
   n <- length(x)

# Calculate indices of nearest neighbors
   dst <- as.matrix(dist(cbind(x, y), diag = TRUE, upper = TRUE))
   idx <- apply(dst, 2, order)[1:nbh,]

# Make neighbors symmetric
   jdx <- rep(1:n, rep(nbh, n))
   kdx <- unique(c((jdx - 1) * n + idx, (idx - 1) * n + jdx)) - 1
   ij <- cbind((kdx) %/% n + 1, (kdx) %% n + 1)
   ij <- ij[order(ij[, 1], dst[ij]),]

# Apply linear taper to the  inverse probability weights

   gct <- tabulate(ij[, 1])
   gwt <- numeric(0)
   for(i in 1:n)
      gwt <- c(gwt, 1 - (1:gct[i] - 1)/(gct[i]))
   gwt <- gwt/prb[ij[, 2]]

   smwt <- sapply(split(gwt, ij[, 1]), sum)
   gwt <- gwt/smwt[ij[, 1]]

   return(list(ij=ij, gwt=gwt))
}
