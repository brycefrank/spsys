library(testthat)
data("hex_pts_small")

# Make a small hex population for testing subsampling and indices
hex_pts_small <- SpatialPointsDataFrame(hex_pts_small[,c('s_1', 's_2')], data=hex_pts_small)
hf <- HexFrame(hex_pts_small, attributes = c('z_1', 'z_50', 'z_100'))
hf@a <- 3
hf@N <- 1000


# Make the population and sample for the variance estimators
hf_pop <- SpatialPointsDataFrame(hex_pts[,c('s_1', 's_2')], hex_pts)
hf_pop <- HexFrame(hf_pop, attributes = c('z_1', 'z_50', 'z_100'))
hf_pop@data$pi <- 1/9
hf_pop@N <- nrow(hf_pop@data)
hf_subsamp <- subsample(hf_pop, c(1,1), 3)
hf_subsamp@N <- nrow(hf_pop@data)

greg_mapping <- list(
  'z_1'   = z_1 ~ x_1 - 1,
  'z_50'  = z_50 ~ x_1 - 1,
  'z_100' = z_100 ~ x_1 - 1
)

greg_hf <- greg(hf_subsamp, greg_mapping, hf_pop@data)
ht_hf <- horvitz_thompson(hf_subsamp)

#' Hexagonal indices should only ever have
#' odd rows aligned with odd columns and vice
#' versa. This expectation checks if this is true
#' for some set of indices.
expect_odd_even <- function(object) {
  act <- quasi_label(rlang::enquo(object), arg='object')
  
  ix <- act$val
  er_oc <- sum(ix$c[ix$r %% 2 == 0] %% 2 != 0)
  or_ec <- sum(ix$c[ix$r %% 2 != 0] %% 2 == 0)
  
  expect(
    (er_oc + or_ec) == 0,
    sprintf('This index had %i even-row odd-row matchings and %i odd-row even-row matchings,
            there should be 0 of each', er_oc, or_ec)
  )
  
}

test_that('HexFrame is properly indexed with default transform.', {
  ix <- hf@data[,c('r', 'c')]
  
  expect_odd_even(ix)
  expect_equal(nrow(ix), 14)
  
  expect_equal(ix[1,1], 4)
  expect_equal(ix[1,2], 2)
})

test_that('First order hexagon subsample returns valid index', {
  # a is even
  hf_subsamp_2 <- subsample(hf, c(1,1), 2)
  ix_2 <- hf_subsamp_2@data[,c('r', 'c')]
  expect_odd_even(ix_2)
  
  
  # a is odd
  hf_subsamp_3 <- subsample(hf, c(1,1), 3)
  ix_3 <- hf_subsamp_3@data[,c('r', 'c')]
  expect_odd_even(ix_3)
})

#' An input set of hexagonal indices can be arbitrary and must be standardized
#' by the instantiation of a HexFrame. The set of indices is first
#' translated to the (1,1) origin. Two important 'modes' could be provided by
#' the user. I refer to these as modes as alpha and beta respectively
#' 
#' An 'alpha' index is one where the minimum column and minimum row
#' index (provided by the user) cohere to an odd-odd index set. For example,  the user
#' inputs a index set that looks like:
#' 
#' (3,5)
#' (5,5)
#' 
#' Here the minimum row index is 3 and the minimum col index is 5. This is a point
#' that exists in our standard index set, and translation to (1,1) is elementary.
#' 
#' However, a user could provide an index set, 'beta', such that this is not the case:
#' 
#' (5,5)
#' (4,6)
#' 
#' These are valid points in our standard index set, but the minimum row index is
#' 4 and the minimum column index is 5 => (4,5). Such a situation needs special consideration
#' for translating to the (1,1) origin that involves adding 1 to the column indices
#' after standard translation.
#' 
#' We test for these two cases below.
test_that('HexFrame is properly indexed with "alpha" input index', {
  # Let's constrain the input points to only those points with
  # r >= 3 and c >= 5 like the above example
  hex_pts_sub <- hex_pts_small[hf$r >= 3 & hf$c >= 5 & hf$c <=8,]
  
  # Make the user input index
  hex_pts_ix <- data.frame(r=c(4,4,3,3) , c=c(6,8,5,7))
  
  hf_sub <- HexFrame(hex_pts_sub, index=hex_pts_ix)
  hf_sub_ix <- hf_sub@data[,c('r', 'c')]
  
  expect_odd_even(hf_sub_ix)
  
  # For alpha-mode indices the minimum row after translation is
  # 1 and its first column index must be 1
  first_row <- min(hf_sub_ix$r)
  first_col_of_first_row <- min(hf_sub_ix$c[hf_sub_ix$r == min(hf_sub_ix$r)])
  
  expect_equal(first_row, 1)
  expect_equal(first_col_of_first_row, 1)
})


test_that('HexFrame is properly indexed with "beta" input index', {
  ix <- hf@data[,c('r', 'c')]
  
  # Let's constrain the input points to only those points with
  # r >= 2 and c >= 5 like the above example
  hex_pts_sub <- hex_pts_small[hf$r >= 2 & hf$c >= 5 & hf$c <=6,]
  
  # Make the user input index
  hex_pts_ix <- data.frame(r=c(4,3,2) , c=c(6,5,6))
  
  hf_sub <- HexFrame(hex_pts_sub, index=hex_pts_ix)
  hf_sub_ix <- hf_sub@data[,c('r', 'c')]
  
  expect_odd_even(hf_sub_ix)
  
  # For beta-mode indices the minimum row after translation is
  # 1 and its first column index must be 3
  first_row <- min(hf_sub_ix$r)
  first_col_of_first_row <- min(hf_sub_ix$c[hf_sub_ix$r == min(hf_sub_ix$r)])
  
  expect_equal(first_row, 1)
  expect_equal(first_col_of_first_row, 3)
})

test_that('HT and GREG point estimators return the correct estimates', {
  # Check greg point estimates
  expect_equal(0.052, round(greg_hf@mu_hat$z_1, 3))
  expect_equal(0.111, round(greg_hf@mu_hat$z_50, 3))
  expect_equal(0.109, round(greg_hf@mu_hat$z_100, 3))
  
  # Check horvitz-thompson point estimates
  expect_equal(-1.467, round(ht_hf@mu_hat$z_1, 3))
  expect_equal(-1.475, round(ht_hf@mu_hat$z_50, 3))
  expect_equal(-5.919, round(ht_hf@mu_hat$z_100, 3))
})


test_that('VarNON works on HexFrame for triangular structure', {
  v_no_fpc <- VarNON(fpc=FALSE, diagnostic=FALSE, nbh='tri')
  v_fpc <- VarNON(fpc=TRUE, diagnostic=FALSE, nbh='tri')
  v_diag <- VarNON(diagnostic=TRUE, nbh='tri')
  
  # Check greg
  est_no_fpc <- v_no_fpc(greg_hf)[[1]]
  expect_equal(0.163, round(est_no_fpc, 3))
  
  est_fpc <- v_fpc(greg_hf)[[1]]
  expect_equal(0.145, round(est_fpc, 3))
  
  # Check horvitz
  
  diagn <- v_diag(greg_hf)
  expect_equal(diagn@n, 114)
  expect_equal(diagn@N, hf_pop@N)
})


test_that('VarNON works on HexFrame for paralellogram structure', {
  v_no_fpc <- VarNON(fpc=FALSE, diagnostic=FALSE, nbh='par')
  v_fpc <- VarNON(fpc=TRUE, diagnostic=FALSE, nbh='par')
  v_diag <- VarNON(diagnostic=TRUE, nbh='tri')
  
  est_no_fpc <- v_no_fpc(greg_hf)[[1]]
  expect_equal(0.194, round(est_no_fpc, 3))
  
  est_fpc <- v_fpc(greg_hf)[[1]]
  expect_equal(0.172, round(est_fpc, 3))
  
  diagn <- v_diag(greg_hf)
  expect_equal(diagn@n, 114)
  expect_equal(diagn@N, hf_pop@N)
})


test_that('VarNON works on HexFrame for hexagonal structure', {
  v_no_fpc <- VarNON(fpc=FALSE, diagnostic=FALSE, nbh='hex')
  v_fpc <- VarNON(fpc=TRUE, diagnostic=FALSE, nbh='hex')
  v_diag <- VarNON(diagnostic=TRUE, nbh='hex')
  
  est_no_fpc <- v_no_fpc(greg_hf)[[1]]
  expect_equal(0.250, round(est_no_fpc, 3))
  
  est_fpc <- v_fpc(greg_hf)[[1]]
  expect_equal(0.222, round(est_fpc, 3))
  
  diagn <- v_diag(greg_hf)
  expect_equal(diagn@n, 114)
  expect_equal(diagn@N, hf_pop@N)
})


test_that('VarMAT works on HexFrame for paralellogram structure', {
  v_no_fpc <- VarMAT(fpc=FALSE, diagnostic=FALSE, nbh='par')
  v_fpc <- VarMAT(fpc=TRUE, diagnostic=FALSE, nbh='par')
  v_diag <- VarMAT(diagnostic=TRUE, nbh='par')
  
  est_no_fpc <- v_no_fpc(greg_hf)[[1]]
  expect_equal(0.217, round(est_no_fpc, 3))
  
  est_fpc <- v_fpc(greg_hf)[[1]]
  expect_equal(0.193, round(est_fpc, 3))
  
  diagn <- v_diag(greg_hf)
  expect_equal(diagn@n, 114)
  expect_equal(diagn@N, hf_pop@N)
})

test_that('VarMAT works on HexFrame for hexagonal structure', {
  v_no_fpc <- VarMAT(fpc=FALSE, diagnostic=FALSE, nbh='hex')
  v_fpc <- VarMAT(fpc=TRUE, diagnostic=FALSE, nbh='hex')
  v_diag <- VarMAT(diagnostic=TRUE, nbh='hex')
  
  est_no_fpc <- v_no_fpc(hf)[[1]]
  expect_equal(0.067, round(est_no_fpc, 3))
  
  est_fpc <- v_fpc(hf)[[1]]
  expect_equal(0.066, round(est_fpc, 3))
  
  diagn <- v_diag(hf)
  expect_equal(diagn@n, 20)
  expect_equal(diagn@N, hf@N)
  expect_equal(round(diagn@mu[[1]], 2), 0.77)
})

test_that('VarDI works on HexFrame', {
  v_no_fpc <- VarDI(fpc=FALSE, diagnostic=FALSE)
  v_fpc <- VarDI(fpc=TRUE, diagnostic=FALSE)
  v_diag <- VarDI(diagnostic=TRUE)
  
  est_no_fpc <- v_no_fpc(hf)[[1]]
  expect_equal(0.089, round(est_no_fpc, 3))
  
  est_fpc <- v_fpc(hf)[[1]]
  expect_equal(0.087, round(est_fpc, 3))
  
  diagn <- v_diag(hf)
  expect_equal(diagn@n, 20)
  expect_equal(diagn@N, hf@N)
  expect_equal(round(diagn@mu[[1]], 2), 0.77)
})

test_that('VarDC works on HexFrame', {
  v_no_fpc <- VarDC(fpc=FALSE, diagnostic=FALSE)
  v_fpc <- VarDC(fpc=TRUE, diagnostic=FALSE)
  v_diag <- VarDC(diagnostic=TRUE)
  
  est_no_fpc <- v_no_fpc(hf)[[1]]
  expect_equal(0.1, round(est_no_fpc, 3))
  
  est_fpc <- v_fpc(hf)[[1]]
  expect_equal(0.098, round(est_fpc, 3))
  
  diagn <- v_diag(hf)
  expect_equal(diagn@n, 20)
  expect_equal(diagn@N, hf@N)
  expect_equal(round(diagn@mu[[1]], 2), 0.77)
})
