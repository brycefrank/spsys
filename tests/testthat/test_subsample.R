library(testthat)

hex_pts_small <- readRDS('data/hex_pts_small.RDS')
hex_pts_small <- SpatialPointsDataFrame(hex_pts_small[,c('s_1', 's_2')], data=hex_pts_small)

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

test_that('First order hexagon subsample returns valid index', {
  hf <- HexFrame(hex_pts_small)
  
  # a is even
  hf_subsamp_2 <- subsample(hf, c(1,1), 2)
  ix_2 <- hf_subsamp_2@data[,c('r', 'c')]
  expect_odd_even(ix_2)
  
  
  # a is odd
  hf_subsamp_3 <- subsample(hf, c(1,1), 3)
  ix_3 <- hf_subsamp_3@data[,c('r', 'c')]
  expect_odd_even(ix_3)
})