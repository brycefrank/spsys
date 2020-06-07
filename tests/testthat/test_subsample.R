library(testthat)
hex_pts_small <- readRDS('data/hex_pts_small.RDS')


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
  hf <- HexFrame(hex_pts_small)
  ix <- hf@data[,c('r', 'c')]
  
  expect_odd_even(ix)
  expect_equal(nrow(ix), 15)
  
  expect_equal(ix[1,1], 5)
  expect_equal(ix[1,2], 1)
})



#' An input set of hexagonal indices can be arbitrary and must be standardized
#' by the instantiation of a HexFrame. The set of indices is first
#' translated to the (1,1) origin. Two important 'modes' could be provided by
#' the user, and spsys has to determine which is which. I refer to these as 
#' alpha and beta respectively
#' 
#' An 'alpha' index is one where the minimum column and minimum row
#' index (provided by the user) cohere to an odd-odd index set. For example,  the user
#' inputs a index set that looks like:
#' 
#' (3,5)
#' (5,5)
#' 
#' Here the minimum row index is 3 and the minimum col index is 5. This is a point
#' that exists in our standard index set.
#' 
#' However, a user could provide an index set, 'beta', such that this is not the case:
#' 
#' (5,5)
#' (6,4)
#' 
#' These are valid points in our standard index set, but the minimum row index is
#' 5 and the minimum column index is 4. Such a situation needs special consideration
#' for translating to the (1,1) origin that involves adding 1 to the column indices
#' after standard translation.
#' 
#' We test for these two cases here.
test_that('HexFrame is properly indexed with "alpha" input index', {
  hf <- HexFrame(hex_pts_small)
  ix <- hf@data[,c('r', 'c')]
  
  expect_odd_even(ix)
  expect_equal(nrow(ix), 15)
  
  expect_equal(ix[1,1], 5)
  expect_equal(ix[1,2], 1)
})


test_that('HexFrame is properly indexed with "beta" input index', {
  hf <- HexFrame(hex_pts_small)
  ix <- hf@data[,c('r', 'c')]
  
  expect_odd_even(ix)
  expect_equal(nrow(ix), 15)
  
  expect_equal(ix[1,1], 5)
  expect_equal(ix[1,2], 1)
})

