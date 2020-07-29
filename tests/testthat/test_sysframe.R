# This testing script tests all variance estimators that do not rely specifically on 
# HexFrame or RectFrame

data("hex_pts_small")
hex_pts_small <- SpatialPointsDataFrame(hex_pts_small[,c('s_1', 's_2')], data=hex_pts_small)
hf <- HexFrame(hex_pts_small, attributes = c('z_1', 'z_50', 'z_100'))
hf@a <- 3
hf@N <- 1000

test_that('VarSRS runs as expected', {
  v_no_fpc <- VarSRS(fpc=FALSE, diagnostic=FALSE)
  v_fpc <- VarSRS(fpc=TRUE, diagnostic=FALSE)
  v_diag <- VarSRS(diagnostic=TRUE)
  
  est_no_fpc <- v_no_fpc(hf)[[1]]
  expect_equal(0.09, round(est_no_fpc, 2))
  
  est_fpc <- v_fpc(hf)[[1]]
  expect_equal(0.09, round(est_fpc, 2))
  
  diagn <- v_diag(hf)
  expect_equal(diagn@n, 20)
  expect_equal(diagn@N, hf@N)
  expect_equal(round(diagn@mu[[1]], 2), 0.77)
})


test_that('VarSO runs as expected', {
  hf@data$pi_i <- 1/16
  v_no_fpc <- VarSO(fpc=FALSE, diagnostic=FALSE)
  v_fpc <- VarSO(fpc=TRUE, diagnostic=FALSE)
  v_diag <- VarSO(diagnostic=TRUE)
  
  est_no_fpc <- v_no_fpc(hf)[[1]]
  expect_equal(0.007, round(est_no_fpc, 3))
  
  est_fpc <- v_fpc(hf)[[1]]
  expect_equal(0.007, round(est_fpc, 3))
  
  diagn <- v_diag(hf)
  expect_equal(diagn@n, 20)
  expect_equal(diagn@N, hf@N)
  expect_equal(round(diagn@mu[[1]], 2), 0.77)
})

test_that('VarSYS runs as expected', {
  v <- VarSYS(hf, a = 2)
  est <- v(hf)[[1]]
  expect_equal(0.92, round(est, 2))
})