data("rect_pts_small")
rect_pts_small <- SpatialPointsDataFrame(rect_pts_small[,c('s_1', 's_2')], data=rect_pts_small)
rf <- RectFrame(rect_pts_small, attributes = c('z_1', 'z_50', 'z_100'))
rf@a <- 3
rf@N <- 1000

test_that('RectFrame is properly indexed with default transform.', {
  ix <- rf@data[,c('r', 'c')]
  
  expect_equal(nrow(ix), 16)
  expect_equal(ix[1,1], 4)
  expect_equal(ix[1,2], 1)
})

test_that('VarNON works on HexFrame for paralellogram structure', {
  v_no_fpc <- VarNON(fpc=FALSE, diagnostic=FALSE)
  v_fpc <- VarNON(fpc=TRUE, diagnostic=FALSE)
  v_diag <- VarNON(diagnostic=TRUE)
  
  est_no_fpc <- v_no_fpc(rf)[[1]]
  expect_equal(0.124, round(est_no_fpc, 3))
  
  est_fpc <- v_fpc(rf)[[1]]
  expect_equal(0.122, round(est_fpc, 3))
  
  diagn <- v_diag(rf)
  expect_equal(diagn@n, 16)
  expect_equal(diagn@N, rf@N)
  expect_equal(round(diagn@mu[[1]], 2), 1.22)
})

test_that('VarMAT works on RectFrame for paralellogram structure', {
  v_no_fpc <- VarMAT(fpc=FALSE, diagnostic=FALSE)
  v_fpc <- VarMAT(fpc=TRUE, diagnostic=FALSE)
  v_diag <- VarMAT(diagnostic=TRUE)
  
  est_no_fpc <- v_no_fpc(rf)[[1]]
  expect_equal(0.381, round(est_no_fpc, 3))
  
  est_fpc <- v_fpc(rf)[[1]]
  expect_equal(0.375, round(est_fpc, 3))
  
  diagn <- v_diag(rf)
  expect_equal(diagn@n, 16)
  expect_equal(diagn@N, rf@N)
  expect_equal(round(diagn@mu[[1]], 2), 1.22)
})

test_that('VarDI works on RectFrame', {
  v_no_fpc <- VarDI(fpc=FALSE, diagnostic=FALSE)
  v_fpc <- VarDI(fpc=TRUE, diagnostic=FALSE)
  v_diag <- VarDI(diagnostic=TRUE)
  
  est_no_fpc <- v_no_fpc(rf)[[1]]
  expect_equal(0.156, round(est_no_fpc, 3))
  
  est_fpc <- v_fpc(rf)[[1]]
  expect_equal(0.153, round(est_fpc, 3))
  
  diagn <- v_diag(rf)
  expect_equal(diagn@n, 16)
  expect_equal(diagn@N, rf@N)
  expect_equal(round(diagn@mu[[1]], 2), 1.22)
})

test_that('VarDC works on HexFrame', {
  v_no_fpc <- VarDC(fpc=FALSE, diagnostic=FALSE)
  v_fpc <- VarDC(fpc=TRUE, diagnostic=FALSE)
  v_diag <- VarDC(diagnostic=TRUE)
  
  est_no_fpc <- v_no_fpc(rf)[[1]]
  expect_equal(0.167, round(est_no_fpc, 3))
  
  est_fpc <- v_fpc(rf)[[1]]
  expect_equal(0.164, round(est_fpc, 3))
  
  diagn <- v_diag(rf)
  expect_equal(diagn@n, 16)
  expect_equal(diagn@N, rf@N)
  expect_equal(round(diagn@mu[[1]], 2), 1.22)
})