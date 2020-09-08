data("rect_pts_small")
rect_pts_small <- SpatialPointsDataFrame(rect_pts_small[,c('s_1', 's_2')], data=rect_pts_small)
rf <- RectFrame(rect_pts_small, attributes = c('z_1', 'z_50', 'z_100'))
rf@a <- 3
rf@N <- 1000

# Make the population and sample for the variance estimators
rf_pop <- SpatialPointsDataFrame(rect_pts[,c('s_1', 's_2')], rect_pts)
rf_pop <- HexFrame(rf_pop, attributes = c('z_1', 'z_50', 'z_100'))
rf_pop@data$pi <- 1/9
rf_pop@N <- nrow(rf_pop@data)
rf_subsamp <- subsample(rf_pop, c(1,1), 3)
rf_subsamp@N <- nrow(rf_pop@data)

greg_mapping <- list(
  'z_1'   = z_1 ~ x_1 - 1,
  'z_50'  = z_50 ~ x_1 - 1,
  'z_100' = z_100 ~ x_1 - 1
)

greg_rf <- greg(rf_subsamp, greg_mapping, rf_pop@data)
ht_rf <- horvitz_thompson(rf_subsamp)

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
  
  est_no_fpc <- v_no_fpc(greg_rf)[[1]]
  expect_equal(0.448, round(est_no_fpc, 3))
  
  est_fpc <- v_fpc(greg_rf)[[1]]
  expect_equal(0.423, round(est_fpc, 3))
  
  diagn <- v_diag(greg_rf)
  expect_equal(diagn@n, 50)
  expect_equal(diagn@N, greg_rf@N)
})

test_that('VarMAT works on RectFrame for paralellogram structure', {
  v_no_fpc <- VarMAT(fpc=FALSE, diagnostic=FALSE)
  v_fpc <- VarMAT(fpc=TRUE, diagnostic=FALSE)
  v_diag <- VarMAT(diagnostic=TRUE)
  
  est_no_fpc <- v_no_fpc(greg_rf)[[1]]
  expect_equal(0.731, round(est_no_fpc, 3))
  
  est_fpc <- v_fpc(greg_rf)[[1]]
  expect_equal(0.690, round(est_fpc, 3))
  
  diagn <- v_diag(greg_rf)
  expect_equal(diagn@n, 50)
  expect_equal(diagn@N, greg_rf@N)
})

test_that('VarDI works on RectFrame', {
  v_no_fpc <- VarDI(fpc=FALSE, diagnostic=FALSE)
  v_fpc <- VarDI(fpc=TRUE, diagnostic=FALSE)
  v_diag <- VarDI(diagnostic=TRUE)
  
  est_no_fpc <- v_no_fpc(greg_rf)[[1]]
  expect_equal(0.219, round(est_no_fpc, 3))
  
  est_fpc <- v_fpc(greg_rf)[[1]]
  expect_equal(0.207, round(est_fpc, 3))
  
  diagn <- v_diag(greg_rf)
  expect_equal(diagn@n, 50)
  expect_equal(diagn@N, greg_rf@N)
})

test_that('VarDC works on HexFrame', {
  v_no_fpc <- VarDC(fpc=FALSE, diagnostic=FALSE)
  v_fpc <- VarDC(fpc=TRUE, diagnostic=FALSE)
  v_diag <- VarDC(diagnostic=TRUE)
  
  est_no_fpc <- v_no_fpc(greg_rf)[[1]]
  expect_equal(0.660, round(est_no_fpc, 3))
  
  est_fpc <- v_fpc(greg_rf)[[1]]
  expect_equal(0.623, round(est_fpc, 3))
  
  diagn <- v_diag(greg_rf)
  expect_equal(diagn@n, 50)
  expect_equal(diagn@N, greg_rf@N)
})