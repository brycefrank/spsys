data("hex_pts")

# Make a small hex population for testing subsampling and indices
hex_pts <- SpatialPointsDataFrame(hex_pts[,c('s_1', 's_2')], data=hex_pts)
hf_pop <- HexFrame(hex_pts, attributes = c('z_1', 'z_50', 'z_100'))
hf_pop@a <- 3
hf_pop@N <- nrow(hf_pop@data)

greg_mapping <- list(
  z_1   ~ 1,
  z_50  ~ 1,
  z_100 ~ 1
)

var_estimators <- list(
  'var_srs' = VarSRS(fpc=TRUE, diagnostic=FALSE),
  'var_non_hex' = VarNON(fpc=TRUE, diagnostic=FALSE, nbh='hex')
)

test_that('Compare estimators returns correct information', {
  a <- 4
  p <- length(hf_pop@attributes)
  assessment <- compare_estimators(hf_pop, c(a), var_estimators, greg_mapping)
  
  expect_equal(nrow(assessment$pop_frame), p)
  expect_equal(nrow(assessment$sys_frame), p)
  expect_equal(nrow(assessment$est_frame), (a^2)*length(var_estimators) * p)
})

