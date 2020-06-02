library(ggplot2)
devtools::load_all()
block_hex <- readRDS('data/hex_pts.RDS')

hex_frame <- HexFrame(block_hex, attributes=c('z_1'))
p <- length(hex_frame@attributes)

estimators <- list(
  'var_srs' = var_srs,
  'var_mat' = var_mat,
  'var_non_overlap'= var_non_overlap
)


a_vec   <- c(4:10)
results <- compare_estimators(hex_frame, a_vec, estimators)