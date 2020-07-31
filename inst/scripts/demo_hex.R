library(ggplot2)
devtools::load_all()

pts <- readRDS('data/hex_pts.RDS')
sys_frame <- HexFrame(pts, attributes=c('z_1'))
p <- length(sys_frame@attributes)

a <- 2
starts <- subsample_starts(sys_frame, a)

N <- nrow(sys_frame)
N_check <- 0
mu <- c()
for(i in 1:nrow(starts)) {
  start <- starts[i,]
  subsamp <- subsample(sys_frame, start, a)
  print(head(subsamp@data))
  N_i <- nrow(subsamp)
  N_check <- N_check + N_i
  mu <- c(mu, mean(subsamp@data$z_1))
}

estimators <- list(
  'var_srs' = var_srs,
  'var_dorazio_i' = var_dorazio_i,
  'var_dorazio_c' = var_dorazio_c,
  'var_non_overlap' = var_non_overlap,
  'var_mat' = var_mat
)

a_vec   <- c(2:12)
results <- compare_estimators(sys_frame, a_vec, estimators)

ggplot(results) +
  geom_point(aes(x=a, y=variance, color=estimator)) +
  geom_line(aes(x=a,  y=variance), data=filter(results, estimator=='var_sys'))
