library(reshape2)
devtools::load_all()

pop_mat <- make_pop(100, 100)
pop_df  <- melt(pop_mat)
pop_df <- data.table(pop_df)
colnames(pop_df) <- c('x', 'y', 'z')

saveRDS(pop_mat, 'data/pop_mat.RDS')
saveRDS(pop_df,  'data/pop_df.RDS')
