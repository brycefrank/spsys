devtools::load_all()


pop_mat <- make_pop(100, 100)
pop_df  <- melt(pop_mat)


grd <- samp_tri_grid_pi(pop_mat, pop_df, strat_mat)

