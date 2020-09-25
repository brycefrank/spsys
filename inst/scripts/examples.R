# Examples given in the paper
devtools::load_all()
data("hex_pts")

hex_spdf <- SpatialPointsDataFrame(
  coords=hex_pts[,c('s_1', 's_2')], 
  data=hex_pts[,c('z_1', 'z_50', 'z_100')]
)

# Create a HexFrame
hf <- HexFrame(hex_spdf, attributes = c('z_1', 'z_50', 'z_100'))

# Prepare the HexFrame for estimation
mapping <- list(
  z_1   ~ x_1 - 1,
  z_50  ~ x_1 - 1,
  z_100 ~ x_1 - 1
)

hf_greg <- greg(hf, mapping, hf@data)
var_estimators <- list(
  'var_srs' = VarSRS(fpc=TRUE, diagnostic=FALSE),
  'var_non_hex' = VarNON(fpc=TRUE, diagnostic=FALSE, nbh='hex'),
  'var_non_par' = VarNON(fpc=TRUE, diagnostic=FALSE, nbh='mat'),
  'var_mat_par' = VarMAT(fpc=TRUE, diagnostic=FALSE, nbh='mat')
)

assessment <- compare_estimators(hf, 4, estimators)
var_ratios <- merge(assessment$est_frame, assessment$sys_frame, 
                    by=c('attribute', 'a'), suffix=c('_est', '_sys'))
var_ratios$ratio <- var_ratios$variance_est / var_ratios$variance_sys

boxplot_fig <- function(var_ratios) {
  var_ratios$attribute <- factor(var_ratios$attribute, levels=c('z_1', 'z_50', 'z_100'))
  levels(var_ratios$attribute) <- c(TeX('$\\phi = 1$'), TeX('$\\phi = 50$'), TeX('$\\phi = 100$'))
  
  ggplot(var_ratios) +
    geom_boxplot(aes(x=estimator_est, y=ratio, color = estimator_est)) +
    scale_color_discrete(name='Estimator') +
    geom_abline(slope=0, intercept=1) +
    ylim(c(0,3)) +
    facet_wrap(~attribute, labeller=label_parsed) +
    theme(strip.text.x = element_text(size = 16),
          axis.text.x = element_blank(), legend.position='bottom',
          axis.ticks.x = element_blank(), axis.title.x=element_blank())
}

fig <- boxplot_fig(var_ratios)
ggsave('vignettes/fig_3.png', plot = fig, units='in', height = 5.5, width=7, dpi=300)

mses <- var_ratios %>%
  group_by(estimator_est, attribute) %>%
  summarize(mse = mean((variance_est - variance_sys)^2)) %>%
  pivot_wider(names_from = 'attribute', values_from='mse')

v_mat <- VarMAT(fpc=TRUE, diagnostic=TRUE, nbh='mat')
samp <- subsample(hf, c(1,1), 5)
v_out <- v_mat(samp)


setGeneric('var_sup', function(sys_frame, ...) {
  standardGeneric('var_sup')
})

setMethod('var_sup', signature(sys_frame = 'SysFrame'), 
          function(sys_frame, z_bar_r, K) {
            att_df  <- sys_frame@data[,sys_frame@attributes, drop=FALSE]
            z_bar_s <- colMeans(att_df)
            
            v_sup <- ((K - 1) / (2*K))^2 * (z_bar_s - z_bar_r)^2
            return(v_sup)
          }
)

not_sampled <- hf@data[!row.names(hf@data) %in% row.names(samp@data), ]
srs_wor_samp_ix <- sample(1:nrow(not_sampled), 20)
srs_wor_samp <- not_sampled[srs_wor_samp_ix,]

z_bar_r <- colMeans(srs_wor_samp[, samp@attributes])
var_sup(samp, z_bar_r, 16)