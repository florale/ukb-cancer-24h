

# https://www.datanovia.com/en/blog/elegant-visualization-of-density-distribution-in-r-using-ridgeline/
(sleep_hist <- 
  ggplot(histd, aes(x = sleep, y = sleep_group)) +
  geom_density(color = "#999999", fill = "#FAF7F3", linetype = "dashed", linewidth = 0.5, alpha = 0.8) +
  # geom_vline(
  #   aes(xintercept = avg_sleep, group = sleep_group, colour = sleep_group),
  #   linetype = "dashed",
  #   linewidth = 1
  # ) +
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = 2, quantile_lines = TRUE
  ) +
  scale_x_continuous(breaks=c(0, 0.25, 0.75, 1)) + scale_fill_brewer(guide="none") +
  scale_colour_manual(values = col_hist) +
  labs(x = "Sleep period") +
  theme_minimal() +
  theme(
    panel.background  = element_blank(),
    panel.border      = element_blank(),
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.title.y      = element_blank(),
    # axis.text.y       = element_blank(),
    legend.position   = "none",
    plot.margin       = unit(c(0,0,0,2), "lines")
  )
)