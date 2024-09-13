source("ukb-cancer-24h-utils.R")
source(paste0(redir, "ukb_utils.R"))
# source("ukb-cancer-24h-data.R")

# main model --------
fit_cancer_adj <- brmcoda(clr_cancer_acc,
                                          mvbind(ilr1, ilr2, ilr3) ~ cancer_before_acc +
                                            s(age_at_acc) + sex + white + working + edu + never_smoked + current_drinker + s(deprivation),
                                          # save_pars = save_pars(all = TRUE),
                                          warmup = 500, chains = 4, cores = 4, backend = "cmdstanr"
)
saveRDS(fit_cancer_adj, paste0(outputdir, "fit_cancer_adj", ".RDS"))

# predicted posteriors ------------
fit_cancer_adj <- readRDS(paste0(outputdir, "fit_cancer_adj", ".RDS"))

# reference grid
d_cancer_adj <- emmeans::ref_grid(fit_cancer_adj$model)@grid

# predict
pred_cancer_adj <- fitted(fit_cancer_adj, newdata = d_cancer_adj, scale = "response", summary = FALSE)

# summarise by cancer types
pred_cancer_adj <- apply(pred_cancer_adj, c(1), function(x)  cbind(d_cancer_adj, x))
pred_cancer_adj <- lapply(pred_cancer_adj, function(d) {
  d <- as.data.table(d)
  
  # estimated means by cancer types
  d[, sleep := mean(V1), by = cancer_before_acc]
  d[, mvpa  := mean(V2), by = cancer_before_acc]
  d[, lpa   := mean(V3), by = cancer_before_acc]
  d[, sb    := mean(V4), by = cancer_before_acc]
  
  # constrast cancer types vs healthy
  d[, sleep_constrast := sleep - d$sleep[1]]
  d[, mvpa_constrast := mvpa - d$mvpa[1]]
  d[, lpa_constrast := lpa - d$lpa[1]]
  d[, sb_constrast := sb - d$sb[1]]
  
  d <- d[, .(cancer_before_acc, 
             sleep, mvpa, lpa, sb,
             sleep_constrast, mvpa_constrast, lpa_constrast, sb_constrast
  )]
  d <- unique(d)
  d
})

# assemble back to summarise posteriors
pred_cancer_adj <- as.data.frame(abind(pred_cancer_adj, along = 1))
pred_cancer_adj <- split(pred_cancer_adj, pred_cancer_adj$cancer_before_acc)

## estimated means  ----------------------
pred_comp_cancer_adj <- lapply(pred_cancer_adj, function(l) {
  l <- as.data.frame(l[, c("sleep", "mvpa", "lpa", "sb")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, part = names(l))
  l <- rbindlist(l)
  l
})
pred_comp_cancer_adj <- Map(cbind, pred_comp_cancer_adj, cancer_before_acc = names(pred_comp_cancer_adj))
pred_comp_cancer_adj <- rbindlist(pred_comp_cancer_adj)

## contrasts --------------------
diff_comp_cancer_adj <- lapply(pred_cancer_adj, function(l) {
  l <- as.data.frame(l[, c("sleep_constrast", "mvpa_constrast", "lpa_constrast", "sb_constrast")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, constrast = names(l))
  l <- rbindlist(l)
  l
})
diff_comp_cancer_adj <- Map(cbind, diff_comp_cancer_adj, diff_from_healthy = names(diff_comp_cancer_adj))
diff_comp_cancer_adj <- rbindlist(diff_comp_cancer_adj)

setnames(diff_comp_cancer_adj, "Mean", "Mean_diff")
setnames(diff_comp_cancer_adj, "CI_low", "CI_low_diff")
setnames(diff_comp_cancer_adj, "CI_high", "CI_high_diff")

# all results  ------------------------
comp_cancer_adj <- cbind(
  pred_comp_cancer_adj[, .(Mean, CI_low, CI_high, part, cancer_before_acc)],
  diff_comp_cancer_adj[, .(Mean_diff, CI_low_diff, CI_high_diff)]
)

# add sig indicators
comp_cancer_adj[, nonsig := between(0, comp_cancer_adj$CI_low_diff, comp_cancer_adj$CI_high_diff)]
comp_cancer_adj[, Sig := ifelse(nonsig == FALSE, paste(intToUtf8(0x2217)), " ")]

comp_cancer_adj[, part := ifelse(part == "sleep", "Sleep", part)]
comp_cancer_adj[, part := ifelse(part == "mvpa", "Moderate-to-vigorous physical activity", part)]
comp_cancer_adj[, part := ifelse(part == "lpa", "Light physical activity", part)]
comp_cancer_adj[, part := ifelse(part == "sb", "Sedentary behaviour", part)]

comp_cancer_adj[, sig_position := min(CI_low), by = part]
comp_cancer_adj[, est_position := max(CI_high), by = part]

comp_cancer_adj[, estimates := paste0(round(Mean, 0), " [", round(CI_low, 0), ", ", round(CI_high, ), "]")]
comp_cancer_adj[, est_sig := paste0(estimates, " ", str_replace_na(Sig, ""))]

# plots -----------------------------
## facet all (olivia's thesis) -----------------------
(plot_comp_cancer_adj <- 
   ggplot(comp_cancer_adj, aes(x = cancer_before_acc, y = Mean, group = part)) +
   # geom_hline(aes(yintercept = yintercept), linewidth = 0.5, linetype = 2, colour = "#a8a8a8") +
   geom_pointrange(aes(ymin = CI_low,
                       ymax = CI_high, colour = cancer_before_acc)) +
   geom_text(aes(y = est_position + 5, label = Sig), 
             colour = "#978787", hjust = 1,
             size = 6, 
             # position = position_dodge2(width = 1),
             show.legend = FALSE) +
   facet_wrap(~part, scales = "free") +
   scale_colour_manual(values = pal) +
   labs(x = "", y = "", colour = "") +
   coord_flip() +      
   theme_ipsum() +
   theme(
     axis.ticks          = element_blank(),
     panel.background    = element_rect(fill = "transparent", colour = "black", linewidth = 0.5),
     plot.background     = element_rect(fill = "transparent", colour = NA),
     panel.grid.major.x  = element_blank(),
     panel.grid.minor    = element_blank(),
     strip.text          = element_text(size = 12, hjust = .5, face = "bold"),
     legend.position     = "none"
   )
)

grDevices::cairo_pdf(
  file = paste0(outputdir, "cancer_adj", ".pdf"),
  width = 8,
  height = 6,
)
plot_comp_cancer_adj
dev.off()

## plot by behaviour -----------------------
(plot_comp_cancer_sleep <- 
   ggplot(comp_cancer_adj[part == "Sleep"], aes(x = cancer_before_acc, y = Mean)) +
   # geom_hline(aes(yintercept = yintercept), linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8") +
   geom_pointrange(aes(ymin = CI_low,
                       ymax = CI_high, colour = cancer_before_acc), size = 0.5, linewidth = 0.75) +
   geom_text(aes(y = 625, label = est_sig),
             hjust = 1, nudge_x = 0, 
             family = "Arial Narrow", size = 4,
             show.legend = FALSE) +
   geom_text(aes(y = 525, label = cancer_before_acc),
             hjust = 1, nudge_x = 0, 
             family = "Arial Narrow", size = 4,
             show.legend = FALSE) +
   geom_segment(aes(x = 0, yend = 525), col = "black", linewidth = 0.5) +
   geom_segment(aes(x = 0, yend = 575), col = "black", linewidth = 0.5) +   
   scale_y_continuous(limits = c(475, 625),
                      breaks = c(525, 550, 575),
                      name = "Sleep") +
   scale_colour_manual(values = pal) +
   
   labs(x = "", y = "", colour = "") +
   coord_flip() +
   theme_ipsum() +
   theme(
     axis.ticks          = element_blank(),
     # panel.background    = element_rect(fill = "transparent", colour = "black", linewidth = 0.5),
     plot.background     = element_rect(fill = "transparent", colour = NA, linewidth = 0.5),
     panel.grid.major    = element_blank(),
     panel.grid.minor    = element_blank(),
     # axis.line.x         = element_line(linewidth = 0.5, colour = "black"),
     axis.title.x        = element_text(size = 13, face = "bold", hjust = .5),
     axis.text.x         = element_text(size = 12),
     axis.text.y         = element_blank(),
     strip.text          = element_text(size = 12, hjust = .5, face = "bold"),
     legend.text         = element_text(size = 13, face = "bold", hjust = .5),
     legend.position     = "none",
     plot.margin         = unit(c(0.5,0,1,0), "lines")
   )
)

(plot_comp_cancer_mvpa <- 
    ggplot(comp_cancer_adj[part == "Moderate-to-vigorous physical activity"], aes(x = cancer_before_acc, y = Mean)) +
    # geom_hline(aes(yintercept = yintercept), linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8") +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_before_acc), size = 0.5, linewidth = 0.75) +
    geom_text(aes(y = 40, label = est_sig),
              hjust = 1, nudge_x = 0, 
              family = "Arial Narrow", size = 4,
              show.legend = FALSE) +
    geom_text(aes(y = 20, label = cancer_before_acc),
              hjust = 1, nudge_x = 0, 
              family = "Arial Narrow", size = 4,
              show.legend = FALSE) +
    geom_segment(aes(x = 0, yend = 20), col = "black", linewidth = 0.5) +
    geom_segment(aes(x = 0, yend = 30), col = "black", linewidth = 0.5) +
    scale_y_continuous(limits = c(10, 40),
                       breaks = c(20, 25, 30),
                       name = "Moderate-to-vigorous physical activity") +
    scale_colour_manual(values = pal) +
    labs(x = "", y = "", colour = "") +
    coord_flip() +
    theme_ipsum() +
    theme(
      axis.ticks          = element_blank(),
      # panel.background    = element_rect(fill = "transparent", colour = "black", linewidth = 0.5),
      plot.background     = element_rect(fill = "transparent", colour = NA, linewidth = 0.5),
      panel.grid.major    = element_blank(),
      panel.grid.minor    = element_blank(),
      # axis.line.x         = element_line(linewidth = 0.5, colour = "black"),
      axis.title.x        = element_text(size = 13, face = "bold", hjust = .5),
      axis.text.x         = element_text(size = 12),
      axis.text.y         = element_blank(),
      strip.text          = element_text(size = 12, hjust = .5, face = "bold"),
      legend.text         = element_text(size = 13, face = "bold", hjust = .5),
      legend.position     = "none",
      plot.margin         = unit(c(0.5,0,1,0), "lines")
    )
)

(plot_comp_cancer_lpa <- 
    ggplot(comp_cancer_adj[part == "Light physical activity"], aes(x = cancer_before_acc, y = Mean)) +
    # geom_hline(aes(yintercept = yintercept), linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8") +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_before_acc), size = 0.5, linewidth = 0.75) +
    geom_text(aes(y = 375, label = est_sig),
              hjust = 1, nudge_x = 0, 
              family = "Arial Narrow", size = 4,
              show.legend = FALSE) +
    geom_text(aes(y = 275, label = cancer_before_acc),
              hjust = 1, nudge_x = 0, 
              family = "Arial Narrow", size = 4,
              show.legend = FALSE) +
    geom_segment(aes(x = 0, yend = 275), col = "black", linewidth = 0.5) +
    geom_segment(aes(x = 0, yend = 325), col = "black", linewidth = 0.5) +
    scale_y_continuous(limits = c(225, 375),
                       breaks = c(275, 300, 325),
                       name = "Light physical activity") +
    scale_colour_manual(values = pal) +
    labs(x = "", y = "", colour = "") +
    coord_flip() +
    theme_ipsum() +
    theme(
      axis.ticks          = element_blank(),
      # panel.background    = element_rect(fill = "transparent", colour = "black", linewidth = 0.5),
      plot.background     = element_rect(fill = "transparent", colour = NA, linewidth = 0.5),
      panel.grid.major    = element_blank(),
      panel.grid.minor    = element_blank(),
      # axis.line.x         = element_line(linewidth = 0.5, colour = "black"),
      axis.title.x        = element_text(size = 13, face = "bold", hjust = .5),
      axis.text.x         = element_text(size = 12),
      axis.text.y         = element_blank(),
      strip.text          = element_text(size = 12, hjust = .5, face = "bold"),
      legend.text         = element_text(size = 13, face = "bold", hjust = .5),
      legend.position     = "none",
      plot.margin         = unit(c(0.5,0,1,0), "lines")
    )
)

(plot_comp_cancer_sb <- 
    ggplot(comp_cancer_adj[part == "Sedentary behaviour"], aes(x = cancer_before_acc, y = Mean)) +
    # geom_hline(aes(yintercept = yintercept), linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8") +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_before_acc), size = 0.5, linewidth = 0.75) +
    geom_text(aes(y = 650, label = est_sig),
              hjust = 1, nudge_x = 0, 
              family = "Arial Narrow", size = 4,
              show.legend = FALSE) +
    geom_text(aes(y = 550, label = cancer_before_acc),
              hjust = 1, nudge_x = 0, 
              family = "Arial Narrow", size = 4,
              show.legend = FALSE) +
    geom_segment(aes(x = 0, yend = 600), col = "black", linewidth = 0.5) +
    geom_segment(aes(x = 0, yend = 550), col = "black", linewidth = 0.5) +
    scale_y_continuous(limits = c(500, 650),
                       breaks = c(550, 575, 600),
                       name = "Sedentary behaviour") +
    scale_colour_manual(values = pal) +
    labs(x = "", y = "", colour = "") +
    coord_flip() +
    theme_ipsum() +
    theme(
      axis.ticks          = element_blank(),
      # panel.background    = element_rect(fill = "transparent", colour = "black", linewidth = 0.5),
      plot.background     = element_rect(fill = "transparent", colour = NA, linewidth = 0.5),
      panel.grid.major    = element_blank(),
      panel.grid.minor    = element_blank(),
      # axis.line.x         = element_line(linewidth = 0.5, colour = "black"),
      axis.title.x        = element_text(size = 13, face = "bold", hjust = .5),
      axis.text.x         = element_text(size = 12),
      axis.text.y         = element_blank(),
      strip.text          = element_text(size = 12, hjust = .5, face = "bold"),
      legend.text         = element_text(size = 13, face = "bold", hjust = .5),
      legend.position     = "none",
      plot.margin         = unit(c(0.5,0,1,0), "lines")
    )
)

grDevices::cairo_pdf(
  file = paste0(outputdir, "cancer_est", ".pdf"),
  width = 8,
  height = 4,
)

ggarrange(
  plot_comp_cancer_mvpa,
  plot_comp_cancer_lpa,
  plot_comp_cancer_sb,
  plot_comp_cancer_sleep,
  nrow = 2, ncol = 2
)
dev.off()

