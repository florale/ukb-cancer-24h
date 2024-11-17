source("ukb-cancer-24h-utils.R")
source(paste0(redir, "ukb_utils.R"))
source("ukb-cancer-24h-data.R")

# MAIN MODEL CANCER VS HEALTHY --------
fit_cancer_adj <- brmcoda(clr_cancer_acc,
                          mvbind(ilr1, ilr2, ilr3) ~ cancer_before_acc +
                            s(age_at_acc) + 
                            # other_conds_at_acc +
                            sex + white + working + edu + never_smoked + current_drinker + s(deprivation),
                          # save_pars = save_pars(all = TRUE),
                          warmup = 500, chains = 4, cores = 4, backend = "cmdstanr"
)
saveRDS(fit_cancer_adj, paste0(outputdir, "fit_cancer_adj", ".RDS"))

# predicted posteriors ------------
# fit_cancer_adj <- readRDS(paste0(outputdir, "fit_cancer_adj", ".RDS"))

# reference grid
d_cancer_adj <- emmeans::ref_grid(fit_cancer_adj$model)@grid

# predict
pred_cancer_adj <- fitted(fit_cancer_adj, newdata = d_cancer_adj, scale = "response", summary = FALSE)

# summarise by cancer types
pred_cancer_adj <- apply(pred_cancer_adj, c(1), function(x)  cbind(d_cancer_adj, x))
pred_cancer_adj <- lapply(pred_cancer_adj, function(d) {
  d <- as.data.table(d)
  
  # estimated means by cancer groups
  d[, sleep := mean(V1), by = cancer_before_acc]
  d[, mvpa  := mean(V2), by = cancer_before_acc]
  d[, lpa   := mean(V3), by = cancer_before_acc]
  d[, sb    := mean(V4), by = cancer_before_acc]
  
  # constrast cancer vs healthy
  d[, sleep_constrast := sleep - d[cancer_other_before_acc == "Healthy"]$sleep[1]]
  d[, mvpa_constrast_ := mvpa - d[cancer_other_before_acc == "Healthy"]$mvpa[1]]
  d[, lpa_constrast := lpa - d[cancer_other_before_acc == "Healthy"]$lpa[1]]
  d[, sb_constrast := sb - d[cancer_other_before_acc == "Healthy"]$sb[1]]
  
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

comp_cancer_adj[, part := ifelse(part == "sleep", "Sleep period", part)]
comp_cancer_adj[, part := ifelse(part == "mvpa", "Moderate-to-vigorous physical activity", part)]
comp_cancer_adj[, part := ifelse(part == "lpa", "Light physical activity", part)]
comp_cancer_adj[, part := ifelse(part == "sb", "Sedentary behaviour", part)]

comp_cancer_adj[, sig_position := min(CI_low), by = part]
comp_cancer_adj[, est_position := max(CI_high), by = part]

comp_cancer_adj[, estimates := paste0(round(Mean, 0), " [", round(CI_low, 0), ", ", round(CI_high, ), "]")]
comp_cancer_adj[, est_sig := paste0(estimates, " ", str_replace_na(Sig, ""))]

# plots -----------------------------
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



# MAIN MODEL CANCER VS OTHERS VS HEALTHY --------
fit_cancer_other_adj <- brmcoda(clr_cancer_acc,
                                mvbind(ilr1, ilr2, ilr3) ~ cancer_other_before_acc +
                                  s(age_at_acc) + sex + white + working + edu + never_smoked + current_drinker + s(deprivation),
                                # save_pars = save_pars(all = TRUE),
                                warmup = 500, chains = 4, cores = 4, backend = "cmdstanr"
)
saveRDS(fit_cancer_other_adj, paste0(outputdir, "fit_cancer_other_adj", ".RDS"))

# predicted posteriors ------------
fit_cancer_other_adj <- readRDS(paste0(outputdir, "fit_cancer_other_adj", ".RDS"))

# reference grid
d_cancer_other_adj <- emmeans::ref_grid(fit_cancer_other_adj$model)@grid

# predict
pred_cancer_other_adj <- fitted(fit_cancer_other_adj, newdata = d_cancer_other_adj, scale = "response", summary = FALSE)

# summarise by cancer types
pred_cancer_other_adj <- apply(pred_cancer_other_adj, c(1), function(x)  cbind(d_cancer_other_adj, x))
pred_cancer_other_adj <- lapply(pred_cancer_other_adj, function(d) {
  d <- as.data.table(d)
  
  # estimated means by cancer groups
  d[, sleep := mean(V1), by = cancer_other_before_acc]
  d[, mvpa  := mean(V2), by = cancer_other_before_acc]
  d[, lpa   := mean(V3), by = cancer_other_before_acc]
  d[, sb    := mean(V4), by = cancer_other_before_acc]
  
  # constrast cancer vs healthy
  d[, sleep_vs_healthy := sleep - d[cancer_other_before_acc == "Healthy"]$sleep[1]]
  d[, mvpa_vs_healthy := mvpa - d[cancer_other_before_acc == "Healthy"]$mvpa[1]]
  d[, lpa_vs_healthy := lpa - d[cancer_other_before_acc == "Healthy"]$lpa[1]]
  d[, sb_vs_healthy := sb - d[cancer_other_before_acc == "Healthy"]$sb[1]]
  
  # vs cancer vs other conds
  d[, sleep_vs_others := sleep - d[cancer_other_before_acc == "Others"]$sleep[1]]
  d[, mvpa_vs_others := mvpa - d[cancer_other_before_acc == "Others"]$mvpa[1]]
  d[, lpa_vs_others := lpa - d[cancer_other_before_acc == "Others"]$lpa[1]]
  d[, sb_vs_others := sb - d[cancer_other_before_acc == "Others"]$sb[1]]
  
  d <- d[, .(cancer_other_before_acc, 
             sleep, mvpa, lpa, sb,
             sleep_vs_healthy, mvpa_vs_healthy, lpa_vs_healthy, sb_vs_healthy,
             sleep_vs_others, mvpa_vs_others, lpa_vs_others, sb_vs_others
  )]
  d <- unique(d)
  d
})

# assemble back to summarise posteriors
pred_cancer_other_adj <- as.data.frame(abind(pred_cancer_other_adj, along = 1))
pred_cancer_other_adj <- split(pred_cancer_other_adj, pred_cancer_other_adj$cancer_other_before_acc)

## estimated means  ----------------------
pred_comp_cancer_other_adj <- lapply(pred_cancer_other_adj, function(l) {
  l <- as.data.frame(l[, c("sleep", "mvpa", "lpa", "sb")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, part = names(l))
  l <- rbindlist(l)
  l
})
pred_comp_cancer_other_adj <- Map(cbind, pred_comp_cancer_other_adj, cancer_other_before_acc = names(pred_comp_cancer_other_adj))
pred_comp_cancer_other_adj <- rbindlist(pred_comp_cancer_other_adj)

## contrasts --------------------
### vs healthy
diff1_comp_cancer_other_adj <- lapply(pred_cancer_other_adj, function(l) {
  l <- as.data.frame(l[, c("sleep_vs_healthy", "mvpa_vs_healthy", "lpa_vs_healthy", "sb_vs_healthy")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, constrast = names(l))
  l <- rbindlist(l)
  l
})
diff1_comp_cancer_other_adj <- Map(cbind, diff1_comp_cancer_other_adj, diff_from_healthy = names(diff1_comp_cancer_other_adj))
diff1_comp_cancer_other_adj <- rbindlist(diff1_comp_cancer_other_adj)

setnames(diff1_comp_cancer_other_adj, "Mean", "Mean_diff_ref_healthy")
setnames(diff1_comp_cancer_other_adj, "CI_low", "CI_low_diff_ref_healthy")
setnames(diff1_comp_cancer_other_adj, "CI_high", "CI_high_diff_ref_healthy")

### vs others
diff2_comp_cancer_other_adj <- lapply(pred_cancer_other_adj, function(l) {
  l <- as.data.frame(l[, c("sleep_vs_others", "mvpa_vs_others", "lpa_vs_others", "sb_vs_others")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, constrast = names(l))
  l <- rbindlist(l)
  l
})
diff2_comp_cancer_other_adj <- Map(cbind, diff2_comp_cancer_other_adj, diff_from_others = names(diff2_comp_cancer_other_adj))
diff2_comp_cancer_other_adj <- rbindlist(diff2_comp_cancer_other_adj)

setnames(diff2_comp_cancer_other_adj, "Mean", "Mean_diff_ref_others")
setnames(diff2_comp_cancer_other_adj, "CI_low", "CI_low_diff_ref_others")
setnames(diff2_comp_cancer_other_adj, "CI_high", "CI_high_diff_ref_others")

# all results  ------------------------
comp_cancer_other_adj <- cbind(
  pred_comp_cancer_other_adj[, .(Mean, CI_low, CI_high, part, cancer_other_before_acc)],
  diff1_comp_cancer_other_adj[, .(Mean_diff_ref_healthy, CI_low_diff_ref_healthy, CI_high_diff_ref_healthy)],
  diff2_comp_cancer_other_adj[, .(Mean_diff_ref_others, CI_low_diff_ref_others, CI_high_diff_ref_others)]
)

# add sig indicators
comp_cancer_other_adj[, nonsig_healthy := between(0, comp_cancer_other_adj$CI_low_diff_ref_healthy, comp_cancer_other_adj$CI_high_diff_ref_healthy)]
comp_cancer_other_adj[, nonsig_others := between(0, comp_cancer_other_adj$CI_low_diff_ref_others, comp_cancer_other_adj$CI_high_diff_ref_others)]

comp_cancer_other_adj[, sig_ref_healthy := ifelse(nonsig_healthy == FALSE, "$^a$", "$\\phantom{a}$")]
comp_cancer_other_adj[, sig_ref_others := ifelse(nonsig_others == FALSE, "$^b$", "$\\phantom{b}$")]

comp_cancer_other_adj[, yintercept_healthy := NA]
comp_cancer_other_adj[, yintercept_healthy := ifelse(part == "sleep", comp_cancer_other_adj[cancer_other_before_acc == "Healthy" & part == "sleep"]$Mean, yintercept_healthy)]
comp_cancer_other_adj[, yintercept_healthy := ifelse(part == "mvpa", comp_cancer_other_adj[cancer_other_before_acc == "Healthy" & part == "mvpa"]$Mean, yintercept_healthy)]
comp_cancer_other_adj[, yintercept_healthy := ifelse(part == "lpa", comp_cancer_other_adj[cancer_other_before_acc == "Healthy" & part == "lpa"]$Mean, yintercept_healthy)]
comp_cancer_other_adj[, yintercept_healthy := ifelse(part == "sb", comp_cancer_other_adj[cancer_other_before_acc == "Healthy" & part == "sb"]$Mean, yintercept_healthy)]

comp_cancer_other_adj[, ci_low_healthy := NA]
comp_cancer_other_adj[, ci_low_healthy := ifelse(part == "sleep", comp_cancer_other_adj[cancer_other_before_acc == "Healthy" & part == "sleep"]$CI_low, ci_low_healthy)]
comp_cancer_other_adj[, ci_low_healthy := ifelse(part == "mvpa", comp_cancer_other_adj[cancer_other_before_acc == "Healthy" & part == "mvpa"]$CI_low, ci_low_healthy)]
comp_cancer_other_adj[, ci_low_healthy := ifelse(part == "lpa", comp_cancer_other_adj[cancer_other_before_acc == "Healthy" & part == "lpa"]$CI_low, ci_low_healthy)]
comp_cancer_other_adj[, ci_low_healthy := ifelse(part == "sb", comp_cancer_other_adj[cancer_other_before_acc == "Healthy" & part == "sb"]$CI_low, ci_low_healthy)]

comp_cancer_other_adj[, ci_high_healthy := NA]
comp_cancer_other_adj[, ci_high_healthy := ifelse(part == "sleep", comp_cancer_other_adj[cancer_other_before_acc == "Healthy" & part == "sleep"]$CI_high, ci_high_healthy)]
comp_cancer_other_adj[, ci_high_healthy := ifelse(part == "mvpa", comp_cancer_other_adj[cancer_other_before_acc == "Healthy" & part == "mvpa"]$CI_high, ci_high_healthy)]
comp_cancer_other_adj[, ci_high_healthy := ifelse(part == "lpa", comp_cancer_other_adj[cancer_other_before_acc == "Healthy" & part == "lpa"]$CI_high, ci_high_healthy)]
comp_cancer_other_adj[, ci_high_healthy := ifelse(part == "sb", comp_cancer_other_adj[cancer_other_before_acc == "Healthy" & part == "sb"]$CI_high, ci_high_healthy)]

comp_cancer_other_adj[, yintercept_others := NA]
comp_cancer_other_adj[, yintercept_others := ifelse(part == "sleep", comp_cancer_other_adj[cancer_other_before_acc == "Others" & part == "sleep"]$Mean, yintercept_others)]
comp_cancer_other_adj[, yintercept_others := ifelse(part == "mvpa", comp_cancer_other_adj[cancer_other_before_acc == "Others" & part == "mvpa"]$Mean, yintercept_others)]
comp_cancer_other_adj[, yintercept_others := ifelse(part == "lpa", comp_cancer_other_adj[cancer_other_before_acc == "Others" & part == "lpa"]$Mean, yintercept_others)]
comp_cancer_other_adj[, yintercept_others := ifelse(part == "sb", comp_cancer_other_adj[cancer_other_before_acc == "Others" & part == "sb"]$Mean, yintercept_others)]

comp_cancer_other_adj[, ci_low_others := NA]
comp_cancer_other_adj[, ci_low_others := ifelse(part == "sleep", comp_cancer_other_adj[cancer_other_before_acc == "Others" & part == "sleep"]$CI_low, ci_low_others)]
comp_cancer_other_adj[, ci_low_others := ifelse(part == "mvpa", comp_cancer_other_adj[cancer_other_before_acc == "Others" & part == "mvpa"]$CI_low, ci_low_others)]
comp_cancer_other_adj[, ci_low_others := ifelse(part == "lpa", comp_cancer_other_adj[cancer_other_before_acc == "Others" & part == "lpa"]$CI_low, ci_low_others)]
comp_cancer_other_adj[, ci_low_others := ifelse(part == "sb", comp_cancer_other_adj[cancer_other_before_acc == "Others" & part == "sb"]$CI_low, ci_low_others)]

comp_cancer_other_adj[, ci_high_others := NA]
comp_cancer_other_adj[, ci_high_others := ifelse(part == "sleep", comp_cancer_other_adj[cancer_other_before_acc == "Others" & part == "sleep"]$CI_high, ci_high_others)]
comp_cancer_other_adj[, ci_high_others := ifelse(part == "mvpa", comp_cancer_other_adj[cancer_other_before_acc == "Others" & part == "mvpa"]$CI_high, ci_high_others)]
comp_cancer_other_adj[, ci_high_others := ifelse(part == "lpa", comp_cancer_other_adj[cancer_other_before_acc == "Others" & part == "lpa"]$CI_high, ci_high_others)]
comp_cancer_other_adj[, ci_high_others := ifelse(part == "sb", comp_cancer_other_adj[cancer_other_before_acc == "Others" & part == "sb"]$CI_high, ci_high_others)]

comp_cancer_other_adj[, cancer_other_before_acc := factor(cancer_other_before_acc, ordered = TRUE,
                                                          levels = c(
                                                            "Healthy",
                                                            "Others",
                                                            "Cancer"
                                                          ))]

comp_cancer_other_adj[, part := ifelse(part == "sleep", "Sleep period", part)]
comp_cancer_other_adj[, part := ifelse(part == "mvpa", "Moderate-to-vigorous physical activity", part)]
comp_cancer_other_adj[, part := ifelse(part == "lpa", "Light physical activity", part)]
comp_cancer_other_adj[, part := ifelse(part == "sb", "Sedentary behaviour", part)]

comp_cancer_other_adj[, sig_position := min(CI_low), by = part]
comp_cancer_other_adj[, est_position := max(CI_high), by = part]

comp_cancer_other_adj[, estimates := paste0(round(Mean, 0), "[", round(CI_low, 0), ", ", round(CI_high, 0), "]")]
comp_cancer_other_adj[, est_sig := paste0(estimates, " ", sig_ref_healthy, sig_ref_others)]
# comp_cancer_other_adj[, est_sig := paste0(estimates, " ", str_replace_na(sig_ref_healthy, " "), str_replace_na(sig_ref_others, " "))]

## plot by behaviour -----------------------
(plot_comp_cancer_sleep <- 
   ggplot(comp_cancer_other_adj[part == "Sleep period"], aes(x = cancer_other_before_acc, y = Mean)) +
   # geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8") +
   # geom_hline(aes(yintercept = yintercept_others), linewidth = 0.5, linetype= "longdash", colour = "#a8a8a8") +
   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#CBD5D0", alpha = 0.1) +
   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#FAF7F4", alpha = 0.2) +
   geom_pointrange(aes(ymin = CI_low,
                       ymax = CI_high, colour = cancer_other_before_acc), size = 0.25, linewidth = 0.5) +
   geom_text(aes(y = 600, label = TeX(est_sig, output = "character")), parse = TRUE,
             hjust = 1, nudge_x = 0, 
             family = "Arial Narrow", size = 3,
             show.legend = FALSE) +
   geom_text(aes(y = 500, label = cancer_other_before_acc),
             hjust = 0, nudge_x = 0, 
             family = "Arial Narrow", size = 3,
             show.legend = FALSE) +
   geom_segment(aes(x = 0, yend = 500), col = "black", linewidth = 0.5) +
   geom_segment(aes(x = 0, yend = 600), col = "black", linewidth = 0.5) +   
   scale_y_continuous(limits = c(500, 600),
                      breaks = c(500, 550, 600),
                      name = "Sleep period") +
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
     axis.title.x        = element_text(size = 10, face = "bold", hjust = .5),
     axis.text.x         = element_text(size = 9),
     axis.text.y         = element_blank(),
     strip.text          = element_text(size = 9, hjust = .5, face = "bold"),
     legend.text         = element_text(size = 10, face = "bold", hjust = .5),
     legend.position     = "none",
     plot.margin         = unit(c(0.5,0,1,0), "lines")
   )
)

(plot_comp_cancer_mvpa <- 
    ggplot(comp_cancer_other_adj[part == "Moderate-to-vigorous physical activity"], aes(x = cancer_other_before_acc, y = Mean)) +
    # geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8") +
    # geom_hline(aes(yintercept = yintercept_others), linewidth = 0.5, linetype= "longdash", colour = "#a8a8a8") +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#CBD5D0", alpha = 0.1) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#FAF7F4", alpha = 0.2) +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_other_before_acc), size = 0.25, linewidth = 0.5) +
    geom_text(aes(y = 35, label = TeX(est_sig, output = "character")), parse = TRUE,
              hjust = 1, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_text(aes(y = 15, label = cancer_other_before_acc),
              hjust = 0, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_segment(aes(x = 0, yend = 15), col = "black", linewidth = 0.5) +
    geom_segment(aes(x = 0, yend = 35), col = "black", linewidth = 0.5) +
    scale_y_continuous(limits = c(15, 35),
                       breaks = c(15, 25, 35),
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
      axis.title.x        = element_text(size = 10, face = "bold", hjust = .5),
      axis.text.x         = element_text(size = 9),
      axis.text.y         = element_blank(),
      strip.text          = element_text(size = 9, hjust = .5, face = "bold"),
      legend.text         = element_text(size = 10, face = "bold", hjust = .5),
      legend.position     = "none",
      plot.margin         = unit(c(0.5,0,1,0), "lines")
    )
)

(plot_comp_cancer_lpa <- 
    ggplot(comp_cancer_other_adj[part == "Light physical activity"], aes(x = cancer_other_before_acc, y = Mean)) +
    # geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8") +
    # geom_hline(aes(yintercept = yintercept_others), linewidth = 0.5, linetype= "longdash", colour = "#a8a8a8") +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#CBD5D0", alpha = 0.1) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#FAF7F4", alpha = 0.2) +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_other_before_acc), size = 0.25, linewidth = 0.5) +
    geom_text(aes(y = 340, label = TeX(est_sig, output = "character")), parse = TRUE,
              hjust = 1, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_text(aes(y = 280, label = cancer_other_before_acc),
              hjust = 0, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_segment(aes(x = 0, yend = 280), col = "black", linewidth = 0.5) +
    geom_segment(aes(x = 0, yend = 340), col = "black", linewidth = 0.5) +
    scale_y_continuous(limits = c(280, 340),
                       breaks = c(280, 310, 340),
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
      axis.title.x        = element_text(size = 10, face = "bold", hjust = .5),
      axis.text.x         = element_text(size = 9),
      axis.text.y         = element_blank(),
      strip.text          = element_text(size = 9, hjust = .5, face = "bold"),
      legend.text         = element_text(size = 10, face = "bold", hjust = .5),
      legend.position     = "none",
      plot.margin         = unit(c(0.5,0,1,0), "lines")
    )
)

(plot_comp_cancer_sb <- 
    ggplot(comp_cancer_other_adj[part == "Sedentary behaviour"], aes(x = cancer_other_before_acc, y = Mean)) +
    # geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8") +
    # geom_hline(aes(yintercept = yintercept_others), linewidth = 0.5, linetype= "longdash", colour = "#a8a8a8") +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#CBD5D0", alpha = 0.1) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#FAF7F4", alpha = 0.2) +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_other_before_acc), size = 0.25, linewidth = 0.5) +
    geom_text(aes(y = 600, label = TeX(est_sig, output = "character")), parse = TRUE,
              hjust = 1, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_text(aes(y = 500, label = cancer_other_before_acc),
              hjust = 0, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_segment(aes(x = 0, yend = 600), col = "black", linewidth = 0.5) +
    geom_segment(aes(x = 0, yend = 500), col = "black", linewidth = 0.5) +
    scale_y_continuous(limits = c(500, 600),
                       breaks = c(500, 550, 600),
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
      axis.title.x        = element_text(size = 10, face = "bold", hjust = .5),
      axis.text.x         = element_text(size = 9),
      axis.text.y         = element_blank(),
      strip.text          = element_text(size = 9, hjust = .5, face = "bold"),
      legend.text         = element_text(size = 10, face = "bold", hjust = .5),
      legend.position     = "none",
      plot.margin         = unit(c(0.5,0,1,0), "lines")
    )
)


grDevices::cairo_pdf(
  file = paste0(outputdir, "cancer_other_est", ".pdf"),
  width = 8,
  height = 5,
)

ggarrange(
  plot_comp_cancer_mvpa,
  plot_comp_cancer_lpa,
  plot_comp_cancer_sb,
  plot_comp_cancer_sleep,
  nrow = 2, ncol = 2
)
dev.off()

grDevices::png(
  file = paste0(outputdir, "cancer_other_est", ".png"),
  width = 6000,
  height = 2000,
  res = 900
)

ggarrange(
  plot_comp_cancer_mvpa,
  plot_comp_cancer_lpa,
  plot_comp_cancer_sb,
  plot_comp_cancer_sleep,
  nrow = 2, ncol = 2
)
dev.off()

# GAM MODEL - CANCER BY TIME -------
# fit_gam_cancer_time_adj <- brmcoda(clr_cancer_acc,
#                                    mvbind(ilr1, ilr2, ilr3) ~ 
#                                      s(age_diff_cancer_acc, by = cancer_before_acc) +
#                                      s(age_at_acc, by = cancer_before_acc) + 
#                                      other_conds_at_acc +
#                                      sex + white + working + edu + never_smoked + current_drinker + s(deprivation),
#                                    warmup = 500, chains = 4, cores = 4, backend = "cmdstanr"
# )
# saveRDS(fit_gam_cancer_time_adj, paste0(outputdir, "fit_gam_cancer_time_adj", ".RDS"))
