source("ukb-cancer-24h-utils.R")
source(paste0(redir, "ukb_utils.R"))
# source("ukb-cancer-24h-data.R")

# CANCER TYPES VS HEALTHY
# main model --------
fit_cancer_type_adj <- brmcoda(clr_cancer_acc,
                               mvbind(ilr1, ilr2, ilr3) ~ cancer_before_acc_type +
                                 s(age_diff_cancer_acc) +
                                 s(age_at_acc) + sex + white + working + edu + never_smoked + current_drinker + s(deprivation),
                               # save_pars = save_pars(all = TRUE),
                               warmup = 500, chains = 4, cores = 4, backend = "cmdstanr"
)
saveRDS(fit_cancer_type_adj, paste0(outputdir, "fit_cancer_type_adj", ".RDS"))

# predicted posteriors ------------
fit_cancer_type_adj <- readRDS(paste0(outputdir, "fit_cancer_type_adj", ".RDS"))

# reference grid
d_cancer_type_adj <- emmeans::ref_grid(fit_cancer_type_adj$model)@grid

# predict
pred_cancer_type_adj <- fitted(fit_cancer_type_adj, newdata = d_cancer_type_adj, scale = "response", summary = FALSE)

# summarise by cancer types
pred_cancer_type_adj <- apply(pred_cancer_type_adj, c(1), function(x)  cbind(d_cancer_type_adj, x))
pred_cancer_type_adj <- lapply(pred_cancer_type_adj, function(d) {
  d <- as.data.table(d)
  
  # estimated means by cancer types
  d[, sleep := mean(V1), by = cancer_before_acc_type]
  d[, mvpa  := mean(V2), by = cancer_before_acc_type]
  d[, lpa   := mean(V3), by = cancer_before_acc_type]
  d[, sb    := mean(V4), by = cancer_before_acc_type]
  
  # constrast cancer types vs healthy
  d[, sleep_vs_healthy := sleep - d$sleep[1]]
  d[, mvpa_vs_healthy := mvpa - d$mvpa[1]]
  d[, lpa_vs_healthy := lpa - d$lpa[1]]
  d[, sb_vs_healthy := sb - d$sb[1]]
  
  d <- d[, .(cancer_before_acc_type, 
             sleep, mvpa, lpa, sb,
             sleep_vs_healthy, mvpa_vs_healthy, lpa_vs_healthy, sb_vs_healthy
  )]
  d <- unique(d)
  d
})

# assemble back to summarise posteriors
pred_cancer_type_adj <- as.data.frame(abind(pred_cancer_type_adj, along = 1))
pred_cancer_type_adj <- split(pred_cancer_type_adj, pred_cancer_type_adj$cancer_before_acc_type)

## estimated means  ----------------------
pred_comp_cancer_type_adj <- lapply(pred_cancer_type_adj, function(l) {
  l <- as.data.frame(l[, c("sleep", "mvpa", "lpa", "sb")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, part = names(l))
  l <- rbindlist(l)
  l
})
pred_comp_cancer_type_adj <- Map(cbind, pred_comp_cancer_type_adj, cancer_before_acc_type = names(pred_comp_cancer_type_adj))
pred_comp_cancer_type_adj <- rbindlist(pred_comp_cancer_type_adj)

## contrasts --------------------
diff_comp_cancer_type_adj <- lapply(pred_cancer_type_adj, function(l) {
  l <- as.data.frame(l[, c("sleep_vs_healthy", "mvpa_vs_healthy", "lpa_vs_healthy", "sb_vs_healthy")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, constrast = names(l))
  l <- rbindlist(l)
  l
})
diff_comp_cancer_type_adj <- Map(cbind, diff_comp_cancer_type_adj, diff_from_healthy = names(diff_comp_cancer_type_adj))
diff_comp_cancer_type_adj <- rbindlist(diff_comp_cancer_type_adj)

setnames(diff_comp_cancer_type_adj, "Mean", "Mean_diff")
setnames(diff_comp_cancer_type_adj, "CI_low", "CI_low_diff")
setnames(diff_comp_cancer_type_adj, "CI_high", "CI_high_diff")

# all results  ------------------------
comp_cancer_type_adj <- cbind(
  pred_comp_cancer_type_adj[, .(Mean, CI_low, CI_high, part, cancer_before_acc_type)],
  diff_comp_cancer_type_adj[, .(Mean_diff, CI_low_diff, CI_high_diff)]
)

# add sig indicators
comp_cancer_type_adj[, nonsig := between(0, comp_cancer_type_adj$CI_low_diff, comp_cancer_type_adj$CI_high_diff)]
comp_cancer_type_adj[, Sig := ifelse(nonsig == FALSE, paste(intToUtf8(0x2217)), " ")]

# sort by MVPA
# healthy as bottom in plot
# comp_cancer_type_adj[, cancer_before_acc_type := factor(cancer_before_acc_type, ordered = TRUE,
#                                                             levels = c(
#                                                               "Healthy",
#                                                               "Other Skin",
#                                                               "Prostate",
#                                                               "Melanoma",
#                                                               "Endocrine Gland",
#                                                               "Breast",
#                                                               "Genitourinary",
#                                                               "Colorectal",
#                                                               "Other Cancer",
#                                                               "Gynaecological",
#                                                               "Head & Neck",
#                                                               "Blood",
#                                                               "Gastrointestinal Tract",
#                                                               "Lung",
#                                                               "Multiple Primary"))]
# healthy as top in plot
comp_cancer_type_adj[, cancer_before_acc_type := factor(cancer_before_acc_type, ordered = TRUE,
                                                        levels = c(
                                                          "Multiple Primary",
                                                          "Lung",
                                                          "Gastrointestinal Tract",
                                                          "Blood",
                                                          "Head & Neck",
                                                          "Colorectal",
                                                          "Gynaecological",
                                                          "Other Cancer",
                                                          "Genitourinary",
                                                          "Endocrine Gland",
                                                          "Breast",
                                                          "Melanoma",
                                                          "Prostate",
                                                          "Other Skin",
                                                          "Healthy"
                                                        ))]

comp_cancer_type_adj[, yintercept := NA]
comp_cancer_type_adj[, yintercept := ifelse(part == "sleep", comp_cancer_type_adj[cancer_before_acc_type == "Healthy" & part == "sleep"]$Mean, yintercept)]
comp_cancer_type_adj[, yintercept := ifelse(part == "mvpa", comp_cancer_type_adj[cancer_before_acc_type == "Healthy" & part == "mvpa"]$Mean, yintercept)]
comp_cancer_type_adj[, yintercept := ifelse(part == "lpa", comp_cancer_type_adj[cancer_before_acc_type == "Healthy" & part == "lpa"]$Mean, yintercept)]
comp_cancer_type_adj[, yintercept := ifelse(part == "sb", comp_cancer_type_adj[cancer_before_acc_type == "Healthy" & part == "sb"]$Mean, yintercept)]

comp_cancer_type_adj[, part := ifelse(part == "sleep", "Sleep", part)]
comp_cancer_type_adj[, part := ifelse(part == "mvpa", "Moderate-to-vigorous physical activity", part)]
comp_cancer_type_adj[, part := ifelse(part == "lpa", "Light physical activity", part)]
comp_cancer_type_adj[, part := ifelse(part == "sb", "Sedentary behaviour", part)]

comp_cancer_type_adj[, sig_position := min(CI_low), by = part]
comp_cancer_type_adj[, est_position := max(CI_high), by = part]

comp_cancer_type_adj[, estimates := paste0(round(Mean, 0), " [", round(CI_low, 0), ", ", round(CI_high, ), "]")]
comp_cancer_type_adj[, est_sig := paste0(estimates, " ", str_replace_na(Sig, " "))]

# plots -----------------------------
## facet all -----------------------
(plot_comp_cancer_type_adj <- 
   ggplot(comp_cancer_type_adj, aes(x = cancer_before_acc_type, y = Mean, group = part)) +
   geom_hline(aes(yintercept = yintercept), linewidth = 0.5, linetype = 2, colour = "#a8a8a8") +
   geom_pointrange(aes(ymin = CI_low,
                       ymax = CI_high, colour = cancer_before_acc_type)) +
   geom_text(aes(y = sig_position + 8, label = Sig, colour = cancer_before_acc_type), 
             size = 5.5, 
             # position = position_dodge2(width = 1),
             show.legend = FALSE) +
   facet_wrap(~part, scales = "free") +
   scale_colour_manual(values = pal_type) +
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
  file = paste0(outputdir, "cancer_type_adj", ".pdf"),
  width = 11,
  height = 10,
)
plot_comp_cancer_type_adj
dev.off()


## plot by behaviour -----------------------
(plot_comp_cancer_type_mvpa <- 
   ggplot(comp_cancer_type_adj[part == "Moderate-to-vigorous physical activity"], aes(x = cancer_before_acc_type, y = Mean)) +
   geom_hline(aes(yintercept = yintercept), linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8") +
   geom_pointrange(aes(ymin = CI_low,
                       ymax = CI_high, colour = cancer_before_acc_type), size = 0.25, linewidth = 0.5) +
   geom_text(aes(y = 40, label = est_sig),
             hjust = 1, nudge_x = 0, 
             family = "Arial Narrow", size = 3,
             show.legend = FALSE) +
   geom_text(aes(y = 0, label = cancer_before_acc_type),
             hjust = 0, nudge_x = 0, 
             family = "Arial Narrow", size = 3,
             show.legend = FALSE) +
   geom_segment(aes(x = 0, yend = 10), col = "black", linewidth = 0.5) +
   geom_segment(aes(x = 0, yend = 30), col = "black", linewidth = 0.5) +
   scale_y_continuous(limits = c(0, 40),
                      breaks = c(10, 20, 30),
                      name = "Moderate-to-vigorous physical activity") +
   scale_colour_manual(values = pal_type) +
   
   labs(x = "", y = "", colour = "") +
   coord_flip() +
   theme_ipsum() +
   theme(
     axis.ticks          = element_blank(),
     # panel.background    = element_rect(fill = "transparent", colour = "black", linewidth = 0.5),
     plot.background     = element_rect(fill = "transparent", colour = NA, linewidth = 0.5),
     panel.grid.major    = element_blank(),
     panel.grid.minor    = element_blank(),
     # axis.line.x         = element_line(linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8"),
     axis.title.x        = element_text(size = 9, face = "bold", hjust = .5),
     axis.text.x         = element_text(size = 8),
     axis.text.y         = element_blank(),
     strip.text          = element_text(size = 8, hjust = .5, face = "bold"),
     legend.text         = element_text(size = 9, face = "bold", hjust = .5),
     legend.position     = "none",
     plot.margin         = unit(c(0.5,0,0.5,0), "lines")
   )
)

(plot_comp_cancer_type_lpa <- 
    ggplot(comp_cancer_type_adj[part == "Light physical activity"], aes(x = cancer_before_acc_type, y = Mean)) +
    geom_hline(aes(yintercept = yintercept), linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8") +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_before_acc_type), size = 0.25, linewidth = 0.5) +
    geom_text(aes(y = 375, label = est_sig),
              hjust = 1, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_text(aes(y = 175, label = cancer_before_acc_type),
              hjust = 0, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_segment(aes(x = 0, yend = 225), col = "black", linewidth = 0.5) +
    geom_segment(aes(x = 0, yend = 325), col = "black", linewidth = 0.5) +    
    scale_y_continuous(limits = c(175, 375),
                       breaks = c(225, 275, 325),
                       name = "Light physical activity") +
    scale_colour_manual(values = pal_type) +
    
    labs(x = "", y = "", colour = "") +
    coord_flip() +
    theme_ipsum() +
    theme(
      axis.ticks          = element_blank(),
      # panel.background    = element_rect(fill = "transparent", colour = "black", linewidth = 0.5),
      plot.background     = element_rect(fill = "transparent", colour = NA, linewidth = 0.5),
      panel.grid.major    = element_blank(),
      panel.grid.minor    = element_blank(),
      # axis.line.x         = element_line(linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8"),
      axis.title.x        = element_text(size = 9, face = "bold", hjust = .5),
      axis.text.x         = element_text(size = 8),
      axis.text.y         = element_blank(),
      strip.text          = element_text(size = 8, hjust = .5, face = "bold"),
      legend.text         = element_text(size = 9, face = "bold", hjust = .5),
      legend.position     = "none",
      plot.margin         = unit(c(0.5,0,0,0), "lines")
    )
)

(plot_comp_cancer_type_sb <- 
    ggplot(comp_cancer_type_adj[part == "Sedentary behaviour"], aes(x = cancer_before_acc_type, y = Mean)) +
    geom_hline(aes(yintercept = yintercept), linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8") +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_before_acc_type), size = 0.25, linewidth = 0.5) +
    geom_text(aes(y = 700, label = est_sig),
              hjust = 1, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_text(aes(y = 500, label = cancer_before_acc_type),
              hjust = 0, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_segment(aes(x = 0, yend = 550), col = "black", linewidth = 0.5) +
    geom_segment(aes(x = 0, yend = 650), col = "black", linewidth = 0.5) +    
    scale_y_continuous(limits = c(500, 700),
                       breaks = c(550, 600, 650),
                       name = "Sedentary behaviour") +
    scale_colour_manual(values = pal_type) +
    
    labs(x = "", y = "", colour = "") +
    coord_flip() +
    theme_ipsum() +
    theme(
      axis.ticks          = element_blank(),
      # panel.background    = element_rect(fill = "transparent", colour = "black", linewidth = 0.5),
      plot.background     = element_rect(fill = "transparent", colour = NA, linewidth = 0.5),
      panel.grid.major    = element_blank(),
      panel.grid.minor    = element_blank(),
      # axis.line.x         = element_line(linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8"),
      axis.title.x        = element_text(size = 9, face = "bold", hjust = .5),
      axis.text.x         = element_text(size = 8),
      axis.text.y         = element_blank(),
      strip.text          = element_text(size = 8, hjust = .5, face = "bold"),
      legend.text         = element_text(size = 9, face = "bold", hjust = .5),
      legend.position     = "none",
      plot.margin         = unit(c(0.5,0,0,0), "lines")
    )
)


(plot_comp_cancer_type_sleep <- 
    ggplot(comp_cancer_type_adj[part == "Sleep"], aes(x = cancer_before_acc_type, y = Mean)) +
    geom_hline(aes(yintercept = yintercept), linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8") +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_before_acc_type), size = 0.25, linewidth = 0.5) +
    geom_text(aes(y = 600, label = est_sig),
              hjust = 1, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_text(aes(y = 500, label = cancer_before_acc_type),
              hjust = 0, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_segment(aes(x = 0, yend = 525), col = "black", linewidth = 0.5) +
    geom_segment(aes(x = 0, yend = 575), col = "black", linewidth = 0.5) +    
    scale_y_continuous(limits = c(500, 600),
                       breaks = c(525, 550, 575),
                       name = "Sleep period") +
    scale_colour_manual(values = pal_type) +
    
    labs(x = "", y = "", colour = "") +
    coord_flip() +
    theme_ipsum() +
    theme(
      axis.ticks          = element_blank(),
      # panel.background    = element_rect(fill = "transparent", colour = "black", linewidth = 0.5),
      plot.background     = element_rect(fill = "transparent", colour = NA, linewidth = 0.5),
      panel.grid.major    = element_blank(),
      panel.grid.minor    = element_blank(),
      # axis.line.x         = element_line(linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8"),
      axis.title.x        = element_text(size = 9, face = "bold", hjust = .5),
      axis.text.x         = element_text(size = 8),
      axis.text.y         = element_blank(),
      strip.text          = element_text(size = 8, hjust = .5, face = "bold"),
      legend.text         = element_text(size = 9, face = "bold", hjust = .5),
      legend.position     = "none",
      plot.margin         = unit(c(0.5,0,0,0), "lines")
    )
)

# save
# grDevices::cairo_pdf(
#   file = paste0(outputdir, "plot_comp_cancer_type_mvpa_est", ".pdf"),
#   width = 6,
#   height = 5,
# )
# plot_comp_cancer_type_mvpa
# dev.off()
# 
# grDevices::cairo_pdf(
#   file = paste0(outputdir, "plot_comp_cancer_type_lpa_est", ".pdf"),
#   width = 6,
#   height = 5,
# )
# plot_comp_cancer_type_lpa
# dev.off()
# 
# grDevices::cairo_pdf(
#   file = paste0(outputdir, "plot_comp_cancer_type_sb_est", ".pdf"),
#   width = 6,
#   height = 5,
# )
# plot_comp_cancer_type_sb
# dev.off()
# 
# grDevices::cairo_pdf(
#   file = paste0(outputdir, "plot_comp_cancer_type_sleep_est", ".pdf"),
#   width = 6,
#   height = 5,
# )
# plot_comp_cancer_type_sleep
# dev.off()

grDevices::cairo_pdf(
  file = paste0(outputdir, "cancer_type_est", ".pdf"),
  width = 6,
  height = 11,
)

ggarrange(
  plot_comp_cancer_type_mvpa,
  plot_comp_cancer_type_lpa,
  plot_comp_cancer_type_sb,
  plot_comp_cancer_type_sleep,
  nrow = 4
)
dev.off()

# CANCER TYPES VS HEALTHY VS OTHERS
# main model --------
fit_cancer_type_other_adj <- brmcoda(clr_cancer_acc,
                               mvbind(ilr1, ilr2, ilr3) ~ cancer_before_acc_type_other +
                                 # s(age_diff_cancer_acc) +
                                 s(age_at_acc) + sex + white + working + edu + never_smoked + current_drinker + s(deprivation),
                               # save_pars = save_pars(all = TRUE),
                               warmup = 500, chains = 4, cores = 4, backend = "cmdstanr"
)
saveRDS(fit_cancer_type_other_adj, paste0(outputdir, "fit_cancer_type_other_adj", ".RDS"))

# predicted posteriors ------------
fit_cancer_type_other_adj <- readRDS(paste0(outputdir, "fit_cancer_type_other_adj", ".RDS"))

# reference grid
d_cancer_type_other_adj <- emmeans::ref_grid(fit_cancer_type_other_adj$model)@grid

# predict
pred_cancer_type_other_adj <- fitted(fit_cancer_type_other_adj, newdata = d_cancer_type_other_adj, scale = "response", summary = FALSE)

# summarise by cancer types
pred_cancer_type_other_adj <- apply(pred_cancer_type_other_adj, c(1), function(x)  cbind(d_cancer_type_other_adj, x))
pred_cancer_type_other_adj <- lapply(pred_cancer_type_other_adj, function(d) {
  d <- as.data.table(d)
  
  # estimated means by cancer types
  d[, sleep := mean(V1), by = cancer_before_acc_type_other]
  d[, mvpa  := mean(V2), by = cancer_before_acc_type_other]
  d[, lpa   := mean(V3), by = cancer_before_acc_type_other]
  d[, sb    := mean(V4), by = cancer_before_acc_type_other]
  
  # constrast cancer types vs healthy
  d[, sleep_vs_healthy := sleep - d[cancer_before_acc_type_other == "Healthy"]$sleep[1]]
  d[, mvpa_vs_healthy := mvpa - d[cancer_before_acc_type_other == "Healthy"]$mvpa[1]]
  d[, lpa_vs_healthy := lpa - d[cancer_before_acc_type_other == "Healthy"]$lpa[1]]
  d[, sb_vs_healthy := sb - d[cancer_before_acc_type_other == "Healthy"]$sb[1]]
  
  # constrast cancer types vs healthy
  d[, sleep_vs_others := sleep - d[cancer_before_acc_type_other == "Others"]$sleep[1]]
  d[, mvpa_vs_others := mvpa - d[cancer_before_acc_type_other == "Others"]$mvpa[1]]
  d[, lpa_vs_others := lpa - d[cancer_before_acc_type_other == "Others"]$lpa[1]]
  d[, sb_vs_others := sb - d[cancer_before_acc_type_other == "Others"]$sb[1]]
  
  d <- d[, .(cancer_before_acc_type_other, 
             sleep, mvpa, lpa, sb,
             
             sleep_vs_healthy, mvpa_vs_healthy, lpa_vs_healthy, sb_vs_healthy,
             sleep_vs_others, mvpa_vs_others, lpa_vs_others, sb_vs_others
             
  )]
  d <- unique(d)
  d
})

# assemble back to summarise posteriors
pred_cancer_type_other_adj <- as.data.frame(abind(pred_cancer_type_other_adj, along = 1))
pred_cancer_type_other_adj <- split(pred_cancer_type_other_adj, pred_cancer_type_other_adj$cancer_before_acc_type_other)

## estimated means  ----------------------
pred_comp_cancer_type_other_adj <- lapply(pred_cancer_type_other_adj, function(l) {
  l <- as.data.frame(l[, c("sleep", "mvpa", "lpa", "sb")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, part = names(l))
  l <- rbindlist(l)
  l
})
pred_comp_cancer_type_other_adj <- Map(cbind, pred_comp_cancer_type_other_adj, cancer_before_acc_type_other = names(pred_comp_cancer_type_other_adj))
pred_comp_cancer_type_other_adj <- rbindlist(pred_comp_cancer_type_other_adj)

## contrasts --------------------
### vs healthy
diff1_comp_cancer_type_other_adj <- lapply(pred_cancer_type_other_adj, function(l) {
  l <- as.data.frame(l[, c("sleep_vs_healthy", "mvpa_vs_healthy", "lpa_vs_healthy", "sb_vs_healthy")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, constrast = names(l))
  l <- rbindlist(l)
  l
})
diff1_comp_cancer_type_other_adj <- Map(cbind, diff1_comp_cancer_type_other_adj, diff_from_healthy = names(diff1_comp_cancer_type_other_adj))
diff1_comp_cancer_type_other_adj <- rbindlist(diff1_comp_cancer_type_other_adj)

setnames(diff1_comp_cancer_type_other_adj, "Mean", "Mean_diff_ref_healthy")
setnames(diff1_comp_cancer_type_other_adj, "CI_low", "CI_low_diff_ref_healthy")
setnames(diff1_comp_cancer_type_other_adj, "CI_high", "CI_high_diff_ref_healthy")

### vs others
diff2_comp_cancer_type_other_adj <- lapply(pred_cancer_type_other_adj, function(l) {
  l <- as.data.frame(l[, c("sleep_vs_others", "mvpa_vs_others", "lpa_vs_others", "sb_vs_others")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, constrast = names(l))
  l <- rbindlist(l)
  l
})
diff2_comp_cancer_type_other_adj <- Map(cbind, diff2_comp_cancer_type_other_adj, diff_from_healthy = names(diff2_comp_cancer_type_other_adj))
diff2_comp_cancer_type_other_adj <- rbindlist(diff2_comp_cancer_type_other_adj)

setnames(diff2_comp_cancer_type_other_adj, "Mean", "Mean_diff_ref_others")
setnames(diff2_comp_cancer_type_other_adj, "CI_low", "CI_low_diff_ref_others")
setnames(diff2_comp_cancer_type_other_adj, "CI_high", "CI_high_diff_ref_others")

# all results  ------------------------
comp_cancer_type_other_adj <- cbind(
  pred_comp_cancer_type_other_adj[, .(Mean, CI_low, CI_high, part, cancer_before_acc_type_other)],
  diff1_comp_cancer_type_other_adj[, .(Mean_diff_ref_healthy, CI_low_diff_ref_healthy, CI_high_diff_ref_healthy)],
  diff2_comp_cancer_type_other_adj[, .(Mean_diff_ref_others, CI_low_diff_ref_others, CI_high_diff_ref_others)]
)

# add sig indicators
comp_cancer_type_other_adj[, nonsig_healthy := between(0, comp_cancer_type_other_adj$CI_low_diff_ref_healthy, comp_cancer_type_other_adj$CI_high_diff_ref_healthy)]
comp_cancer_type_other_adj[, nonsig_others := between(0, comp_cancer_type_other_adj$CI_low_diff_ref_others, comp_cancer_type_other_adj$CI_high_diff_ref_others)]

comp_cancer_type_other_adj[, sig_ref_healthy := ifelse(nonsig_healthy == FALSE, "$^a$", "$\\phantom{a}$")]
comp_cancer_type_other_adj[, sig_ref_others := ifelse(nonsig_others == FALSE, "$^b$", "$\\phantom{b}$")]

comp_cancer_type_other_adj[, yintercept_healthy := NA]
comp_cancer_type_other_adj[, yintercept_healthy := ifelse(part == "sleep", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Healthy" & part == "sleep"]$Mean, yintercept_healthy)]
comp_cancer_type_other_adj[, yintercept_healthy := ifelse(part == "mvpa", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Healthy" & part == "mvpa"]$Mean, yintercept_healthy)]
comp_cancer_type_other_adj[, yintercept_healthy := ifelse(part == "lpa", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Healthy" & part == "lpa"]$Mean, yintercept_healthy)]
comp_cancer_type_other_adj[, yintercept_healthy := ifelse(part == "sb", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Healthy" & part == "sb"]$Mean, yintercept_healthy)]

comp_cancer_type_other_adj[, ci_low_healthy := NA]
comp_cancer_type_other_adj[, ci_low_healthy := ifelse(part == "sleep", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Healthy" & part == "sleep"]$CI_low, ci_low_healthy)]
comp_cancer_type_other_adj[, ci_low_healthy := ifelse(part == "mvpa", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Healthy" & part == "mvpa"]$CI_low, ci_low_healthy)]
comp_cancer_type_other_adj[, ci_low_healthy := ifelse(part == "lpa", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Healthy" & part == "lpa"]$CI_low, ci_low_healthy)]
comp_cancer_type_other_adj[, ci_low_healthy := ifelse(part == "sb", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Healthy" & part == "sb"]$CI_low, ci_low_healthy)]

comp_cancer_type_other_adj[, ci_high_healthy := NA]
comp_cancer_type_other_adj[, ci_high_healthy := ifelse(part == "sleep", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Healthy" & part == "sleep"]$CI_high, ci_high_healthy)]
comp_cancer_type_other_adj[, ci_high_healthy := ifelse(part == "mvpa", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Healthy" & part == "mvpa"]$CI_high, ci_high_healthy)]
comp_cancer_type_other_adj[, ci_high_healthy := ifelse(part == "lpa", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Healthy" & part == "lpa"]$CI_high, ci_high_healthy)]
comp_cancer_type_other_adj[, ci_high_healthy := ifelse(part == "sb", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Healthy" & part == "sb"]$CI_high, ci_high_healthy)]

comp_cancer_type_other_adj[, yintercept_others := NA]
comp_cancer_type_other_adj[, yintercept_others := ifelse(part == "sleep", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Others" & part == "sleep"]$Mean, yintercept_others)]
comp_cancer_type_other_adj[, yintercept_others := ifelse(part == "mvpa", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Others" & part == "mvpa"]$Mean, yintercept_others)]
comp_cancer_type_other_adj[, yintercept_others := ifelse(part == "lpa", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Others" & part == "lpa"]$Mean, yintercept_others)]
comp_cancer_type_other_adj[, yintercept_others := ifelse(part == "sb", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Others" & part == "sb"]$Mean, yintercept_others)]

comp_cancer_type_other_adj[, ci_low_others := NA]
comp_cancer_type_other_adj[, ci_low_others := ifelse(part == "sleep", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Others" & part == "sleep"]$CI_low, ci_low_others)]
comp_cancer_type_other_adj[, ci_low_others := ifelse(part == "mvpa", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Others" & part == "mvpa"]$CI_low, ci_low_others)]
comp_cancer_type_other_adj[, ci_low_others := ifelse(part == "lpa", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Others" & part == "lpa"]$CI_low, ci_low_others)]
comp_cancer_type_other_adj[, ci_low_others := ifelse(part == "sb", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Others" & part == "sb"]$CI_low, ci_low_others)]

comp_cancer_type_other_adj[, ci_high_others := NA]
comp_cancer_type_other_adj[, ci_high_others := ifelse(part == "sleep", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Others" & part == "sleep"]$CI_high, ci_high_others)]
comp_cancer_type_other_adj[, ci_high_others := ifelse(part == "mvpa", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Others" & part == "mvpa"]$CI_high, ci_high_others)]
comp_cancer_type_other_adj[, ci_high_others := ifelse(part == "lpa", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Others" & part == "lpa"]$CI_high, ci_high_others)]
comp_cancer_type_other_adj[, ci_high_others := ifelse(part == "sb", comp_cancer_type_other_adj[cancer_before_acc_type_other == "Others" & part == "sb"]$CI_high, ci_high_others)]

## sort by MVPA
## healthy as bottom in plot
comp_cancer_type_other_adj[, cancer_before_acc_type_other := factor(cancer_before_acc_type_other, ordered = TRUE,
                                                            levels = c(
                                                              "Healthy",
                                                              "Others",
                                                              "Other Skin",
                                                              "Prostate",
                                                              "Melanoma",
                                                              "Endocrine Gland",
                                                              "Breast",
                                                              "Other Cancer",
                                                              "Colorectal",
                                                              "Gynaecological",
                                                              "Genitourinary",
                                                              "Head & Neck",
                                                              "Blood",
                                                              "Gastrointestinal Tract",
                                                              "Lung",
                                                              "Multiple Primary"))]
# # healthy as top in plot
# comp_cancer_type_other_adj[, cancer_before_acc_type_other := factor(cancer_before_acc_type_other, ordered = TRUE,
#                                                         levels = c(
#                                                           "Multiple Primary",
#                                                           "Lung",
#                                                           "Gastrointestinal Tract",
#                                                           "Blood",
#                                                           "Head & Neck",
#                                                           "Colorectal",
#                                                           "Gynaecological",
#                                                           "Other Cancer",
#                                                           "Genitourinary",
#                                                           "Endocrine Gland",
#                                                           "Breast",
#                                                           "Melanoma",
#                                                           "Prostate",
#                                                           "Other Skin",
#                                                           "Healthy"
#                                                         ))]

comp_cancer_type_other_adj[, part := ifelse(part == "sleep", "Sleep period", part)]
comp_cancer_type_other_adj[, part := ifelse(part == "mvpa", "Moderate-to-vigorous physical activity", part)]
comp_cancer_type_other_adj[, part := ifelse(part == "lpa", "Light physical activity", part)]
comp_cancer_type_other_adj[, part := ifelse(part == "sb", "Sedentary behaviour", part)]

comp_cancer_type_other_adj[, sig_position := min(CI_low), by = part]
comp_cancer_type_other_adj[, est_position := max(CI_high), by = part]

comp_cancer_type_other_adj[, estimates := paste0(round(Mean, 0), "[", round(CI_low, 0), ", ", round(CI_high, 0), "]")]
comp_cancer_type_other_adj[, est_sig := paste0(estimates, " ", sig_ref_healthy, sig_ref_others)]
# comp_cancer_type_other_adj[, est_sig := paste0(estimates, " ", str_replace_na(sig_ref_healthy, " "), str_replace_na(sig_ref_others, " "))]

## plot -----------------------
(plot_comp_cancer_type_other_sleep <- 
   ggplot(comp_cancer_type_other_adj[part == "Sleep period"], aes(x = cancer_before_acc_type_other, y = Mean)) +
   # geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8") +
   # geom_hline(aes(yintercept = yintercept_others), linewidth = 0.5, linetype= "longdash", colour = "#a8a8a8") +
   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#CBD5D0", alpha = 0.1) +
   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#FAF7F4", alpha = 0.2) +
   geom_pointrange(aes(ymin = CI_low,
                       ymax = CI_high, colour = cancer_before_acc_type_other), size = 0.25, linewidth = 0.5) +
   geom_text(aes(y = 600, label = TeX(est_sig, output = "character")), parse = TRUE,
             hjust = 1, nudge_x = 0, 
             family = "Arial Narrow", size = 3,
             show.legend = FALSE) +
   geom_text(aes(y = 500, label = cancer_before_acc_type_other),
             hjust = 0, nudge_x = 0, 
             family = "Arial Narrow", size = 3,
             show.legend = FALSE) +
   geom_segment(aes(x = 0, yend = 500), col = "black", linewidth = 0.5) +
   geom_segment(aes(x = 0, yend = 600), col = "black", linewidth = 0.5) +   
   scale_y_continuous(limits = c(500, 600),
                      breaks = c(500, 550, 600),
                      name = "Sleep period") +
   scale_colour_manual(values = pal_type) +
   labs(x = "", y = "", colour = "") +
   coord_flip() +
   theme_ipsum() +
   theme(
     axis.ticks          = element_blank(),
     # panel.background    = element_rect(fill = "transparent", colour = "black", linewidth = 0.5),
     plot.background     = element_rect(fill = "transparent", colour = NA, linewidth = 0.5),
     panel.grid.major    = element_blank(),
     panel.grid.minor    = element_blank(),
     # axis.line.x         = element_line(linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8"),
     axis.title.x        = element_text(size = 9, face = "bold", hjust = .5),
     axis.text.x         = element_text(size = 8),
     axis.text.y         = element_blank(),
     strip.text          = element_text(size = 8, hjust = .5, face = "bold"),
     legend.text         = element_text(size = 9, face = "bold", hjust = .5),
     legend.position     = "none",
     plot.margin         = unit(c(0.5,0,0,0), "lines")
   )
)
(plot_comp_cancer_type_other_mvpa <- 
   ggplot(comp_cancer_type_other_adj[part == "Moderate-to-vigorous physical activity"], aes(x = cancer_before_acc_type_other, y = Mean)) +
    # geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8") +
    # geom_hline(aes(yintercept = yintercept_others), linewidth = 0.5, linetype= "longdash", colour = "#a8a8a8") +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#CBD5D0", alpha = 0.1) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#FAF7F4", alpha = 0.2) +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_before_acc_type_other), size = 0.25, linewidth = 0.5) +
    geom_text(aes(y = 35, label = TeX(est_sig, output = "character")), parse = TRUE,
              hjust = 1, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_text(aes(y = 5, label = cancer_before_acc_type_other),
              hjust = 0, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_segment(aes(x = 0, yend = 5), col = "black", linewidth = 0.5) +
    geom_segment(aes(x = 0, yend = 35), col = "black", linewidth = 0.5) +
    scale_y_continuous(limits = c(5, 35),
                       breaks = c(5, 20, 35),
                       name = "Moderate-to-vigorous physical activity") +
    scale_colour_manual(values = pal_type) +
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
      plot.margin         = unit(c(0.5,0,0,0), "lines")
    )
)

(plot_comp_cancer_type_other_lpa <- 
    ggplot(comp_cancer_type_other_adj[part == "Light physical activity"], aes(x = cancer_before_acc_type_other, y = Mean)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#CBD5D0", alpha = 0.1) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#FAF7F4", alpha = 0.2) +    
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_before_acc_type_other), size = 0.25, linewidth = 0.5) +
    geom_text(aes(y = 350, label = TeX(est_sig, output = "character")), parse = TRUE,
              hjust = 1, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_text(aes(y = 200, label = cancer_before_acc_type_other),
              hjust = 0, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_segment(aes(x = 0, yend = 200), col = "black", linewidth = 0.5) +
    geom_segment(aes(x = 0, yend = 350), col = "black", linewidth = 0.5) +
    scale_y_continuous(limits = c(200, 350),
                       breaks = c(200, 275, 350),
                       name = "Light physical activity") +
    scale_colour_manual(values = pal_type) +
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
      plot.margin         = unit(c(0.5,0,0,0), "lines")
    )
)

(plot_comp_cancer_type_other_sb <- 
    ggplot(comp_cancer_type_other_adj[part == "Sedentary behaviour"], aes(x = cancer_before_acc_type_other, y = Mean)) +
    # geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8") +
    # geom_hline(aes(yintercept = yintercept_others), linewidth = 0.5, linetype= "longdash", colour = "#a8a8a8") +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#CBD5D0", alpha = 0.1) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#FAF7F4", alpha = 0.2) + 
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_before_acc_type_other), size = 0.25, linewidth = 0.5) +
    geom_text(aes(y = 650, label = TeX(est_sig, output = "character")), parse = TRUE,
              hjust = 1, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_text(aes(y = 500, label = cancer_before_acc_type_other),
              hjust = 0, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_segment(aes(x = 0, yend = 650), col = "black", linewidth = 0.5) +
    geom_segment(aes(x = 0, yend = 500), col = "black", linewidth = 0.5) +
    scale_y_continuous(limits = c(500, 650),
                       breaks = c(500, 575, 650),
                       name = "Sedentary behaviour") +
    scale_colour_manual(values = pal_type) +
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
      plot.margin         = unit(c(0.5,0,0,0), "lines")
    )
)

# save
# grDevices::cairo_pdf(
#   file = paste0(outputdir, "plot_comp_cancer_type_other_mvpa_est", ".pdf"),
#   width = 6,
#   height = 5,
# )
# plot_comp_cancer_type_other_mvpa
# dev.off()
# 
# grDevices::cairo_pdf(
#   file = paste0(outputdir, "plot_comp_cancer_type_other_lpa_est", ".pdf"),
#   width = 6,
#   height = 5,
# )
# plot_comp_cancer_type_other_lpa
# dev.off()
# 
# grDevices::cairo_pdf(
#   file = paste0(outputdir, "plot_comp_cancer_type_other_sb_est", ".pdf"),
#   width = 6,
#   height = 5,
# )
# plot_comp_cancer_type_other_sb
# dev.off()
# 
# grDevices::cairo_pdf(
#   file = paste0(outputdir, "plot_comp_cancer_type_other_sleep_est", ".pdf"),
#   width = 6,
#   height = 5,
# )
# plot_comp_cancer_type_other_sleep
# dev.off()

grDevices::cairo_pdf(
  file = paste0(outputdir, "cancer_type_other_est", ".pdf"),
  width = 8,
  height = 14,
)

ggarrange(
  plot_comp_cancer_type_other_mvpa,
  plot_comp_cancer_type_other_lpa,
  plot_comp_cancer_type_other_sb,
  plot_comp_cancer_type_other_sleep,
  nrow = 4
)
dev.off()

grDevices::png(
  file = paste0(outputdir, "cancer_type_other_est", ".png"),
  width = 6500,
  height = 12000,
  res = 900
)

ggarrange(
  plot_comp_cancer_type_other_mvpa,
  plot_comp_cancer_type_other_lpa,
  plot_comp_cancer_type_other_sb,
  plot_comp_cancer_type_other_sleep,
  nrow = 4
)
dev.off()
