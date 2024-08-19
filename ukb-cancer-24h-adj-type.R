source("ukb-cancer-24h-utils.R")
source(paste0(redir, "ukb_utils.R"))
# source("ukb-cancer-24h-data.R")

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
  d[, sleep_constrast := sleep - d$sleep[1]]
  d[, mvpa_constrast := mvpa - d$mvpa[1]]
  d[, lpa_constrast := lpa - d$lpa[1]]
  d[, sb_constrast := sb - d$sb[1]]
  
  d <- d[, .(cancer_before_acc_type, 
             sleep, mvpa, lpa, sb,
             sleep_constrast, mvpa_constrast, lpa_constrast, sb_constrast
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
  l <- as.data.frame(l[, c("sleep_constrast", "mvpa_constrast", "lpa_constrast", "sb_constrast")])
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
comp_cancer_type_adj[, est_sig := paste0(estimates, " ", str_replace_na(Sig, ""))]

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
                        ymax = CI_high, colour = cancer_before_acc_type), size = 0.5, linewidth = 0.75) +
    # geom_text(aes(y = 35 - 1.5, label = sig_ref_healthy, colour = cancer_before_acc_type), 
    #           size = 6, nudge_x = 0, 
    #           show.legend = FALSE) +
    # geom_text(aes(y = 35 - 1.25, label = sig_ref_cancer, colour = cancer_before_acc_type), 
    #           size = 4, nudge_x = 0,
    #           show.legend = FALSE) +
    geom_text(aes(y = 40, label = est_sig),
              vjust = "outward", hjust = 1, nudge_x = 0, 
              family = "Arial Narrow", size = 4,
              show.legend = FALSE) +
    geom_text(aes(y = 0, label = cancer_before_acc_type),
              vjust = "outward", hjust = 0, nudge_x = 0, 
              family = "Arial Narrow", size = 4,
              show.legend = FALSE) +
    # facet_wrap(~part, scales = "free", nrow = 4) +
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
      axis.line.x         = element_line(linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8"),
      axis.title.x        = element_text(size = 13, face = "bold", hjust = .5),
      axis.text.x         = element_text(size = 12),
      axis.text.y         = element_blank(),
      strip.text          = element_text(size = 12, hjust = .5, face = "bold"),
      legend.text         = element_text(size = 13, face = "bold", hjust = .5),
      legend.position     = "none",
      plot.margin         = unit(c(0.5,0,0.5,0), "lines")
    )
)

(plot_comp_cancer_type_lpa <- 
    ggplot(comp_cancer_type_adj[part == "Light physical activity"], aes(x = cancer_before_acc_type, y = Mean)) +
    geom_hline(aes(yintercept = yintercept), linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8") +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_before_acc_type), size = 0.5, linewidth = 0.75) +
    # geom_text(aes(y = 350 - 8, label = sig_ref_healthy, colour = cancer_before_acc_type), 
    #           size = 6, nudge_x = 0, 
    #           show.legend = FALSE) +
    # geom_text(aes(y = 350 - 7, label = sig_ref_cancer, colour = cancer_before_acc_type), 
    #           size = 4, nudge_x = 0,
    #           show.legend = FALSE) +
    geom_text(aes(y = 375, label = est_sig),
              vjust = "outward", hjust = 1, nudge_x = 0, 
              family = "Arial Narrow", size = 4,
              show.legend = FALSE) +
    geom_text(aes(y = 175, label = cancer_before_acc_type),
              vjust = "outward", hjust = 0, nudge_x = 0, 
              family = "Arial Narrow", size = 4,
              show.legend = FALSE) +
    # facet_wrap(~part, scales = "free", nrow = 4) +
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
      axis.line.x         = element_line(linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8"),
      axis.title.x        = element_text(size = 13, face = "bold", hjust = .5),
      axis.text.x         = element_text(size = 12),
      axis.text.y         = element_blank(),
      strip.text          = element_text(size = 12, hjust = .5, face = "bold"),
      legend.text         = element_text(size = 13, face = "bold", hjust = .5),
      legend.position     = "none",
      plot.margin         = unit(c(2,0,0,0), "lines")
    )
)

(plot_comp_cancer_type_sb <- 
    ggplot(comp_cancer_type_adj[part == "Sedentary behaviour"], aes(x = cancer_before_acc_type, y = Mean)) +
    geom_hline(aes(yintercept = yintercept), linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8") +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_before_acc_type), size = 0.5, linewidth = 0.75) +
    # geom_text(aes(y = 600 - 8, label = sig_ref_healthy, colour = cancer_before_acc_type), 
    #           size = 6, nudge_x = 0, 
    #           show.legend = FALSE) +
    # geom_text(aes(y = 600 - 7, label = sig_ref_cancer, colour = cancer_before_acc_type), 
    #           size = 4, nudge_x = 0,
    #           show.legend = FALSE) +
    geom_text(aes(y = 700, label = est_sig),
              vjust = "outward", hjust = 1, nudge_x = 0, 
              family = "Arial Narrow", size = 4,
              show.legend = FALSE) +
    geom_text(aes(y = 500, label = cancer_before_acc_type),
              vjust = "outward", hjust = 0, nudge_x = 0, 
              family = "Arial Narrow", size = 4,
              show.legend = FALSE) +
    # facet_wrap(~part, scales = "free", nrow = 4) +
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
      axis.line.x         = element_line(linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8"),
      axis.title.x        = element_text(size = 13, face = "bold", hjust = .5),
      axis.text.x         = element_text(size = 12),
      axis.text.y         = element_blank(),
      strip.text          = element_text(size = 12, hjust = .5, face = "bold"),
      legend.text         = element_text(size = 13, face = "bold", hjust = .5),
      legend.position     = "none",
      plot.margin         = unit(c(2,0,0,0), "lines")
    )
)


(plot_comp_cancer_type_sleep <- 
    ggplot(comp_cancer_type_adj[part == "Sleep"], aes(x = cancer_before_acc_type, y = Mean)) +
    geom_hline(aes(yintercept = yintercept), linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8") +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_before_acc_type), size = 0.5, linewidth = 0.75) +
    geom_text(aes(y = 600, label = est_sig),
              vjust = "outward", hjust = 1, nudge_x = 0, 
              family = "Arial Narrow", size = 4,
              show.legend = FALSE) +
    geom_text(aes(y = 500, label = cancer_before_acc_type),
              vjust = "outward", hjust = 0, nudge_x = 0, 
              family = "Arial Narrow", size = 4,
              show.legend = FALSE) +
    # facet_wrap(~part, scales = "free", nrow = 4) +
    scale_y_continuous(limits = c(500, 600),
                       breaks = c(525, 550, 575),
                       name = "Sleep") +
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
      axis.line.x         = element_line(linewidth = 0.5, linetype= "dashed", colour = "#a8a8a8"),
      axis.title.x        = element_text(size = 13, face = "bold", hjust = .5),
      axis.text.x         = element_text(size = 12),
      axis.text.y         = element_blank(),
      strip.text          = element_text(size = 12, hjust = .5, face = "bold"),
      legend.text         = element_text(size = 13, face = "bold", hjust = .5),
      legend.position     = "none",
      plot.margin         = unit(c(0.5,0,0.5,0), "lines")
    )
)

# save
grDevices::cairo_pdf(
  file = paste0(outputdir, "cancer_type_est_mvpa_lpa", ".pdf"),
  width = 7,
  height = 12,
)
ggarrange(
  plot_comp_cancer_type_mvpa,
  plot_comp_cancer_type_lpa,
  nrow = 2
)
dev.off()

grDevices::cairo_pdf(
  file = paste0(outputdir, "cancer_type_est_sleep_sb", ".pdf"),
  width = 7,
  height = 12,
)
ggarrange(
  plot_comp_cancer_type_sleep,
  plot_comp_cancer_type_sb,
  nrow = 2
)
dev.off()

# olivia's thesis -----------
# plot separate behaviours ------------
(plot_comp_cancer_type_sleep <- 
   ggplot(comp_cancer_type_adj[part == "Sleep Period"], aes(x = cancer_before_acc_type, y = Mean)) +
   geom_hline(aes(yintercept = yintercept), linewidth = 0.5, linetype = 2, colour = "#a8a8a8") +
   geom_pointrange(aes(ymin = CI_low,
                       ymax = CI_high, colour = cancer_before_acc_type)) +
   geom_text(aes(y = text_position + 8, label = Sig, colour = cancer_before_acc_type), 
             size = 5.5, 
             # position = position_dodge2(width = 1),
             show.legend = FALSE) +
   # facet_wrap(~part, scales = "free") +
   scale_colour_manual(values = pal_type) +
   labs(x = "", title = "Sleep Period", colour = "") +
   coord_flip() +      
   theme_ipsum() +
   theme(
     axis.ticks          = element_blank(),
     panel.background    = element_rect(fill = "transparent", colour = "black", linewidth = 0.5),
     plot.background     = element_rect(fill = "transparent", colour = NA),
     panel.grid.major.x  = element_blank(),
     panel.grid.minor    = element_blank(),
     # strip.text          = element_text(size = 12, hjust = .5, face = "bold"),
     axis.title.x        = element_blank(),
     plot.title          = element_text(size = 12, hjust = .5, face = "bold"),
     legend.position     = "none"
   )
)

grDevices::cairo_pdf(
  file = paste0(outputdir, "cancer_type_adj_sleep", ".pdf"),
  width = 6,
  height = 5,
)
plot_comp_cancer_type_sleep
dev.off()

(plot_comp_cancer_type_mvpa <- 
    ggplot(comp_cancer_type_adj[part == "Moderate-to-vigorous physical activity"], aes(x = cancer_before_acc_type, y = Mean)) +
    geom_hline(aes(yintercept = yintercept), linewidth = 0.5, linetype = 2, colour = "#a8a8a8") +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_before_acc_type)) +
    geom_text(aes(y = text_position + 8, label = Sig, colour = cancer_before_acc_type), 
              size = 5.5, 
              # position = position_dodge2(width = 1),
              show.legend = FALSE) +
    # facet_wrap(~part, scales = "free") +
    scale_colour_manual(values = pal_type) +
    labs(x = "", title = "Moderate-to-vigorous physical activity", colour = "") +
    coord_flip() +      
    theme_ipsum() +
    theme(
      axis.ticks          = element_blank(),
      panel.background    = element_rect(fill = "transparent", colour = "black", linewidth = 0.5),
      plot.background     = element_rect(fill = "transparent", colour = NA),
      panel.grid.major.x  = element_blank(),
      panel.grid.minor    = element_blank(),
      # strip.text          = element_text(size = 12, hjust = .5, face = "bold"),
      axis.title.x        = element_blank(),
      plot.title          = element_text(size = 12, hjust = .5, face = "bold"),
      legend.position     = "none"
    )
)

grDevices::cairo_pdf(
  file = paste0(outputdir, "cancer_type_adj_mvpa", ".pdf"),
  width = 6,
  height = 5,
)
plot_comp_cancer_type_mvpa
dev.off()

(plot_comp_cancer_type_lpa <- 
    ggplot(comp_cancer_type_adj[part == "Light physical activity"], aes(x = cancer_before_acc_type, y = Mean)) +
    geom_hline(aes(yintercept = yintercept), linewidth = 0.5, linetype = 2, colour = "#a8a8a8") +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_before_acc_type)) +
    geom_text(aes(y = text_position + 8, label = Sig, colour = cancer_before_acc_type), 
              size = 5.5, 
              # position = position_dodge2(width = 1),
              show.legend = FALSE) +
    # facet_wrap(~part, scales = "free") +
    scale_colour_manual(values = pal_type) +
    labs(x = "", title = "Light physical activity", colour = "") +
    coord_flip() +      
    theme_ipsum() +
    theme(
      axis.ticks          = element_blank(),
      panel.background    = element_rect(fill = "transparent", colour = "black", linewidth = 0.5),
      plot.background     = element_rect(fill = "transparent", colour = NA),
      panel.grid.major.x  = element_blank(),
      panel.grid.minor    = element_blank(),
      # strip.text          = element_text(size = 12, hjust = .5, face = "bold"),
      axis.title.x        = element_blank(),
      plot.title          = element_text(size = 12, hjust = .5, face = "bold"),
      legend.position     = "none"
    )
)

grDevices::cairo_pdf(
  file = paste0(outputdir, "cancer_type_adj_lpa", ".pdf"),
  width = 6,
  height = 5,
)
plot_comp_cancer_type_lpa
dev.off()

(plot_comp_cancer_type_sb <- 
    ggplot(comp_cancer_type_adj[part == "Sedentary behaviour"], aes(x = cancer_before_acc_type, y = Mean)) +
    geom_hline(aes(yintercept = yintercept), linewidth = 0.5, linetype = 2, colour = "#a8a8a8") +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_before_acc_type)) +
    geom_text(aes(y = text_position + 8, label = Sig, colour = cancer_before_acc_type), 
              size = 5.5, 
              # position = position_dodge2(width = 1),
              show.legend = FALSE) +
    # facet_wrap(~part, scales = "free") +
    scale_colour_manual(values = pal_type) +
    labs(x = "", title = "Sedentary behaviour", colour = "") +
    coord_flip() +      
    theme_ipsum() +
    theme(
      axis.ticks          = element_blank(),
      panel.background    = element_rect(fill = "transparent", colour = "black", linewidth = 0.5),
      plot.background     = element_rect(fill = "transparent", colour = NA),
      panel.grid.major.x  = element_blank(),
      panel.grid.minor    = element_blank(),
      # strip.text          = element_text(size = 12, hjust = .5, face = "bold"),
      axis.title.x        = element_blank(),
      plot.title          = element_text(size = 12, hjust = .5, face = "bold"),
      legend.position     = "none"
    )
)

grDevices::cairo_pdf(
  file = paste0(outputdir, "cancer_type_adj_sb", ".pdf"),
  width = 6,
  height = 5,
)
plot_comp_cancer_type_sb
dev.off()

# scatter plot lpa and sb --------
wided <- reshape(comp_cancer_type_adj, direction = "wide", idvar = "cancer_before_acc_type", timevar = "part")
names(wided)

wided[, cancer_before_acc_type := factor(cancer_before_acc_type, ordered = TRUE,
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

(lpa_sb_cor <- ggplot(wided, aes(x = `Mean.Light physical activity`, y = `Mean.Sedentary behaviour`, colour = cancer_before_acc_type)) +
    geom_point(size = 3) +
    geom_text(aes(x = `Mean.Light physical activity` + 1, label = cancer_before_acc_type, colour = cancer_before_acc_type),
              size = 3
              # position = position_dodge2(width = 3, preserve = "single", reverse = TRUE)
    ) +
    scale_colour_manual(values = pal_type) +
    labs(x = "Sedentary Behaviour", y = "Light Physical Activity", colour = "") +
    coord_flip() +
    theme_ipsum() +
    theme(
      axis.ticks          = element_blank(),
      panel.background    = element_rect(fill = "transparent", colour = "black", linewidth = 0.5),
      plot.background     = element_rect(fill = "transparent", colour = NA),
      axis.title.x        = element_text(size = 13, hjust = .5, face = "bold"),
      axis.title.y        = element_text(size = 13, hjust = .5, face = "bold"),
      panel.grid.major    = element_blank(),
      panel.grid.minor    = element_blank(),
      legend.position     = "none"
    )
)

grDevices::cairo_pdf(
  file = paste0(outputdir, "cancer_type_adj_lpa_sb", ".pdf"),
  width = 9,
  height = 10,
)
lpa_sb_cor
dev.off()

