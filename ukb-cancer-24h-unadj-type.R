source("ukb-cancer-24h-utils.R")
source(paste0(redir, "ukb_utils.R"))
source("ukb-cancer-24h-data.R")

# main model --------
fit_cancer_type_unadj <- brmcoda(clr_cancer_acc,
                               mvbind(ilr1, ilr2, ilr3) ~ cancer_before_acc_type,
                               # save_pars = save_pars(all = TRUE),
                               warmup = 500, chains = 4, cores = 4, backend = "cmdstanr"
)
saveRDS(fit_cancer_type_unadj, paste0("fit_cancer_type_unadj", ".RDS"))
# saveRDS(fit_cancer_type_unadj, paste0(outputdir, "fit_cancer_type_unadj", ".RDS"))

# Predicted posteriors ------------
fit_cancer_type_unadj <- readRDS(paste0(outputdir, "fit_cancer_type_unadj", ".RDS"))

# reference grid
d_cancer_type_unadj <- emmeans::ref_grid(fit_cancer_type_unadj$model)@grid

# predict
pred_cancer_type_unadj <- fitted(fit_cancer_type_unadj, newdata = d_cancer_type_unadj, scale = "response", summary = FALSE)

# summarise by cancer types
pred_cancer_type_unadj <- apply(pred_cancer_type_unadj, c(1), function(x)  cbind(d_cancer_type_unadj, x))
pred_cancer_type_unadj <- lapply(pred_cancer_type_unadj, function(d) {
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
pred_cancer_type_unadj <- as.data.frame(abind(pred_cancer_type_unadj, along = 1))
pred_cancer_type_unadj <- split(pred_cancer_type_unadj, pred_cancer_type_unadj$cancer_before_acc_type)

## Est means  ----------------------
pred_comp_cancer_type_unadj <- lapply(pred_cancer_type_unadj, function(l) {
  l <- as.data.frame(l[, c("sleep", "mvpa", "lpa", "sb")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, part = names(l))
  l <- rbindlist(l)
  l
})
pred_comp_cancer_type_unadj <- Map(cbind, pred_comp_cancer_type_unadj, cancer_before_acc_type = names(pred_comp_cancer_type_unadj))
pred_comp_cancer_type_unadj <- rbindlist(pred_comp_cancer_type_unadj)

## Contrasts --------------------
diff_comp_cancer_type_unadj <- lapply(pred_cancer_type_unadj, function(l) {
  l <- as.data.frame(l[, c("sleep_constrast", "mvpa_constrast", "lpa_constrast", "sb_constrast")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, constrast = names(l))
  l <- rbindlist(l)
  l
})
diff_comp_cancer_type_unadj <- Map(cbind, diff_comp_cancer_type_unadj, diff_from_healthy = names(diff_comp_cancer_type_unadj))
diff_comp_cancer_type_unadj <- rbindlist(diff_comp_cancer_type_unadj)

setnames(diff_comp_cancer_type_unadj, "Mean", "Mean_diff")
setnames(diff_comp_cancer_type_unadj, "CI_low", "CI_low_diff")
setnames(diff_comp_cancer_type_unadj, "CI_high", "CI_high_diff")

# All results  ------------------------
comp_cancer_type_unadj <- cbind(
  pred_comp_cancer_type_unadj[, .(Mean, CI_low, CI_high, part, cancer_before_acc_type)],
  diff_comp_cancer_type_unadj[, .(Mean_diff, CI_low_diff, CI_high_diff)]
)

# add sig indicators
comp_cancer_type_unadj[, nonsig := between(0, comp_cancer_type_unadj$CI_low_diff, comp_cancer_type_unadj$CI_high_diff)]
comp_cancer_type_unadj[, Sig := ifelse(nonsig == FALSE, paste(intToUtf8(0x2217)), "")]

# sort by MVPA
# healthy as bottom in plot
# comp_cancer_type_unadj[, cancer_before_acc_type := factor(cancer_before_acc_type, ordered = TRUE,
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
comp_cancer_type_unadj[, cancer_before_acc_type := factor(cancer_before_acc_type, ordered = TRUE,
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

comp_cancer_type_unadj[, yintercept := NA]
comp_cancer_type_unadj[, yintercept := ifelse(part == "sleep", comp_cancer_type_unadj[cancer_before_acc_type == "Healthy" & part == "sleep"]$Mean, yintercept)]
comp_cancer_type_unadj[, yintercept := ifelse(part == "mvpa", comp_cancer_type_unadj[cancer_before_acc_type == "Healthy" & part == "mvpa"]$Mean, yintercept)]
comp_cancer_type_unadj[, yintercept := ifelse(part == "lpa", comp_cancer_type_unadj[cancer_before_acc_type == "Healthy" & part == "lpa"]$Mean, yintercept)]
comp_cancer_type_unadj[, yintercept := ifelse(part == "sb", comp_cancer_type_unadj[cancer_before_acc_type == "Healthy" & part == "sb"]$Mean, yintercept)]

comp_cancer_type_unadj[, part := ifelse(part == "sleep", "Sleep", part)]
comp_cancer_type_unadj[, part := ifelse(part == "mvpa", "Moderate-to-vigorous physical activity", part)]
comp_cancer_type_unadj[, part := ifelse(part == "lpa", "Light physical activity", part)]
comp_cancer_type_unadj[, part := ifelse(part == "sb", "Sedentary behaviour", part)]

comp_cancer_type_unadj[, text_position := max(CI_high), by = part]

(plot_comp_cancer_type_unadj <- 
    ggplot(comp_cancer_type_unadj, aes(x = cancer_before_acc_type, y = Mean, group = part)) +
    geom_hline(aes(yintercept = yintercept), linewidth = 0.5, linetype = 2, colour = "#a8a8a8") +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_before_acc_type)) +
    geom_text(aes(y = text_position + 8, label = Sig, colour = cancer_before_acc_type), 
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
  file = paste0(outputdir, "cancer_type_unadj", ".pdf"),
  width = 11,
  height = 10,
)
plot_comp_cancer_type_unadj
dev.off()





