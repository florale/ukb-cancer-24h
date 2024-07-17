source("ukb-cancer-24h-utils.R")
source(paste0(redir, "ukb_utils.R"))
# source("ukb-cancer-24h-data.R")

# main model --------
fit_cancer_time_since_diag_unadj <- brmcoda(clr_cancer_acc,
                                          mvbind(ilr1, ilr2, ilr3) ~ cancer_time_since_diag,
                                          # save_pars = save_pars(all = TRUE),
                                          warmup = 500, chains = 4, cores = 4, backend = "cmdstanr"
)
saveRDS(fit_cancer_time_since_diag_unadj, paste0( "fit_cancer_time_since_diag_unadj", ".RDS"))
# saveRDS(fit_cancer_time_since_diag_unadj, paste0(outputdir, "fit_cancer_time_since_diag_unadj", ".RDS"))

# estimates ------
fit_cancer_time_since_diag_unadj <- readRDS(paste0(outputdir, "fit_cancer_time_since_diag_unadj", ".RDS"))

# reference grid
d_cancer_time_since_diag_unadj <- emmeans::ref_grid(fit_cancer_time_since_diag_unadj$model)@grid

# predict
pred_cancer_time_since_diag_unadj <- fitted(fit_cancer_time_since_diag_unadj, newdata = d_cancer_time_since_diag_unadj, scale = "response", summary = FALSE)

# summarise by cancer group
pred_cancer_time_since_diag_unadj <- apply(pred_cancer_time_since_diag_unadj, c(1), function(x)  cbind(d_cancer_time_since_diag_unadj, x))
pred_cancer_time_since_diag_unadj <- lapply(pred_cancer_time_since_diag_unadj, function(d) {
  d <- as.data.table(d)
  
  # estimated means by time
  d[, sleep := mean(V1), by = cancer_time_since_diag]
  d[, mvpa  := mean(V2), by = cancer_time_since_diag]
  d[, lpa   := mean(V3), by = cancer_time_since_diag]
  d[, sb    := mean(V4), by = cancer_time_since_diag]
  
  # contrast vs healthy
  d[, sleep_vs_healthy := sleep - d[cancer_time_since_diag == "Healthy"]$sleep]
  d[, mvpa_vs_healthy := mvpa - d[cancer_time_since_diag == "Healthy"]$mvpa]
  d[, lpa_vs_healthy := lpa - d[cancer_time_since_diag == "Healthy"]$lpa]
  d[, sb_vs_healthy := sb - d[cancer_time_since_diag == "Healthy"]$sb]
  
  # pairwise contrast within cancer 
  d[, sleep_lessthan1_vs_1to5 := d[cancer_time_since_diag == "1-5 years since diagnosis"]$sleep[1] - d[cancer_time_since_diag == "Less than 1 year since diagnosis"]$sleep[1]]
  d[, sleep_lessthan1_vs_morethan5 := d[cancer_time_since_diag == "More than 5 years since diagnosis"]$sleep[1] - d[cancer_time_since_diag == "Less than 1 year since diagnosis"]$sleep[1]]
  d[, sleep_1to5_vs_morethan5 := d[cancer_time_since_diag == "More than 5 years since diagnosis"]$sleep[1] - d[cancer_time_since_diag == "1-5 years since diagnosis"]$sleep[1]]
  
  d[, mvpa_lessthan1_vs_1to5 := d[cancer_time_since_diag == "1-5 years since diagnosis"]$mvpa[1] - d[cancer_time_since_diag == "Less than 1 year since diagnosis"]$mvpa[1]]
  d[, mvpa_lessthan1_vs_morethan5 := d[cancer_time_since_diag == "More than 5 years since diagnosis"]$mvpa[1] - d[cancer_time_since_diag == "Less than 1 year since diagnosis"]$mvpa[1]]
  d[, mvpa_1to5_vs_morethan5 := d[cancer_time_since_diag == "More than 5 years since diagnosis"]$mvpa[1] - d[cancer_time_since_diag == "1-5 years since diagnosis"]$mvpa[1]]
  
  d[, lpa_lessthan1_vs_1to5 := d[cancer_time_since_diag == "1-5 years since diagnosis"]$lpa[1] - d[cancer_time_since_diag == "Less than 1 year since diagnosis"]$lpa[1]]
  d[, lpa_lessthan1_vs_morethan5 := d[cancer_time_since_diag == "More than 5 years since diagnosis"]$lpa[1] - d[cancer_time_since_diag == "Less than 1 year since diagnosis"]$lpa[1]]
  d[, lpa_1to5_vs_morethan5 := d[cancer_time_since_diag == "More than 5 years since diagnosis"]$lpa[1] - d[cancer_time_since_diag == "1-5 years since diagnosis"]$lpa[1]]
  
  d[, sb_lessthan1_vs_1to5 := d[cancer_time_since_diag == "1-5 years since diagnosis"]$sb[1] - d[cancer_time_since_diag == "Less than 1 year since diagnosis"]$sb[1]]
  d[, sb_lessthan1_vs_morethan5 := d[cancer_time_since_diag == "More than 5 years since diagnosis"]$sb[1] - d[cancer_time_since_diag == "Less than 1 year since diagnosis"]$sb[1]]
  d[, sb_1to5_vs_morethan5 := d[cancer_time_since_diag == "More than 5 years since diagnosis"]$sb[1] - d[cancer_time_since_diag == "1-5 years since diagnosis"]$sb[1]]
  
  d <- d[, .(cancer_time_since_diag, 
             sleep, mvpa, lpa, sb,
             
             sleep_vs_healthy, mvpa_vs_healthy, lpa_vs_healthy, sb_vs_healthy,
             
             sleep_lessthan1_vs_1to5, sleep_lessthan1_vs_morethan5, sleep_1to5_vs_morethan5,
             mvpa_lessthan1_vs_1to5, mvpa_lessthan1_vs_morethan5, mvpa_1to5_vs_morethan5,
             lpa_lessthan1_vs_1to5, lpa_lessthan1_vs_morethan5, lpa_1to5_vs_morethan5,
             sb_lessthan1_vs_1to5, sb_lessthan1_vs_morethan5, sb_1to5_vs_morethan5
             
  )]
  d <- unique(d)
  d
})

# assemble back to summarise posteriors
pred_cancer_time_since_diag_unadj <- as.data.frame(abind(pred_cancer_time_since_diag_unadj, along = 1))
pred_cancer_time_since_diag_unadj <- split(pred_cancer_time_since_diag_unadj, pred_cancer_time_since_diag_unadj$cancer_time_since_diag)
# est means ----
pred_comp_cancer_time_since_diag_unadj <- lapply(pred_comp_cancer_time_since_diag_unadj, function(l) {
  l <- as.data.frame(l[, c("sleep", "mvpa", "lpa", "sb")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, part = names(l))
  l <- rbindlist(l)
  l
})
pred_comp_cancer_time_since_diag_unadj <- Map(cbind, pred_comp_cancer_time_since_diag_unadj, cancer_time_since_diag = names(pred_comp_cancer_time_since_diag_unadj))
pred_comp_cancer_time_since_diag_unadj <- rbindlist(pred_comp_cancer_time_since_diag_unadj)
## Contrasts --------------------
### vs healthy
diff1_comp_cancer_time_since_diag_unadj <- lapply(pred_cancer_time_since_diag_unadj, function(l) {
  l <- as.data.frame(l[, c("sleep_vs_healthy", "mvpa_vs_healthy", "lpa_vs_healthy", "sb_vs_healthy")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, contrast = names(l))
  l <- rbindlist(l)
  l
})
diff1_comp_cancer_time_since_diag_unadj <- Map(cbind, diff1_comp_cancer_time_since_diag_unadj, cancer_time_since_diag = names(diff1_comp_cancer_time_since_diag_unadj))
diff1_comp_cancer_time_since_diag_unadj <- rbindlist(diff1_comp_cancer_time_since_diag_unadj)

setnames(diff1_comp_cancer_time_since_diag_unadj, "Mean", "Mean_diff_ref_healthy")
setnames(diff1_comp_cancer_time_since_diag_unadj, "CI_low", "CI_low_diff_ref_healthy")
setnames(diff1_comp_cancer_time_since_diag_unadj, "CI_high", "CI_high_diff_ref_healthy")

### pairwise cancer
diff2_comp_cancer_time_since_diag_unadj <- lapply(pred_cancer_time_since_diag_unadj, function(l) {
  l <- as.data.frame(l[, c(
    "sleep_lessthan1_vs_1to5", "sleep_lessthan1_vs_morethan5", "sleep_1to5_vs_morethan5",
    "mvpa_lessthan1_vs_1to5", "mvpa_lessthan1_vs_morethan5", "mvpa_1to5_vs_morethan5",
    "lpa_lessthan1_vs_1to5", "lpa_lessthan1_vs_morethan5", "lpa_1to5_vs_morethan5",
    "sb_lessthan1_vs_1to5", "sb_lessthan1_vs_morethan5", "sb_1to5_vs_morethan5"
  )])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, contrast = names(l))
  l <- rbindlist(l)
  l
})
# diff2_comp_cancer_time_since_diag_unadj <- Map(cbind, diff2_comp_cancer_time_since_diag_unadj, cancer_time_since_diag = names(diff2_comp_cancer_time_since_diag_unadj))
diff2_comp_cancer_time_since_diag_unadj <- rbindlist(diff2_comp_cancer_time_since_diag_unadj)

diff2_comp_cancer_time_since_diag_unadj[, cancer_time_since_diag := NA]
diff2_comp_cancer_time_since_diag_unadj[, cancer_time_since_diag := ifelse(grepl("lessthan1_vs_1to5", diff2_comp_cancer_time_since_diag_unadj$contrast), "Less than 1 year since diagnosis", cancer_time_since_diag)]
diff2_comp_cancer_time_since_diag_unadj[, cancer_time_since_diag := ifelse(grepl("1to5_vs_morethan5", diff2_comp_cancer_time_since_diag_unadj$contrast), "1-5 years since diagnosis", cancer_time_since_diag)]
diff2_comp_cancer_time_since_diag_unadj[, cancer_time_since_diag := ifelse(grepl("lessthan1_vs_morethan5", diff2_comp_cancer_time_since_diag_unadj$contrast), "More than 5 years since diagnosis", cancer_time_since_diag)]

setnames(diff2_comp_cancer_time_since_diag_unadj, "Mean", "Mean_diff_ref_cancer")
setnames(diff2_comp_cancer_time_since_diag_unadj, "CI_low", "CI_low_diff_ref_cancer")
setnames(diff2_comp_cancer_time_since_diag_unadj, "CI_high", "CI_high_diff_ref_cancer")

# take only unique values
diff2_comp_cancer_time_since_diag_unadj <- unique(diff2_comp_cancer_time_since_diag_unadj)

# add id to merge later
diff2_comp_cancer_time_since_diag_unadj[, id := 1:.N, by = cancer_time_since_diag]

# All results  ------------------------
ç <- cbind(
  pred_comp_cancer_time_since_diag_unadj[, .(Mean, CI_low, CI_high, part, cancer_time_since_diag)],
  pred_comp_cancer_time_since_diag_unadj[, .(Mean_diff_ref_healthy, CI_low_diff_ref_healthy, CI_high_diff_ref_healthy)]
)
comp_cancer_time_since_diag_unadj[, id := 1:.N, by = cancer_time_since_diag]

comp_cancer_time_since_diag_unadj <- merge(
  comp_cancer_time_since_diag_unadj,
  diff2_comp_cancer_time_since_diag_unadj[, .(Mean_diff_ref_cancer,
                                            CI_low_diff_ref_cancer,
                                            CI_high_diff_ref_cancer,
                                            cancer_time_since_diag,
                                            id
  )],
  by = c("cancer_time_since_diag", "id"),
  all.x = T
)


# add sig indicators
comp_cancer_time_since_diag_unadj[, nonsig_vs_healthy := between(0, comp_cancer_time_since_diag_unadj$CI_low_diff_ref_healthy, comp_cancer_time_since_diag_unadj$CI_high_diff_ref_healthy)]
comp_cancer_time_since_diag_unadj[, sig_ref_healthy := ifelse(nonsig_vs_healthy == FALSE & Mean_diff_ref_healthy != 0, paste(intToUtf8(0x2217)), "")]

comp_cancer_time_since_diag_unadj[, nonsig_vs_cancer := between(0, comp_cancer_time_since_diag_unadj$CI_low_diff_ref_cancer, comp_cancer_time_since_diag_unadj$CI_high_diff_ref_cancer)]

comp_cancer_time_since_diag_unadj[, sig_ref_cancer := NA]
comp_cancer_time_since_diag_unadj[, sig_ref_cancer := ifelse(nonsig_vs_cancer == FALSE & !is.na(Mean_diff_ref_cancer) & cancer_time_since_diag == "Less than 1 year since diagnosis",
                                                           paste(intToUtf8(8224)), sig_ref_cancer)]
comp_cancer_time_since_diag_unadj[, sig_ref_cancer := ifelse(nonsig_vs_cancer == FALSE & !is.na(Mean_diff_ref_cancer) & cancer_time_since_diag == "1-5 years since diagnosis",
                                                           paste(intToUtf8(8225)), sig_ref_cancer)] 
comp_cancer_time_since_diag_unadj[, sig_ref_cancer := ifelse(nonsig_vs_cancer == FALSE & !is.na(Mean_diff_ref_cancer) & cancer_time_since_diag == "More than 5 years since diagnosis",
                                                           paste(intToUtf8(0x00A7)), sig_ref_cancer)] 

# sort by time since diagnoses
# comp_cancer_time_since_diag_unadj[, cancer_time_since_diag := factor(cancer_time_since_diag, ordered = TRUE,
#                                              levels = c(
#                                                "Healthy",
#                                                "Less than 1 year since diagnosis",
#                                                "1-5 years since diagnosis",
#                                                "More than 5 years since diagnosis"))]

# healthy in first row of plot
comp_cancer_time_since_diag_unadj[, cancer_time_since_diag := ifelse(cancer_time_since_diag == "More than 5 years since diagnosis", "More than 5 years since Cancer diagnosis", cancer_time_since_diag)]
comp_cancer_time_since_diag_unadj[, cancer_time_since_diag := ifelse(cancer_time_since_diag == "1-5 years since diagnosis", "1-5 years since Cancer diagnosis", cancer_time_since_diag)]
comp_cancer_time_since_diag_unadj[, cancer_time_since_diag := ifelse(cancer_time_since_diag == "Less than 1 year since diagnosis", "Less than 1 year since Cancer diagnosis", cancer_time_since_diag)]

comp_cancer_time_since_diag_unadj[, cancer_time_since_diag := factor(cancer_time_since_diag, ordered = TRUE,
                                                                   levels = c(
                                                                     "More than 5 years since Cancer diagnosis",
                                                                     "1-5 years since Cancer diagnosis",
                                                                     "Less than 1 year since Cancer diagnosis",
                                                                     "Healthy"))]

comp_cancer_time_since_diag_unadj[, yintercept := NA]
comp_cancer_time_since_diag_unadj[, yintercept := ifelse(part == "sleep", comp_cancer_time_since_diag_unadj[cancer_time_since_diag == "Healthy" & part == "sleep"]$Mean, yintercept)]
comp_cancer_time_since_diag_unadj[, yintercept := ifelse(part == "mvpa", comp_cancer_time_since_diag_unadj[cancer_time_since_diag == "Healthy" & part == "mvpa"]$Mean, yintercept)]
comp_cancer_time_since_diag_unadj[, yintercept := ifelse(part == "lpa", comp_cancer_time_since_diag_unadj[cancer_time_since_diag == "Healthy" & part == "lpa"]$Mean, yintercept)]
comp_cancer_time_since_diag_unadj[, yintercept := ifelse(part == "sb", comp_cancer_time_since_diag_unadj[cancer_time_since_diag == "Healthy" & part == "sb"]$Mean, yintercept)]

comp_cancer_time_since_diag_unadj[, part := ifelse(part == "sleep", "Sleep", part)]
comp_cancer_time_since_diag_unadj[, part := ifelse(part == "mvpa", "Moderate-to-vigorous physical activity", part)]
comp_cancer_time_since_diag_unadj[, part := ifelse(part == "lpa", "Light physical activity", part)]
comp_cancer_time_since_diag_unadj[, part := ifelse(part == "sb", "Sedentary behaviour", part)]

comp_cancer_time_since_diag_unadj[, text_position := max(CI_high), by = part]

(plot_comp_cancer_time_since_diag_unadj <- 
    ggplot(comp_cancer_time_since_diag_unadj, aes(x = cancer_time_since_diag, y = Mean, group = part)) +
    geom_hline(aes(yintercept = yintercept), linewidth = 0.5, linetype = 2, colour = "#a8a8a8") +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_time_since_diag), size = 0.75, linewidth = 0.75) +
    geom_text(aes(y = text_position + 3, label = sig_ref_healthy, colour = cancer_time_since_diag), 
              size = 6, nudge_x = 0.2, 
              show.legend = FALSE) +
    geom_text(aes(y = text_position + 4, label = sig_ref_cancer, colour = cancer_time_since_diag), 
              size = 4, nudge_x = 0.2,
              show.legend = FALSE) +
    facet_wrap(~part, scales = "free") +
    scale_colour_manual(values = pal_time) +
    # scale_colour_jco() +
    labs(x = "", y = "", colour = "") +
    coord_flip() +
    theme_ipsum() +
    theme(
      axis.ticks          = element_blank(),
      panel.background    = element_rect(fill = "transparent", colour = "black", linewidth = 0.5),
      plot.background     = element_rect(fill = "transparent", colour = NA),
      panel.grid.major.x  = element_blank(),
      panel.grid.minor    = element_blank(),
      strip.text          = element_text(size = 13, hjust = .5, face = "bold"),
      legend.text         = element_text(size = 14, face = "bold", hjust = .5),
      legend.position     = "none"
    )
)

grDevices::cairo_pdf(
  file = paste0(outputdir, "cancer_time_since_diag_unadj", ".pdf"),
  width = 11.5,
  height = 8,
)
plot_comp_cancer_time_since_diag_unadj
dev.off()

