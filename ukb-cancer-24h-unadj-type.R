source("ukb-cancer-24h-setup.R")
source(paste0(redir, "ukb_utils.R"))
# source("ukb-cancer-24h-data.R")

# main model --------
# fit_cancer_type_other_unadj <- brmcoda(clr_cancer_acc,
#                                  mvbind(ilr1, ilr2, ilr3) ~ cancer_before_acc_type_other,
#                                  warmup = 500, chains = 4, cores = 4, backend = "cmdstanr"
# )
# saveRDS(fit_cancer_type_other_unadj, paste0(outputdir, "fit_cancer_type_other_unadj", ".RDS"))

# predicted posteriors ------------
fit_cancer_type_other_unadj <- readRDS(paste0(outputdir, "fit_cancer_type_other_unadj", ".RDS"))

# reference grid
d_cancer_type_other_unadj <- emmeans::ref_grid(fit_cancer_type_other_unadj$model)@grid

# predict
pred_cancer_type_other_unadj <- fitted(fit_cancer_type_other_unadj, newdata = d_cancer_type_other_unadj, scale = "response", summary = FALSE)

# summarise by cancer types
pred_cancer_type_other_unadj <- apply(pred_cancer_type_other_unadj, c(1), function(x)  cbind(d_cancer_type_other_unadj, x))
pred_cancer_type_other_unadj <- lapply(pred_cancer_type_other_unadj, function(d) {
  
  parts <- c("sleep", "mvpa", "lpa", "sb")
  parts0 <- c("V1", "V2", "V3", "V4")
  
  d <- as.data.table(d)
  d[, cancer_weight := sum(.wgt.), by = cancer_before_acc_type_other]
  d[, denominator := sum(cancer_weight)]
  d[, cancer_denominator := denominator - 
      (d[cancer_before_acc_type_other == "Others"]$cancer_weight[[1]])*nrow(d[cancer_before_acc_type_other == "Others"]) - 
      (d[cancer_before_acc_type_other == "Healthy"]$cancer_weight[[1]])*nrow(d[cancer_before_acc_type_other == "Healthy"])]
  
  d[, (parts) := lapply(.SD, mean), by = cancer_before_acc_type_other, .SDcols =  parts0]
  
  d[, sleep_weighted := ifelse(cancer_before_acc_type_other %nin% c("Others", "Healthy"), sleep*(cancer_weight/cancer_denominator), NA)]
  d[, mvpa_weighted := ifelse(cancer_before_acc_type_other %nin% c("Others", "Healthy"), mvpa*(cancer_weight/cancer_denominator), NA)]
  d[, lpa_weighted := ifelse(cancer_before_acc_type_other %nin% c("Others", "Healthy"), lpa*(cancer_weight/cancer_denominator), NA)]
  d[, sb_weighted := ifelse(cancer_before_acc_type_other %nin% c("Others", "Healthy"), sb*(cancer_weight/cancer_denominator), NA)]
  
  d[, (paste0(parts, "_cancer")) := lapply(.SD, sum, na.rm = TRUE), .SDcols = paste0(parts, "_weighted")]
  
  d <- rbind(d,
             data.table(cancer_before_acc_type_other = "Cancer"),
             fill = TRUE
  )
  
  d[cancer_before_acc_type_other == "Cancer", sleep := d$sleep_cancer[[1]]]
  d[cancer_before_acc_type_other == "Cancer", mvpa := d$mvpa_cancer[[1]]]
  d[cancer_before_acc_type_other == "Cancer", lpa := d$lpa_cancer[[1]]]
  d[cancer_before_acc_type_other == "Cancer", sb := d$sb_cancer[[1]]]
  
  # constrast cancer types vs healthy
  d[, sleep_vs_healthy := sleep - d[cancer_before_acc_type_other == "Healthy"]$sleep[1]]
  d[, mvpa_vs_healthy := mvpa - d[cancer_before_acc_type_other == "Healthy"]$mvpa[1]]
  d[, lpa_vs_healthy := lpa - d[cancer_before_acc_type_other == "Healthy"]$lpa[1]]
  d[, sb_vs_healthy := sb - d[cancer_before_acc_type_other == "Healthy"]$sb[1]]
  
  # constrast cancer types vs others
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
pred_cancer_type_other_unadj <- as.data.frame(abind(pred_cancer_type_other_unadj, along = 1))
pred_cancer_type_other_unadj <- split(pred_cancer_type_other_unadj, pred_cancer_type_other_unadj$cancer_before_acc_type_other)

## estimated means  ----------------------
pred_comp_cancer_type_other_unadj <- lapply(pred_cancer_type_other_unadj, function(l) {
  l <- as.data.frame(l[, c("sleep", "mvpa", "lpa", "sb")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, function(p) describe_posterior(p, centrality = "mean", ci = 0.99))
  l <- Map(cbind, l, part = names(l))
  l <- rbindlist(l)
  l
})
pred_comp_cancer_type_other_unadj <- Map(cbind, pred_comp_cancer_type_other_unadj, cancer_before_acc_type_other = names(pred_comp_cancer_type_other_unadj))
pred_comp_cancer_type_other_unadj <- rbindlist(pred_comp_cancer_type_other_unadj)

## contrasts --------------------
### vs healthy
diff1_comp_cancer_type_other_unadj <- lapply(pred_cancer_type_other_unadj, function(l) {
  l <- as.data.frame(l[, c("sleep_vs_healthy", "mvpa_vs_healthy", "lpa_vs_healthy", "sb_vs_healthy")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, function(p) describe_posterior(p, centrality = "mean", ci = 0.99))
  l <- Map(cbind, l, constrast = names(l))
  l <- rbindlist(l)
  l
})
diff1_comp_cancer_type_other_unadj <- Map(cbind, diff1_comp_cancer_type_other_unadj, diff_from_healthy = names(diff1_comp_cancer_type_other_unadj))
diff1_comp_cancer_type_other_unadj <- rbindlist(diff1_comp_cancer_type_other_unadj)

setnames(diff1_comp_cancer_type_other_unadj, "Mean", "Mean_diff_ref_healthy")
setnames(diff1_comp_cancer_type_other_unadj, "CI_low", "CI_low_diff_ref_healthy")
setnames(diff1_comp_cancer_type_other_unadj, "CI_high", "CI_high_diff_ref_healthy")

### vs others
diff2_comp_cancer_type_other_unadj <- lapply(pred_cancer_type_other_unadj, function(l) {
  l <- as.data.frame(l[, c("sleep_vs_others", "mvpa_vs_others", "lpa_vs_others", "sb_vs_others")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, function(p) describe_posterior(p, centrality = "mean", ci = 0.99))
  l <- Map(cbind, l, constrast = names(l))
  l <- rbindlist(l)
  l
})
diff2_comp_cancer_type_other_unadj <- Map(cbind, diff2_comp_cancer_type_other_unadj, diff_from_healthy = names(diff2_comp_cancer_type_other_unadj))
diff2_comp_cancer_type_other_unadj <- rbindlist(diff2_comp_cancer_type_other_unadj)

setnames(diff2_comp_cancer_type_other_unadj, "Mean", "Mean_diff_ref_others")
setnames(diff2_comp_cancer_type_other_unadj, "CI_low", "CI_low_diff_ref_others")
setnames(diff2_comp_cancer_type_other_unadj, "CI_high", "CI_high_diff_ref_others")

# all results  ------------------------
comp_cancer_type_other_unadj <- cbind(
  pred_comp_cancer_type_other_unadj[, .(Mean, CI_low, CI_high, part, cancer_before_acc_type_other)],
  diff1_comp_cancer_type_other_unadj[, .(Mean_diff_ref_healthy, CI_low_diff_ref_healthy, CI_high_diff_ref_healthy)],
  diff2_comp_cancer_type_other_unadj[, .(Mean_diff_ref_others, CI_low_diff_ref_others, CI_high_diff_ref_others)]
)

# add sig indicators
comp_cancer_type_other_unadj[, nonsig_healthy := between(0, comp_cancer_type_other_unadj$CI_low_diff_ref_healthy, comp_cancer_type_other_unadj$CI_high_diff_ref_healthy)]
comp_cancer_type_other_unadj[, nonsig_others := between(0, comp_cancer_type_other_unadj$CI_low_diff_ref_others, comp_cancer_type_other_unadj$CI_high_diff_ref_others)]

comp_cancer_type_other_unadj[, sig_ref_healthy := ifelse(nonsig_healthy == FALSE, "$^a$", "$\\phantom{^a}$")]
comp_cancer_type_other_unadj[, sig_ref_others := ifelse(nonsig_others == FALSE, "$^b$", "$\\phantom{^b}$")]

# leave healthy empty
comp_cancer_type_other_unadj[, sig_ref_healthy := ifelse(cancer_before_acc_type_other == "Healthy", "$\\phantom{^a}$", sig_ref_healthy)]
comp_cancer_type_other_unadj[, sig_ref_others := ifelse(cancer_before_acc_type_other == "Healthy", "$\\phantom{^b}$", sig_ref_others)]

comp_cancer_type_other_unadj[, yintercept_healthy := NA]
comp_cancer_type_other_unadj[, yintercept_healthy := ifelse(part == "sleep", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Healthy" & part == "sleep"]$Mean, yintercept_healthy)]
comp_cancer_type_other_unadj[, yintercept_healthy := ifelse(part == "mvpa", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Healthy" & part == "mvpa"]$Mean, yintercept_healthy)]
comp_cancer_type_other_unadj[, yintercept_healthy := ifelse(part == "lpa", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Healthy" & part == "lpa"]$Mean, yintercept_healthy)]
comp_cancer_type_other_unadj[, yintercept_healthy := ifelse(part == "sb", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Healthy" & part == "sb"]$Mean, yintercept_healthy)]

comp_cancer_type_other_unadj[, ci_low_healthy := NA]
comp_cancer_type_other_unadj[, ci_low_healthy := ifelse(part == "sleep", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Healthy" & part == "sleep"]$CI_low, ci_low_healthy)]
comp_cancer_type_other_unadj[, ci_low_healthy := ifelse(part == "mvpa", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Healthy" & part == "mvpa"]$CI_low, ci_low_healthy)]
comp_cancer_type_other_unadj[, ci_low_healthy := ifelse(part == "lpa", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Healthy" & part == "lpa"]$CI_low, ci_low_healthy)]
comp_cancer_type_other_unadj[, ci_low_healthy := ifelse(part == "sb", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Healthy" & part == "sb"]$CI_low, ci_low_healthy)]

comp_cancer_type_other_unadj[, ci_high_healthy := NA]
comp_cancer_type_other_unadj[, ci_high_healthy := ifelse(part == "sleep", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Healthy" & part == "sleep"]$CI_high, ci_high_healthy)]
comp_cancer_type_other_unadj[, ci_high_healthy := ifelse(part == "mvpa", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Healthy" & part == "mvpa"]$CI_high, ci_high_healthy)]
comp_cancer_type_other_unadj[, ci_high_healthy := ifelse(part == "lpa", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Healthy" & part == "lpa"]$CI_high, ci_high_healthy)]
comp_cancer_type_other_unadj[, ci_high_healthy := ifelse(part == "sb", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Healthy" & part == "sb"]$CI_high, ci_high_healthy)]

comp_cancer_type_other_unadj[, yintercept_others := NA]
comp_cancer_type_other_unadj[, yintercept_others := ifelse(part == "sleep", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Others" & part == "sleep"]$Mean, yintercept_others)]
comp_cancer_type_other_unadj[, yintercept_others := ifelse(part == "mvpa", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Others" & part == "mvpa"]$Mean, yintercept_others)]
comp_cancer_type_other_unadj[, yintercept_others := ifelse(part == "lpa", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Others" & part == "lpa"]$Mean, yintercept_others)]
comp_cancer_type_other_unadj[, yintercept_others := ifelse(part == "sb", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Others" & part == "sb"]$Mean, yintercept_others)]

comp_cancer_type_other_unadj[, ci_low_others := NA]
comp_cancer_type_other_unadj[, ci_low_others := ifelse(part == "sleep", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Others" & part == "sleep"]$CI_low, ci_low_others)]
comp_cancer_type_other_unadj[, ci_low_others := ifelse(part == "mvpa", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Others" & part == "mvpa"]$CI_low, ci_low_others)]
comp_cancer_type_other_unadj[, ci_low_others := ifelse(part == "lpa", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Others" & part == "lpa"]$CI_low, ci_low_others)]
comp_cancer_type_other_unadj[, ci_low_others := ifelse(part == "sb", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Others" & part == "sb"]$CI_low, ci_low_others)]

comp_cancer_type_other_unadj[, ci_high_others := NA]
comp_cancer_type_other_unadj[, ci_high_others := ifelse(part == "sleep", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Others" & part == "sleep"]$CI_high, ci_high_others)]
comp_cancer_type_other_unadj[, ci_high_others := ifelse(part == "mvpa", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Others" & part == "mvpa"]$CI_high, ci_high_others)]
comp_cancer_type_other_unadj[, ci_high_others := ifelse(part == "lpa", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Others" & part == "lpa"]$CI_high, ci_high_others)]
comp_cancer_type_other_unadj[, ci_high_others := ifelse(part == "sb", comp_cancer_type_other_unadj[cancer_before_acc_type_other == "Others" & part == "sb"]$CI_high, ci_high_others)]

# n
comp_cancer_type_other_unadj[, Cases := NA]
comp_cancer_type_other_unadj[, Cases := ifelse(cancer_before_acc_type_other == "Healthy", "13 722", Cases)]
comp_cancer_type_other_unadj[, Cases := ifelse(cancer_before_acc_type_other == "Cancer", "10 152", Cases)]
comp_cancer_type_other_unadj[, Cases := ifelse(cancer_before_acc_type_other == "Others", "67 478", Cases)]

comp_cancer_type_other_unadj[, Cases := ifelse(cancer_before_acc_type_other == "Blood", "412", Cases)]
comp_cancer_type_other_unadj[, Cases := ifelse(cancer_before_acc_type_other == "Breast", "1 911", Cases)]
comp_cancer_type_other_unadj[, Cases := ifelse(cancer_before_acc_type_other == "Colorectal", "391", Cases)]
comp_cancer_type_other_unadj[, Cases := ifelse(cancer_before_acc_type_other == "Endocrine Gland", "79", Cases)]
comp_cancer_type_other_unadj[, Cases := ifelse(cancer_before_acc_type_other == "Gastrointestinal Tract", "100", Cases)]
comp_cancer_type_other_unadj[, Cases := ifelse(cancer_before_acc_type_other == "Genitourinary", "296", Cases)]
comp_cancer_type_other_unadj[, Cases := ifelse(cancer_before_acc_type_other == "Gynaecological", "401", Cases)]
comp_cancer_type_other_unadj[, Cases := ifelse(cancer_before_acc_type_other == "Head & Neck", "115", Cases)]
comp_cancer_type_other_unadj[, Cases := ifelse(cancer_before_acc_type_other == "Lung", "57", Cases)]
comp_cancer_type_other_unadj[, Cases := ifelse(cancer_before_acc_type_other == "Melanoma", "395", Cases)]
comp_cancer_type_other_unadj[, Cases := ifelse(cancer_before_acc_type_other == "Multiple Primary", "1 714", Cases)]
comp_cancer_type_other_unadj[, Cases := ifelse(cancer_before_acc_type_other == "Other Cancer", "148", Cases)]
comp_cancer_type_other_unadj[, Cases := ifelse(cancer_before_acc_type_other == "Other Skin", "3 017", Cases)]
comp_cancer_type_other_unadj[, Cases := ifelse(cancer_before_acc_type_other == "Prostate", "1 116", Cases)]

comp_cancer_type_other_unadj[, cancer_before_acc_type_other := ifelse(cancer_before_acc_type_other == "Others", "Other Conditions", cancer_before_acc_type_other)]
comp_cancer_type_other_unadj[, cancer_before_acc_type_other := ifelse(cancer_before_acc_type_other == "Other Skin", "   Skin (non-melanoma)", cancer_before_acc_type_other)]
comp_cancer_type_other_unadj[, cancer_before_acc_type_other := ifelse(cancer_before_acc_type_other == "Lung", "   Lung", cancer_before_acc_type_other)]
comp_cancer_type_other_unadj[, cancer_before_acc_type_other := ifelse(cancer_before_acc_type_other == "Gastrointestinal Tract", "   Gastrointestinal Tract", cancer_before_acc_type_other)]
comp_cancer_type_other_unadj[, cancer_before_acc_type_other := ifelse(cancer_before_acc_type_other == "Blood", "   Blood", cancer_before_acc_type_other)]
comp_cancer_type_other_unadj[, cancer_before_acc_type_other := ifelse(cancer_before_acc_type_other == "Head & Neck", "   Head and Neck", cancer_before_acc_type_other)]
comp_cancer_type_other_unadj[, cancer_before_acc_type_other := ifelse(cancer_before_acc_type_other == "Colorectal", "   Colorectal", cancer_before_acc_type_other)]
comp_cancer_type_other_unadj[, cancer_before_acc_type_other := ifelse(cancer_before_acc_type_other == "Gynaecological", "   Gynaecological", cancer_before_acc_type_other)]
comp_cancer_type_other_unadj[, cancer_before_acc_type_other := ifelse(cancer_before_acc_type_other == "Other Cancer", "   Other Cancer", cancer_before_acc_type_other)]
comp_cancer_type_other_unadj[, cancer_before_acc_type_other := ifelse(cancer_before_acc_type_other == "Genitourinary", "   Genitourinary", cancer_before_acc_type_other)]
comp_cancer_type_other_unadj[, cancer_before_acc_type_other := ifelse(cancer_before_acc_type_other == "Endocrine Gland", "   Endocrine Gland", cancer_before_acc_type_other)]
comp_cancer_type_other_unadj[, cancer_before_acc_type_other := ifelse(cancer_before_acc_type_other == "Breast", "   Breast", cancer_before_acc_type_other)]
comp_cancer_type_other_unadj[, cancer_before_acc_type_other := ifelse(cancer_before_acc_type_other == "Melanoma", "   Melanoma", cancer_before_acc_type_other)]
comp_cancer_type_other_unadj[, cancer_before_acc_type_other := ifelse(cancer_before_acc_type_other == "Prostate", "   Prostate", cancer_before_acc_type_other)]
comp_cancer_type_other_unadj[, cancer_before_acc_type_other := ifelse(cancer_before_acc_type_other == "Multiple Primary", "   Multiple Primary", cancer_before_acc_type_other)]

## sort by MVPA
## healthy as bottom in plot
# comp_cancer_type_other_unadj[, cancer_before_acc_type_other := factor(cancer_before_acc_type_other, ordered = TRUE,
#                                                             levels = c(
#                                                               "Healthy",
#                                                               "Non-cancer Conditions",
#                                                               "Skin (non-melanoma)",
#                                                               "Prostate",
#                                                               "Melanoma",
#                                                               "Endocrine Gland",
#                                                               "Breast",
#                                                               "Other Cancer",
#                                                               "Colorectal",
#                                                               "Gynaecological",
#                                                               "Genitourinary",
#                                                               "Head & Neck",
#                                                               "Blood",
#                                                               "Gastrointestinal Tract",
#                                                               "Lung",
#                                                               "Multiple Primary"))]
# # healthy as top in plot
comp_cancer_type_other_unadj[, cancer_before_acc_type_other := factor(cancer_before_acc_type_other, ordered = TRUE,
                                                                    levels = c(
                                                                      "   Multiple Primary",
                                                                      "   Lung",
                                                                      "   Gastrointestinal Tract",
                                                                      "   Blood",
                                                                      "   Head and Neck",
                                                                      "   Genitourinary",
                                                                      "   Gynaecological",
                                                                      "   Colorectal",
                                                                      "   Other Cancer",
                                                                      "   Breast",
                                                                      "   Endocrine Gland",
                                                                      "   Melanoma",
                                                                      "   Prostate",
                                                                      "   Skin (non-melanoma)",
                                                                      "Cancer",
                                                                      "Other Conditions",
                                                                      "Healthy"
                                                                    ))]

comp_cancer_type_other_unadj[, part := ifelse(part == "sleep", "Sleep period", part)]
comp_cancer_type_other_unadj[, part := ifelse(part == "mvpa", "Moderate-to-vigorous physical activity", part)]
comp_cancer_type_other_unadj[, part := ifelse(part == "lpa", "Light physical activity", part)]
comp_cancer_type_other_unadj[, part := ifelse(part == "sb", "Sedentary behaviour", part)]

comp_cancer_type_other_unadj[, sig_position := min(CI_low), by = part]
comp_cancer_type_other_unadj[, est_position := max(CI_high), by = part]

comp_cancer_type_other_unadj[, estimates := paste0(round(Mean, 0), "[", round(CI_low, 0), ", ", round(CI_high, 0), "]")]
comp_cancer_type_other_unadj[, est_sig := paste0(estimates, " ", sig_ref_healthy, sig_ref_others)]
# comp_cancer_type_other_unadj[, est_sig := paste0(estimates, " ", str_replace_na(sig_ref_healthy, " "), str_replace_na(sig_ref_others, " "))]

# for tables
comp_cancer_type_other_unadj[, estimates_contrast_healthy := paste0(round(Mean_diff_ref_healthy, 2), "[", round(CI_low_diff_ref_healthy, 2), ", ", round(CI_high_diff_ref_healthy, 2), "]")]
comp_cancer_type_other_unadj[, estimates_contrast_others := paste0(round(Mean_diff_ref_others, 2), "[", round(CI_low_diff_ref_others, 2), ", ", round(CI_high_diff_ref_others, 2), "]")]

## plot -----------------------
(plot_comp_cancer_type_other_sleep_unadj <- 
   ggplot(comp_cancer_type_other_unadj[part == "Sleep period"], aes(x = cancer_before_acc_type_other, y = Mean)) +
   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#D9E0DD", alpha = 0.075) +
   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#F2F2F2", alpha = 0.2) +
   geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#708885") +
   geom_hline(aes(yintercept = yintercept_others), linewidth = 0.5, linetype= "dashed", colour = "#A9A9A9") +
   geom_pointrange(aes(ymin = CI_low,
                       ymax = CI_high, colour = cancer_before_acc_type_other), size = 0.25, linewidth = 0.5) +
   geom_text(aes(y = 500, label = cancer_before_acc_type_other),
             hjust = 0, nudge_x = 0, 
             family = "Arial Narrow", size = 3,
             show.legend = FALSE) +
   geom_text(aes(y = 632.5, label = TeX(est_sig, output = "character")), parse = TRUE,
             hjust = 0.5, nudge_x = 0, 
             family = "Arial Narrow", size = 3,
             show.legend = FALSE) +
   geom_text(aes(y = 650, label = Cases),
             hjust = 1, nudge_x = 0, 
             family = "Arial Narrow", size = 3,
             show.legend = FALSE) +
   geom_segment(aes(x = 0, yend = 650), col = "black", linewidth = 0.5) +
   geom_segment(aes(x = 0, yend = 500), col = "black", linewidth = 0.5) +
   scale_y_continuous(limits = c(500, 650),
                      breaks = c(500, 650),
                      name = "Sleep period (mins/day)") +
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
     # axis.line.x         = element_line(linewidth = 0.5, linetype= "dashed", colour = "#708885"),
     axis.title.x        = element_text(size = 10, face = "bold", hjust = .5),
     axis.text.x         = element_text(size = 9),
     axis.text.y         = element_blank(),
     strip.text          = element_text(size = , hjust = .5, face = "bold"),
     legend.text         = element_text(size = 10, face = "bold", hjust = .5),
     legend.position     = "none",
     plot.margin         = unit(c(0.5,0,0,0), "lines")
   )
)
(plot_comp_cancer_type_other_mvpa_unadj <- 
    ggplot(comp_cancer_type_other_unadj[part == "Moderate-to-vigorous physical activity"], aes(x = cancer_before_acc_type_other, y = Mean)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#D9E0DD", alpha = 0.075) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#F2F2F2", alpha = 0.2) +
    geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#708885") +
    geom_hline(aes(yintercept = yintercept_others), linewidth = 0.5, linetype= "dashed", colour = "#A9A9A9") +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_before_acc_type_other), size = 0.25, linewidth = 0.5) +
    geom_text(aes(y = 0, label = cancer_before_acc_type_other),
              hjust = 0, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_text(aes(y = 44, label = TeX(est_sig, output = "character")), parse = TRUE,
              hjust = 0.5, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_text(aes(y = 50, label = Cases),
              hjust = 1, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_segment(aes(x = 0, yend = 0), col = "black", linewidth = 0.5) +
    geom_segment(aes(x = 0, yend = 50), col = "black", linewidth = 0.5) +
    scale_y_continuous(limits = c(0, 50),
                       breaks = c(0, 50),
                       name = "Moderate-to-vigorous physical activity (mins/day)") +
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

(plot_comp_cancer_type_other_lpa_unadj <- 
    ggplot(comp_cancer_type_other_unadj[part == "Light physical activity"], aes(x = cancer_before_acc_type_other, y = Mean)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#D9E0DD", alpha = 0.075) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#F2F2F2", alpha = 0.2) +
    geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#708885") +
    geom_hline(aes(yintercept = yintercept_others), linewidth = 0.5, linetype= "dashed", colour = "#A9A9A9") +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_before_acc_type_other), size = 0.25, linewidth = 0.5) +
    geom_text(aes(y = 200, label = cancer_before_acc_type_other),
              hjust = 0, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_text(aes(y = 377.5, label = TeX(est_sig, output = "character")), parse = TRUE,
              hjust = 0.5, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_text(aes(y = 400, label = Cases),
              hjust = 1, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_segment(aes(x = 0, yend = 200), col = "black", linewidth = 0.5) +
    geom_segment(aes(x = 0, yend = 400), col = "black", linewidth = 0.5) +
    scale_y_continuous(limits = c(200, 400),
                       breaks = c(200, 400),
                       name = "Light physical activity (mins/day)") +
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

(plot_comp_cancer_type_other_sb_unadj <- 
    ggplot(comp_cancer_type_other_unadj[part == "Sedentary behaviour"], aes(x = cancer_before_acc_type_other, y = Mean)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#D9E0DD", alpha = 0.075) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#F2F2F2", alpha = 0.2) +
    geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#708885") +
    geom_hline(aes(yintercept = yintercept_others), linewidth = 0.5, linetype= "dashed", colour = "#A9A9A9") +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_before_acc_type_other), size = 0.25, linewidth = 0.5) +
    geom_text(aes(y = 500, label = cancer_before_acc_type_other),
              hjust = 0, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_text(aes(y = 677.5, label = TeX(est_sig, output = "character")), parse = TRUE,
              hjust = 0.5, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_text(aes(y = 700, label = Cases),
              hjust = 1, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_segment(aes(x = 0, yend = 700), col = "black", linewidth = 0.5) +
    geom_segment(aes(x = 0, yend = 500), col = "black", linewidth = 0.5) +
    scale_y_continuous(limits = c(500, 700),
                       breaks = c(500, 700),
                       name = "Sedentary (mins/day)") +
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

grDevices::cairo_pdf(
  file = paste0(outputdir, "cancer_type_other_est_unadj", ".pdf"),
  width = 8,
  height = 14,
)

ggarrange(
  plot_comp_cancer_type_other_mvpa_unadj,
  plot_comp_cancer_type_other_lpa_unadj,
  plot_comp_cancer_type_other_sb_unadj,
  plot_comp_cancer_type_other_sleep_unadj,
  nrow = 4
)
dev.off()

grDevices::png(
  file = paste0(outputdir, "cancer_type_other_est_unadj", ".png"),
  width = 8500,
  height = 12000,
  res = 900
)

ggarrange(
  plot_comp_cancer_type_other_mvpa_unadj,
  plot_comp_cancer_type_other_lpa_unadj,
  plot_comp_cancer_type_other_sb_unadj,
  plot_comp_cancer_type_other_sleep_unadj,
  nrow = 4
)
dev.off()

# estimates for tables --------------
comp_cancer_type_other_unadj[part == "Moderate-to-vigorous physical activity", .(cancer_before_acc_type_other, estimates_contrast_healthy)]
comp_cancer_type_other_unadj[part == "Light physical activity", .(cancer_before_acc_type_other, estimates_contrast_healthy)]
comp_cancer_type_other_unadj[part == "Sedentary behaviour", .(cancer_before_acc_type_other, estimates_contrast_healthy)]
comp_cancer_type_other_unadj[part == "Sleep period", .(cancer_before_acc_type_other, estimates_contrast_healthy)]

comp_cancer_type_other_unadj[part == "Moderate-to-vigorous physical activity", .(cancer_before_acc_type_other, estimates_contrast_others)]
comp_cancer_type_other_unadj[part == "Light physical activity", .(cancer_before_acc_type_other, estimates_contrast_others)]
comp_cancer_type_other_unadj[part == "Sedentary behaviour", .(cancer_before_acc_type_other, estimates_contrast_others)]
comp_cancer_type_other_unadj[part == "Sleep period", .(cancer_before_acc_type_other, estimates_contrast_others)]

