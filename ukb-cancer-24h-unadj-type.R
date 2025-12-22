source("ukb-cancer-24h-setup.R")
source(paste0(redir, "ukb_utils.R"))
# source("ukb-cancer-24h-data.R")

# main model --------
# fit_cancer_type_other_unadj <- brmcoda(clr_cancer_acc,
#   mvbind(z1_1, z2_1, z3_1) ~ cancer_before_acc_type_other,
#   warmup = 500, chains = 4, cores = 4, backend = "cmdstanr"
# )
# saveRDS(fit_cancer_type_other_unadj, paste0(outputdir, "fit_cancer_type_other_unadj", ".RDS"))

# predicted posteriors ------------
fit_cancer_type_other_unadj <- readRDS(paste0(outputdir, "fit_cancer_type_other_unadj", ".RDS"))

# reference grid
d_cancer_type_other_unadj <- emmeans::ref_grid(fit_cancer_type_other_unadj$model)@grid

# predict
pred_cancer_type_other_unadj <- fitted(fit_cancer_type_other_unadj, newdata = d_cancer_type_other_unadj, scale = "response", summary = FALSE)

# summarise by cancer types
pred_cancer_type_other_unadj <- apply(pred_cancer_type_other_unadj, c(1), function(x) cbind(d_cancer_type_other_unadj, x))
pred_cancer_type_other_unadj <- lapply(pred_cancer_type_other_unadj, function(d) {
  parts <- c("sleep", "mvpa", "lpa", "sb")
  parts0 <- c("tsleep_comp", "tmvpa_comp", "tlpa_comp", "tsb_comp")

  d <- as.data.table(d)

  d[, (paste0(parts)) :=
    lapply(.SD, function(x) {
      weighted.mean(x, .wgt.)
    }),
  .SDcols = parts0, by = cancer_before_acc_type_other
  ]

  d[, cancer_wgt := sum(.wgt.), by = cancer_before_acc_type_other]

  d[cancer_before_acc_type_other %nin% c("Healthy", "Others"), (paste0(parts, "_cancer")) := lapply(.SD, function(x) weighted.mean(x, cancer_wgt)), .SDcols = paste0(parts)]

  d <- rbind(d,
    data.table(cancer_before_acc_type_other = "Cancer"),
    fill = TRUE
  )

  d[cancer_before_acc_type_other == "Cancer", (parts) := d[cancer_before_acc_type_other %nin% c("Healthy", "Others"),
    lapply(.SD, function(x) unique(na.omit(x))),
    .SDcols = paste0(parts, "_cancer")
  ]]

  # contrast vs healthy
  d[, (paste0(parts, "_vs_healthy")) :=
    Map(`-`, .SD, d[cancer_before_acc_type_other == "Healthy",
      .SD,
      .SDcols = parts
    ][1]),
  .SDcols = parts
  ]

  # contrast vs others
  d[, (paste0(parts, "_vs_others")) :=
    Map(`-`, .SD, d[cancer_before_acc_type_other == "Others",
      .SD,
      .SDcols = parts
    ][1]),
  .SDcols = parts
  ]

  d <- d[, .(
    cancer_before_acc_type_other,
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
comp_cancer_type_other_unadj[, cancer_before_acc_type_other := factor(cancer_before_acc_type_other,
  ordered = TRUE,
  levels = c(
    "   Multiple Primary",
    "   Lung",
    "   Gastrointestinal Tract",
    "   Gynaecological",
    "   Breast",
    "   Blood",
    "   Endocrine Gland",
    "   Head and Neck",
    "   Colorectal",
    "   Other Cancer",
    "   Genitourinary",
    "   Melanoma",
    "   Skin (non-melanoma)",
    "   Prostate",
    "Cancer",
    "Other Conditions",
    "Healthy"
  )
)]

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
comp_cancer_type_other_unadj[, estimates_contrast_healthy := paste0(format(round(Mean_diff_ref_healthy, 1)), " [", format(round(CI_low_diff_ref_healthy, 1)), ", ", format(round(CI_high_diff_ref_healthy, 1)), "]")]
comp_cancer_type_other_unadj[, estimates_contrast_others := paste0(format(round(Mean_diff_ref_others, 1)), " [", format(round(CI_low_diff_ref_others, 1)), ", ", format(round(CI_high_diff_ref_others, 1)), "]")]

## plot -----------------------
plot_specs <- data.table(
  part = c(
    "Sleep period",
    "Moderate-to-vigorous physical activity",
    "Light physical activity",
    "Sedentary behaviour"
  ),
  y_limits = list(c(450, 650), c(0, 55), c(200, 400), c(500, 700)),
  text_y = c(450, 0, 200, 500),
  sig_y = c(615, 45, 365.5, 665),
  cases_y = c(635, 51, 385, 685),
  surv_y = c(650, 55, 400, 700),
  seg_y = list(c(450, 650), c(0, 55), c(200, 400), c(500, 700)),
  ylab = c(
    "Sleep period (min/day)",
    "Moderate-to-vigorous physical activity (min/day)",
    "Light physical activity (min/day)",
    "Sedentary behaviour (min/day)"
  )
)

plot_list <- lapply(seq_len(nrow(plot_specs)), function(i) {
  spec <- plot_specs[i]
  seg_vals <- spec$seg_y[[1]]
  limits <- spec$y_limits[[1]]
  data_i <- comp_cancer_type_other_unadj[part == spec$part]

  ggplot(data_i, aes(x = cancer_before_acc_type_other, y = Mean)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others),
      fill = "#F2F2F2", alpha = 0.2
    ) +
    geom_hline(aes(yintercept = yintercept_others),
      linewidth = 0.5, linetype = "dashed", colour = "#A9A9A9"
    ) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy),
      fill = "#CAD4CF", alpha = 0.075
    ) + ## D9E0DD
    geom_hline(aes(yintercept = yintercept_healthy),
      linewidth = 0.5, linetype = "dashed", colour = "#708885"
    ) +
    geom_pointrange(aes(ymin = CI_low, ymax = CI_high, colour = cancer_before_acc_type_other),
      size = 0.15, linewidth = 0.5
    ) +
    geom_text(aes(y = spec$text_y, label = cancer_before_acc_type_other),
      hjust = 0, family = "Arial Narrow", size = 2.25, show.legend = FALSE
    ) +
    geom_text(aes(y = spec$sig_y, label = TeX(est_sig, output = "character")),
      parse = TRUE, hjust = 0.5, family = "Arial Narrow", size = 2.25, show.legend = FALSE
    ) +
    geom_text(aes(y = spec$cases_y, label = Cases),
      hjust = 1, family = "Arial Narrow", size = 2.25, show.legend = FALSE
    ) +
    # geom_text(aes(y = spec$surv_y, label = Survival),
    #   hjust = 1, family = "Arial Narrow", size = 2.25, show.legend = FALSE
    # ) +
    # annotate("segment",
    #          x = 0.5,
    #          xend = length(unique(data_i$cancer_before_acc_type_other)) + 0.5,
    #          y = seg_vals[1], yend = seg_vals[1], linewidth = 0.5) +
    # annotate("segment",
    #          x = 0.5,
    #          xend = length(unique(data_i$cancer_before_acc_type_other)) + 0.5,
    #          y = seg_vals[2], yend = seg_vals[2], linewidth = 0.5) +
    geom_segment(aes(x = 0, yend = seg_vals[1]), col = "black", linewidth = 0.5) +
    geom_segment(aes(x = 0, yend = seg_vals[2]), col = "black", linewidth = 0.5) +
    scale_y_continuous(limits = limits, breaks = limits, name = spec$ylab) +
    scale_colour_manual(values = pal_type) +
    labs(x = "", y = "", colour = "") +
    coord_flip() +
    theme_ipsum() +
    theme(
      axis.ticks = element_blank(),
      plot.background = element_rect(fill = "transparent", colour = NA, linewidth = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = 8, face = "bold", hjust = .5, margin = margin(t = -9)),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_blank(),
      strip.text = element_text(size = 8, hjust = .5, face = "bold"),
      legend.text = element_text(size = 8, face = "bold", hjust = .5),
      legend.position = "none",
      plot.margin = unit(c(0.5, 0, 1, 0), "lines")
    )
})

names(plot_list) <- tolower(gsub(" ", "_", plot_specs$part))
plot_comp_cancer_type_other_sleep_unadj <- plot_list[["sleep_period"]]
plot_comp_cancer_type_other_mvpa_unadj <- plot_list[["moderate-to-vigorous_physical_activity"]]
plot_comp_cancer_type_other_lpa_unadj <- plot_list[["light_physical_activity"]]
plot_comp_cancer_type_other_sb_unadj <- plot_list[["sedentary_behaviour"]]

grDevices::cairo_pdf(
  file = paste0(outputdir, "cancer_type_other_est_unadj", ".pdf"),
  width = 6,
  height = 8,
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
  width = 6000,
  height = 8000,
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
comp_cancer_type_other_unadj[part == "Moderate-to-vigorous physical activity", .(cancer_before_acc_type_other, estimates_contrast_healthy)][order(-cancer_before_acc_type_other)]
comp_cancer_type_other_unadj[part == "Light physical activity", .(cancer_before_acc_type_other, estimates_contrast_healthy)][order(-cancer_before_acc_type_other)]
comp_cancer_type_other_unadj[part == "Sedentary behaviour", .(cancer_before_acc_type_other, estimates_contrast_healthy)][order(-cancer_before_acc_type_other)]
comp_cancer_type_other_unadj[part == "Sleep period", .(cancer_before_acc_type_other, estimates_contrast_healthy)][order(-cancer_before_acc_type_other)]

comp_cancer_type_other_unadj[part == "Moderate-to-vigorous physical activity", .(cancer_before_acc_type_other, estimates_contrast_others)][order(-cancer_before_acc_type_other)]
comp_cancer_type_other_unadj[part == "Light physical activity", .(cancer_before_acc_type_other, estimates_contrast_others)][order(-cancer_before_acc_type_other)]
comp_cancer_type_other_unadj[part == "Sedentary behaviour", .(cancer_before_acc_type_other, estimates_contrast_others)][order(-cancer_before_acc_type_other)]
comp_cancer_type_other_unadj[part == "Sleep period", .(cancer_before_acc_type_other, estimates_contrast_others)][order(-cancer_before_acc_type_other)]
