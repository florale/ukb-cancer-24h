source("ukb-cancer-24h-setup.R")
source(paste0(redir, "ukb_utils.R"))
source("ukb-cancer-24h-data.R")

# main model --------
# fit_cancer_time_since_diag_other_f_adj <- brmcoda(clr_cancer_acc,
#                                           mvbind(z1_1, z2_1, z3_1) ~ cancer_time_since_diag_other_f +
#                                             # + other_conds_at_acc +
#                                             s(age_at_acc) + sex + white + working + edu + never_smoked + current_drinker + s(deprivation) + season,
#                                           # save_pars = save_pars(all = TRUE),
#                                           warmup = 500, chains = 4, cores = 4, backend = "cmdstanr"
# )
# saveRDS(fit_cancer_time_since_diag_other_f_adj, paste0(outputdir, "fit_cancer_time_since_diag_other_f_adj", ".RDS"))
# summary(fit_cancer_time_since_diag_other_f_adj)

# fit_cancer_time_since_diag_f_adj <- brmcoda(clr_cancer_acc,
#  mvbind(z1_1, z2_1, z3_1) ~ cancer_time_since_diag_f +
#   # + other_conds_at_acc +
#   s(age_at_acc) + sex + white + working + edu + never_smoked + current_drinker + s(deprivation),
#  # save_pars = save_pars(all = TRUE),
#  warmup = 500, chains = 4, cores = 4, backend = "cmdstanr"
# )
# saveRDS(fit_cancer_time_since_diag_f_adj, paste0(outputdir, "fit_cancer_time_since_diag_f_adj", ".RDS"))
# summary(fit_cancer_time_since_diag_f_adj)

# predicted posteriors ------------
fit_cancer_time_since_diag_f_adj <- readRDS(paste0(outputdir, "fit_cancer_time_since_diag_f_adj", ".RDS"))

# reference grid
d_cancer_time_since_diag_f_adj <- emmeans::ref_grid(fit_cancer_time_since_diag_f_adj$model)@grid

# predict
pred_cancer_time_since_diag_f_adj <- fitted(fit_cancer_time_since_diag_f_adj, newdata = d_cancer_time_since_diag_f_adj, scale = "response", summary = FALSE)

# weight by length equal to the number of observations in the dataset
# summarise by cancer group
pred_cancer_time_since_diag_f_adj <- apply(pred_cancer_time_since_diag_f_adj, c(1), function(x) cbind(d_cancer_time_since_diag_f_adj, x))
pred_cancer_time_since_diag_f_adj <- lapply(pred_cancer_time_since_diag_f_adj, function(d) {
 parts <- c("sleep", "mvpa", "lpa", "sb")
 parts0 <- c("tsleep_comp", "tmvpa_comp", "tlpa_comp", "tsb_comp")

 d <- as.data.table(d)

 # weighted mean by group
 d[, (paste0(parts)) :=
  lapply(.SD, function(x) {
   weighted.mean(x, .wgt.)
  }),
 .SDcols = parts0, by = cancer_time_since_diag_f
 ]
  d[, cancer_wgt := sum(.wgt.), by = cancer_time_since_diag_f]

  d[cancer_time_since_diag_f %nin% c("Healthy", "Others"), (paste0(parts, "_cancer")) := lapply(.SD, function(x) weighted.mean(x, cancer_wgt)), .SDcols = paste0(parts)]

 d <- rbind(d,
  data.table(cancer_time_since_diag_f = "Cancer"),
  fill = TRUE
 )

 d[cancer_time_since_diag_f == "Cancer", (parts) := d[cancer_time_since_diag_f != "Healthy",
  lapply(.SD, function(x) unique(na.omit(x))),
  .SDcols = paste0(parts, "_cancer")
 ]]

 d[, (paste0(parts, "_vs_healthy")) :=
  Map(`-`, .SD, d[cancer_time_since_diag_f == "Healthy",
   .SD,
   .SDcols = parts
  ][1]),
 .SDcols = parts
 ]

 d <- d[, .(
  cancer_time_since_diag_f,
  sleep, mvpa, lpa, sb,

  # sleep_cancer, mvpa_cancer, lpa_cancer, sb_cancer,

  sleep_vs_healthy, mvpa_vs_healthy, lpa_vs_healthy, sb_vs_healthy
  # sleep_vs_others, mvpa_vs_others, lpa_vs_others, sb_vs_others

  # sleep_lessthan1_vs_1to5, sleep_lessthan1_vs_morethan5, sleep_1to5_vs_morethan5,
  # mvpa_lessthan1_vs_1to5, mvpa_lessthan1_vs_morethan5, mvpa_1to5_vs_morethan5,
  # lpa_lessthan1_vs_1to5, lpa_lessthan1_vs_morethan5, lpa_1to5_vs_morethan5,
  # sb_lessthan1_vs_1to5, sb_lessthan1_vs_morethan5, sb_1to5_vs_morethan5
 )]
 d <- unique(d)
 d
})
saveRDS(pred_cancer_time_since_diag_f_adj, paste0(outputdir, "pred_cancer_time_since_diag_f_adj", ".RDS"))

# assemble back to summarise posteriors
pred_cancer_time_since_diag_f_adj <- as.data.frame(abind::abind(pred_cancer_time_since_diag_f_adj, along = 1))
pred_cancer_time_since_diag_f_adj <- split(pred_cancer_time_since_diag_f_adj, pred_cancer_time_since_diag_f_adj$cancer_time_since_diag_f)

## estimated means  ----------------------
pred_comp_cancer_time_since_diag_f_adj <- lapply(pred_cancer_time_since_diag_f_adj, function(l) {
 l <- as.data.frame(l[, c("sleep", "mvpa", "lpa", "sb")])
 l <- apply(l, 2, as.numeric)
 l <- apply(l, 2, function(p) describe_posterior(p, centrality = "mean", ci = 0.99))
 l <- Map(cbind, l, part = names(l))
 l <- rbindlist(l)
 l
})
pred_comp_cancer_time_since_diag_f_adj <- Map(cbind, pred_comp_cancer_time_since_diag_f_adj, cancer_time_since_diag_f = names(pred_comp_cancer_time_since_diag_f_adj))
pred_comp_cancer_time_since_diag_f_adj <- rbindlist(pred_comp_cancer_time_since_diag_f_adj)

## contrasts --------------------
### vs healthy
diff1_comp_cancer_time_since_diag_f_adj <- lapply(pred_cancer_time_since_diag_f_adj, function(l) {
 l <- as.data.frame(l[, c("sleep_vs_healthy", "mvpa_vs_healthy", "lpa_vs_healthy", "sb_vs_healthy")])
 l <- apply(l, 2, as.numeric)
 l <- apply(l, 2, function(p) describe_posterior(p, centrality = "mean", ci = 0.99))
 l <- Map(cbind, l, contrast = names(l))
 l <- rbindlist(l)
 l
})
diff1_comp_cancer_time_since_diag_f_adj <- Map(cbind, diff1_comp_cancer_time_since_diag_f_adj, cancer_time_since_diag_f = names(diff1_comp_cancer_time_since_diag_f_adj))
diff1_comp_cancer_time_since_diag_f_adj <- rbindlist(diff1_comp_cancer_time_since_diag_f_adj)

setnames(diff1_comp_cancer_time_since_diag_f_adj, "Mean", "Mean_diff_ref_healthy")
setnames(diff1_comp_cancer_time_since_diag_f_adj, "CI_low", "CI_low_diff_ref_healthy")
setnames(diff1_comp_cancer_time_since_diag_f_adj, "CI_high", "CI_high_diff_ref_healthy")

# all results  ------------------------
comp_cancer_time_since_diag_f_adj <- cbind(
 pred_comp_cancer_time_since_diag_f_adj[, .(Mean, CI_low, CI_high, part, cancer_time_since_diag_f)],
 diff1_comp_cancer_time_since_diag_f_adj[, .(Mean_diff_ref_healthy, CI_low_diff_ref_healthy, CI_high_diff_ref_healthy)]
 # diff3_comp_cancer_time_since_diag_f_adj[, .(Mean_diff_ref_others, CI_low_diff_ref_others, CI_high_diff_ref_others)]
)
comp_cancer_time_since_diag_f_adj[, id := 1:.N, by = cancer_time_since_diag_f]

# add sig indicators
comp_cancer_time_since_diag_f_adj[, nonsig_vs_healthy := between(0, comp_cancer_time_since_diag_f_adj$CI_low_diff_ref_healthy, comp_cancer_time_since_diag_f_adj$CI_high_diff_ref_healthy)]
# comp_cancer_time_since_diag_f_adj[, nonsig_vs_others  := between(0, comp_cancer_time_since_diag_f_adj$CI_low_diff_ref_others, comp_cancer_time_since_diag_f_adj$CI_high_diff_ref_others)]
comp_cancer_time_since_diag_f_adj[, nonsig_vs_cancer := between(0, comp_cancer_time_since_diag_f_adj$CI_low_diff_ref_cancer, comp_cancer_time_since_diag_f_adj$CI_high_diff_ref_cancer)]

comp_cancer_time_since_diag_f_adj[, sig_ref_healthy := ifelse(nonsig_vs_healthy == FALSE & Mean_diff_ref_healthy != 0,
 "$^a$", "$\\phantom{^a}$"
)]

# leave healthy empty
comp_cancer_time_since_diag_f_adj[, sig_ref_healthy := ifelse(cancer_time_since_diag_f == "Healthy", "$\\phantom{^a}$", sig_ref_healthy)]
# comp_cancer_time_since_diag_f_adj[, sig_ref_others := ifelse(cancer_time_since_diag_f == "Healthy", "$\\phantom{^b}$", sig_ref_others)]

comp_cancer_time_since_diag_f_adj[, yintercept_healthy := NA]
comp_cancer_time_since_diag_f_adj[, yintercept_healthy := ifelse(part == "sleep", comp_cancer_time_since_diag_f_adj[cancer_time_since_diag_f == "Healthy" & part == "sleep"]$Mean, yintercept_healthy)]
comp_cancer_time_since_diag_f_adj[, yintercept_healthy := ifelse(part == "mvpa", comp_cancer_time_since_diag_f_adj[cancer_time_since_diag_f == "Healthy" & part == "mvpa"]$Mean, yintercept_healthy)]
comp_cancer_time_since_diag_f_adj[, yintercept_healthy := ifelse(part == "lpa", comp_cancer_time_since_diag_f_adj[cancer_time_since_diag_f == "Healthy" & part == "lpa"]$Mean, yintercept_healthy)]
comp_cancer_time_since_diag_f_adj[, yintercept_healthy := ifelse(part == "sb", comp_cancer_time_since_diag_f_adj[cancer_time_since_diag_f == "Healthy" & part == "sb"]$Mean, yintercept_healthy)]

comp_cancer_time_since_diag_f_adj[, ci_low_healthy := NA]
comp_cancer_time_since_diag_f_adj[, ci_low_healthy := ifelse(part == "sleep", comp_cancer_time_since_diag_f_adj[cancer_time_since_diag_f == "Healthy" & part == "sleep"]$CI_low, ci_low_healthy)]
comp_cancer_time_since_diag_f_adj[, ci_low_healthy := ifelse(part == "mvpa", comp_cancer_time_since_diag_f_adj[cancer_time_since_diag_f == "Healthy" & part == "mvpa"]$CI_low, ci_low_healthy)]
comp_cancer_time_since_diag_f_adj[, ci_low_healthy := ifelse(part == "lpa", comp_cancer_time_since_diag_f_adj[cancer_time_since_diag_f == "Healthy" & part == "lpa"]$CI_low, ci_low_healthy)]
comp_cancer_time_since_diag_f_adj[, ci_low_healthy := ifelse(part == "sb", comp_cancer_time_since_diag_f_adj[cancer_time_since_diag_f == "Healthy" & part == "sb"]$CI_low, ci_low_healthy)]

comp_cancer_time_since_diag_f_adj[, ci_high_healthy := NA]
comp_cancer_time_since_diag_f_adj[, ci_high_healthy := ifelse(part == "sleep", comp_cancer_time_since_diag_f_adj[cancer_time_since_diag_f == "Healthy" & part == "sleep"]$CI_high, ci_high_healthy)]
comp_cancer_time_since_diag_f_adj[, ci_high_healthy := ifelse(part == "mvpa", comp_cancer_time_since_diag_f_adj[cancer_time_since_diag_f == "Healthy" & part == "mvpa"]$CI_high, ci_high_healthy)]
comp_cancer_time_since_diag_f_adj[, ci_high_healthy := ifelse(part == "lpa", comp_cancer_time_since_diag_f_adj[cancer_time_since_diag_f == "Healthy" & part == "lpa"]$CI_high, ci_high_healthy)]
comp_cancer_time_since_diag_f_adj[, ci_high_healthy := ifelse(part == "sb", comp_cancer_time_since_diag_f_adj[cancer_time_since_diag_f == "Healthy" & part == "sb"]$CI_high, ci_high_healthy)]

# n
table(model.frame(fit_cancer_time_since_diag_f_adj)$cancer_time_since_diag_f)
comp_cancer_time_since_diag_f_adj[, Cases := NA]
comp_cancer_time_since_diag_f_adj[, Cases := ifelse(cancer_time_since_diag_f == "Healthy", "13 722", Cases)]
comp_cancer_time_since_diag_f_adj[, Cases := ifelse(cancer_time_since_diag_f == "Cancer", "10 152", Cases)]
comp_cancer_time_since_diag_f_adj[, Cases := ifelse(cancer_time_since_diag_f == "Less than 1 year since diagnosis", "971", Cases)]
comp_cancer_time_since_diag_f_adj[, Cases := ifelse(cancer_time_since_diag_f == "1-2 years since diagnosis", "1 014", Cases)]
comp_cancer_time_since_diag_f_adj[, Cases := ifelse(cancer_time_since_diag_f == "2-3 years since diagnosis", "911", Cases)]
comp_cancer_time_since_diag_f_adj[, Cases := ifelse(cancer_time_since_diag_f == "3-4 years since diagnosis", "823", Cases)]
comp_cancer_time_since_diag_f_adj[, Cases := ifelse(cancer_time_since_diag_f == "4-5 years since diagnosis", "703", Cases)]
comp_cancer_time_since_diag_f_adj[, Cases := ifelse(cancer_time_since_diag_f == "5-6 years since diagnosis", "709", Cases)]
comp_cancer_time_since_diag_f_adj[, Cases := ifelse(cancer_time_since_diag_f == "6-7 years since diagnosis", "636", Cases)]
comp_cancer_time_since_diag_f_adj[, Cases := ifelse(cancer_time_since_diag_f == "7-8 years since diagnosis", "597", Cases)]
comp_cancer_time_since_diag_f_adj[, Cases := ifelse(cancer_time_since_diag_f == "8-9 years since diagnosis", "496", Cases)]
comp_cancer_time_since_diag_f_adj[, Cases := ifelse(cancer_time_since_diag_f == "9-10 years since diagnosis", "465", Cases)]
comp_cancer_time_since_diag_f_adj[, Cases := ifelse(cancer_time_since_diag_f == "More than 10 years since diagnosis", "2 827", Cases)]

## healthy in first row of plot
comp_cancer_time_since_diag_f_adj[, cancer_time_since_diag_f := ifelse(cancer_time_since_diag_f == "More than 10 years since diagnosis", "> 10 years since diagnosis", cancer_time_since_diag_f)]
# comp_cancer_time_since_diag_f_adj[, cancer_time_since_diag_f := ifelse(cancer_time_since_diag_f == "1-5 years since diagnosis",         "    1-5 years since diagnosis", cancer_time_since_diag_f)]
comp_cancer_time_since_diag_f_adj[, cancer_time_since_diag_f := ifelse(cancer_time_since_diag_f == "Less than 1 year since diagnosis",  "< 1 year since diagnosis", cancer_time_since_diag_f)]
# comp_cancer_time_since_diag_f_adj[, cancer_time_since_diag_f := ifelse(cancer_time_since_diag_f == "Others", "Other Conditions", cancer_time_since_diag_f)]

comp_cancer_time_since_diag_f_adj[, cancer_time_since_diag_f := factor(cancer_time_since_diag_f,
 ordered = TRUE,
 levels = c(
  "> 10 years since diagnosis",
  "9-10 years since diagnosis",
  "8-9 years since diagnosis",
  "7-8 years since diagnosis",
  "6-7 years since diagnosis",
  "5-6 years since diagnosis",
  "4-5 years since diagnosis",
  "3-4 years since diagnosis",
  "2-3 years since diagnosis",
  "1-2 years since diagnosis",
  "< 1 year since diagnosis",
  "Cancer",
  "Healthy"
 )
)]

comp_cancer_time_since_diag_f_adj[, part := ifelse(part == "sleep", "Sleep period", part)]
comp_cancer_time_since_diag_f_adj[, part := ifelse(part == "mvpa", "Moderate-to-vigorous physical activity", part)]
comp_cancer_time_since_diag_f_adj[, part := ifelse(part == "lpa", "Light physical activity", part)]
comp_cancer_time_since_diag_f_adj[, part := ifelse(part == "sb", "Sedentary behaviour", part)]

comp_cancer_time_since_diag_f_adj[, sig_position := min(CI_low), by = part]
comp_cancer_time_since_diag_f_adj[, est_position := max(CI_high), by = part]

comp_cancer_time_since_diag_f_adj[, estimates := paste0(round(Mean, 0), "[", round(CI_low, 0), ", ", round(CI_high, 0), "]")]
comp_cancer_time_since_diag_f_adj[, est_sig := paste0(estimates, " ", sig_ref_healthy)]
comp_cancer_time_since_diag_f_adj[, est_sig := paste0(estimates, " ", str_replace_na(sig_ref_healthy, " "))]

# for tables
comp_cancer_time_since_diag_f_adj[, estimates_contrast_healthy := paste0(round(Mean_diff_ref_healthy, 2), "[", round(CI_low_diff_ref_healthy, 2), ", ", round(CI_high_diff_ref_healthy, 2), "]")]
# comp_cancer_time_since_diag_f_adj[, estimates_contrast_others := paste0(round(Mean_diff_ref_others, 2), "[", round(CI_low_diff_ref_others, 2), ", ", round(CI_high_diff_ref_others, 2), "]")]
# comp_cancer_time_since_diag_f_adj[, estimates_contrast_cancer := paste0(round(Mean_diff_ref_cancer, 2), "[", round(CI_low_diff_ref_cancer, 2), ", ", round(CI_high_diff_ref_cancer, 2), "]")]

## plot -----------------------
plot_specs <- data.table(
 part = c(
  "Sleep period",
  "Moderate-to-vigorous physical activity",
  "Light physical activity",
  "Sedentary behaviour"
 ),
 ymin = c(500, 0, 200, 500),
 ymax = c(650, 50, 400, 650),
 ysig = c(627.5, 43.5, 371.5, 627.5),
 pal = list(pal_type, pal_type, pal_type, pal_type),
 text_y = c(500, 0, 200, 500),
 seg_y = list(c(500, 650), c(0, 50), c(200, 400), c(500, 650)),
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

 ggplot(comp_cancer_time_since_diag_f_adj[part == spec$part], aes(x = cancer_time_since_diag_f, y = Mean)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy),
   fill = "#CBD5D0", alpha = 0.1
  ) +
  geom_hline(aes(yintercept = yintercept_healthy),
   linewidth = 0.5,
   linetype = "dashed", colour = "#708885"
  ) +
  geom_pointrange(aes(ymin = CI_low, ymax = CI_high, colour = cancer_time_since_diag_f),
   size = 0.25, linewidth = 0.5
  ) +
  geom_text(aes(y = spec$text_y, label = cancer_time_since_diag_f),
   hjust = 0, nudge_x = 0, family = "Arial Narrow", size = 2.5,
   show.legend = FALSE
  ) +
  geom_text(aes(y = spec$ysig, label = TeX(est_sig, output = "character")),
   parse = TRUE,
   hjust = 0.5, nudge_x = 0, family = "Arial Narrow", size = 2.5,
   show.legend = FALSE
  ) +
  geom_text(aes(y = spec$ymax, label = Cases),
   hjust = 1, nudge_x = 0,
   family = "Arial Narrow", size = 2.5,
   show.legend = FALSE
  ) +
  geom_segment(aes(x = 0, yend = seg_vals[1]), col = "black", linewidth = 0.5) +
  geom_segment(aes(x = 0, yend = seg_vals[2]), col = "black", linewidth = 0.5) +
  scale_y_continuous(
   limits = c(spec$ymin, spec$ymax),
   breaks = c(spec$ymin, spec$ymax),
   name = spec$ylab
  ) +
  scale_colour_manual(values = pal_sen) +
  labs(x = "", y = "", colour = "") +
  coord_flip() +
  theme_ipsum() +
  theme(
   axis.ticks = element_blank(),
   plot.background = element_rect(fill = "transparent", colour = NA, linewidth = 0.5),
   panel.grid.major = element_blank(),
   panel.grid.minor = element_blank(),
   axis.title.x = element_text(size = 9, face = "bold", hjust = .5, margin = margin(t = -10)),
   axis.text.x = element_text(size = 9),
   axis.text.y = element_blank(),
   strip.text = element_text(size = 9, hjust = .5, face = "bold"),
   legend.text = element_text(size = 9, face = "bold", hjust = .5),
   legend.position = "none",
   plot.margin = unit(c(0.5, 0, 1, 0), "lines")
  )
})

names(plot_list) <- tolower(gsub(" ", "_", plot_specs$part))

plot_comp_cancer_time_since_diag_f_sleep <- plot_list[["sleep_period"]]
plot_comp_cancer_time_since_diag_f_mvpa <- plot_list[["moderate-to-vigorous_physical_activity"]]
plot_comp_cancer_time_since_diag_f_lpa <- plot_list[["light_physical_activity"]]
plot_comp_cancer_time_since_diag_f_sb <- plot_list[["sedentary_behaviour"]]

grDevices::cairo_pdf(
 file = paste0(outputdir, "cancer_time_since_diag_sen_est", ".pdf"),
 width = 6,
 height = 8,
)

ggarrange(
 plot_comp_cancer_time_since_diag_f_mvpa,
 plot_comp_cancer_time_since_diag_f_lpa,
 plot_comp_cancer_time_since_diag_f_sb,
 plot_comp_cancer_time_since_diag_f_sleep,
 nrow = 4
)
dev.off()

grDevices::png(
 file = paste0(outputdir, "cancer_time_since_diag_sen_est", ".png"),
 width = 6000,
 height = 8000,
 res = 900
)

ggarrange(
 plot_comp_cancer_time_since_diag_f_mvpa,
 plot_comp_cancer_time_since_diag_f_lpa,
 plot_comp_cancer_time_since_diag_f_sb,
 plot_comp_cancer_time_since_diag_f_sleep,
 nrow = 4
)
dev.off()

# estimates for tables --------------
comp_cancer_time_since_diag_f_adj[part == "Moderate-to-vigorous physical activity", .(cancer_time_since_diag_f, estimates_contrast_healthy)]
comp_cancer_time_since_diag_f_adj[part == "Light physical activity", .(cancer_time_since_diag_f, estimates_contrast_healthy)]
comp_cancer_time_since_diag_f_adj[part == "Sedentary behaviour", .(cancer_time_since_diag_f, estimates_contrast_healthy)]
comp_cancer_time_since_diag_f_adj[part == "Sleep period", .(cancer_time_since_diag_f, estimates_contrast_healthy)]

comp_cancer_time_since_diag_f_adj[part == "Moderate-to-vigorous physical activity", .(cancer_time_since_diag_f, estimates_contrast_others)]
comp_cancer_time_since_diag_f_adj[part == "Light physical activity", .(cancer_time_since_diag_f, estimates_contrast_others)]
comp_cancer_time_since_diag_f_adj[part == "Sedentary behaviour", .(cancer_time_since_diag_f, estimates_contrast_others)]
comp_cancer_time_since_diag_f_adj[part == "Sleep period", .(cancer_time_since_diag_f, estimates_contrast_others)]

comp_cancer_time_since_diag_f_adj[part == "Moderate-to-vigorous physical activity", .(cancer_time_since_diag_f, estimates_contrast_cancer)]
comp_cancer_time_since_diag_f_adj[part == "Light physical activity", .(cancer_time_since_diag_f, estimates_contrast_cancer)]
comp_cancer_time_since_diag_f_adj[part == "Sedentary behaviour", .(cancer_time_since_diag_f, estimates_contrast_cancer)]
comp_cancer_time_since_diag_f_adj[part == "Sleep period", .(cancer_time_since_diag_f, estimates_contrast_cancer)]

### Time since diagnosis x type --------------------------------------------------
fit_cancer_time_since_diag_type_adj <- brmcoda(clr_cancer_acc_sub,
                                          mvbind(z1_1, z2_1, z3_1) ~ cancer_time_since_diag_other * cancer_before_acc_type_other +
                                            s(age_at_acc) + sex + white + working + edu + never_smoked + current_drinker + s(deprivation) + season,
                                          # save_pars = save_pars(all = TRUE),
                                          warmup = 500, chains = 4, cores = 4, backend = "cmdstanr"
)
saveRDS(fit_cancer_time_since_diag_type_adj, paste0(outputdir, "fit_cancer_time_since_diag_type_adj", ".RDS"))

fit_cancer_time_since_diag_type_adj <- readRDS(paste0(outputdir, "fit_cancer_time_since_diag_type_adj", ".RDS"))

# reference grid
d_cancer_time_since_diag_type_adj <- emmeans::ref_grid(fit_cancer_time_since_diag_type_adj$model)@grid
d_cancer_time_since_diag_type_adj <- d_cancer_time_since_diag_type_adj[
  !(
    (d_cancer_time_since_diag_type_adj$cancer_time_since_diag_other == "Healthy" &
     d_cancer_time_since_diag_type_adj$cancer_before_acc_type_other != "Healthy") |
    (d_cancer_time_since_diag_type_adj$cancer_time_since_diag_other != "Healthy" &
     d_cancer_time_since_diag_type_adj$cancer_before_acc_type_other == "Healthy")
  ),
  , drop = FALSE
]

# predict
pred_cancer_time_since_diag_type_adj <- fitted(fit_cancer_time_since_diag_type_adj, newdata = d_cancer_time_since_diag_type_adj, scale = "response", summary = FALSE)

# weight by length equal to the number of observations in the dataset
# summarise by cancer group
pred_cancer_time_since_diag_type_adj <- apply(pred_cancer_time_since_diag_type_adj, c(1), function(x) cbind(d_cancer_time_since_diag_type_adj, x))
pred_cancer_time_since_diag_type_adj <- lapply(pred_cancer_time_since_diag_type_adj, function(d) {
 parts <- c("sleep", "mvpa", "lpa", "sb")
 parts0 <- c("tsleep_comp", "tmvpa_comp", "tlpa_comp", "tsb_comp")

 d <- as.data.table(d)

 # weighted mean by group
 d[, (paste0(parts)) :=
  lapply(.SD, function(x) {
   weighted.mean(x, .wgt.)
  }),
 .SDcols = parts0, by = .(cancer_time_since_diag_other, cancer_before_acc_type_other)
 ]
  d[, cancer_wgt := sum(.wgt.), by = .(cancer_time_since_diag_other, cancer_before_acc_type_other)]

  d[cancer_before_acc_type_other %nin% c("Healthy", "Others"), (paste0(parts, "_cancer")) := lapply(.SD, function(x) weighted.mean(x, cancer_wgt)), .SDcols = paste0(parts)]

 d <- rbind(d,
  data.table(cancer_before_acc_type_other = "Cancer"),
  fill = TRUE
 )

 d[cancer_before_acc_type_other == "Cancer", (parts) := d[cancer_before_acc_type_other != "Healthy",
  lapply(.SD, function(x) unique(na.omit(x))),
  .SDcols = paste0(parts, "_cancer")
 ]]

 d[, (paste0(parts, "_vs_healthy")) :=
  Map(`-`, .SD, d[cancer_before_acc_type_other == "Healthy",
   .SD,
   .SDcols = parts
  ][1]),
 .SDcols = parts
 ]

 d <- d[, .(
  cancer_before_acc_type_other, cancer_time_since_diag_other,
  sleep, mvpa, lpa, sb,

  sleep_vs_healthy, mvpa_vs_healthy, lpa_vs_healthy, sb_vs_healthy
 )]
 d <- unique(d)
 d
})
saveRDS(pred_cancer_time_since_diag_type_adj, paste0(outputdir, "pred_cancer_time_since_diag_type_adj", ".RDS"))
pred_cancer_time_since_diag_type_adj <- readRDS(paste0(outputdir, "pred_cancer_time_since_diag_type_adj", ".RDS"))

# assemble back to summarise posteriors
pred_cancer_time_since_diag_type_adj <- as.data.frame(abind::abind(pred_cancer_time_since_diag_type_adj, along = 1))
pred_cancer_time_since_diag_type_adj <- split(
  pred_cancer_time_since_diag_type_adj,
  interaction(
        pred_cancer_time_since_diag_type_adj$cancer_time_since_diag_other,
        pred_cancer_time_since_diag_type_adj$cancer_before_acc_type_other,
    drop = TRUE
  )
)

## estimated means  ----------------------
pred_comp_cancer_time_since_diag_type_adj <- lapply(pred_cancer_time_since_diag_type_adj, function(l) {
 l <- as.data.frame(l[, c("sleep", "mvpa", "lpa", "sb")])
 l <- apply(l, 2, as.numeric)
 l <- apply(l, 2, function(p) describe_posterior(p, centrality = "mean", ci = 0.99))
 l <- Map(cbind, l, part = names(l))
 l <- rbindlist(l)
 l
})
pred_comp_cancer_time_since_diag_type_adj <- Map(cbind, pred_comp_cancer_time_since_diag_type_adj, cancer_time_since_diag_type_other = names(pred_comp_cancer_time_since_diag_type_adj))
pred_comp_cancer_time_since_diag_type_adj <- rbindlist(pred_comp_cancer_time_since_diag_type_adj)

## contrasts --------------------
### vs healthy
diff1_comp_cancer_time_since_diag_type_adj <- lapply(pred_cancer_time_since_diag_type_adj, function(l) {
 l <- as.data.frame(l[, c("sleep_vs_healthy", "mvpa_vs_healthy", "lpa_vs_healthy", "sb_vs_healthy")])
 l <- apply(l, 2, as.numeric)
 l <- apply(l, 2, function(p) describe_posterior(p, centrality = "mean", ci = 0.99))
 l <- Map(cbind, l, contrast = names(l))
 l <- rbindlist(l)
 l
})
diff1_comp_cancer_time_since_diag_type_adj <- Map(cbind, diff1_comp_cancer_time_since_diag_type_adj, cancer_time_since_diag_type_other = names(diff1_comp_cancer_time_since_diag_type_adj))
diff1_comp_cancer_time_since_diag_type_adj <- rbindlist(diff1_comp_cancer_time_since_diag_type_adj)

setnames(diff1_comp_cancer_time_since_diag_type_adj, "Mean", "Mean_diff_ref_healthy")
setnames(diff1_comp_cancer_time_since_diag_type_adj, "CI_low", "CI_low_diff_ref_healthy")
setnames(diff1_comp_cancer_time_since_diag_type_adj, "CI_high", "CI_high_diff_ref_healthy")

# all results  ------------------------
comp_cancer_time_since_diag_type_adj <- cbind(
 pred_comp_cancer_time_since_diag_type_adj[, .(Mean, CI_low, CI_high, part, cancer_time_since_diag_type_other)],
 diff1_comp_cancer_time_since_diag_type_adj[, .(Mean_diff_ref_healthy, CI_low_diff_ref_healthy, CI_high_diff_ref_healthy)]
 # diff3_comp_cancer_time_since_diag_f_adj[, .(Mean_diff_ref_others, CI_low_diff_ref_others, CI_high_diff_ref_others)]
)
comp_cancer_time_since_diag_type_adj[, id := 1:.N, by = cancer_time_since_diag_type_other]

# add sig indicators
comp_cancer_time_since_diag_type_adj[, nonsig_vs_healthy := between(0, comp_cancer_time_since_diag_type_adj$CI_low_diff_ref_healthy, comp_cancer_time_since_diag_type_adj$CI_high_diff_ref_healthy)]
# comp_cancer_time_since_diag_type_adj[, nonsig_vs_others  := between(0, comp_cancer_time_since_diag_type_adj$CI_low_diff_ref_others, comp_cancer_time_since_diag_type_adj$CI_high_diff_ref_others)]
comp_cancer_time_since_diag_type_adj[, nonsig_vs_cancer := between(0, comp_cancer_time_since_diag_type_adj$CI_low_diff_ref_cancer, comp_cancer_time_since_diag_type_adj$CI_high_diff_ref_cancer)]

comp_cancer_time_since_diag_type_adj[, sig_ref_healthy := ifelse(nonsig_vs_healthy == FALSE & Mean_diff_ref_healthy != 0,
 "$^a$", "$\\phantom{^a}$"
)]

# leave healthy empty
comp_cancer_time_since_diag_type_adj[, sig_ref_healthy := ifelse(cancer_time_since_diag_type_other == "Healthy.Healthy", "$\\phantom{^a}$", sig_ref_healthy)]
# comp_cancer_time_since_diag_f_adj[, sig_ref_others := ifelse(cancer_time_since_diag_f == "Healthy", "$\\phantom{^b}$", sig_ref_others)]

# healthy in first row of plot
table(comp_cancer_time_since_diag_type_adj$cancer_time_since_diag_type_other)

comp_cancer_time_since_diag_type_adj[, cancer_time_since_diag_type_other := ifelse(cancer_time_since_diag_type_other == "Healthy.Healthy", "Healthy", cancer_time_since_diag_type_other)]

comp_cancer_time_since_diag_type_adj[, cancer_time_since_diag_type_other := ifelse(cancer_time_since_diag_type_other == "Less than 1 year since diagnosis.Blood",  "Blood - <1y", cancer_time_since_diag_type_other)]
comp_cancer_time_since_diag_type_adj[, cancer_time_since_diag_type_other := ifelse(cancer_time_since_diag_type_other == "Less than 1 year since diagnosis.Breast",  "Breast - <1y", cancer_time_since_diag_type_other)]
comp_cancer_time_since_diag_type_adj[, cancer_time_since_diag_type_other := ifelse(cancer_time_since_diag_type_other == "Less than 1 year since diagnosis.Other Skin",  "Skin (non-melanoma) - <1y", cancer_time_since_diag_type_other)]
comp_cancer_time_since_diag_type_adj[, cancer_time_since_diag_type_other := ifelse(cancer_time_since_diag_type_other == "Less than 1 year since diagnosis.Prostate",  "Prostate - <1y", cancer_time_since_diag_type_other)]
comp_cancer_time_since_diag_type_adj[, cancer_time_since_diag_type_other := ifelse(cancer_time_since_diag_type_other == "Less than 1 year since diagnosis.Gynaecological",  "Gynaecological - <1y", cancer_time_since_diag_type_other)]
comp_cancer_time_since_diag_type_adj[, cancer_time_since_diag_type_other := ifelse(cancer_time_since_diag_type_other == "Less than 1 year since diagnosis.Colorectal",  "Colorectal - <1y", cancer_time_since_diag_type_other)]
comp_cancer_time_since_diag_type_adj[, cancer_time_since_diag_type_other := ifelse(cancer_time_since_diag_type_other == "Less than 1 year since diagnosis.Cancer",  "Cancer - <1y", cancer_time_since_diag_type_other)]

comp_cancer_time_since_diag_type_adj[, cancer_time_since_diag_type_other := ifelse(cancer_time_since_diag_type_other == "1-5 years since diagnosis.Blood",  "Blood - 1-5y", cancer_time_since_diag_type_other)]
comp_cancer_time_since_diag_type_adj[, cancer_time_since_diag_type_other := ifelse(cancer_time_since_diag_type_other == "1-5 years since diagnosis.Breast",  "Breast - 1-5y", cancer_time_since_diag_type_other)]
comp_cancer_time_since_diag_type_adj[, cancer_time_since_diag_type_other := ifelse(cancer_time_since_diag_type_other == "1-5 years since diagnosis.Other Skin",  "Skin (non-melanoma) - 1-5y", cancer_time_since_diag_type_other)]
comp_cancer_time_since_diag_type_adj[, cancer_time_since_diag_type_other := ifelse(cancer_time_since_diag_type_other == "1-5 years since diagnosis.Prostate",  "Prostate - 1-5y", cancer_time_since_diag_type_other)]
comp_cancer_time_since_diag_type_adj[, cancer_time_since_diag_type_other := ifelse(cancer_time_since_diag_type_other == "1-5 years since diagnosis.Gynaecological",  "Gynaecological - 1-5y", cancer_time_since_diag_type_other)]
comp_cancer_time_since_diag_type_adj[, cancer_time_since_diag_type_other := ifelse(cancer_time_since_diag_type_other == "1-5 years since diagnosis.Colorectal",  "Colorectal - 1-5y", cancer_time_since_diag_type_other)]
comp_cancer_time_since_diag_type_adj[, cancer_time_since_diag_type_other := ifelse(cancer_time_since_diag_type_other == "1-5 years since diagnosis.Cancer",  "Cancer - 1-5y", cancer_time_since_diag_type_other)]

comp_cancer_time_since_diag_type_adj[, cancer_time_since_diag_type_other := ifelse(cancer_time_since_diag_type_other == "More than 5 years since diagnosis.Blood",  "Blood - >5y", cancer_time_since_diag_type_other)]
comp_cancer_time_since_diag_type_adj[, cancer_time_since_diag_type_other := ifelse(cancer_time_since_diag_type_other == "More than 5 years since diagnosis.Breast",  "Breast - >5y", cancer_time_since_diag_type_other)]
comp_cancer_time_since_diag_type_adj[, cancer_time_since_diag_type_other := ifelse(cancer_time_since_diag_type_other == "More than 5 years since diagnosis.Other Skin",  "Skin (non-melanoma) - >5y", cancer_time_since_diag_type_other)]
comp_cancer_time_since_diag_type_adj[, cancer_time_since_diag_type_other := ifelse(cancer_time_since_diag_type_other == "More than 5 years since diagnosis.Prostate",  "Prostate - >5y", cancer_time_since_diag_type_other)]
comp_cancer_time_since_diag_type_adj[, cancer_time_since_diag_type_other := ifelse(cancer_time_since_diag_type_other == "More than 5 years since diagnosis.Gynaecological",  "Gynaecological - >5y", cancer_time_since_diag_type_other)]
comp_cancer_time_since_diag_type_adj[, cancer_time_since_diag_type_other := ifelse(cancer_time_since_diag_type_other == "More than 5 years since diagnosis.Colorectal",  "Colorectal - >5y", cancer_time_since_diag_type_other)]
comp_cancer_time_since_diag_type_adj[, cancer_time_since_diag_type_other := ifelse(cancer_time_since_diag_type_other == "More than 5 years since diagnosis.Cancer",  "Cancer - >5y", cancer_time_since_diag_type_other)]


comp_cancer_time_since_diag_type_adj[, cancer_time_since_diag_type_other := factor(cancer_time_since_diag_type_other,
 ordered = TRUE,
 levels = c(
  "Gynaecological - >5y",
  "Gynaecological - 1-5y",
  "Gynaecological - <1y",

  "Breast - >5y",
  "Breast - 1-5y",
  "Breast - <1y",

  "Blood - >5y",
  "Blood - 1-5y",
  "Blood - <1y",

  "Colorectal - >5y",
  "Colorectal - 1-5y",
  "Colorectal - <1y",

  "Skin (non-melanoma) - >5y",
  "Skin (non-melanoma) - 1-5y",
  "Skin (non-melanoma) - <1y",

  "Prostate - >5y",
  "Prostate - 1-5y",
  "Prostate - <1y",

  "Cancer - >5y",
  "Cancer - 1-5y",
  "Cancer - <1y",
  "Healthy"
 )
)]

comp_cancer_time_since_diag_type_adj[, yintercept_healthy := NA]
comp_cancer_time_since_diag_type_adj[, yintercept_healthy := ifelse(part == "sleep", comp_cancer_time_since_diag_type_adj[cancer_time_since_diag_type_other == "Healthy" & part == "sleep"]$Mean, yintercept_healthy)]
comp_cancer_time_since_diag_type_adj[, yintercept_healthy := ifelse(part == "mvpa", comp_cancer_time_since_diag_type_adj[cancer_time_since_diag_type_other == "Healthy" & part == "mvpa"]$Mean, yintercept_healthy)]
comp_cancer_time_since_diag_type_adj[, yintercept_healthy := ifelse(part == "lpa", comp_cancer_time_since_diag_type_adj[cancer_time_since_diag_type_other == "Healthy" & part == "lpa"]$Mean, yintercept_healthy)]
comp_cancer_time_since_diag_type_adj[, yintercept_healthy := ifelse(part == "sb", comp_cancer_time_since_diag_type_adj[cancer_time_since_diag_type_other == "Healthy" & part == "sb"]$Mean, yintercept_healthy)]

comp_cancer_time_since_diag_type_adj[, ci_low_healthy := NA]
comp_cancer_time_since_diag_type_adj[, ci_low_healthy := ifelse(part == "sleep", comp_cancer_time_since_diag_type_adj[cancer_time_since_diag_type_other == "Healthy" & part == "sleep"]$CI_low, ci_low_healthy)]
comp_cancer_time_since_diag_type_adj[, ci_low_healthy := ifelse(part == "mvpa", comp_cancer_time_since_diag_type_adj[cancer_time_since_diag_type_other == "Healthy" & part == "mvpa"]$CI_low, ci_low_healthy)]
comp_cancer_time_since_diag_type_adj[, ci_low_healthy := ifelse(part == "lpa", comp_cancer_time_since_diag_type_adj[cancer_time_since_diag_type_other == "Healthy" & part == "lpa"]$CI_low, ci_low_healthy)]
comp_cancer_time_since_diag_type_adj[, ci_low_healthy := ifelse(part == "sb", comp_cancer_time_since_diag_type_adj[cancer_time_since_diag_type_other == "Healthy" & part == "sb"]$CI_low, ci_low_healthy)]

comp_cancer_time_since_diag_type_adj[, ci_high_healthy := NA]
comp_cancer_time_since_diag_type_adj[, ci_high_healthy := ifelse(part == "sleep", comp_cancer_time_since_diag_type_adj[cancer_time_since_diag_type_other == "Healthy" & part == "sleep"]$CI_high, ci_high_healthy)]
comp_cancer_time_since_diag_type_adj[, ci_high_healthy := ifelse(part == "mvpa", comp_cancer_time_since_diag_type_adj[cancer_time_since_diag_type_other == "Healthy" & part == "mvpa"]$CI_high, ci_high_healthy)]
comp_cancer_time_since_diag_type_adj[, ci_high_healthy := ifelse(part == "lpa", comp_cancer_time_since_diag_type_adj[cancer_time_since_diag_type_other == "Healthy" & part == "lpa"]$CI_high, ci_high_healthy)]
comp_cancer_time_since_diag_type_adj[, ci_high_healthy := ifelse(part == "sb", comp_cancer_time_since_diag_type_adj[cancer_time_since_diag_type_other == "Healthy" & part == "sb"]$CI_high, ci_high_healthy)]

# n
table(comp_cancer_time_since_diag_type_adj$cancer_time_since_diag_type_other)
as.data.table(model.frame(fit_cancer_time_since_diag_type_adj))[ , .N, by = .(cancer_time_since_diag_other, cancer_before_acc_type_other)][
  order(cancer_time_since_diag_other, cancer_before_acc_type_other)]

comp_cancer_time_since_diag_type_adj[, Cases := NA]
comp_cancer_time_since_diag_type_adj[, Cases := ifelse(cancer_time_since_diag_type_other == "Healthy", "13 722", Cases)]

comp_cancer_time_since_diag_type_adj[, Cases := ifelse(cancer_time_since_diag_type_other == "Cancer - <1y", "635", Cases)]
comp_cancer_time_since_diag_type_adj[, Cases := ifelse(cancer_time_since_diag_type_other == "Cancer - 1-5y", "2482", Cases)]
comp_cancer_time_since_diag_type_adj[, Cases := ifelse(cancer_time_since_diag_type_other == "Cancer - >5y", "4131", Cases)]

comp_cancer_time_since_diag_type_adj[, Cases := ifelse(cancer_time_since_diag_type_other == "Prostate - <1y", "116", Cases)]
comp_cancer_time_since_diag_type_adj[, Cases := ifelse(cancer_time_since_diag_type_other == "Prostate - 1-5y", "514", Cases)]
comp_cancer_time_since_diag_type_adj[, Cases := ifelse(cancer_time_since_diag_type_other == "Prostate - >5y", "486", Cases)]

comp_cancer_time_since_diag_type_adj[, Cases := ifelse(cancer_time_since_diag_type_other == "Skin (non-melanoma) - <1y", "319", Cases)]
comp_cancer_time_since_diag_type_adj[, Cases := ifelse(cancer_time_since_diag_type_other == "Skin (non-melanoma) - 1-5y", "1097", Cases)]
comp_cancer_time_since_diag_type_adj[, Cases := ifelse(cancer_time_since_diag_type_other == "Skin (non-melanoma) - >5y", "1601", Cases)]

comp_cancer_time_since_diag_type_adj[, Cases := ifelse(cancer_time_since_diag_type_other == "Breast - <1y", "85", Cases)]
comp_cancer_time_since_diag_type_adj[, Cases := ifelse(cancer_time_since_diag_type_other == "Breast - 1-5y", "491", Cases)]
comp_cancer_time_since_diag_type_adj[, Cases := ifelse(cancer_time_since_diag_type_other == "Breast - >5y", "1335", Cases)]

comp_cancer_time_since_diag_type_adj[, Cases := ifelse(cancer_time_since_diag_type_other == "Blood - <1y", "39", Cases)]
comp_cancer_time_since_diag_type_adj[, Cases := ifelse(cancer_time_since_diag_type_other == "Blood - 1-5y", "126", Cases)]
comp_cancer_time_since_diag_type_adj[, Cases := ifelse(cancer_time_since_diag_type_other == "Blood - >5y", "247", Cases)]

comp_cancer_time_since_diag_type_adj[, Cases := ifelse(cancer_time_since_diag_type_other == "Colorectal - <1y", "44", Cases)]
comp_cancer_time_since_diag_type_adj[, Cases := ifelse(cancer_time_since_diag_type_other == "Colorectal - 1-5y", "152", Cases)]
comp_cancer_time_since_diag_type_adj[, Cases := ifelse(cancer_time_since_diag_type_other == "Colorectal - >5y", "195", Cases)]

comp_cancer_time_since_diag_type_adj[, Cases := ifelse(cancer_time_since_diag_type_other == "Gynaecological - <1y", "32", Cases)]
comp_cancer_time_since_diag_type_adj[, Cases := ifelse(cancer_time_since_diag_type_other == "Gynaecological - 1-5y", "102", Cases)]
comp_cancer_time_since_diag_type_adj[, Cases := ifelse(cancer_time_since_diag_type_other == "Gynaecological - >5y", "267", Cases)]

comp_cancer_time_since_diag_type_adj[, part := ifelse(part == "sleep", "Sleep period", part)]
comp_cancer_time_since_diag_type_adj[, part := ifelse(part == "mvpa", "Moderate-to-vigorous physical activity", part)]
comp_cancer_time_since_diag_type_adj[, part := ifelse(part == "lpa", "Light physical activity", part)]
comp_cancer_time_since_diag_type_adj[, part := ifelse(part == "sb", "Sedentary behaviour", part)]

comp_cancer_time_since_diag_type_adj[, sig_position := min(CI_low), by = part]
comp_cancer_time_since_diag_type_adj[, est_position := max(CI_high), by = part]

comp_cancer_time_since_diag_type_adj[, estimates := paste0(round(Mean, 0), "[", round(CI_low, 0), ", ", round(CI_high, 0), "]")]
comp_cancer_time_since_diag_type_adj[, est_sig := paste0(estimates, " ", sig_ref_healthy)]
comp_cancer_time_since_diag_type_adj[, est_sig := paste0(estimates, " ", str_replace_na(sig_ref_healthy, " "))]

# for tables
comp_cancer_time_since_diag_type_adj[, estimates_contrast_healthy := paste0(round(Mean_diff_ref_healthy, 2), "[", round(CI_low_diff_ref_healthy, 2), ", ", round(CI_high_diff_ref_healthy, 2), "]")]
# comp_cancer_time_since_diag_type_adj[, estimates_contrast_others := paste0(round(Mean_diff_ref_others, 2), "[", round(CI_low_diff_ref_others, 2), ", ", round(CI_high_diff_ref_others, 2), "]")]
# comp_cancer_time_since_diag_type_adj[, estimates_contrast_cancer := paste0(round(Mean_diff_ref_cancer, 2), "[", round(CI_low_diff_ref_cancer, 2), ", ", round(CI_high_diff_ref_cancer, 2), "]")]

## plot -----------------------
plot_specs <- data.table(
  part = c(
    "Sleep period",
    "Moderate-to-vigorous physical activity",
    "Light physical activity",
    "Sedentary behaviour"
  ),
  y_limits = list(c(450, 650), c(0, 60), c(200, 400), c(500, 700)),
  text_y = c(450, 0, 200, 500),
  sig_y = c(615, 49, 365.5, 665),
  cases_y = c(635, 55.5, 385, 685),
  surv_y = c(650, 60, 400, 700),
  seg_y = list(c(450, 650), c(0, 60), c(200, 400), c(500, 700)),
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

 ggplot(comp_cancer_time_since_diag_type_adj[part == spec$part], aes(x = cancer_time_since_diag_type_other, y = Mean, group = cancer_time_since_diag_type_other)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy),
   fill = "#CBD5D0", alpha = 0.1
  ) +
  geom_hline(aes(yintercept = yintercept_healthy),
   linewidth = 0.5,
   linetype = "dashed", colour = "#708885"
  ) +
  geom_pointrange(aes(ymin = CI_low, ymax = CI_high, colour = cancer_time_since_diag_type_other),
   size = 0.2, linewidth = 0.5
  ) +
  geom_text(aes(y = spec$text_y, label = cancer_time_since_diag_type_other),
   hjust = 0, nudge_x = 0, family = "Arial Narrow", size = 2.5,
   show.legend = FALSE
  ) +
  geom_text(aes(y = spec$sig_y, label = TeX(est_sig, output = "character")),
   parse = TRUE,
   hjust = 0.5, nudge_x = 0, family = "Arial Narrow", size = 2.5,
   show.legend = FALSE
  ) +
  geom_text(aes(y = spec$cases_y, label = Cases),
   hjust = 1, nudge_x = 0,
   family = "Arial Narrow", size = 2.5,
   show.legend = FALSE
  ) +
  geom_segment(aes(x = 0, yend = seg_vals[1]), col = "black", linewidth = 0.5) +
  geom_segment(aes(x = 0, yend = seg_vals[2]), col = "black", linewidth = 0.5) +
  scale_y_continuous(limits = limits, breaks = limits, name = spec$ylab) +
  scale_colour_manual(values = pal_type_quantile) +
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

plot_comp_cancer_time_since_diag_type_sleep <- plot_list[["sleep_period"]]
plot_comp_cancer_time_since_diag_type_mvpa <- plot_list[["moderate-to-vigorous_physical_activity"]]
plot_comp_cancer_time_since_diag_type_lpa <- plot_list[["light_physical_activity"]]
plot_comp_cancer_time_since_diag_type_sb <- plot_list[["sedentary_behaviour"]]

grDevices::cairo_pdf(
 file = paste0(outputdir, "cancer_time_since_diag_type_est", ".pdf"),
 width = 6,
 height = 8,
)

ggarrange(
 plot_comp_cancer_time_since_diag_type_mvpa,
 plot_comp_cancer_time_since_diag_type_lpa,
 plot_comp_cancer_time_since_diag_type_sb,
 plot_comp_cancer_time_since_diag_type_sleep,
 nrow = 4
)
dev.off()

grDevices::png(
 file = paste0(outputdir, "cancer_time_since_diag_type_est", ".png"),
 width = 6000,
 height = 8000,
 res = 900
)

ggarrange(
 plot_comp_cancer_time_since_diag_type_mvpa,
 plot_comp_cancer_time_since_diag_type_lpa,
 plot_comp_cancer_time_since_diag_type_sb,
 plot_comp_cancer_time_since_diag_type_sleep,
 nrow = 4
)
dev.off()

### Excluding multiple cancers --------------------------------------------------

fit_cancer_time_since_diag_one_primary <- brmcoda(clr_cancer_acc_one_primary,
                                          mvbind(z1_1, z2_1, z3_1) ~ cancer_time_since_diag_other +
                                            # + other_conds_at_acc +
                                            s(age_at_acc) + sex + white + working + edu + never_smoked + current_drinker + s(deprivation) + season,
                                          # save_pars = save_pars(all = TRUE),
                                          warmup = 500, chains = 4, cores = 4, backend = "cmdstanr",
                                          seed = 2025
)
saveRDS(fit_cancer_time_since_diag_one_primary, paste0(outputdir, "fit_cancer_time_since_diag_one_primary", ".RDS"))
fit_cancer_time_since_diag_one_primary <- readRDS(paste0(outputdir, "fit_cancer_time_since_diag_one_primary", ".RDS"))

table(model.frame(fit_cancer_time_since_diag_one_primary)$cancer_time_since_diag_other)

# predicted posteriors ------------

# reference grid
d_cancer_time_since_diag_one_primary_adj <- emmeans::ref_grid(fit_cancer_time_since_diag_one_primary$model)@grid

# predict
pred_cancer_time_since_diag_one_primary_adj <- fitted(fit_cancer_time_since_diag_one_primary, newdata = d_cancer_time_since_diag_one_primary_adj, scale = "response", summary = FALSE)

# weight by length equal to the number of observations in the dataset
# summarise by cancer group
pred_cancer_time_since_diag_one_primary_adj <- apply(pred_cancer_time_since_diag_one_primary_adj, c(1), function(x)  cbind(d_cancer_time_since_diag_one_primary_adj, x))
pred_cancer_time_since_diag_one_primary_adj <- lapply(pred_cancer_time_since_diag_one_primary_adj, function(d) {
  
  parts <- c("sleep", "mvpa", "lpa", "sb")
  parts0 <- c("tsleep_comp", "tmvpa_comp", "tlpa_comp", "tsb_comp")
  
  d <- as.data.table(d)
  
  d[, (paste0(parts)) :=
        lapply(.SD, function(x)
          weighted.mean(x, .wgt.)),
        .SDcols = parts0, by = cancer_time_since_diag_other]

  d[, cancer_wgt := sum(.wgt.), by = cancer_time_since_diag_other]

  d[cancer_time_since_diag_other %nin% c("Healthy", "Others"), (paste0(parts, "_cancer")) := lapply(.SD, function(x) weighted.mean(x, cancer_wgt)), .SDcols = paste0(parts)]
  
  d <- rbind(d, 
             data.table(cancer_time_since_diag_other = "Cancer"),
             fill = TRUE
  )

  d[cancer_time_since_diag_other == "Cancer", (parts) := d[cancer_time_since_diag_other %nin% c("Healthy", "Others"),
   lapply(.SD, function(x) unique(na.omit(x))),
   .SDcols = paste0(parts, "_cancer")
  ]]

  d[, (paste0(parts, "_vs_healthy")) :=
   Map(`-`, .SD, d[cancer_time_since_diag_other == "Healthy",
    .SD,
    .SDcols = parts
   ][1]),
  .SDcols = parts
  ]

  # pairwise contrast within cancer 
  d[, sleep_lessthan1_vs_1to5 := d[cancer_time_since_diag_other == "1-5 years since diagnosis"]$sleep[1] - d[cancer_time_since_diag_other == "Less than 1 year since diagnosis"]$sleep[1]]
  d[, sleep_lessthan1_vs_morethan5 := d[cancer_time_since_diag_other == "More than 5 years since diagnosis"]$sleep[1] - d[cancer_time_since_diag_other == "Less than 1 year since diagnosis"]$sleep[1]]
  d[, sleep_1to5_vs_morethan5 := d[cancer_time_since_diag_other == "More than 5 years since diagnosis"]$sleep[1] - d[cancer_time_since_diag_other == "1-5 years since diagnosis"]$sleep[1]]
  
  d[, mvpa_lessthan1_vs_1to5 := d[cancer_time_since_diag_other == "1-5 years since diagnosis"]$mvpa[1] - d[cancer_time_since_diag_other == "Less than 1 year since diagnosis"]$mvpa[1]]
  d[, mvpa_lessthan1_vs_morethan5 := d[cancer_time_since_diag_other == "More than 5 years since diagnosis"]$mvpa[1] - d[cancer_time_since_diag_other == "Less than 1 year since diagnosis"]$mvpa[1]]
  d[, mvpa_1to5_vs_morethan5 := d[cancer_time_since_diag_other == "More than 5 years since diagnosis"]$mvpa[1] - d[cancer_time_since_diag_other == "1-5 years since diagnosis"]$mvpa[1]]
  
  d[, lpa_lessthan1_vs_1to5 := d[cancer_time_since_diag_other == "1-5 years since diagnosis"]$lpa[1] - d[cancer_time_since_diag_other == "Less than 1 year since diagnosis"]$lpa[1]]
  d[, lpa_lessthan1_vs_morethan5 := d[cancer_time_since_diag_other == "More than 5 years since diagnosis"]$lpa[1] - d[cancer_time_since_diag_other == "Less than 1 year since diagnosis"]$lpa[1]]
  d[, lpa_1to5_vs_morethan5 := d[cancer_time_since_diag_other == "More than 5 years since diagnosis"]$lpa[1] - d[cancer_time_since_diag_other == "1-5 years since diagnosis"]$lpa[1]]
  
  d[, sb_lessthan1_vs_1to5 := d[cancer_time_since_diag_other == "1-5 years since diagnosis"]$sb[1] - d[cancer_time_since_diag_other == "Less than 1 year since diagnosis"]$sb[1]]
  d[, sb_lessthan1_vs_morethan5 := d[cancer_time_since_diag_other == "More than 5 years since diagnosis"]$sb[1] - d[cancer_time_since_diag_other == "Less than 1 year since diagnosis"]$sb[1]]
  d[, sb_1to5_vs_morethan5 := d[cancer_time_since_diag_other == "More than 5 years since diagnosis"]$sb[1] - d[cancer_time_since_diag_other == "1-5 years since diagnosis"]$sb[1]]
  
  d <- d[, .(cancer_time_since_diag_other, 
             sleep, mvpa, lpa, sb,
             
             # sleep_cancer, mvpa_cancer, lpa_cancer, sb_cancer,
             
             sleep_vs_healthy, mvpa_vs_healthy, lpa_vs_healthy, sb_vs_healthy,
            #  sleep_vs_others, mvpa_vs_others, lpa_vs_others, sb_vs_others,
             
             sleep_lessthan1_vs_1to5, sleep_lessthan1_vs_morethan5, sleep_1to5_vs_morethan5,
             mvpa_lessthan1_vs_1to5, mvpa_lessthan1_vs_morethan5, mvpa_1to5_vs_morethan5,
             lpa_lessthan1_vs_1to5, lpa_lessthan1_vs_morethan5, lpa_1to5_vs_morethan5,
             sb_lessthan1_vs_1to5, sb_lessthan1_vs_morethan5, sb_1to5_vs_morethan5
             
  )]
  d <- unique(d)
  d
})

saveRDS(pred_cancer_time_since_diag_one_primary_adj, paste0(outputdir, "pred_cancer_time_since_diag_one_primary_adj", ".RDS"))
pred_cancer_time_since_diag_one_primary_adj <- readRDS(paste0(outputdir, "pred_cancer_time_since_diag_one_primary_adj", ".RDS"))

# assemble back to summarise posteriors
pred_cancer_time_since_diag_one_primary_adj <- as.data.frame(abind::abind(pred_cancer_time_since_diag_one_primary_adj, along = 1))
pred_cancer_time_since_diag_one_primary_adj <- split(pred_cancer_time_since_diag_one_primary_adj, pred_cancer_time_since_diag_one_primary_adj$cancer_time_since_diag_other)

## estimated means  ----------------------
pred_comp_cancer_time_since_diag_one_primary_adj <- lapply(pred_cancer_time_since_diag_one_primary_adj, function(l) {
  l <- as.data.frame(l[, c("sleep", "mvpa", "lpa", "sb")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, function(p) describe_posterior(p, centrality = "mean", ci = 0.99))
  l <- Map(cbind, l, part = names(l))
  l <- rbindlist(l)
  l
})
pred_comp_cancer_time_since_diag_one_primary_adj <- Map(cbind, pred_comp_cancer_time_since_diag_one_primary_adj, cancer_time_since_diag_other = names(pred_comp_cancer_time_since_diag_one_primary_adj))
pred_comp_cancer_time_since_diag_one_primary_adj <- rbindlist(pred_comp_cancer_time_since_diag_one_primary_adj)

## contrasts --------------------
### vs healthy
diff1_comp_cancer_time_since_diag_one_primary_adj <- lapply(pred_cancer_time_since_diag_one_primary_adj, function(l) {
  l <- as.data.frame(l[, c("sleep_vs_healthy", "mvpa_vs_healthy", "lpa_vs_healthy", "sb_vs_healthy")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, function(p) describe_posterior(p, centrality = "mean", ci = 0.99))
  l <- Map(cbind, l, contrast = names(l))
  l <- rbindlist(l)
  l
})
diff1_comp_cancer_time_since_diag_one_primary_adj <- Map(cbind, diff1_comp_cancer_time_since_diag_one_primary_adj, cancer_time_since_diag_other = names(diff1_comp_cancer_time_since_diag_one_primary_adj))
diff1_comp_cancer_time_since_diag_one_primary_adj <- rbindlist(diff1_comp_cancer_time_since_diag_one_primary_adj)

setnames(diff1_comp_cancer_time_since_diag_one_primary_adj, "Mean", "Mean_diff_ref_healthy")
setnames(diff1_comp_cancer_time_since_diag_one_primary_adj, "CI_low", "CI_low_diff_ref_healthy")
setnames(diff1_comp_cancer_time_since_diag_one_primary_adj, "CI_high", "CI_high_diff_ref_healthy")

# ### vs others
# diff3_comp_cancer_time_since_diag_one_primary_adj <- lapply(pred_cancer_time_since_diag_one_primary_adj, function(l) {
#   l <- as.data.frame(l[, c("sleep_vs_others", "mvpa_vs_others", "lpa_vs_others", "sb_vs_others")])
#   l <- apply(l, 2, as.numeric)
#   l <- apply(l, 2, function(p) describe_posterior(p, centrality = "mean", ci = 0.99))
#   l <- Map(cbind, l, contrast = names(l))
#   l <- rbindlist(l)
#   l
# })
# diff3_comp_cancer_time_since_diag_one_primary_adj <- Map(cbind, diff3_comp_cancer_time_since_diag_one_primary_adj, cancer_time_since_diag_other = names(diff3_comp_cancer_time_since_diag_one_primary_adj))
# diff3_comp_cancer_time_since_diag_one_primary_adj <- rbindlist(diff3_comp_cancer_time_since_diag_one_primary_adj)

# setnames(diff3_comp_cancer_time_since_diag_one_primary_adj, "Mean", "Mean_diff_ref_others")
# setnames(diff3_comp_cancer_time_since_diag_one_primary_adj, "CI_low", "CI_low_diff_ref_others")
# setnames(diff3_comp_cancer_time_since_diag_one_primary_adj, "CI_high", "CI_high_diff_ref_others")

### pairwise cancer
diff2_comp_cancer_time_since_diag_one_primary_adj <- lapply(pred_cancer_time_since_diag_one_primary_adj, function(l) {
  l <- as.data.frame(l[, c(
    "sleep_lessthan1_vs_1to5", "sleep_lessthan1_vs_morethan5", "sleep_1to5_vs_morethan5",
    "mvpa_lessthan1_vs_1to5", "mvpa_lessthan1_vs_morethan5", "mvpa_1to5_vs_morethan5",
    "lpa_lessthan1_vs_1to5", "lpa_lessthan1_vs_morethan5", "lpa_1to5_vs_morethan5",
    "sb_lessthan1_vs_1to5", "sb_lessthan1_vs_morethan5", "sb_1to5_vs_morethan5"
  )])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, function(p) describe_posterior(p, centrality = "mean", ci = 0.99))
  l <- Map(cbind, l, contrast = names(l))
  l <- rbindlist(l)
  l
})
# diff2_comp_cancer_time_since_diag_one_primary_adj <- Map(cbind, diff2_comp_cancer_time_since_diag_one_primary_adj, cancer_time_since_diag_other = names(diff2_comp_cancer_time_since_diag_one_primary_adj))
diff2_comp_cancer_time_since_diag_one_primary_adj <- rbindlist(diff2_comp_cancer_time_since_diag_one_primary_adj)

diff2_comp_cancer_time_since_diag_one_primary_adj[, cancer_time_since_diag_other := NA]
diff2_comp_cancer_time_since_diag_one_primary_adj[, cancer_time_since_diag_other := ifelse(grepl("lessthan1_vs_1to5", diff2_comp_cancer_time_since_diag_one_primary_adj$contrast), "Less than 1 year since diagnosis", cancer_time_since_diag_other)]
diff2_comp_cancer_time_since_diag_one_primary_adj[, cancer_time_since_diag_other := ifelse(grepl("1to5_vs_morethan5", diff2_comp_cancer_time_since_diag_one_primary_adj$contrast), "1-5 years since diagnosis", cancer_time_since_diag_other)]
diff2_comp_cancer_time_since_diag_one_primary_adj[, cancer_time_since_diag_other := ifelse(grepl("lessthan1_vs_morethan5", diff2_comp_cancer_time_since_diag_one_primary_adj$contrast), "More than 5 years since diagnosis", cancer_time_since_diag_other)]

setnames(diff2_comp_cancer_time_since_diag_one_primary_adj, "Mean", "Mean_diff_ref_cancer")
setnames(diff2_comp_cancer_time_since_diag_one_primary_adj, "CI_low", "CI_low_diff_ref_cancer")
setnames(diff2_comp_cancer_time_since_diag_one_primary_adj, "CI_high", "CI_high_diff_ref_cancer")

# take only unique values
diff2_comp_cancer_time_since_diag_one_primary_adj <- unique(diff2_comp_cancer_time_since_diag_one_primary_adj)

# add id to merge later
diff2_comp_cancer_time_since_diag_one_primary_adj[, id := 1:.N, by = cancer_time_since_diag_other]

# all results  ------------------------
comp_cancer_time_since_diag_one_primary_adj <- cbind(
  pred_comp_cancer_time_since_diag_one_primary_adj[, .(Mean, CI_low, CI_high, part, cancer_time_since_diag_other)],
  diff1_comp_cancer_time_since_diag_one_primary_adj[, .(Mean_diff_ref_healthy, CI_low_diff_ref_healthy, CI_high_diff_ref_healthy)]
    # diff3_comp_cancer_time_since_diag_one_primary_adj[, .(Mean_diff_ref_others, CI_low_diff_ref_others, CI_high_diff_ref_others)]
)
comp_cancer_time_since_diag_one_primary_adj[, id := 1:.N, by = cancer_time_since_diag_other]

comp_cancer_time_since_diag_one_primary_adj <- merge(
  comp_cancer_time_since_diag_one_primary_adj,
  diff2_comp_cancer_time_since_diag_one_primary_adj[, .(
    Mean_diff_ref_cancer,
    CI_low_diff_ref_cancer,
    CI_high_diff_ref_cancer,
    cancer_time_since_diag_other,
    id
  )],
  by = c("cancer_time_since_diag_other", "id"),
  all.x = T
)

# add sig indicators
comp_cancer_time_since_diag_one_primary_adj[, nonsig_vs_healthy := between(0, comp_cancer_time_since_diag_one_primary_adj$CI_low_diff_ref_healthy, comp_cancer_time_since_diag_one_primary_adj$CI_high_diff_ref_healthy)]
# comp_cancer_time_since_diag_one_primary_adj[, nonsig_vs_others  := between(0, comp_cancer_time_since_diag_one_primary_adj$CI_low_diff_ref_others, comp_cancer_time_since_diag_one_primary_adj$CI_high_diff_ref_others)]
comp_cancer_time_since_diag_one_primary_adj[, nonsig_vs_cancer  := between(0, comp_cancer_time_since_diag_one_primary_adj$CI_low_diff_ref_cancer, comp_cancer_time_since_diag_one_primary_adj$CI_high_diff_ref_cancer)]

comp_cancer_time_since_diag_one_primary_adj[, sig_ref_healthy := ifelse(nonsig_vs_healthy == FALSE & Mean_diff_ref_healthy != 0, 
                                                                  "$^a$", "$\\phantom{^a}$")]
# comp_cancer_time_since_diag_one_primary_adj[, sig_ref_others  := ifelse(nonsig_vs_others == FALSE & Mean_diff_ref_others != 0, 
#                                                                   "$^b$", "$\\phantom{^b}$")]

comp_cancer_time_since_diag_one_primary_adj[, sig_ref_cancer_15vs1  := ifelse(nonsig_vs_cancer == FALSE & !is.na(Mean_diff_ref_cancer) & cancer_time_since_diag_other == "Less than 1 year since diagnosis",
                                                                          "$^c$", "$\\phantom{^c}$")]
comp_cancer_time_since_diag_one_primary_adj[, sig_ref_cancer_5vs15  := ifelse(nonsig_vs_cancer == FALSE & !is.na(Mean_diff_ref_cancer) & cancer_time_since_diag_other == "1-5 years since diagnosis",
                                                                          "$^d$", "$\\phantom{^d}$")]
comp_cancer_time_since_diag_one_primary_adj[, sig_ref_cancer_5vs1  := ifelse(nonsig_vs_cancer == FALSE & !is.na(Mean_diff_ref_cancer) & cancer_time_since_diag_other == "More than 5 years since diagnosis",
                                                                         "$^e$", "$\\phantom{^e}$")]
comp_cancer_time_since_diag_one_primary_adj[, sig_ref_cancer  := paste0(sig_ref_cancer_15vs1, sig_ref_cancer_5vs15, sig_ref_cancer_5vs1)]

# leave healthy empty
comp_cancer_time_since_diag_one_primary_adj[, sig_ref_healthy := ifelse(cancer_time_since_diag_other == "Healthy", "$\\phantom{^a}$", sig_ref_healthy)]
# comp_cancer_time_since_diag_one_primary_adj[, sig_ref_others := ifelse(cancer_time_since_diag_other == "Healthy", "$\\phantom{^b}$", sig_ref_others)]

comp_cancer_time_since_diag_one_primary_adj[, yintercept_healthy := NA]
comp_cancer_time_since_diag_one_primary_adj[, yintercept_healthy := ifelse(part == "sleep", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Healthy" & part == "sleep"]$Mean, yintercept_healthy)]
comp_cancer_time_since_diag_one_primary_adj[, yintercept_healthy := ifelse(part == "mvpa", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Healthy" & part == "mvpa"]$Mean, yintercept_healthy)]
comp_cancer_time_since_diag_one_primary_adj[, yintercept_healthy := ifelse(part == "lpa", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Healthy" & part == "lpa"]$Mean, yintercept_healthy)]
comp_cancer_time_since_diag_one_primary_adj[, yintercept_healthy := ifelse(part == "sb", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Healthy" & part == "sb"]$Mean, yintercept_healthy)]

comp_cancer_time_since_diag_one_primary_adj[, ci_low_healthy := NA]
comp_cancer_time_since_diag_one_primary_adj[, ci_low_healthy := ifelse(part == "sleep", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Healthy" & part == "sleep"]$CI_low, ci_low_healthy)]
comp_cancer_time_since_diag_one_primary_adj[, ci_low_healthy := ifelse(part == "mvpa", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Healthy" & part == "mvpa"]$CI_low, ci_low_healthy)]
comp_cancer_time_since_diag_one_primary_adj[, ci_low_healthy := ifelse(part == "lpa", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Healthy" & part == "lpa"]$CI_low, ci_low_healthy)]
comp_cancer_time_since_diag_one_primary_adj[, ci_low_healthy := ifelse(part == "sb", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Healthy" & part == "sb"]$CI_low, ci_low_healthy)]

comp_cancer_time_since_diag_one_primary_adj[, ci_high_healthy := NA]
comp_cancer_time_since_diag_one_primary_adj[, ci_high_healthy := ifelse(part == "sleep", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Healthy" & part == "sleep"]$CI_high, ci_high_healthy)]
comp_cancer_time_since_diag_one_primary_adj[, ci_high_healthy := ifelse(part == "mvpa", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Healthy" & part == "mvpa"]$CI_high, ci_high_healthy)]
comp_cancer_time_since_diag_one_primary_adj[, ci_high_healthy := ifelse(part == "lpa", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Healthy" & part == "lpa"]$CI_high, ci_high_healthy)]
comp_cancer_time_since_diag_one_primary_adj[, ci_high_healthy := ifelse(part == "sb", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Healthy" & part == "sb"]$CI_high, ci_high_healthy)]

comp_cancer_time_since_diag_one_primary_adj[, yintercept_others := NA]
comp_cancer_time_since_diag_one_primary_adj[, yintercept_others := ifelse(part == "sleep", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Others" & part == "sleep"]$Mean, yintercept_others)]
comp_cancer_time_since_diag_one_primary_adj[, yintercept_others := ifelse(part == "mvpa", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Others" & part == "mvpa"]$Mean, yintercept_others)]
comp_cancer_time_since_diag_one_primary_adj[, yintercept_others := ifelse(part == "lpa", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Others" & part == "lpa"]$Mean, yintercept_others)]
comp_cancer_time_since_diag_one_primary_adj[, yintercept_others := ifelse(part == "sb", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Others" & part == "sb"]$Mean, yintercept_others)]

comp_cancer_time_since_diag_one_primary_adj[, ci_low_others := NA]
comp_cancer_time_since_diag_one_primary_adj[, ci_low_others := ifelse(part == "sleep", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Others" & part == "sleep"]$CI_low, ci_low_others)]
comp_cancer_time_since_diag_one_primary_adj[, ci_low_others := ifelse(part == "mvpa", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Others" & part == "mvpa"]$CI_low, ci_low_others)]
comp_cancer_time_since_diag_one_primary_adj[, ci_low_others := ifelse(part == "lpa", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Others" & part == "lpa"]$CI_low, ci_low_others)]
comp_cancer_time_since_diag_one_primary_adj[, ci_low_others := ifelse(part == "sb", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Others" & part == "sb"]$CI_low, ci_low_others)]

comp_cancer_time_since_diag_one_primary_adj[, ci_high_others := NA]
comp_cancer_time_since_diag_one_primary_adj[, ci_high_others := ifelse(part == "sleep", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Others" & part == "sleep"]$CI_high, ci_high_others)]
comp_cancer_time_since_diag_one_primary_adj[, ci_high_others := ifelse(part == "mvpa", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Others" & part == "mvpa"]$CI_high, ci_high_others)]
comp_cancer_time_since_diag_one_primary_adj[, ci_high_others := ifelse(part == "lpa", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Others" & part == "lpa"]$CI_high, ci_high_others)]
comp_cancer_time_since_diag_one_primary_adj[, ci_high_others := ifelse(part == "sb", comp_cancer_time_since_diag_one_primary_adj[cancer_time_since_diag_other == "Others" & part == "sb"]$CI_high, ci_high_others)]

# n
comp_cancer_time_since_diag_one_primary_adj[, Cases := NA]
comp_cancer_time_since_diag_one_primary_adj[, Cases := ifelse(cancer_time_since_diag_other == "Healthy", "13 722", Cases)]
comp_cancer_time_since_diag_one_primary_adj[, Cases := ifelse(cancer_time_since_diag_other == "Cancer", "8 438", Cases)]
comp_cancer_time_since_diag_one_primary_adj[, Cases := ifelse(cancer_time_since_diag_other == "Others", "67 478", Cases)]
comp_cancer_time_since_diag_one_primary_adj[, Cases := ifelse(cancer_time_since_diag_other == "More than 5 years since diagnosis", "4 884", Cases)]
comp_cancer_time_since_diag_one_primary_adj[, Cases := ifelse(cancer_time_since_diag_other == "1-5 years since diagnosis", "2 803", Cases)]
comp_cancer_time_since_diag_one_primary_adj[, Cases := ifelse(cancer_time_since_diag_other == "Less than 1 year since diagnosis", "751", Cases)]

## healthy in first row of plot
comp_cancer_time_since_diag_one_primary_adj[, cancer_time_since_diag_other := ifelse(cancer_time_since_diag_other == "More than 5 years since diagnosis", "    >5 years since diagnosis", cancer_time_since_diag_other)]
comp_cancer_time_since_diag_one_primary_adj[, cancer_time_since_diag_other := ifelse(cancer_time_since_diag_other == "1-5 years since diagnosis",         "    1-5 years since diagnosis", cancer_time_since_diag_other)]
comp_cancer_time_since_diag_one_primary_adj[, cancer_time_since_diag_other := ifelse(cancer_time_since_diag_other == "Less than 1 year since diagnosis",  "    <1 year since diagnosis", cancer_time_since_diag_other)]
comp_cancer_time_since_diag_one_primary_adj[, cancer_time_since_diag_other := ifelse(cancer_time_since_diag_other == "Others", "Other Conditions", cancer_time_since_diag_other)]

comp_cancer_time_since_diag_one_primary_adj[, cancer_time_since_diag_other := factor(cancer_time_since_diag_other, ordered = TRUE,
                                                                               levels = c(
                                                                                 "    >5 years since diagnosis",
                                                                                 "    1-5 years since diagnosis",
                                                                                 "    <1 year since diagnosis",
                                                                                 "Cancer",
                                                                                 "Other Conditions",
                                                                                 "Healthy"))]

comp_cancer_time_since_diag_one_primary_adj[, part := ifelse(part == "sleep", "Sleep period", part)]
comp_cancer_time_since_diag_one_primary_adj[, part := ifelse(part == "mvpa", "Moderate-to-vigorous physical activity", part)]
comp_cancer_time_since_diag_one_primary_adj[, part := ifelse(part == "lpa", "Light physical activity", part)]
comp_cancer_time_since_diag_one_primary_adj[, part := ifelse(part == "sb", "Sedentary behaviour", part)]

comp_cancer_time_since_diag_one_primary_adj[, sig_position := min(CI_low), by = part]
comp_cancer_time_since_diag_one_primary_adj[, est_position := max(CI_high), by = part]

comp_cancer_time_since_diag_one_primary_adj[, estimates := paste0(round(Mean, 0), "[", round(CI_low, 0), ", ", round(CI_high, 0), "]")]
comp_cancer_time_since_diag_one_primary_adj[, est_sig := paste0(estimates, " ", sig_ref_healthy, sig_ref_cancer)]
# comp_cancer_time_since_diag_one_primary_adj[, est_sig := paste0(estimates, " ", str_replace_na(sig_ref_healthy, " "), str_replace_na(sig_ref_others, " "), str_replace_na(sig_ref_cancer, " "))]

# for tables
comp_cancer_time_since_diag_one_primary_adj[, estimates_contrast_healthy := paste0(format(round(Mean_diff_ref_healthy, 1)), " [", format(round(CI_low_diff_ref_healthy, 1)), ", ", format(round(CI_high_diff_ref_healthy, 1)), "]")]
# comp_cancer_time_since_diag_one_primary_adj[, estimates_contrast_others := paste0(format(round(Mean_diff_ref_others, 1)), " [", format(round(CI_low_diff_ref_others, 1)), ", ", format(round(CI_high_diff_ref_others, 1)), "]")]
comp_cancer_time_since_diag_one_primary_adj[, estimates_contrast_cancer := paste0(format(round(Mean_diff_ref_cancer, 1)), " [", format(round(CI_low_diff_ref_cancer, 1)), ", ", format(round(CI_high_diff_ref_cancer, 1)), "]")]

## plot -----------------------
plot_specs <- data.table(
 part = c(
  "Sleep period",
  "Moderate-to-vigorous physical activity",
  "Light physical activity",
  "Sedentary behaviour"
 ),
 ymin = c(500, 0, 200, 500),
 ymax = c(650, 50, 400, 650),
 ysig = c(627.5, 43, 370, 627.5),
 pal = list(pal_combined[-5], pal_combined[-5], pal_combined[-5], pal_combined[-5]),
 text_y = c(500, 0, 200, 500),
 seg_y = list(c(500, 650), c(0, 50), c(200, 400), c(500, 650)),
 ycase  = c(650, 50, 400, 650),
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

  data_i   <- comp_cancer_time_since_diag_one_primary_adj[part == spec$part]

  ggplot(data_i, aes(x = cancer_time_since_diag_other, y = Mean)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others),
              fill = "#F2F2F2", alpha = 0.2) +
    geom_hline(aes(yintercept = yintercept_others),
               linewidth = 0.5, linetype = "dashed", colour = "#A9A9A9") +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy),
              fill = "#A9B4AE", alpha = 0.075) + ## CAD4CF
    geom_hline(aes(yintercept = yintercept_healthy),
               linewidth = 0.5, linetype = "dashed", colour = "#708885") +

    geom_pointrange(aes(ymin = CI_low, ymax = CI_high, colour = cancer_time_since_diag_other),
                    size = 0.25, linewidth = 0.5) +
    geom_text(aes(y = spec$text_y, label = cancer_time_since_diag_other),
              hjust = 0, nudge_x = 0, family = "Arial Narrow", size = 3, show.legend = FALSE) +
    geom_text(aes(y = spec$ysig, label = TeX(est_sig, output = "character")),
              parse = TRUE, hjust = 0.5, family = "Arial Narrow", size = 3, show.legend = FALSE) +
    geom_text(aes(y = spec$ycase, label = Cases),
              hjust = 1, nudge_x = 0, family = "Arial Narrow", size = 3, show.legend = FALSE) +
    # annotate("segment",
    #          x = 0.5,
    #          xend = length(unique(data_i$cancer_time_since_diag_other)) + 0.5,
    #          y = seg_vals[1],
    #          yend = seg_vals[1],
    #          linewidth = 0.5) +
    # annotate("segment",
    #          x = 0.5,
    #          xend = length(unique(data_i$cancer_time_since_diag_other)) + 0.5,
    #          y = seg_vals[2],
    #          yend = seg_vals[2],
    #          linewidth = 0.5) +
    geom_segment(aes(x = 0, yend = seg_vals[1]), col = "black", linewidth = 0.5) +
    geom_segment(aes(x = 0, yend = seg_vals[2]), col = "black", linewidth = 0.5) +
    scale_y_continuous(
   limits = c(spec$ymin, spec$ymax),
   breaks = c(spec$ymin, spec$ymax),
   name = spec$ylab
  ) +
    scale_colour_manual(values = spec$pal[[1]]) +
    labs(x = "", y = "", colour = "") +
    coord_flip() +
    theme_ipsum() +
    theme(
      axis.ticks      = element_blank(),
      plot.background = element_rect(fill = "transparent", colour = NA, linewidth = 0.5),
      panel.grid.major= element_blank(),
      panel.grid.minor= element_blank(),
      axis.title.x    = element_text(size = 9, face = "bold", hjust = .5, margin = margin(t = -9)),
      axis.text.x     = element_text(size = 9),
      axis.text.y     = element_blank(),
      strip.text      = element_text(size = 9, hjust = .5, face = "bold"),
      legend.text     = element_text(size = 9, face = "bold", hjust = .5),
      legend.position = "none",
      plot.margin     = unit(c(0.5, 0, 1, 0), "lines")
    )
})

names(plot_list) <- tolower(gsub(" ", "_", plot_specs$part))
plot_comp_cancer_time_since_diag_one_primary_sleep <- plot_list[["sleep_period"]]
plot_comp_cancer_time_since_diag_one_primary_mvpa  <- plot_list[["moderate-to-vigorous_physical_activity"]]
plot_comp_cancer_time_since_diag_one_primary_lpa   <- plot_list[["light_physical_activity"]]
plot_comp_cancer_time_since_diag_one_primary_sb    <- plot_list[["sedentary_behaviour"]]

grDevices::cairo_pdf(
  file = paste0(outputdir, "cancer_time_since_diag_one_primary_est", ".pdf"),
  width = 6,
  height = 8,
)

ggarrange(
  plot_comp_cancer_time_since_diag_one_primary_mvpa,
  plot_comp_cancer_time_since_diag_one_primary_lpa,
  plot_comp_cancer_time_since_diag_one_primary_sb,
  plot_comp_cancer_time_since_diag_one_primary_sleep,
  nrow = 4
)
dev.off()

grDevices::png(
  file = paste0(outputdir, "cancer_time_since_diag_one_primary_est", ".png"),
  width = 6000,
  height = 8000,
  res = 900
)

ggarrange(
  plot_comp_cancer_time_since_diag_one_primary_mvpa,
  plot_comp_cancer_time_since_diag_one_primary_lpa,
  plot_comp_cancer_time_since_diag_one_primary_sb,
  plot_comp_cancer_time_since_diag_one_primary_sleep,
  nrow = 4
)
dev.off()
