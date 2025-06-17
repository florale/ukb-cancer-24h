source("ukb-cancer-24h-setup.R")
source(paste0(redir, "ukb_utils.R"))
# source("ukb-cancer-24h-data.R")

# main model --------
# fit_cancer_time_since_diag_other_adj <- brmcoda(clr_cancer_acc,
#                                           mvbind(ilr1, ilr2, ilr3) ~ cancer_time_since_diag_other +
#                                             # + other_conds_at_acc +
#                                             s(age_at_acc) + sex + white + working + edu + never_smoked + current_drinker + s(deprivation),
#                                           # save_pars = save_pars(all = TRUE),
#                                           warmup = 500, chains = 4, cores = 4, backend = "cmdstanr"
# )
# saveRDS(fit_cancer_time_since_diag_other_adj, paste0(outputdir, "fit_cancer_time_since_diag_other_adj", ".RDS"))

# predicted posteriors ------------
fit_cancer_time_since_diag_other_adj <- readRDS(paste0(outputdir, "fit_cancer_time_since_diag_other_adj", ".RDS"))

# reference grid
d_cancer_time_since_diag_other_adj <- emmeans::ref_grid(fit_cancer_time_since_diag_other_adj$model)@grid

# predict
pred_cancer_time_since_diag_other_adj <- fitted(fit_cancer_time_since_diag_other_adj, newdata = d_cancer_time_since_diag_other_adj, scale = "response", summary = FALSE)

# weight by length equal to the number of observations in the dataset
# summarise by cancer group
pred_cancer_time_since_diag_other_adj <- apply(pred_cancer_time_since_diag_other_adj, c(1), function(x)  cbind(d_cancer_time_since_diag_other_adj, x))
pred_cancer_time_since_diag_other_adj <- lapply(pred_cancer_time_since_diag_other_adj, function(d) {
  
  parts <- c("sleep", "mvpa", "lpa", "sb")
  parts0 <- c("V1", "V2", "V3", "V4")
  
  d <- as.data.table(d)
  
  d[, cancer_weight := sum(.wgt.), by = cancer_time_since_diag_other]
  d[, denominator := sum(cancer_weight)]
  d[, cancer_denominator := denominator - 
      (d[cancer_time_since_diag_other == "Others"]$cancer_weight[[1]])*nrow(d[cancer_time_since_diag_other == "Others"]) - 
      (d[cancer_time_since_diag_other == "Healthy"]$cancer_weight[[1]])*nrow(d[cancer_time_since_diag_other == "Healthy"])]
  
  d[, (parts) := lapply(.SD, mean), by = cancer_time_since_diag_other, .SDcols =  parts0]

  d[, sleep_weighted := ifelse(cancer_time_since_diag_other %nin% c("Others", "Healthy"), sleep*(cancer_weight/cancer_denominator), NA)]
  d[, mvpa_weighted := ifelse(cancer_time_since_diag_other %nin% c("Others", "Healthy"), mvpa*(cancer_weight/cancer_denominator), NA)]
  d[, lpa_weighted := ifelse(cancer_time_since_diag_other %nin% c("Others", "Healthy"), lpa*(cancer_weight/cancer_denominator), NA)]
  d[, sb_weighted := ifelse(cancer_time_since_diag_other %nin% c("Others", "Healthy"), sb*(cancer_weight/cancer_denominator), NA)]
  
  d[, (paste0(parts, "_cancer")) := lapply(.SD, sum, na.rm = TRUE), .SDcols = paste0(parts, "_weighted")]
  
  d <- rbind(d, 
             data.table(cancer_time_since_diag_other = "Cancer"),
             fill = TRUE
  )
  
  d[cancer_time_since_diag_other == "Cancer", sleep := d$sleep_cancer[[1]]]
  d[cancer_time_since_diag_other == "Cancer", mvpa := d$mvpa_cancer[[1]]]
  d[cancer_time_since_diag_other == "Cancer", lpa := d$lpa_cancer[[1]]]
  d[cancer_time_since_diag_other == "Cancer", sb := d$sb_cancer[[1]]]
  
  # contrast vs healthy
  d[, sleep_vs_healthy := sleep - d[cancer_time_since_diag_other == "Healthy"]$sleep]
  d[, mvpa_vs_healthy := mvpa - d[cancer_time_since_diag_other == "Healthy"]$mvpa]
  d[, lpa_vs_healthy := lpa - d[cancer_time_since_diag_other == "Healthy"]$lpa]
  d[, sb_vs_healthy := sb - d[cancer_time_since_diag_other == "Healthy"]$sb]
  
  # contrast vs others
  d[, sleep_vs_others := sleep - d[cancer_time_since_diag_other == "Others"]$sleep]
  d[, mvpa_vs_others := mvpa - d[cancer_time_since_diag_other == "Others"]$mvpa]
  d[, lpa_vs_others := lpa - d[cancer_time_since_diag_other == "Others"]$lpa]
  d[, sb_vs_others := sb - d[cancer_time_since_diag_other == "Others"]$sb]
  
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
             sleep_vs_others, mvpa_vs_others, lpa_vs_others, sb_vs_others,
             
             sleep_lessthan1_vs_1to5, sleep_lessthan1_vs_morethan5, sleep_1to5_vs_morethan5,
             mvpa_lessthan1_vs_1to5, mvpa_lessthan1_vs_morethan5, mvpa_1to5_vs_morethan5,
             lpa_lessthan1_vs_1to5, lpa_lessthan1_vs_morethan5, lpa_1to5_vs_morethan5,
             sb_lessthan1_vs_1to5, sb_lessthan1_vs_morethan5, sb_1to5_vs_morethan5
             
  )]
  d <- unique(d)
  d
})

# assemble back to summarise posteriors
pred_cancer_time_since_diag_other_adj <- as.data.frame(abind::abind(pred_cancer_time_since_diag_other_adj, along = 1))
pred_cancer_time_since_diag_other_adj <- split(pred_cancer_time_since_diag_other_adj, pred_cancer_time_since_diag_other_adj$cancer_time_since_diag_other)

## estimated means  ----------------------
pred_comp_cancer_time_since_diag_other_adj <- lapply(pred_cancer_time_since_diag_other_adj, function(l) {
  l <- as.data.frame(l[, c("sleep", "mvpa", "lpa", "sb")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, function(p) describe_posterior(p, centrality = "mean", ci = 0.99))
  l <- Map(cbind, l, part = names(l))
  l <- rbindlist(l)
  l
})
pred_comp_cancer_time_since_diag_other_adj <- Map(cbind, pred_comp_cancer_time_since_diag_other_adj, cancer_time_since_diag_other = names(pred_comp_cancer_time_since_diag_other_adj))
pred_comp_cancer_time_since_diag_other_adj <- rbindlist(pred_comp_cancer_time_since_diag_other_adj)

## contrasts --------------------
### vs healthy
diff1_comp_cancer_time_since_diag_other_adj <- lapply(pred_cancer_time_since_diag_other_adj, function(l) {
  l <- as.data.frame(l[, c("sleep_vs_healthy", "mvpa_vs_healthy", "lpa_vs_healthy", "sb_vs_healthy")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, function(p) describe_posterior(p, centrality = "mean", ci = 0.99))
  l <- Map(cbind, l, contrast = names(l))
  l <- rbindlist(l)
  l
})
diff1_comp_cancer_time_since_diag_other_adj <- Map(cbind, diff1_comp_cancer_time_since_diag_other_adj, cancer_time_since_diag_other = names(diff1_comp_cancer_time_since_diag_other_adj))
diff1_comp_cancer_time_since_diag_other_adj <- rbindlist(diff1_comp_cancer_time_since_diag_other_adj)

setnames(diff1_comp_cancer_time_since_diag_other_adj, "Mean", "Mean_diff_ref_healthy")
setnames(diff1_comp_cancer_time_since_diag_other_adj, "CI_low", "CI_low_diff_ref_healthy")
setnames(diff1_comp_cancer_time_since_diag_other_adj, "CI_high", "CI_high_diff_ref_healthy")

### vs others
diff3_comp_cancer_time_since_diag_other_adj <- lapply(pred_cancer_time_since_diag_other_adj, function(l) {
  l <- as.data.frame(l[, c("sleep_vs_others", "mvpa_vs_others", "lpa_vs_others", "sb_vs_others")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, function(p) describe_posterior(p, centrality = "mean", ci = 0.99))
  l <- Map(cbind, l, contrast = names(l))
  l <- rbindlist(l)
  l
})
diff3_comp_cancer_time_since_diag_other_adj <- Map(cbind, diff3_comp_cancer_time_since_diag_other_adj, cancer_time_since_diag_other = names(diff3_comp_cancer_time_since_diag_other_adj))
diff3_comp_cancer_time_since_diag_other_adj <- rbindlist(diff3_comp_cancer_time_since_diag_other_adj)

setnames(diff3_comp_cancer_time_since_diag_other_adj, "Mean", "Mean_diff_ref_others")
setnames(diff3_comp_cancer_time_since_diag_other_adj, "CI_low", "CI_low_diff_ref_others")
setnames(diff3_comp_cancer_time_since_diag_other_adj, "CI_high", "CI_high_diff_ref_others")

### pairwise cancer
diff2_comp_cancer_time_since_diag_other_adj <- lapply(pred_cancer_time_since_diag_other_adj, function(l) {
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
# diff2_comp_cancer_time_since_diag_other_adj <- Map(cbind, diff2_comp_cancer_time_since_diag_other_adj, cancer_time_since_diag_other = names(diff2_comp_cancer_time_since_diag_other_adj))
diff2_comp_cancer_time_since_diag_other_adj <- rbindlist(diff2_comp_cancer_time_since_diag_other_adj)

diff2_comp_cancer_time_since_diag_other_adj[, cancer_time_since_diag_other := NA]
diff2_comp_cancer_time_since_diag_other_adj[, cancer_time_since_diag_other := ifelse(grepl("lessthan1_vs_1to5", diff2_comp_cancer_time_since_diag_other_adj$contrast), "Less than 1 year since diagnosis", cancer_time_since_diag_other)]
diff2_comp_cancer_time_since_diag_other_adj[, cancer_time_since_diag_other := ifelse(grepl("1to5_vs_morethan5", diff2_comp_cancer_time_since_diag_other_adj$contrast), "1-5 years since diagnosis", cancer_time_since_diag_other)]
diff2_comp_cancer_time_since_diag_other_adj[, cancer_time_since_diag_other := ifelse(grepl("lessthan1_vs_morethan5", diff2_comp_cancer_time_since_diag_other_adj$contrast), "More than 5 years since diagnosis", cancer_time_since_diag_other)]

setnames(diff2_comp_cancer_time_since_diag_other_adj, "Mean", "Mean_diff_ref_cancer")
setnames(diff2_comp_cancer_time_since_diag_other_adj, "CI_low", "CI_low_diff_ref_cancer")
setnames(diff2_comp_cancer_time_since_diag_other_adj, "CI_high", "CI_high_diff_ref_cancer")

# take only unique values
diff2_comp_cancer_time_since_diag_other_adj <- unique(diff2_comp_cancer_time_since_diag_other_adj)

# add id to merge later
diff2_comp_cancer_time_since_diag_other_adj[, id := 1:.N, by = cancer_time_since_diag_other]

# all results  ------------------------
comp_cancer_time_since_diag_other_adj <- cbind(
  pred_comp_cancer_time_since_diag_other_adj[, .(Mean, CI_low, CI_high, part, cancer_time_since_diag_other)],
  diff1_comp_cancer_time_since_diag_other_adj[, .(Mean_diff_ref_healthy, CI_low_diff_ref_healthy, CI_high_diff_ref_healthy)],
  diff3_comp_cancer_time_since_diag_other_adj[, .(Mean_diff_ref_others, CI_low_diff_ref_others, CI_high_diff_ref_others)]
  
)
comp_cancer_time_since_diag_other_adj[, id := 1:.N, by = cancer_time_since_diag_other]

comp_cancer_time_since_diag_other_adj <- merge(
  comp_cancer_time_since_diag_other_adj,
  diff2_comp_cancer_time_since_diag_other_adj[, .(Mean_diff_ref_cancer,
                                            CI_low_diff_ref_cancer,
                                            CI_high_diff_ref_cancer,
                                            cancer_time_since_diag_other,
                                            id
  )],
  by = c("cancer_time_since_diag_other", "id"),
  all.x = T
)


# add sig indicators
comp_cancer_time_since_diag_other_adj[, nonsig_vs_healthy := between(0, comp_cancer_time_since_diag_other_adj$CI_low_diff_ref_healthy, comp_cancer_time_since_diag_other_adj$CI_high_diff_ref_healthy)]
comp_cancer_time_since_diag_other_adj[, nonsig_vs_others  := between(0, comp_cancer_time_since_diag_other_adj$CI_low_diff_ref_others, comp_cancer_time_since_diag_other_adj$CI_high_diff_ref_others)]
comp_cancer_time_since_diag_other_adj[, nonsig_vs_cancer  := between(0, comp_cancer_time_since_diag_other_adj$CI_low_diff_ref_cancer, comp_cancer_time_since_diag_other_adj$CI_high_diff_ref_cancer)]

comp_cancer_time_since_diag_other_adj[, sig_ref_healthy := ifelse(nonsig_vs_healthy == FALSE & Mean_diff_ref_healthy != 0, 
                                                                  "$^a$", "$\\phantom{^a}$")]
comp_cancer_time_since_diag_other_adj[, sig_ref_others  := ifelse(nonsig_vs_others == FALSE & Mean_diff_ref_others != 0, 
                                                                  "$^b$", "$\\phantom{^b}$")]

comp_cancer_time_since_diag_other_adj[, sig_ref_cancer_15vs1  := ifelse(nonsig_vs_cancer == FALSE & !is.na(Mean_diff_ref_cancer) & cancer_time_since_diag_other == "Less than 1 year since diagnosis",
                                                                          "$^c$", "$\\phantom{^c}$")]
comp_cancer_time_since_diag_other_adj[, sig_ref_cancer_5vs15  := ifelse(nonsig_vs_cancer == FALSE & !is.na(Mean_diff_ref_cancer) & cancer_time_since_diag_other == "1-5 years since diagnosis",
                                                                          "$^d$", "$\\phantom{^d}$")]
comp_cancer_time_since_diag_other_adj[, sig_ref_cancer_5vs1  := ifelse(nonsig_vs_cancer == FALSE & !is.na(Mean_diff_ref_cancer) & cancer_time_since_diag_other == "More than 5 years since diagnosis",
                                                                         "$^e$", "$\\phantom{^e}$")]
comp_cancer_time_since_diag_other_adj[, sig_ref_cancer  := paste0(sig_ref_cancer_15vs1, sig_ref_cancer_5vs15, sig_ref_cancer_5vs1)]

# leave healthy empty
comp_cancer_time_since_diag_other_adj[, sig_ref_healthy := ifelse(cancer_time_since_diag_other == "Healthy", "$\\phantom{^a}$", sig_ref_healthy)]
comp_cancer_time_since_diag_other_adj[, sig_ref_others := ifelse(cancer_time_since_diag_other == "Healthy", "$\\phantom{^b}$", sig_ref_others)]

comp_cancer_time_since_diag_other_adj[, yintercept_healthy := NA]
comp_cancer_time_since_diag_other_adj[, yintercept_healthy := ifelse(part == "sleep", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Healthy" & part == "sleep"]$Mean, yintercept_healthy)]
comp_cancer_time_since_diag_other_adj[, yintercept_healthy := ifelse(part == "mvpa", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Healthy" & part == "mvpa"]$Mean, yintercept_healthy)]
comp_cancer_time_since_diag_other_adj[, yintercept_healthy := ifelse(part == "lpa", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Healthy" & part == "lpa"]$Mean, yintercept_healthy)]
comp_cancer_time_since_diag_other_adj[, yintercept_healthy := ifelse(part == "sb", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Healthy" & part == "sb"]$Mean, yintercept_healthy)]

comp_cancer_time_since_diag_other_adj[, ci_low_healthy := NA]
comp_cancer_time_since_diag_other_adj[, ci_low_healthy := ifelse(part == "sleep", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Healthy" & part == "sleep"]$CI_low, ci_low_healthy)]
comp_cancer_time_since_diag_other_adj[, ci_low_healthy := ifelse(part == "mvpa", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Healthy" & part == "mvpa"]$CI_low, ci_low_healthy)]
comp_cancer_time_since_diag_other_adj[, ci_low_healthy := ifelse(part == "lpa", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Healthy" & part == "lpa"]$CI_low, ci_low_healthy)]
comp_cancer_time_since_diag_other_adj[, ci_low_healthy := ifelse(part == "sb", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Healthy" & part == "sb"]$CI_low, ci_low_healthy)]

comp_cancer_time_since_diag_other_adj[, ci_high_healthy := NA]
comp_cancer_time_since_diag_other_adj[, ci_high_healthy := ifelse(part == "sleep", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Healthy" & part == "sleep"]$CI_high, ci_high_healthy)]
comp_cancer_time_since_diag_other_adj[, ci_high_healthy := ifelse(part == "mvpa", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Healthy" & part == "mvpa"]$CI_high, ci_high_healthy)]
comp_cancer_time_since_diag_other_adj[, ci_high_healthy := ifelse(part == "lpa", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Healthy" & part == "lpa"]$CI_high, ci_high_healthy)]
comp_cancer_time_since_diag_other_adj[, ci_high_healthy := ifelse(part == "sb", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Healthy" & part == "sb"]$CI_high, ci_high_healthy)]

comp_cancer_time_since_diag_other_adj[, yintercept_others := NA]
comp_cancer_time_since_diag_other_adj[, yintercept_others := ifelse(part == "sleep", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Others" & part == "sleep"]$Mean, yintercept_others)]
comp_cancer_time_since_diag_other_adj[, yintercept_others := ifelse(part == "mvpa", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Others" & part == "mvpa"]$Mean, yintercept_others)]
comp_cancer_time_since_diag_other_adj[, yintercept_others := ifelse(part == "lpa", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Others" & part == "lpa"]$Mean, yintercept_others)]
comp_cancer_time_since_diag_other_adj[, yintercept_others := ifelse(part == "sb", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Others" & part == "sb"]$Mean, yintercept_others)]

comp_cancer_time_since_diag_other_adj[, ci_low_others := NA]
comp_cancer_time_since_diag_other_adj[, ci_low_others := ifelse(part == "sleep", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Others" & part == "sleep"]$CI_low, ci_low_others)]
comp_cancer_time_since_diag_other_adj[, ci_low_others := ifelse(part == "mvpa", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Others" & part == "mvpa"]$CI_low, ci_low_others)]
comp_cancer_time_since_diag_other_adj[, ci_low_others := ifelse(part == "lpa", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Others" & part == "lpa"]$CI_low, ci_low_others)]
comp_cancer_time_since_diag_other_adj[, ci_low_others := ifelse(part == "sb", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Others" & part == "sb"]$CI_low, ci_low_others)]

comp_cancer_time_since_diag_other_adj[, ci_high_others := NA]
comp_cancer_time_since_diag_other_adj[, ci_high_others := ifelse(part == "sleep", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Others" & part == "sleep"]$CI_high, ci_high_others)]
comp_cancer_time_since_diag_other_adj[, ci_high_others := ifelse(part == "mvpa", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Others" & part == "mvpa"]$CI_high, ci_high_others)]
comp_cancer_time_since_diag_other_adj[, ci_high_others := ifelse(part == "lpa", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Others" & part == "lpa"]$CI_high, ci_high_others)]
comp_cancer_time_since_diag_other_adj[, ci_high_others := ifelse(part == "sb", comp_cancer_time_since_diag_other_adj[cancer_time_since_diag_other == "Others" & part == "sb"]$CI_high, ci_high_others)]

# n
comp_cancer_time_since_diag_other_adj[, Cases := NA]
comp_cancer_time_since_diag_other_adj[, Cases := ifelse(cancer_time_since_diag_other == "Healthy", "13 722", Cases)]
comp_cancer_time_since_diag_other_adj[, Cases := ifelse(cancer_time_since_diag_other == "Cancer", "10 152", Cases)]
comp_cancer_time_since_diag_other_adj[, Cases := ifelse(cancer_time_since_diag_other == "Others", "67 478", Cases)]
comp_cancer_time_since_diag_other_adj[, Cases := ifelse(cancer_time_since_diag_other == "More than 5 years since diagnosis", "5 730", Cases)]
comp_cancer_time_since_diag_other_adj[, Cases := ifelse(cancer_time_since_diag_other == "1-5 years since diagnosis", "3 451", Cases)]
comp_cancer_time_since_diag_other_adj[, Cases := ifelse(cancer_time_since_diag_other == "Less than 1 year since diagnosis", "971", Cases)]

## healthy in first row of plot
comp_cancer_time_since_diag_other_adj[, cancer_time_since_diag_other := ifelse(cancer_time_since_diag_other == "More than 5 years since diagnosis", "    >5 years since diagnosis", cancer_time_since_diag_other)]
comp_cancer_time_since_diag_other_adj[, cancer_time_since_diag_other := ifelse(cancer_time_since_diag_other == "1-5 years since diagnosis",         "    1-5 years since diagnosis", cancer_time_since_diag_other)]
comp_cancer_time_since_diag_other_adj[, cancer_time_since_diag_other := ifelse(cancer_time_since_diag_other == "Less than 1 year since diagnosis",  "    <1 year since diagnosis", cancer_time_since_diag_other)]
comp_cancer_time_since_diag_other_adj[, cancer_time_since_diag_other := ifelse(cancer_time_since_diag_other == "Others", "Other Conditions", cancer_time_since_diag_other)]

comp_cancer_time_since_diag_other_adj[, cancer_time_since_diag_other := factor(cancer_time_since_diag_other, ordered = TRUE,
                                                                               levels = c(
                                                                                 "    >5 years since diagnosis",
                                                                                 "    1-5 years since diagnosis",
                                                                                 "    <1 year since diagnosis",
                                                                                 "Cancer",
                                                                                 "Other Conditions",
                                                                                 "Healthy"))]

comp_cancer_time_since_diag_other_adj[, part := ifelse(part == "sleep", "Sleep period", part)]
comp_cancer_time_since_diag_other_adj[, part := ifelse(part == "mvpa", "Moderate-to-vigorous physical activity", part)]
comp_cancer_time_since_diag_other_adj[, part := ifelse(part == "lpa", "Light physical activity", part)]
comp_cancer_time_since_diag_other_adj[, part := ifelse(part == "sb", "Sedentary behaviour", part)]

comp_cancer_time_since_diag_other_adj[, sig_position := min(CI_low), by = part]
comp_cancer_time_since_diag_other_adj[, est_position := max(CI_high), by = part]

comp_cancer_time_since_diag_other_adj[, estimates := paste0(round(Mean, 0), "[", round(CI_low, 0), ", ", round(CI_high, 0), "]")]
comp_cancer_time_since_diag_other_adj[, est_sig := paste0(estimates, " ", sig_ref_healthy, sig_ref_others, sig_ref_cancer)]
# comp_cancer_time_since_diag_other_adj[, est_sig := paste0(estimates, " ", str_replace_na(sig_ref_healthy, " "), str_replace_na(sig_ref_others, " "), str_replace_na(sig_ref_cancer, " "))]

# for tables
comp_cancer_time_since_diag_other_adj[, estimates_contrast_healthy := paste0(round(Mean_diff_ref_healthy, 2), "[", round(CI_low_diff_ref_healthy, 2), ", ", round(CI_high_diff_ref_healthy, 2), "]")]
comp_cancer_time_since_diag_other_adj[, estimates_contrast_others := paste0(round(Mean_diff_ref_others, 2), "[", round(CI_low_diff_ref_others, 2), ", ", round(CI_high_diff_ref_others, 2), "]")]
comp_cancer_time_since_diag_other_adj[, estimates_contrast_cancer := paste0(round(Mean_diff_ref_cancer, 2), "[", round(CI_low_diff_ref_cancer, 2), ", ", round(CI_high_diff_ref_cancer, 2), "]")]

## plot -----------------------
(plot_comp_cancer_time_since_diag_other_sleep <- 
   ggplot(comp_cancer_time_since_diag_other_adj[part == "Sleep period"], aes(x = cancer_time_since_diag_other, y = Mean)) +
   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#CBD5D0", alpha = 0.1) +
   geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#F2F2F2", alpha = 0.2) + 
   geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#708885") +
   geom_hline(aes(yintercept = yintercept_others), linewidth = 0.5, linetype= "dashed", colour = "#A9A9A9") +
   geom_pointrange(aes(ymin = CI_low,
                       ymax = CI_high, colour = cancer_time_since_diag_other), size = 0.25, linewidth = 0.5) +
   geom_text(aes(y = 500, label = cancer_time_since_diag_other),
             hjust = 0, nudge_x = 0, 
             family = "Arial Narrow", size = 3,
             show.legend = FALSE) +
   geom_text(aes(y = 627.5, label = TeX(est_sig, output = "character")), parse = TRUE,
             hjust = 0.5, nudge_x = 0, 
             family = "Arial Narrow", size = 3,
             show.legend = FALSE) +
   geom_text(aes(y = 650, label = Cases),
             hjust = 1, nudge_x = 0, 
             family = "Arial Narrow", size = 3,
             show.legend = FALSE) +
   geom_segment(aes(x = 0, yend = 500), col = "black", linewidth = 0.5) +
   geom_segment(aes(x = 0, yend = 650), col = "black", linewidth = 0.5) +   
   scale_y_continuous(limits = c(500, 650),
                      breaks = c(500,  650),
                      name = "Sleep period (mins/day)") +
   scale_colour_manual(values = pal_combined) +
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

(plot_comp_cancer_time_since_diag_other_mvpa <- 
    ggplot(comp_cancer_time_since_diag_other_adj[part == "Moderate-to-vigorous physical activity"], aes(x = cancer_time_since_diag_other, y = Mean)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#CBD5D0", alpha = 0.1) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#F2F2F2", alpha = 0.2) + 
    geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#708885") +
    geom_hline(aes(yintercept = yintercept_others), linewidth = 0.5, linetype= "dashed", colour = "#A9A9A9") +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_time_since_diag_other), size = 0.25, linewidth = 0.5) +
    geom_text(aes(y = 0, label = cancer_time_since_diag_other),
              hjust = 0, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_text(aes(y = 33.75, label = TeX(est_sig, output = "character")), parse = TRUE,
              hjust = 0.5, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_text(aes(y = 40, label = Cases),
              hjust = 1, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_segment(aes(x = 0, yend = 0), col = "black", linewidth = 0.5) +
    geom_segment(aes(x = 0, yend = 40), col = "black", linewidth = 0.5) +
    scale_y_continuous(limits = c(0, 40),
                       breaks = c(0, 40),
                       name = "Moderate-to-vigorous physical activity (mins/day)") +
    scale_colour_manual(values = pal_combined) +
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

(plot_comp_cancer_time_since_diag_other_lpa <- 
    ggplot(comp_cancer_time_since_diag_other_adj[part == "Light physical activity"], aes(x = cancer_time_since_diag_other, y = Mean)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#CBD5D0", alpha = 0.1) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#F2F2F2", alpha = 0.2) + 
    geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#708885") +
    geom_hline(aes(yintercept = yintercept_others), linewidth = 0.5, linetype= "dashed", colour = "#A9A9A9") +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_time_since_diag_other), size = 0.25, linewidth = 0.5) +
    geom_text(aes(y = 200, label = cancer_time_since_diag_other),
              hjust = 0, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_text(aes(y = 371.5, label = TeX(est_sig, output = "character")), parse = TRUE,
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
    scale_colour_manual(values = pal_combined) +
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

(plot_comp_cancer_time_since_diag_other_sb <- 
    ggplot(comp_cancer_time_since_diag_other_adj[part == "Sedentary behaviour"], aes(x = cancer_time_since_diag_other, y = Mean)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_healthy, ymax = ci_high_healthy), fill = "#CBD5D0", alpha = 0.1) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = ci_low_others, ymax = ci_high_others), fill = "#F2F2F2", alpha = 0.2) +
    geom_hline(aes(yintercept = yintercept_healthy), linewidth = 0.5, linetype= "dashed", colour = "#708885") +
    geom_hline(aes(yintercept = yintercept_others), linewidth = 0.5, linetype= "dashed", colour = "#A9A9A9") +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_time_since_diag_other), size = 0.25, linewidth = 0.5) +
    geom_text(aes(y = 500, label = cancer_time_since_diag_other),
              hjust = 0, nudge_x = 0, 
              family = "Arial Narrow", size = 3,
              show.legend = FALSE) +
    geom_text(aes(y = 627.5, label = TeX(est_sig, output = "character")), parse = TRUE,
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
                       name = "Sedentary (mins/day)") +
    scale_colour_manual(values = pal_combined) +
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
  file = paste0(outputdir, "cancer_time_since_diag_other_est", ".pdf"),
  width = 6,
  height = 8,
)

ggarrange(
  plot_comp_cancer_time_since_diag_other_mvpa,
  plot_comp_cancer_time_since_diag_other_lpa,
  plot_comp_cancer_time_since_diag_other_sb,
  plot_comp_cancer_time_since_diag_other_sleep,
  nrow = 4
)
dev.off()

grDevices::png(
  file = paste0(outputdir, "cancer_time_since_diag_other_est", ".png"),
  width = 6000,
  height = 8000,
  res = 900
)

ggarrange(
  plot_comp_cancer_time_since_diag_other_mvpa,
  plot_comp_cancer_time_since_diag_other_lpa,
  plot_comp_cancer_time_since_diag_other_sb,
  plot_comp_cancer_time_since_diag_other_sleep,
  nrow = 4
)
dev.off()

# estimates for tables --------------
comp_cancer_time_since_diag_other_adj[part == "Moderate-to-vigorous physical activity", .(cancer_time_since_diag_other, estimates_contrast_healthy)]
comp_cancer_time_since_diag_other_adj[part == "Light physical activity", .(cancer_time_since_diag_other, estimates_contrast_healthy)]
comp_cancer_time_since_diag_other_adj[part == "Sedentary behaviour", .(cancer_time_since_diag_other, estimates_contrast_healthy)]
comp_cancer_time_since_diag_other_adj[part == "Sleep period", .(cancer_time_since_diag_other, estimates_contrast_healthy)]

comp_cancer_time_since_diag_other_adj[part == "Moderate-to-vigorous physical activity", .(cancer_time_since_diag_other, estimates_contrast_others)]
comp_cancer_time_since_diag_other_adj[part == "Light physical activity", .(cancer_time_since_diag_other, estimates_contrast_others)]
comp_cancer_time_since_diag_other_adj[part == "Sedentary behaviour", .(cancer_time_since_diag_other, estimates_contrast_others)]
comp_cancer_time_since_diag_other_adj[part == "Sleep period", .(cancer_time_since_diag_other, estimates_contrast_others)]

comp_cancer_time_since_diag_other_adj[part == "Moderate-to-vigorous physical activity", .(cancer_time_since_diag_other, estimates_contrast_cancer)]
comp_cancer_time_since_diag_other_adj[part == "Light physical activity", .(cancer_time_since_diag_other, estimates_contrast_cancer)]
comp_cancer_time_since_diag_other_adj[part == "Sedentary behaviour", .(cancer_time_since_diag_other, estimates_contrast_cancer)]
comp_cancer_time_since_diag_other_adj[part == "Sleep period", .(cancer_time_since_diag_other, estimates_contrast_cancer)]
