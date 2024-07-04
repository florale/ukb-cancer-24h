source("ukb-cancer-24h-utils.R")
source(paste0(redir, "ukb_utils.R"))
# source("ukb-cancer-24h-data.R")

# main model --------
# fit_cancer_time_adj <- brmcoda(clr_cancer_acc,
#                                mvbind(ilr1, ilr2, ilr3) ~ cancer_time +
#                                  s(age) + sex + white + working + edu + never_smoked + current_drinker + deprivation,
#                                # save_pars = save_pars(all = TRUE),
#                                warmup = 500, chains = 4, cores = 4, backend = "cmdstanr",
# )
# saveRDS(fit_cancer_time_adj, paste0(outputdir, "fit_cancer_time_adj", ".RDS"))

# estimates ------------
fit_cancer_time_adj <- readRDS(paste0(outputdir, "fit_cancer_time_adj", ".RDS"))

# reference grid
d_cancer_time_adj <- emmeans::ref_grid(fit_cancer_time_adj$model)@grid

# predict
pred_cancer_time_adj <- fitted(fit_cancer_time_adj, newdata = d_cancer_time_adj, scale = "response", summary = FALSE)

# summarise by cancer types
pred <- apply(pred_cancer_time_adj, c(1), function(x)  cbind(d_cancer_time_adj, x))
pred <- lapply(pred, function(d) {
  d <- as.data.table(d)
  
  # estimated means by time
  d[, sleep := mean(V1), by = cancer_time]
  d[, mvpa  := mean(V2), by = cancer_time]
  d[, lpa   := mean(V3), by = cancer_time]
  d[, sb    := mean(V4), by = cancer_time]
  
  # constrastvs healthy
  d[, sleep_vs_healthy := sleep - d$sleep[1]]
  d[, mvpa_vs_healthy := mvpa - d$mvpa[1]]
  d[, lpa_vs_healthy := lpa - d$lpa[1]]
  d[, sb_vs_healthy := sb - d$sb[1]]
  
  # constrast vs lag group
  d[, sleep_vs_lag := sleep - shift(sleep, type = "lag")]
  d[, mvpa_vs_lag := mvpa - shift(mvpa, type = "lag")]
  d[, lpa_vs_lag := lpa - shift(lpa, type = "lag")]
  d[, sb_vs_lag := sb - shift(sb, type = "lag")]
  
  d <- d[, .(cancer_time, 
             sleep, mvpa, lpa, sb,
             sleep_vs_healthy, mvpa_vs_healthy, lpa_vs_healthy, sb_vs_healthy,
             sleep_vs_lag, mvpa_vs_lag, lpa_vs_lag, sb_vs_lag
  )]
  d <- unique(d)
})

# assemble back to summarise posteriors
pred <- as.data.frame(abind(pred, along = 1))
pred <- split(pred, pred$cancer_time)

## Est means  ----------------------
pred_comp_cancer_time_adj <- lapply(pred, function(l) {
  l <- as.data.frame(l[, c("sleep", "mvpa", "lpa", "sb")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, part = names(l))
  l <- rbindlist(l)
  l
})
pred_comp_cancer_time_adj <- Map(cbind, pred_comp_cancer_time_adj, cancer_time = names(pred_comp_cancer_time_adj))
pred_comp_cancer_time_adj <- rbindlist(pred_comp_cancer_time_adj)

## Contrasts --------------------
### vs healthy
diff1_comp_cancer_time_adj <- lapply(pred, function(l) {
  l <- as.data.frame(l[, c("sleep_vs_healthy", "mvpa_vs_healthy", "lpa_vs_healthy", "sb_vs_healthy")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, constrast = names(l))
  l <- rbindlist(l)
  l
})
diff1_comp_cancer_time_adj <- Map(cbind, diff1_comp_cancer_time_adj, cancer_time = names(diff1_comp_cancer_time_adj))
diff1_comp_cancer_time_adj <- rbindlist(diff1_comp_cancer_time_adj)

setnames(diff1_comp_cancer_time_adj, "Mean", "Mean_diff_ref_healthy")
setnames(diff1_comp_cancer_time_adj, "CI_low", "CI_low_diff_ref_healthy")
setnames(diff1_comp_cancer_time_adj, "CI_high", "CI_high_diff_ref_healthy")

### vs lag
diff2_comp_cancer_time_adj <- lapply(pred, function(l) {
  l <- as.data.frame(l[, c("sleep_vs_lag", "mvpa_vs_lag", "lpa_vs_lag", "sb_vs_lag" )])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, constrast = names(l))
  l <- rbindlist(l)
  l
})
diff2_comp_cancer_time_adj <- Map(cbind, diff2_comp_cancer_time_adj, cancer_time = names(diff2_comp_cancer_time_adj))
diff2_comp_cancer_time_adj <- rbindlist(diff2_comp_cancer_time_adj)

setnames(diff2_comp_cancer_time_adj, "Mean", "Mean_diff_ref_lag")
setnames(diff2_comp_cancer_time_adj, "CI_low", "CI_low_diff_ref_lag")
setnames(diff2_comp_cancer_time_adj, "CI_high", "CI_high_diff_ref_lag")

# All results  ------------------------
comp_cancer_time_adj <- cbind(
  pred_comp_cancer_time_adj[, .(Mean, CI_low, CI_high, part, cancer_time)],
  diff1_comp_cancer_time_adj[, .(Mean_diff_ref_healthy, CI_low_diff_ref_healthy, CI_high_diff_ref_healthy)],
  diff2_comp_cancer_time_adj[, .(Mean_diff_ref_lag, CI_low_diff_ref_lag, CI_high_diff_ref_lag)]
)

# add sig indicators
comp_cancer_time_adj[, nonsig_vs_healthy := between(0, comp_cancer_time_adj$CI_low_diff_ref_healthy, comp_cancer_time_adj$CI_high_diff_ref_healthy)]
comp_cancer_time_adj[, sig_ref_healthy := ifelse(nonsig_vs_healthy == FALSE & Mean_diff_ref_healthy != 0, paste(intToUtf8(8224)), "")]

comp_cancer_time_adj[, nonsig_vs_lag := between(0, comp_cancer_time_adj$CI_low_diff_ref_lag, comp_cancer_time_adj$CI_high_diff_ref_lag)]
comp_cancer_time_adj[, sig_ref_lag := ifelse(nonsig_vs_lag == FALSE & !is.na(Mean_diff_ref_lag), paste(intToUtf8(8225)), "")] ## double dagger 


# sort by time since diagnoses
# comp_cancer_time_adj[, cancer_time := factor(cancer_time, ordered = TRUE,
#                                              levels = c(
#                                                "Healthy",
#                                                "Less than 1 year since diagnosis",
#                                                "1-5 years since diagnosis",
#                                                "More than 5 years since diagnosis"))]

# healthy in first row of plot
comp_cancer_time_adj[, cancer_time := factor(cancer_time, ordered = TRUE,
                                             levels = c(
                                               "More than 5 years since diagnosis",
                                               "1-5 years since diagnosis",
                                               "Less than 1 year since diagnosis",
                                               "Healthy"))]

comp_cancer_time_adj[, yintercept := NA]
comp_cancer_time_adj[, yintercept := ifelse(part == "sleep", comp_cancer_time_adj[cancer_time == "Healthy" & part == "sleep"]$Mean, yintercept)]
comp_cancer_time_adj[, yintercept := ifelse(part == "mvpa", comp_cancer_time_adj[cancer_time == "Healthy" & part == "mvpa"]$Mean, yintercept)]
comp_cancer_time_adj[, yintercept := ifelse(part == "lpa", comp_cancer_time_adj[cancer_time == "Healthy" & part == "lpa"]$Mean, yintercept)]
comp_cancer_time_adj[, yintercept := ifelse(part == "sb", comp_cancer_time_adj[cancer_time == "Healthy" & part == "sb"]$Mean, yintercept)]

comp_cancer_time_adj[, part := ifelse(part == "sleep", "Minutes in Sleep", part)]
comp_cancer_time_adj[, part := ifelse(part == "mvpa", "Minutes in Moderate-to-vigorous physical activity", part)]
comp_cancer_time_adj[, part := ifelse(part == "lpa", "Minutes in Light physical activity", part)]
comp_cancer_time_adj[, part := ifelse(part == "sb", "Minutes in Sedentary behaviour", part)]

comp_cancer_time_adj[, text_position := max(CI_high), by = part]

(plot_comp_cancer_time_adj <- 
    ggplot(comp_cancer_time_adj, aes(x = cancer_time, y = Mean, group = part)) +
    geom_hline(aes(yintercept = yintercept), linewidth = 0.5, linetype = 2, colour = "#a8a8a8") +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_time), size = 0.75, linewidth = 0.75) +
    geom_text(aes(y = text_position + 3, label = sig_ref_healthy, colour = cancer_time), 
              size = 4, nudge_x = 0.2, 
              show.legend = FALSE) +
    geom_text(aes(y = text_position + 4, label = sig_ref_lag, colour = cancer_time), 
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
      strip.text          = element_text(size = 12, hjust = .5, face = "bold"),
      legend.position     = "none"
    )
)
