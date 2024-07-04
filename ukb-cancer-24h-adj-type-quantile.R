
source("ukb-cancer-24h-utils.R")
source(paste0(redir, "ukb_utils.R"))
source("ukb-cancer-24h-data.R")

# main model --------
# fit_cancer_24h_type_mdn_adj <- brmcoda(clr_cancer_acc,
#                                        bf(mvbind(ilr1, ilr2, ilr3) ~ cancer_before_acc_type +
#                                             age_diff_cancer_acc +
#                                             s(age) + sex + white + working + edu + never_smoked + current_drinker + deprivation,
#                                           quantile = 0.5),
#                                        family = asym_laplace(),
#                                        # save_pars = save_pars(all = TRUE),
#                                        warmup = 500, chains = 4, cores = 4, backend = "cmdstanr",
# )
# saveRDS(fit_cancer_24h_type_mdn_adj, paste0(outputdir, "fit_cancer_24h_type_mdn_adj", ".RDS"))
# 
# fit_cancer_24h_type_low_adj <- brmcoda(clr_cancer_acc,
#                                        bf(mvbind(ilr1, ilr2, ilr3) ~ cancer_before_acc_type +
#                                             age_diff_cancer_acc +
#                                             s(age) + sex + white + working + edu + never_smoked + current_drinker + deprivation,
#                                           quantile = 0.25),
#                                        family = asym_laplace(),
#                                        # save_pars = save_pars(all = TRUE),
#                                        warmup = 500, chains = 4, cores = 4, backend = "cmdstanr",
# )
# saveRDS(fit_cancer_24h_type_low_adj, paste0(outputdir, "fit_cancer_24h_type_low_adj", ".RDS"))
# 
# fit_cancer_24h_type_high_adj <- brmcoda(clr_cancer_acc,
#                                        bf(mvbind(ilr1, ilr2, ilr3) ~ cancer_before_acc_type +
#                                             age_diff_cancer_acc +
#                                             s(age) + sex + white + working + edu + never_smoked + current_drinker + deprivation,
#                                           quantile = 0.75),
#                                        family = asym_laplace(),
#                                        # save_pars = save_pars(all = TRUE),
#                                        warmup = 500, chains = 4, cores = 4, backend = "cmdstanr",
# )
# saveRDS(fit_cancer_24h_type_high_adj, paste0(outputdir, "fit_cancer_24h_type_high_adj", ".RDS"))


# Q1 -----------------------
# Predicted posteriors ------------
fit_cancer_24h_type_low_adj <- readRDS(paste0(outputdir, "fit_cancer_24h_type_low_adj", ".RDS"))

# reference grid
d_cancer_24h_type_low_adj <- emmeans::ref_grid(fit_cancer_24h_type_low_adj$model)@grid

# predict
pred_cancer_24h_type_low_adj <- fitted(fit_cancer_24h_type_low_adj, newdata = d_cancer_24h_type_low_adj, scale = "response", summary = FALSE)

# summarise by cancer types
pred_cancer_24h_type_low_adj <- apply(pred_cancer_24h_type_low_adj, c(1), function(x)  cbind(d_cancer_24h_type_low_adj, x))
pred_cancer_24h_type_low_adj <- lapply(pred_cancer_24h_type_low_adj, function(d) {
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
})

# assemble back to summarise posteriors
pred_cancer_24h_type_low_adj <- as.data.frame(abind(pred_cancer_24h_type_low_adj, along = 1))
pred_cancer_24h_type_low_adj <- split(pred_cancer_24h_type_low_adj, pred_cancer_24h_type_low_adj$cancer_before_acc_type)

## Est means  ----------------------
pred_comp_cancer_24h_type_low_adj <- lapply(pred_cancer_24h_type_low_adj, function(l) {
  l <- as.data.frame(l[, c("sleep", "mvpa", "lpa", "sb")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, part = names(l))
  l <- rbindlist(l)
  l
})
pred_comp_cancer_24h_type_low_adj <- Map(cbind, pred_comp_cancer_24h_type_low_adj, cancer_before_acc_type = names(pred_comp_cancer_24h_type_low_adj))
pred_comp_cancer_24h_type_low_adj <- rbindlist(pred_comp_cancer_24h_type_low_adj)

## Contrasts --------------------
diff_comp_cancer_24h_type_low_adj <- lapply(pred_cancer_24h_type_low_adj, function(l) {
  l <- as.data.frame(l[, c("sleep_constrast", "mvpa_constrast", "lpa_constrast", "sb_constrast")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, constrast = names(l))
  l <- rbindlist(l)
  l
})
diff_comp_cancer_24h_type_low_adj <- Map(cbind, diff_comp_cancer_24h_type_low_adj, diff_from_healthy = names(diff_comp_cancer_24h_type_low_adj))
diff_comp_cancer_24h_type_low_adj <- rbindlist(diff_comp_cancer_24h_type_low_adj)

setnames(diff_comp_cancer_24h_type_low_adj, "Mean", "Mean_diff")
setnames(diff_comp_cancer_24h_type_low_adj, "CI_low", "CI_low_diff")
setnames(diff_comp_cancer_24h_type_low_adj, "CI_high", "CI_high_diff")

# All results  ------------------------
comp_cancer_24h_type_low_adj <- cbind(
  pred_comp_cancer_24h_type_low_adj[, .(Mean, CI_low, CI_high, part, cancer_before_acc_type)],
  diff_comp_cancer_24h_type_low_adj[, .(Mean_diff, CI_low_diff, CI_high_diff)]
)

# add sig indicators
comp_cancer_24h_type_low_adj[, nonsig := between(0, comp_cancer_24h_type_low_adj$CI_low_diff, comp_cancer_24h_type_low_adj$CI_high_diff)]
comp_cancer_24h_type_low_adj[cancer_before_acc_type != "Healthy", Sig := ifelse(nonsig == FALSE, paste(intToUtf8(0x2217)), "")]
# comp_cancer_24h_type_low_adj[cancer_before_acc_type != "Healthy", Sig := ifelse(nonsig == FALSE, paste(intToUtf8(0x2217)), paste(intToUtf8(0x2010)))]

# sort by MVPA
# comp_cancer_24h_type_low_adj[, cancer_before_acc_type := factor(cancer_before_acc_type, ordered = TRUE,
#                                                                 levels = c(
#                                                                   "Healthy",
#                                                                   "Other Skin",
#                                                                   "Prostate",
#                                                                   "Melanoma",
#                                                                   "Endocrine Gland",
#                                                                   "Breast",
#                                                                   "Genitourinary",
#                                                                   "Colorectal",
#                                                                   "Other Cancer",
#                                                                   "Gynaecological",
#                                                                   "Head & Neck",
#                                                                   "Blood",
#                                                                   "Gastrointestinal Tract",
#                                                                   "Lung",
#                                                                   "Multiple Primary"))]

comp_cancer_24h_type_low_adj[, cancer_before_acc_type := factor(cancer_before_acc_type, ordered = TRUE,
                                                                levels = c(
                                                                  "Multiple Primary",
                                                                  "Lung",
                                                                  "Gastrointestinal Tract",
                                                                  "Blood",
                                                                  "Head & Neck",
                                                                  "Gynaecological",
                                                                  "Other Cancer",
                                                                  "Colorectal",
                                                                  "Genitourinary",
                                                                  "Breast",
                                                                  "Endocrine Gland",
                                                                  "Melanoma",
                                                                  "Prostate",
                                                                  "Other Skin",
                                                                  "Healthy"
                                                                ))]

comp_cancer_24h_type_low_adj[, yintercept := NA]
comp_cancer_24h_type_low_adj[, yintercept := ifelse(part == "sleep", comp_cancer_24h_type_low_adj[cancer_before_acc_type == "Healthy" & part == "sleep"]$Mean, yintercept)]
comp_cancer_24h_type_low_adj[, yintercept := ifelse(part == "mvpa", comp_cancer_24h_type_low_adj[cancer_before_acc_type == "Healthy" & part == "mvpa"]$Mean, yintercept)]
comp_cancer_24h_type_low_adj[, yintercept := ifelse(part == "lpa", comp_cancer_24h_type_low_adj[cancer_before_acc_type == "Healthy" & part == "lpa"]$Mean, yintercept)]
comp_cancer_24h_type_low_adj[, yintercept := ifelse(part == "sb", comp_cancer_24h_type_low_adj[cancer_before_acc_type == "Healthy" & part == "sb"]$Mean, yintercept)]

comp_cancer_24h_type_low_adj[, part := ifelse(part == "sleep", "Minutes in Sleep", part)]
comp_cancer_24h_type_low_adj[, part := ifelse(part == "mvpa", "Minutes in Moderate-to-vigorous physical activity", part)]
comp_cancer_24h_type_low_adj[, part := ifelse(part == "lpa", "Minutes in Light physical activity", part)]
comp_cancer_24h_type_low_adj[, part := ifelse(part == "sb", "Minutes in Sedentary behaviour", part)]

(plot_comp_cancer_24h_type_low_adj <- 
    ggplot(comp_cancer_24h_type_low_adj, aes(x = cancer_before_acc_type, y = Mean, group = part)) +
    geom_hline(aes(yintercept = yintercept), linewidth = 0.5, linetype = 2, colour = pal_type[1]) +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_before_acc_type)) +
    geom_text(aes(label = Sig, colour = cancer_before_acc_type), 
              size = 5, nudge_x = 0.15, nudge_y = 1.25, 
              show.legend = FALSE) +
    facet_wrap(~part, scales = "free") +
    scale_colour_manual(values = pal_type) +
    labs(x = "", y = "", colour = "") +
    coord_flip() +      
    theme_ipsum() +
    theme(
      axis.ticks        = element_blank(),
      panel.background  = element_rect(fill = "transparent", colour = "black", linewidth = 0.5),
      plot.background   = element_rect(fill = "transparent", colour = NA),
      # panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      strip.text        = element_text(size = 12, hjust = .5, face = "bold"),
      legend.position   = "none"
    )
)

# Q2 and Q3 -----------------------
# Predicted posteriors ------------
fit_cancer_24h_type_mdn_adj <- readRDS(paste0(outputdir, "fit_cancer_24h_type_mdn_adj", ".RDS"))

# reference grid
d_cancer_24h_type_mdn_adj <- emmeans::ref_grid(fit_cancer_24h_type_mdn_adj$model)@grid

# predict
pred_cancer_24h_type_mdn_adj <- fitted(fit_cancer_24h_type_mdn_adj, newdata = d_cancer_24h_type_mdn_adj, scale = "response", summary = FALSE)

# summarise by cancer types
pred_cancer_24h_type_mdn_adj <- apply(pred_cancer_24h_type_mdn_adj, c(1), function(x)  cbind(d_cancer_24h_type_mdn_adj, x))
pred_cancer_24h_type_mdn_adj <- lapply(pred_cancer_24h_type_mdn_adj, function(d) {
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
})

# assemble back to summarise posteriors
pred_cancer_24h_type_mdn_adj <- as.data.frame(abind(pred_cancer_24h_type_mdn_adj, along = 1))
pred_cancer_24h_type_mdn_adj <- split(pred_cancer_24h_type_mdn_adj, pred_cancer_24h_type_mdn_adj$cancer_before_acc_type)

## Est means  ----------------------
pred_comp_cancer_24h_type_mdn_adj <- lapply(pred_cancer_24h_type_mdn_adj, function(l) {
  l <- as.data.frame(l[, c("sleep", "mvpa", "lpa", "sb")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, part = names(l))
  l <- rbindlist(l)
  l
})
pred_comp_cancer_24h_type_mdn_adj <- Map(cbind, pred_comp_cancer_24h_type_mdn_adj, cancer_before_acc_type = names(pred_comp_cancer_24h_type_mdn_adj))
pred_comp_cancer_24h_type_mdn_adj <- rbindlist(pred_comp_cancer_24h_type_mdn_adj)

## Contrasts --------------------
diff_comp_cancer_24h_type_mdn_adj <- lapply(pred_cancer_24h_type_mdn_adj, function(l) {
  l <- as.data.frame(l[, c("sleep_constrast", "mvpa_constrast", "lpa_constrast", "sb_constrast")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, constrast = names(l))
  l <- rbindlist(l)
  l
})
diff_comp_cancer_24h_type_mdn_adj <- Map(cbind, diff_comp_cancer_24h_type_mdn_adj, diff_from_healthy = names(diff_comp_cancer_24h_type_mdn_adj))
diff_comp_cancer_24h_type_mdn_adj <- rbindlist(diff_comp_cancer_24h_type_mdn_adj)

setnames(diff_comp_cancer_24h_type_mdn_adj, "Mean", "Mean_diff")
setnames(diff_comp_cancer_24h_type_mdn_adj, "CI_low", "CI_low_diff")
setnames(diff_comp_cancer_24h_type_mdn_adj, "CI_high", "CI_high_diff")

# All results Q2 and Q3  ------------------------
comp_cancer_24h_type_mdn_adj <- cbind(
  pred_comp_cancer_24h_type_mdn_adj[, .(Mean, CI_low, CI_high, part, cancer_before_acc_type)],
  diff_comp_cancer_24h_type_mdn_adj[, .(Mean_diff, CI_low_diff, CI_high_diff)]
)

# add sig indicators
comp_cancer_24h_type_mdn_adj[, nonsig := between(0, comp_cancer_24h_type_mdn_adj$CI_low_diff, comp_cancer_24h_type_mdn_adj$CI_high_diff)]
comp_cancer_24h_type_mdn_adj[cancer_before_acc_type != "Healthy", Sig := ifelse(nonsig == FALSE, paste(intToUtf8(0x2217)), "")]
# comp_cancer_24h_type_mdn_adj[cancer_before_acc_type != "Healthy", Sig := ifelse(nonsig == FALSE, paste(intToUtf8(0x2217)), paste(intToUtf8(0x2010)))]

# sort by MVPA
# comp_cancer_24h_type_mdn_adj[, cancer_before_acc_type := factor(cancer_before_acc_type, ordered = TRUE,
#                                                                 levels = c(
#                                                                   "Healthy",
#                                                                   "Other Skin",
#                                                                   "Prostate",
#                                                                   "Melanoma",
#                                                                   "Endocrine Gland",
#                                                                   "Breast",
#                                                                   "Genitourinary",
#                                                                   "Colorectal",
#                                                                   "Other Cancer",
#                                                                   "Gynaecological",
#                                                                   "Head & Neck",
#                                                                   "Blood",
#                                                                   "Gastrointestinal Tract",
#                                                                   "Lung",
#                                                                   "Multiple Primary"))]

comp_cancer_24h_type_mdn_adj[, cancer_before_acc_type := factor(cancer_before_acc_type, ordered = TRUE,
                                                                levels = c(
                                                                  "Multiple Primary",
                                                                  "Lung",
                                                                  "Gastrointestinal Tract",
                                                                  "Blood",
                                                                  "Head & Neck",
                                                                  "Gynaecological",
                                                                  "Other Cancer",
                                                                  "Colorectal",
                                                                  "Genitourinary",
                                                                  "Breast",
                                                                  "Endocrine Gland",
                                                                  "Melanoma",
                                                                  "Prostate",
                                                                  "Other Skin",
                                                                  "Healthy"
                                                                ))]

comp_cancer_24h_type_mdn_adj[, yintercept := NA]
comp_cancer_24h_type_mdn_adj[, yintercept := ifelse(part == "sleep", comp_cancer_24h_type_mdn_adj[cancer_before_acc_type == "Healthy" & part == "sleep"]$Mean, yintercept)]
comp_cancer_24h_type_mdn_adj[, yintercept := ifelse(part == "mvpa", comp_cancer_24h_type_mdn_adj[cancer_before_acc_type == "Healthy" & part == "mvpa"]$Mean, yintercept)]
comp_cancer_24h_type_mdn_adj[, yintercept := ifelse(part == "lpa", comp_cancer_24h_type_mdn_adj[cancer_before_acc_type == "Healthy" & part == "lpa"]$Mean, yintercept)]
comp_cancer_24h_type_mdn_adj[, yintercept := ifelse(part == "sb", comp_cancer_24h_type_mdn_adj[cancer_before_acc_type == "Healthy" & part == "sb"]$Mean, yintercept)]

comp_cancer_24h_type_mdn_adj[, part := ifelse(part == "sleep", "Minutes in Sleep", part)]
comp_cancer_24h_type_mdn_adj[, part := ifelse(part == "mvpa", "Minutes in Moderate-to-vigorous physical activity", part)]
comp_cancer_24h_type_mdn_adj[, part := ifelse(part == "lpa", "Minutes in Light physical activity", part)]
comp_cancer_24h_type_mdn_adj[, part := ifelse(part == "sb", "Minutes in Sedentary behaviour", part)]

(plot_comp_cancer_24h_type_mdn_adj <- 
    ggplot(comp_cancer_24h_type_mdn_adj, aes(x = cancer_before_acc_type, y = Mean, group = part)) +
    geom_hline(aes(yintercept = yintercept), linewidth = 0.5, linetype = 2, colour = pal_type[1]) +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_before_acc_type)) +
    geom_text(aes(label = Sig, colour = cancer_before_acc_type), 
              size = 5, nudge_x = 0.15, nudge_y = 1.25, 
              show.legend = FALSE) +
    facet_wrap(~part, scales = "free") +
    scale_colour_manual(values = pal_type) +
    labs(x = "", y = "", colour = "") +
    coord_flip() +      
    theme_ipsum() +
    theme(
      axis.ticks        = element_blank(),
      panel.background  = element_rect(fill = "transparent", colour = "black", linewidth = 0.5),
      plot.background   = element_rect(fill = "transparent", colour = NA),
      # panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      strip.text        = element_text(size = 12, hjust = .5, face = "bold"),
      legend.position   = "none"
    )
)

# Q4 -----------------------
# Predicted posteriors ------------
fit_cancer_24h_type_high_adj <- readRDS(paste0(outputdir, "fit_cancer_24h_type_high_adj", ".RDS"))

# reference grid
d_cancer_24h_type_high_adj <- emmeans::ref_grid(fit_cancer_24h_type_high_adj$model)@grid

# predict
pred_cancer_24h_type_high_adj <- fitted(fit_cancer_24h_type_high_adj, newdata = d_cancer_24h_type_high_adj, scale = "response", summary = FALSE)

# summarise by cancer types
pred_cancer_24h_type_high_adj <- apply(pred_cancer_24h_type_high_adj, c(1), function(x)  cbind(d_cancer_24h_type_high_adj, x))
pred_cancer_24h_type_high_adj <- lapply(pred_cancer_24h_type_high_adj, function(d) {
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
})

# assemble back to summarise posteriors
pred_cancer_24h_type_high_adj <- as.data.frame(abind(pred_cancer_24h_type_high_adj, along = 1))
pred_cancer_24h_type_high_adj <- split(pred_cancer_24h_type_high_adj, pred_cancer_24h_type_high_adj$cancer_before_acc_type)

## Est means  ----------------------
pred_comp_cancer_24h_type_high_adj <- lapply(pred_cancer_24h_type_high_adj, function(l) {
  l <- as.data.frame(l[, c("sleep", "mvpa", "lpa", "sb")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, part = names(l))
  l <- rbindlist(l)
  l
})
pred_comp_cancer_24h_type_high_adj <- Map(cbind, pred_comp_cancer_24h_type_high_adj, cancer_before_acc_type = names(pred_comp_cancer_24h_type_high_adj))
pred_comp_cancer_24h_type_high_adj <- rbindlist(pred_comp_cancer_24h_type_high_adj)

## Contrasts --------------------
diff_comp_cancer_24h_type_high_adj <- lapply(pred_cancer_24h_type_high_adj, function(l) {
  l <- as.data.frame(l[, c("sleep_constrast", "mvpa_constrast", "lpa_constrast", "sb_constrast")])
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, describe_posterior, centrality = "mean")
  l <- Map(cbind, l, constrast = names(l))
  l <- rbindlist(l)
  l
})
diff_comp_cancer_24h_type_high_adj <- Map(cbind, diff_comp_cancer_24h_type_high_adj, diff_from_healthy = names(diff_comp_cancer_24h_type_high_adj))
diff_comp_cancer_24h_type_high_adj <- rbindlist(diff_comp_cancer_24h_type_high_adj)

setnames(diff_comp_cancer_24h_type_high_adj, "Mean", "Mean_diff")
setnames(diff_comp_cancer_24h_type_high_adj, "CI_low", "CI_low_diff")
setnames(diff_comp_cancer_24h_type_high_adj, "CI_high", "CI_high_diff")

# All results Q4 ------------------------
comp_cancer_24h_type_high_adj <- cbind(
  pred_comp_cancer_24h_type_high_adj[, .(Mean, CI_low, CI_high, part, cancer_before_acc_type)],
  diff_comp_cancer_24h_type_high_adj[, .(Mean_diff, CI_low_diff, CI_high_diff)]
)

# add sig indicators
comp_cancer_24h_type_high_adj[, nonsig := between(0, comp_cancer_24h_type_high_adj$CI_low_diff, comp_cancer_24h_type_high_adj$CI_high_diff)]
comp_cancer_24h_type_high_adj[cancer_before_acc_type != "Healthy", Sig := ifelse(nonsig == FALSE, paste(intToUtf8(0x2217)), "")]
# comp_cancer_24h_type_high_adj[cancer_before_acc_type != "Healthy", Sig := ifelse(nonsig == FALSE, paste(intToUtf8(0x2217)), paste(intToUtf8(0x2010)))]

# sort by MVPA
# comp_cancer_24h_type_high_adj[, cancer_before_acc_type := factor(cancer_before_acc_type, ordered = TRUE,
#                                                                  levels = c(
#                                                                    "Healthy",
#                                                                    "Other Skin",
#                                                                    "Prostate",
#                                                                    "Melanoma",
#                                                                    "Endocrine Gland",
#                                                                    "Breast",
#                                                                    "Genitourinary",
#                                                                    "Colorectal",
#                                                                    "Other Cancer",
#                                                                    "Gynaecological",
#                                                                    "Head & Neck",
#                                                                    "Blood",
#                                                                    "Gastrointestinal Tract",
#                                                                    "Lung",
#                                                                    "Multiple Primary"))]

comp_cancer_24h_type_high_adj[, cancer_before_acc_type := factor(cancer_before_acc_type, ordered = TRUE,
                                                                 levels = c(
                                                                   "Multiple Primary",
                                                                   "Lung",
                                                                   "Gastrointestinal Tract",
                                                                   "Blood",
                                                                   "Head & Neck",
                                                                   "Gynaecological",
                                                                   "Other Cancer",
                                                                   "Colorectal",
                                                                   "Genitourinary",
                                                                   "Breast",
                                                                   "Endocrine Gland",
                                                                   "Melanoma",
                                                                   "Prostate",
                                                                   "Other Skin",
                                                                   "Healthy"
                                                                 ))]

comp_cancer_24h_type_high_adj[, yintercept := NA]
comp_cancer_24h_type_high_adj[, yintercept := ifelse(part == "sleep", comp_cancer_24h_type_high_adj[cancer_before_acc_type == "Healthy" & part == "sleep"]$Mean, yintercept)]
comp_cancer_24h_type_high_adj[, yintercept := ifelse(part == "mvpa", comp_cancer_24h_type_high_adj[cancer_before_acc_type == "Healthy" & part == "mvpa"]$Mean, yintercept)]
comp_cancer_24h_type_high_adj[, yintercept := ifelse(part == "lpa", comp_cancer_24h_type_high_adj[cancer_before_acc_type == "Healthy" & part == "lpa"]$Mean, yintercept)]
comp_cancer_24h_type_high_adj[, yintercept := ifelse(part == "sb", comp_cancer_24h_type_high_adj[cancer_before_acc_type == "Healthy" & part == "sb"]$Mean, yintercept)]

comp_cancer_24h_type_high_adj[, part := ifelse(part == "sleep", "Minutes in Sleep", part)]
comp_cancer_24h_type_high_adj[, part := ifelse(part == "mvpa", "Minutes in Moderate-to-vigorous physical activity", part)]
comp_cancer_24h_type_high_adj[, part := ifelse(part == "lpa", "Minutes in Light physical activity", part)]
comp_cancer_24h_type_high_adj[, part := ifelse(part == "sb", "Minutes in Sedentary behaviour", part)]

(plot_comp_cancer_24h_type_high_adj <- 
    ggplot(comp_cancer_24h_type_high_adj, aes(x = cancer_before_acc_type, y = Mean, group = part)) +
    geom_hline(aes(yintercept = yintercept), linewidth = 0.5, linetype = 2, colour = pal_type[1]) +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, colour = cancer_before_acc_type)) +
    geom_text(aes(label = Sig, colour = cancer_before_acc_type), 
              size = 5, nudge_x = 0.15, nudge_y = 1.25, 
              show.legend = FALSE) +
    facet_wrap(~part, scales = "free") +
    scale_colour_manual(values = pal_type) +
    labs(x = "", y = "", colour = "") +
    coord_flip() +      
    theme_ipsum() +
    theme(
      axis.ticks        = element_blank(),
      panel.background  = element_rect(fill = "transparent", colour = "black", linewidth = 0.5),
      plot.background   = element_rect(fill = "transparent", colour = NA),
      # panel.grid.major  = element_blank(),
      panel.grid.minor  = element_blank(),
      strip.text        = element_text(size = 12, hjust = .5, face = "bold"),
      legend.position   = "none"
    )
)

# ALL Qs -----------------------
comp_cancer_24h_type_low_adj[, quantile := "Q1"]
comp_cancer_24h_type_mdn_adj[, quantile := "Q2 and Q3"]
comp_cancer_24h_type_high_adj[, quantile := "Q4"]

comp_cancer_24h_type_quantile_adj <- rbind(comp_cancer_24h_type_low_adj,
                                           comp_cancer_24h_type_mdn_adj,
                                           comp_cancer_24h_type_high_adj
)
comp_cancer_24h_type_quantile_adj[, text_position := max(CI_high), by = part]

(plot_comp_cancer_24h_type_quantile_adj <- 
    ggplot(comp_cancer_24h_type_quantile_adj, aes(x = cancer_before_acc_type, y = Mean, colour = interaction(cancer_before_acc_type, quantile))) +
    geom_hline(aes(yintercept = yintercept, linetype = quantile), linewidth = 0.5, colour = "#a8a8a8") +
    # geom_point(aes(y = Mean, shape = quantile),
    #           position = position_dodge(width = 1)) +
    # geom_line( position = position_dodge(width = .75), colour = "#c1c1c1") +
    geom_pointrange(aes(ymin = CI_low,
                        ymax = CI_high, shape = quantile, linetype = quantile, size = quantile, linewidth = quantile),
                    position = position_dodge(width = .75)) +
    geom_text(aes(y = text_position + 5, label = Sig),
              size = 5, 
              position = position_dodge2(width = .75),
              show.legend = FALSE) +
    # geom_label(aes(y = text_position + 5, label = round(Mean_diff, 0)),
    #           size = 3, 
    #           # colour = "#000000",
    #           position = position_dodge2(width = .75),
    #           show.legend = FALSE) +
    # geom_label(aes(y = Mean, label = Sig)) + 
    facet_wrap(~part, scales = "free") +
    scale_colour_manual(values = pal_type_quantile) +
    scale_shape_manual(values = c(13, 16, 4)) +
    scale_size_manual(values = c(.75, .75, .75)) +
    scale_linewidth_manual(values = c(.75, .75, .75)) +
    scale_linetype_manual(values = c("dashed", "solid", "longdash")) +
    # scale_x_discrete(sec.axis = sec_axis(labels = comp_cancer_24h_type_quantile_adj$Sig)) +
    labs(x = "", y = "", colour = "") +
    coord_flip() +
    theme_ipsum() +
    theme(
      axis.ticks        = element_blank(),
      panel.background  = element_rect(fill = "transparent", colour = "black", linewidth = 0.5),
      plot.background   = element_rect(fill = "transparent", colour = NA),
      panel.grid.major.x  = element_blank(),
      panel.grid.minor  = element_blank(),
      strip.text        = element_text(size = 13, hjust = .5, face = "bold"),
      legend.position   = "none"
    )
)









