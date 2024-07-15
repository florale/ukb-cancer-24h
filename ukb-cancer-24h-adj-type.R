source("ukb-cancer-24h-utils.R")
source(paste0(redir, "ukb_utils.R"))
source("ukb-cancer-24h-data.R")

# main model --------
# fit_cancer_type_adj <- brmcoda(clr_cancer_acc,
#                                mvbind(ilr1, ilr2, ilr3) ~ cancer_before_acc_type +
#                                  age_diff_cancer_acc +
#                                  s(age) + sex + white + working + edu + never_smoked + current_drinker + deprivation,
#                                # save_pars = save_pars(all = TRUE),
#                                warmup = 500, chains = 4, cores = 4, backend = "cmdstanr",
# )
# saveRDS(fit_cancer_type_adj, paste0(outputdir, "fit_cancer_type_adj", ".RDS"))

# Predicted posteriors ------------
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

## Est means  ----------------------
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

## Contrasts --------------------
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

# All results  ------------------------
comp_cancer_type_adj <- cbind(
  pred_comp_cancer_type_adj[, .(Mean, CI_low, CI_high, part, cancer_before_acc_type)],
  diff_comp_cancer_type_adj[, .(Mean_diff, CI_low_diff, CI_high_diff)]
)

# add sig indicators
comp_cancer_type_adj[, nonsig := between(0, comp_cancer_type_adj$CI_low_diff, comp_cancer_type_adj$CI_high_diff)]
comp_cancer_type_adj[, Sig := ifelse(nonsig == FALSE, paste(intToUtf8(0x2217)), "")]

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
                                                          "Gynaecological",
                                                          "Colorectal",
                                                          "Breast",
                                                          "Genitourinary",
                                                          "Endocrine Gland",
                                                          "Other Cancer",
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

comp_cancer_type_adj[, text_position := max(CI_high), by = part]

(plot_comp_cancer_type_adj <- 
    ggplot(comp_cancer_type_adj, aes(x = cancer_before_acc_type, y = Mean, group = part)) +
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
  file = paste0(outputdir, "cancer_type_adj", ".pdf"),
  width = 11,
  height = 10,
)
plot_comp_cancer_type_adj
dev.off()





