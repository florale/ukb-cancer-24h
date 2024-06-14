source("ukb-cancer-24h-data.R")

# cancer vs healthy ------------------------
fit_cancer <- brmcoda(clr_cancer_acc,
                      bf(mvbind(ilr1, ilr2, ilr3) ~ cancer_before_acc +
                           age + sex + white + working + edu + never_smoked + current_drinker + deprivation,
                         quantile = 0.5),
                      # save_pars = save_pars(all = TRUE),
                      warmup = 500, chains = 4, cores = 4, backend = "cmdstanr",
                      family = asym_laplace()
)
saveRDS(fit_cancer, paste0(outputdir, "fit_cancer", ".RDS"))
fit_cancer <- readRDS(paste0(outputdir, "fit_cancer", ".RDS"))

# estimates
d_cancer <- build.rg(fit_cancer)
d_cancer <- unique(d_cancer)
pred_cancer_24h <- fitted(fit_cancer, newdata = d_cancer, scale = "response")

comp_cancer <- as.data.table(pred_cancer_24h)

comp_cancer[, cancer_before_acc := NA]
comp_cancer[, cancer_before_acc := ifelse(V1 == 1, "Healthy", cancer_before_acc)]
comp_cancer[, cancer_before_acc := ifelse(V1 == 2, "Cancer", cancer_before_acc)]
comp_cancer <- comp_cancer[!is.na(cancer_before_acc)]
comp_cancer <- dcast(comp_cancer, cancer_before_acc + V3 ~ V2, value.var = "value")

ggplot(comp_cancer, aes(x = cancer_before_acc, y = Estimate, group = V3)) +
  geom_pointrange(aes(ymin = Q2.5,
                      ymax = Q97.5, colour = V3)) +
  facet_wrap(~V3, scales = "free") +
  scale_colour_jco() +
  scale_fill_jco() +
  coord_flip() +
  theme_bw()

# cancer time ------------------------
fit_cancer_time <- brmcoda(clr_cancer_acc,
                           mvbind(ilr1, ilr2, ilr3) ~ cancer_time +
                             age + sex + white + working + edu + never_smoked + current_drinker + deprivation,
                           # save_pars = save_pars(all = TRUE),
                           warmup = 500, chains = 4, cores = 4, backend = "cmdstanr",
)
saveRDS(fit_cancer_time, paste0(outputdir, "fit_cancer_time", ".RDS"))
fit_cancer_time <- readRDS(paste0(outputdir, "fit_cancer_time", ".RDS"))

# estimates
d_cancer_time <- model.frame(fit_cancer_time)
d_cancer_time <- insight::get_datagrid(d_cancer_time)

d_cancer_time <- unique(d_cancer_time)
pred_cancer_time <- fitted(fit_cancer_time, newdata = d_cancer_time, scale = "response")
pred_comp_cancer_time <- as.data.table(pred_cancer_time)

pred_comp_cancer_time[, cancer_time := NA]
pred_comp_cancer_time[, cancer_time := ifelse(V1 == 1, "Healthy", cancer_time)]
pred_comp_cancer_time[, cancer_time := ifelse(V1 == 2, "1-5 years since cancer diagnosis", cancer_time)]
pred_comp_cancer_time[, cancer_time := ifelse(V1 == 3, "More than 5 years cancer since diagnosis", cancer_time)]
pred_comp_cancer_time[, cancer_time := factor(cancer_time, ordered = TRUE,
                                              levels = c("Healthy",
                                                         "1-5 years since cancer diagnosis",
                                                         "More than 5 years cancer since diagnosis"))]

pred_comp_cancer_time <- pred_comp_cancer_time[!is.na(cancer_time)]

pred_comp_cancer_time <- dcast(pred_comp_cancer_time, cancer_time + V3 ~ V2, value.var = "value")

ggplot(pred_comp_cancer_time, aes(x = cancer_time, y = Estimate, group = V3)) +
  geom_pointrange(aes(ymin = Q2.5,
                      ymax = Q97.5, colour = V3)) +
  facet_wrap(~V3, scales = "free") +
  scale_colour_jco() +
  scale_fill_jco() +
  coord_flip() +      
  theme_ipsum() +
  theme(
    axis.ticks        = element_blank(),
    legend.position   = "bottom",
    panel.background  = element_rect(
      fill = "transparent",
      colour = "black",
      linewidth = 0.5
    ))

# estimated differences between group
pred_cancer_24h_draws <- fitted(fit_quantile_cancer_time, newdata = d0, scale = "response", summary = FALSE)

# cancer time 5g ------------------------
fit_cancer_time_5g <- brmcoda(clr_cancer_acc,
                              mvbind(ilr1, ilr2, ilr3) ~ cancer_time_5g +
                                age + sex + white + working + edu + never_smoked + current_drinker + deprivation,
                              # save_pars = save_pars(all = TRUE),
                              warmup = 500, chains = 4, cores = 4, backend = "cmdstanr",
)
saveRDS(fit_cancer_time_5g, paste0(outputdir, "fit_cancer_time_5g", ".RDS"))
fit_cancer_time_5g <- readRDS(paste0(outputdir, "fit_cancer_time_5g", ".RDS"))

# estimates
d_cancer_time_5g <- build.rg(fit_cancer_time_5g)
d_cancer_time_5g <- unique(d_cancer_time_5g)
pred_cancer_time_5g <- fitted(fit_cancer_time_5g, newdata = d_cancer_time_5g, scale = "response")
pred_comp_cancer_time_5g <- as.data.table(pred_cancer_time_5g)

pred_comp_cancer_time_5g[, cancer_time_5g := NA]
pred_comp_cancer_time_5g[, cancer_time_5g := ifelse(V1 == 1, "Healthy", cancer_time_5g)]
pred_comp_cancer_time_5g[, cancer_time_5g := ifelse(V1 == 2, "Less than 1", cancer_time_5g)]
pred_comp_cancer_time_5g[, cancer_time_5g := ifelse(V1 == 3, "1-5 years since cancer diagnosis", cancer_time_5g)]
pred_comp_cancer_time_5g[, cancer_time_5g := ifelse(V1 == 4, "5-10 years since cancer diagnosis", cancer_time_5g)]
pred_comp_cancer_time_5g[, cancer_time_5g := ifelse(V1 == 5, "More than 10 years cancer since diagnosis", cancer_time_5g)]
pred_comp_cancer_time_5g[, cancer_time_5g := factor(cancer_time_5g, ordered = TRUE,
                                                    levels = c("Healthy",
                                                               "Less than 1",
                                                               "1-5 years since cancer diagnosis",
                                                               "5-10 years since cancer diagnosis",
                                                               "More than 10 years cancer since diagnosis"))]

pred_comp_cancer_time_5g <- pred_comp_cancer_time_5g[!is.na(cancer_time_5g)]
pred_comp_cancer_time_5g <- dcast(pred_comp_cancer_time_5g, cancer_time_5g + V3 ~ V2, value.var = "value")

ggplot(pred_comp_cancer_time_5g, aes(x = cancer_time_5g, y = Estimate, group = V3)) +
  geom_pointrange(aes(ymin = Q2.5,
                      ymax = Q97.5, colour = V3)) +
  facet_wrap(~V3, scales = "free") +
  scale_colour_jco() +
  scale_fill_jco() +
  coord_flip() +      
  theme_ipsum() +
  theme(
    axis.ticks        = element_blank(),
    legend.position   = "bottom",
    panel.background  = element_rect(
      fill = "transparent",
      colour = "black",
      linewidth = 0.5
    ))

# cancer type ------------------------
fit_cancer_24h_type <- brmcoda(clr_cancer_acc,
                               mvbind(ilr1, ilr2, ilr3) ~ cancer_before_acc_type +
                                 sex + white + working + edu + never_smoked + current_drinker + deprivation,
                               # save_pars = save_pars(all = TRUE),
                               warmup = 500, chains = 4, cores = 4, backend = "cmdstanr",
)
saveRDS(fit_cancer_24h_type, paste0(outputdir, "fit_cancer_24h_type", ".RDS"))
fit_cancer_24h_type <- readRDS(paste0(outputdir, "fit_cancer_24h_type", ".RDS"))

# estimates
d_cancer_24h_type <- build.rg(fit_cancer_24h_type)
d_cancer_24h_type <- unique(d_cancer_24h_type)
pred_cancer_24h_type <- fitted(fit_cancer_24h_type, newdata = d_cancer_24h_type, scale = "response")

pred_comp_cancer_24h_type <- as.data.table(pred_cancer_24h_type)

pred_comp_cancer_24h_type[, cancer_before_acc_type := NA]
pred_comp_cancer_24h_type[, cancer_before_acc_type := ifelse(V1 == 1, "Healthy", cancer_before_acc_type)]
pred_comp_cancer_24h_type[, cancer_before_acc_type := ifelse(V1 == 2, "Blood", cancer_before_acc_type)]
pred_comp_cancer_24h_type[, cancer_before_acc_type := ifelse(V1 == 3, "Breast", cancer_before_acc_type)]
pred_comp_cancer_24h_type[, cancer_before_acc_type := ifelse(V1 == 4, "Colorectal", cancer_before_acc_type)]
pred_comp_cancer_24h_type[, cancer_before_acc_type := ifelse(V1 == 5, "Endocrine Gland", cancer_before_acc_type)]
pred_comp_cancer_24h_type[, cancer_before_acc_type := ifelse(V1 == 6, "Gastrointestinal Tract", cancer_before_acc_type)]
pred_comp_cancer_24h_type[, cancer_before_acc_type := ifelse(V1 == 7, "Genitourinary", cancer_before_acc_type)]
pred_comp_cancer_24h_type[, cancer_before_acc_type := ifelse(V1 == 8, "Gynaecological", cancer_before_acc_type)]
pred_comp_cancer_24h_type[, cancer_before_acc_type := ifelse(V1 == 9, "Head & Neck", cancer_before_acc_type)]
pred_comp_cancer_24h_type[, cancer_before_acc_type := ifelse(V1 == 10, "Lung", cancer_before_acc_type)]
pred_comp_cancer_24h_type[, cancer_before_acc_type := ifelse(V1 == 11, "Melanoma", cancer_before_acc_type)]
pred_comp_cancer_24h_type[, cancer_before_acc_type := ifelse(V1 == 12, "Multiple Primary", cancer_before_acc_type)]
pred_comp_cancer_24h_type[, cancer_before_acc_type := ifelse(V1 == 13, "Other Cancer", cancer_before_acc_type)]
pred_comp_cancer_24h_type[, cancer_before_acc_type := ifelse(V1 == 14, "Other Skin", cancer_before_acc_type)]
pred_comp_cancer_24h_type[, cancer_before_acc_type := ifelse(V1 == 15, "Prostate", cancer_before_acc_type)]

pred_comp_cancer_24h_type <- pred_comp_cancer_24h_type[!is.na(cancer_before_acc_type)]
pred_comp_cancer_24h_type <- dcast(pred_comp_cancer_24h_type, cancer_before_acc_type + V3 ~ V2, value.var = "value")

ggplot(pred_comp_cancer_24h_type, aes(x = cancer_before_acc_type, y = Estimate, group = V3)) +
  geom_pointrange(aes(ymin = Q2.5,
                      ymax = Q97.5, colour = V3)) +
  facet_wrap(~V3, scales = "free") +
  scale_colour_jco() +
  scale_fill_jco() +
  coord_flip() +      
  theme_ipsum() +
  theme(
    axis.ticks        = element_blank(),
    legend.position   = "bottom",
    panel.background  = element_rect(
      fill = "transparent",
      colour = "black",
      linewidth = 0.5
    ))

# cancer by time since diagnosis ------------------------
fit_cancer_24h_timex <- brmcoda(clr_cancer_acc,
                                mvbind(ilr1, ilr2, ilr3) ~ cancer_before_acc*age_diff_cancer_acc +
                                  sex + white + working + edu + never_smoked + current_drinker + deprivation,
                                # save_pars = save_pars(all = TRUE),
                                warmup = 500, chains = 4, cores = 4, backend = "cmdstanr",
)
saveRDS(fit_cancer_24h_timex, paste0(outputdir, "fit_cancer_24h_timex", ".RDS"))
fit_cancer_24h_timex <- readRDS(paste0(outputdir, "fit_cancer_24h_timex", ".RDS"))

# estimates
age_diff_cancer_acc <- data.table(
  age_diff_cancer_acc = c(0:10))

# d_cancer_24h_timex <- insight::get_datagrid(model.frame(fit_cancer_24h_timex), factors = "mode")
d_cancer_24h_timex <- build.rg(fit_cancer_24h_timex, ref = "grandmean", level = "aggregate", factors = "mode", length = NA)
d_cancer_24h_timex <- model.frame(fit_cancer_24h_timex)
d_cancer_24h_timex <- insight::get_datagrid(d_cancer_24h_timex, 
                                            at = "cancer_before_acc",
                                            factors = "mode",
                                            length = NA)

d_cancer_24h_timex <- expand.grid.df(as.data.table(d_cancer_24h_timex)[, -c("age_diff_cancer_acc")], age_diff_cancer_acc)
d_cancer_24h_timex <- unique(d_cancer_24h_timex)
pred_cancer_24h_timex <- fitted(fit_cancer_24h_timex, newdata = d_cancer_24h_timex, scale = "response")

pred_comp_cancer_24h_timex <- as.data.table(pred_cancer_24h_timex)
pred_comp_cancer_24h_timex <- dcast(pred_comp_cancer_24h_timex, V3 + V1 ~ V2, value.var = "value")
pred_comp_cancer_24h_timex <- cbind(pred_comp_cancer_24h_timex, d_cancer_24h_timex)
pred_comp_cancer_24h_timex <- pred_comp_cancer_24h_timex[!is.na(cancer_before_acc)]

pred_comp_cancer_24h_timex[, cancer_before_acc := ifelse(cancer_before_acc == 0, "Healthy", "Cancer")]

ggplot(pred_comp_cancer_24h_timex, aes(x = age_diff_cancer_acc, y = Estimate, group = cancer_before_acc)) +
  geom_smooth()+
  # geom_ribbon(aes(ymin = Q2.5,
  #                 ymax = Q97.5),
  #             alpha = 1/10, show.legend = TRUE) +
  facet_wrap(~V3, scales = "free") +
  theme_bw()
