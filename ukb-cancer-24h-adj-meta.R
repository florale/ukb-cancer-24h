library(meta)
library(metafor)
source("ukb-cancer-24h-setup.R")
source(paste0(redir, "ukb_utils.R"))

# pred_cancer_type_other_adj <- fitted(fit_cancer_type_other_adj, scale = "response")

pred_cancer_type_other_adj <- readRDS(paste0(outputdir, "pred_cancer_type_other_adj", ".RDS"))

# take estimate yi and variance vi
## https://cran.r-project.org/web/packages/metafor/metafor.pdf
# est_cancer_type_other_adj <- asplit(pred_cancer_type_other_adj, 2)$`Estimate`
# se_cancer_type_other_adj <- asplit(pred_cancer_type_other_adj, 2)$`Est.Error`
# se_cancer_type_other_adj <- sqrt(se_cancer_type_other_adj)
# 
# colnames(est_cancer_type_other_adj) <- paste0("est_", colnames(est_cancer_type_other_adj))
# colnames(se_cancer_type_other_adj) <- paste0("se_", colnames(se_cancer_type_other_adj))

# d_meta <- as.data.table(cbind(est_cancer_type_other_adj, se_cancer_type_other_adj))
# 
# # take only estimates for cancer types
# est_cancer_type_other_adj <- fixef(fit_cancer_type_other_adj)
# est_cancer_type_other_adj <- as.data.table(est_cancer_type_other_adj, keep.rownames = TRUE)
# est_cancer_type_other_adj <- est_cancer_type_other_adj[rn %in% grep("cancer_before_acc_type_other", est_cancer_type_other_adj$rn, value = TRUE)]
# est_cancer_type_other_adj <- est_cancer_type_other_adj[rn %nin% grep("Healthy|Others|MultiplePrimary|OtherCancer", est_cancer_type_other_adj$rn, value = TRUE)]
# est_cancer_type_other_adj[, cancer_type := gsub("ilr[1-3]_cancer_before_acc_type_other", "", rn)]
# names(est_cancer_type_other_adj)

# meta regression
meta_comp_cancer_type_other_adj <- lapply(pred_cancer_type_other_adj, function(l) {
  l <- as.data.frame(l[, c("sleep", "mvpa", "lpa", "sb")])
  
  # browser()
  l <- apply(l, 2, as.numeric)
  l <- apply(l, 2, function(p) posterior_summary(p, probs = c(0.005, 0.995)))
  l <- as.data.frame(l)
  l$par <- c("est", "se", "ci_low", "ci_high")
  # l <- Map(cbind, l, part = names(l))
  # l <- rbindlist(l)
  l
})
meta_comp_cancer_type_other_adj <- Map(cbind, meta_comp_cancer_type_other_adj, cancer_before_acc_type_other = names(meta_comp_cancer_type_other_adj))
meta_comp_cancer_type_other_adj <- rbindlist(meta_comp_cancer_type_other_adj)
meta_comp_cancer_type_other_adj <- meta_comp_cancer_type_other_adj[par %in% c("est", "se")]
meta_comp_cancer_type_other_adj <- meta_comp_cancer_type_other_adj[cancer_before_acc_type_other %nin% grep("Cancer|Healthy|Others|Multiple Primary|Other Cancer", meta_comp_cancer_type_other_adj$cancer_before_acc_type_other, value = TRUE)]

meta_comp_cancer_type_other_adj <- dcast(meta_comp_cancer_type_other_adj, 
                                         cancer_before_acc_type_other ~ par, 
                                         value.var = c("sleep", "mvpa", "lpa", "sb"))


# https://crukcancerintelligence.shinyapps.io/CancerStatsDataHub/
d_survival_england <- fread(paste0(inputdir, "CancerSurvival_England2025-04-24.csv"))
d_survival_scotland <- fread(paste0(inputdir, "CancerSurvival_Scotland2025-04-24.csv"))
d_survival_wales <- fread(paste0(inputdir, "CancerSurvival_Wales2025-04-24.csv"))

# check if colnames are identical across three datasets
all(colnames(d_survival_england) == colnames(d_survival_scotland) & colnames(d_survival_england) == colnames(d_survival_wales))

d_survival <- rbind(d_survival_england, d_survival_scotland, d_survival_wales)

# d_survival <- d_survival[Gender == "Persons" & measure == "Five-Year" & Standardisation == "Age-standardised"]
d_survival <- d_survival[measure == "Five-Year" & Standardisation == "Age-standardised" & CancerSite != "All Cancers Combined"]
d_survival[CancerSite == "Prostate", Gender := "Persons"]

table(d_survival$CancerSite)
table(d_survival$Age)

d_survival[, cancertype := NA]
d_survival[, cancertype := ifelse(grepl("C8[1-6]|C88|C9[0-6]", ICD10Code), "Blood", cancertype)]
d_survival[, cancertype := ifelse(grepl("C50", ICD10Code), "Breast", cancertype)]
d_survival[, cancertype := ifelse(grepl("C18|C19|C21", ICD10Code), "Colorectal", cancertype)]
d_survival[, cancertype := ifelse(grepl("C7[3-5]", ICD10Code), "Endocrine Gland", cancertype)]
d_survival[, cancertype := ifelse(grepl("C1[5-7]|C2[2-6]", ICD10Code), "Gastrointestinal Tract", cancertype)]
d_survival[, cancertype := ifelse(grepl("C62|C6[4-8]", ICD10Code), "Genitourinary", cancertype)]
d_survival[, cancertype := ifelse(grepl("C5[3-6]", ICD10Code), "Gynaecological", cancertype)]
d_survival[, cancertype := ifelse(grepl("C0[1-9]|C1[0-4]", ICD10Code), "Head & Neck", cancertype)]
d_survival[, cancertype := ifelse(grepl("C34", ICD10Code), "Lung", cancertype)]
d_survival[, cancertype := ifelse(grepl("C43", ICD10Code), "Melanoma", cancertype)]
d_survival[, cancertype := ifelse(grepl("C44", ICD10Code), "Other Skin", cancertype)]
d_survival[, cancertype := ifelse(grepl("C61", ICD10Code), "Prostate", cancertype)]
d_survival[, cancertype := as.factor(cancertype)]

d_survival <- d_survival[!is.na(cancertype)]
table(d_survival$cancertype)

d_survival <- d_survival[order(cancertype, CountryKey)]

# d_survival[, survival_rate := mean(as.numeric(NetSurvival), na.rm = TRUE), by = c("cancertype", "CountryKey")]
# d_survival[, survival_rate := survival_rate / 100]

d_survival[, NetSurvivalSiteWeighted := stats::weighted.mean(as.numeric(NetSurvival), Count, na.rm = TRUE), by = c("CancerSite")]

d_survival[, CountCountry := sum(Count, na.rm = TRUE), by = c("cancertype", "CountryKey")]

d_survival[, NetSurvivalTypeWeighted := stats::weighted.mean(as.numeric(NetSurvival), CountCountry, na.rm = TRUE), by = c("cancertype")]

d_survival_clean <- d_survival[!is.na(cancertype) & !duplicated(cancertype), .(cancertype, NetSurvivalSiteWeighted, CountCountry, NetSurvivalTypeWeighted)]

meta_comp_cancer_type_other_adj <- merge(meta_comp_cancer_type_other_adj, d_survival_clean, by.x = "cancer_before_acc_type_other", by.y = "cancertype", all.x = TRUE) 
meta_comp_cancer_type_other_adj <- meta_comp_cancer_type_other_adj[order(mvpa_est, decreasing = TRUE)]
meta_comp_cancer_type_other_adj[, survival := NetSurvivalTypeWeighted / 100]

fit_meta_cancer_sleep <- metafor::rma(yi = sleep_est, sei = sleep_se, mods = ~ survival, data = meta_comp_cancer_type_other_adj)
fit_meta_cancer_mvpa <- metafor::rma(yi = mvpa_est, sei = mvpa_se, mods = ~ survival, data = meta_comp_cancer_type_other_adj)
fit_meta_cancer_lpa <- metafor::rma(yi = lpa_est, sei = lpa_se, mods = ~ survival, data = meta_comp_cancer_type_other_adj)
fit_meta_cancer_sb <- metafor::rma(yi = sb_est, sei = sb_se, mods = ~ survival, data = meta_comp_cancer_type_other_adj)

saveRDS(fit_meta_cancer_sleep, paste0(outputdir, "fit_meta_cancer_sleep", ".RDS"))
saveRDS(fit_meta_cancer_mvpa, paste0(outputdir, "fit_meta_cancer_mvpa", ".RDS"))
saveRDS(fit_meta_cancer_lpa, paste0(outputdir, "fit_meta_cancer_mvpa", ".RDS"))
saveRDS(fit_meta_cancer_sb, paste0(outputdir, "fit_meta_cancer_mvpa", ".RDS"))

summary(fit_meta_cancer_mvpa)
summary(fit_meta_cancer_lpa)
summary(fit_meta_cancer_sb)
summary(fit_meta_cancer_sleep)

predict(fit_meta_cancer_mvpa, digits=2)

plot(fit_meta_cancer_sleep)
plot(fit_meta_cancer_mvpa)
plot(fit_meta_cancer_lpa)
plot(fit_meta_cancer_sb)

set.seed(13)
result <- list()
for(i in 1:1000){
  boot.d <- meta_comp_cancer_type_other_adj[sample(1:nrow(meta_comp_cancer_type_other_adj), nrow(meta_comp_cancer_type_other_adj), replace = TRUE),]
  boot.m.mvpa <- try(suppressWarnings(metafor::rma(yi = mvpa_est, sei = mvpa_se, mods = ~ survival, data = boot.d)), silent=TRUE)
  boot.m.lpa <- try(suppressWarnings(metafor::rma(yi = lpa_est, sei = lpa_se, mods = ~ survival, data = boot.d)), silent=TRUE)
  boot.m.sb <- try(suppressWarnings(metafor::rma(yi = sb_est, sei = sb_se, mods = ~ survival, data = boot.d)), silent=TRUE)
  boot.m.sleep <- try(suppressWarnings(metafor::rma(yi = sleep_est, sei = sleep_se, mods = ~ survival, data = boot.d)), silent=TRUE)
  
  if (inherits(boot.m.mvpa, "try-error") | 
      inherits(boot.m.lpa, "try-error") | 
      inherits(boot.m.sb, "try-error") | 
      inherits(boot.m.sleep, "try-error")) {
    result[[i]] <- c(NA, NA, NA, NA)
  } else {
    result[[i]] <- c(boot.m.mvpa$R2, boot.m.lpa$R2, boot.m.sb$R2, boot.m.sleep$R2)
  }
}
result <- do.call(rbind, result)
result <- result[complete.cases(result), ]
colnames(result) <- c("mvpa", "lpa", "sb", "sleep")

# apply quantile to each column
apply(result, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
apply(result, 2, mean)

