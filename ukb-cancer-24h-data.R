# olivia's
# d_acc_icd <- readRDS("~/Desktop/GitHub/ukb-cancer-24h/d_acc_icd.RDS")
# source("/Users/oliviadelia/Desktop/GitHub/ukbiobank/ukb_utils.R")
# source("ukb-cancer-24h-utils.R")
# d_acc_icd <- readRDS(paste0("d_acc_icd", ".RDS"))

# flora's
source("ukb-cancer-24h-utils.R")
source(paste0(redir, "ukb_utils.R"))
d_acc_icd <- readRDS(paste0(inputdir, "d_acc_icd", ".RDS"))

# var names --------------
icd_not_cancer_any_vars <- c(
  "icd_i_any", "icd_ix_any", "icd_v_any",
  "icd_vii_any", "icd_viii_any", "icd_i_iii_any",
  "icd_xi_any", "icd_xiv_any", "icd_iv_any", "icd_xiii_any", "icd_vi_any",
  "icd_x_any"
)
icd_ii_fo_vars <- grep("icd_ii_.*_fo", names(d_acc_icd), value = TRUE)
icd_ii_gynaecological_vars <- c("icd_ii_c53", "icd_ii_c56", "icd_ii_c54_c55")
icd_ii_genitourinary_vars <- c("icd_ii_c64", "icd_ii_c65_c66","icd_ii_c67", "icd_ii_c68", "icd_ii_c62")
icd_ii_headneck_vars <- c("icd_ii_c01_c14", "icd_ii_c31_c33")
icd_ii_lung_vars <- c("icd_ii_c34")
icd_ii_colorectal_vars <- c("icd_ii_c18_c21")
icd_ii_gastrointestinal_vars <- c("icd_ii_c15", "icd_ii_c16", "icd_ii_c17", "icd_ii_c22", "icd_ii_c23_c24", "icd_ii_c25", "icd_ii_c26")
icd_ii_breast_vars <- c("icd_ii_c50")
icd_ii_prostate_vars <- c("icd_ii_c61")
icd_ii_blood_vars <- c("icd_ii_lymphoid")
icd_ii_melanoma_vars <- c("icd_ii_c43")
icd_ii_skin_vars <- c("icd_ii_c44")
icd_ii_endocrine_vars <- c("icd_ii_c73_c75")
icd_ii_other_vars <- c("icd_ii_c30_c39", "icd_ii_c40_c41", "icd_ii_various", "icd_ii_c45", "icd_ii_c69_c72")
icd_ii_subtype_vars <- c(
  "icd_ii_c53", "icd_ii_c56", "icd_ii_c62", "icd_ii_c01_c14", "icd_ii_c15", 
  "icd_ii_c16", "icd_ii_c17", "icd_ii_c18_c21", "icd_ii_c22", "icd_ii_c23_c24", 
  "icd_ii_c25", "icd_ii_c26", "icd_ii_c31_c33", "icd_ii_c34", "icd_ii_c30_c39", 
  "icd_ii_c40_c41", "icd_ii_c43", "icd_ii_c44", "icd_ii_c45", "icd_ii_various", "icd_ii_c50", 
  "icd_ii_c54_c55", "icd_ii_c61", "icd_ii_c64", "icd_ii_c65_c66", "icd_ii_c67", 
  "icd_ii_c68", "icd_ii_c69_c72", "icd_ii_c73_c75", "icd_ii_lymphoid" 
)
icd_any_vars <- c(
  "icd_i_any", "icd_ii_any", "icd_ix_any", "icd_v_any",
  "icd_vii_any", "icd_viii_any", "icd_i_iii_any",
  "icd_xi_any", "icd_xiv_any", "icd_iv_any", "icd_xiii_any", "icd_vi_any",
  "icd_x_any"
)

# cancer data management  ---------------
# notes for Flora: 13 cancer composite related to PA
# (bladder, breast, colon, endometrial, oesophageal adenocarcinoma, gastric cardia, head and neck, kidney, liver, lung, myeloid leukaemia, myeloma, and rectum)
d_acc_icd[, icd_not_cancer := ifelse((Reduce(`|`, lapply(icd_not_cancer_any_vars, function(v) f1(get(v), 1)))),
                                     1, 0)]

d_acc_icd[, icd_any := ifelse((Reduce(`|`, lapply(icd_any_vars, function(v) f1(get(v), 1)))),
                              1, 0)]

table(d_acc_icd$icd_not_cancer, useNA = "always")
table(d_acc_icd$icd_ii_any, useNA = "always")
table(d_acc_icd$icd_ii_group, useNA = "always")
table(d_acc_icd$icd_any, useNA = "always")

d_acc_icd[, icd_ii_comorbid := NA]
d_acc_icd[, icd_ii_comorbid := ifelse(icd_ii_group %in% c("One Primary", "Multiple Primary") & icd_not_cancer == 0, "Only cancer", icd_ii_comorbid)]
d_acc_icd[, icd_ii_comorbid := ifelse(icd_ii_group %in% c("One Primary", "Multiple Primary") & icd_not_cancer == 1, "Cancer comorbid", icd_ii_comorbid)]
d_acc_icd[, icd_ii_comorbid := ifelse(icd_ii_group == "No" & icd_not_cancer == 1, "Only other conditions", icd_ii_comorbid)]
d_acc_icd[, icd_ii_comorbid := ifelse(icd_ii_group == "No" & icd_not_cancer == 0, "Healthy", icd_ii_comorbid)]

table(d_acc_icd$icd_ii_comorbid, useNA = "always")

d_acc_icd[, health_condition := NA]
d_acc_icd[, health_condition := ifelse(icd_ii_group %in% c("One Primary", "Multiple Primary"), "Cancer", health_condition)]
d_acc_icd[, health_condition := ifelse(icd_ii_group == "No" & icd_not_cancer == 1, "Only other conditions", health_condition)]
d_acc_icd[, health_condition := ifelse(icd_ii_group == "No" & icd_not_cancer == 0, "Healthy", health_condition)]

d_acc_icd[, health_condition := as.factor(health_condition)]
d_acc_icd[, health_condition := relevel(health_condition, ref = "Healthy")]

table(d_acc_icd$health_condition, useNA = "always")

# censor 1 years to exclude from healthy --------
nrow(d_acc_icd[!is.na(icd_ii_sub_lo)])
nrow(d_acc_icd[!is.na(icd_not_ii_sub_lo)])

# time since first cancer diagnoses
## negative value means acc is before diagnosis
# d_acc_icd[, age_diff_cancer_acc := year(acc_startdate) - year(icd_ii_sub_fo)]
d_acc_icd[, age_diff_cancer_acc := (acc_startdate - icd_ii_sub_fo)/365.25]
# table(d_acc_icd$age_diff_cancer_acc, useNA = "always")
table(round(d_acc_icd$age_diff_cancer_acc), useNA = "always")

# time since most recent other diagnoses
# d_acc_icd[, age_diff_other_cond_acc := year(acc_startdate) - year(icd_not_ii_sub_fo)]
d_acc_icd[, age_diff_other_cond_acc := (acc_startdate - icd_not_ii_sub_fo)/365.25]
table(round(d_acc_icd$age_diff_other_cond_acc), useNA = "always")
# table(round(as.integer(as.Date("2015-03-02") - as.Date("1938-03-16"))/365.25))

# number of cancer diags after acc
nrow(d_acc_icd[((acc_startdate - icd_ii_sub_fo)/365.25) <= 0])

# number of cancer diags within 1y
nrow(d_acc_icd[((acc_startdate - icd_ii_sub_fo)/365.25) %gele% c(-1, 0)])

nrow(d_acc_icd[((acc_startdate - icd_ii_sub_fo)/365.25) < -1])

# number of other diags 1y after acc
nrow(d_acc_icd[((acc_startdate - icd_not_ii_sub_fo)/365.25) < -1])

# number of other diags before acc
nrow(d_acc_icd[((acc_startdate - icd_not_ii_sub_fo)/365.25) >= 0])

# number of other diags within 1y
nrow(d_acc_icd[((acc_startdate - icd_not_ii_sub_fo)/365.25) %gele% c(-1, 0)])

# number of other diags up to 1y after acc
nrow(d_acc_icd[((acc_startdate - icd_not_ii_sub_fo)/365.25) >= -1])


# time since most recent any diagnoses
d_acc_icd[, time_diff_acc_icd_any := (acc_startdate - icd_sub_fo)/365.25]
d_acc_icd[time_diff_acc_icd_any %gele% c(-1, 0), time_diff_acc_icd_any := NA]
d_acc_icd[, time_diff_acc_icd_any := ifelse(icd_any == 0 | time_diff_acc_icd_any < - 1, 0, time_diff_acc_icd_any)]

# as.integer(as.Date("2015-03-02") - as.Date("1938-03-16"))/365.25
# table(round(as.integer(as.Date("2015-03-02") - as.Date("1938-03-16"))/365.25))

# censor 1 years to exclude from healthy - set to NA to exclude in time since diagnoses
# d_acc_icd <- d_acc_icd[age_diff_cancer_acc %gele% c(-1, 0) | age_diff_other_cond_acc %gele% c(-1, 0)]

d_acc_icd[age_diff_cancer_acc %gele% c(-1, 0), age_diff_cancer_acc := NA]
d_acc_icd[age_diff_other_cond_acc %gele% c(-1, 0), age_diff_other_cond_acc := NA]

# d_acc_icd[between(age_diff_cancer_acc, -1, 0, incbounds = TRUE), age_diff_cancer_acc := NA]
# d_acc_icd[between(age_diff_other_cond_acc, -1, 0, incbounds = TRUE), age_diff_cancer_acc := NA]

# add cancer diag after 1y and always healthy  to healthy
d_acc_icd[, age_diff_cancer_acc := ifelse(age_diff_cancer_acc < -1 | icd_any == 0, 0, age_diff_cancer_acc)]

# add and any conds after 1y since acc to healthy (0)
d_acc_icd[, age_diff_cancer_acc := ifelse(time_diff_acc_icd_any == 0 & (age_diff_cancer_acc == 0 | is.na(age_diff_cancer_acc)), 0, age_diff_cancer_acc)]

# remove the ones with health diag up to 1y after acc from healthy
d_acc_icd[, age_diff_cancer_acc := ifelse(age_diff_other_cond_acc > - 1 & age_diff_cancer_acc == 0 & icd_not_cancer == 1, NA, age_diff_cancer_acc)]

table(round(d_acc_icd$age_diff_cancer_acc), useNA = "always")
nrow(d_acc_icd[age_diff_cancer_acc == 0])


# subset -------------
d_acc_icd <- d_acc_icd[!is.na(age_diff_cancer_acc)]

# main cancer predictor variables ----------------
# time since cancer diagnosis - healthy should be 22599
d_acc_icd[, cancer_time_since_diag := cut(age_diff_cancer_acc,
                                          breaks = c(-Inf, 0, 1, 5, Inf),
                                          labels = c("Healthy",
                                                     "Less than 1 year since diagnosis",
                                                     "1-5 years since diagnosis",
                                                     "More than 5 years since diagnosis"))]
d_acc_icd[, cancer_time_since_diag := as.factor(cancer_time_since_diag)]
table(d_acc_icd$cancer_time_since_diag, useNA = "always")

d_acc_icd[, cancer_time_since_diag_5g := cut(age_diff_cancer_acc, 
                                  breaks = c(-Inf, 0, 1, 5, 10, Inf),
                                  labels = c("Healthy",
                                             "Less than 1 year since diagnosis",
                                             "1-5 years since diagnosis", 
                                             "5-10 years since diagnosis", 
                                             "More than 10 years since diagnosis"))]
d_acc_icd[, cancer_time_since_diag_5g := as.factor(cancer_time_since_diag_5g)]
table(d_acc_icd$cancer_time_since_diag_5g, useNA = "always")

# factor cancer vs healthy
d_acc_icd[, cancer_before_acc := NA]
d_acc_icd[, cancer_before_acc := ifelse(cancer_time_since_diag == "Healthy", 0, cancer_before_acc)]
d_acc_icd[, cancer_before_acc := ifelse(cancer_time_since_diag != "Healthy", 1, cancer_before_acc)]
d_acc_icd[, cancer_before_acc := factor(cancer_before_acc, levels = c(0, 1), labels = c("Healthy", "Cancer"))]
table(d_acc_icd$cancer_before_acc, useNA = "always")

# type
d_acc_icd[, cancer_before_acc_type := NA]
d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_gynaecological_vars, function(v) f1(get(v), 1)))) &
                                               (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_gynaecological_vars, function(v) f1(get(v), 0)))) &
                                               cancer_before_acc == "Cancer",
                                             "Gynaecological", cancer_before_acc_type
)]
d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_genitourinary_vars, function(v) f1(get(v), 1)))) &
                                               (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_genitourinary_vars, function(v) f1(get(v), 0)))) &
                                               cancer_before_acc == "Cancer",
                                             "Genitourinary", cancer_before_acc_type
)]

d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_headneck_vars, function(v) f1(get(v), 1)))) &
                                               (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_headneck_vars, function(v) f1(get(v), 0))))  &
                                               cancer_before_acc == "Cancer",
                                             "Head & Neck", cancer_before_acc_type
)]

d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_lung_vars, function(v) f1(get(v), 1)))) &
                                               (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_lung_vars, function(v) f1(get(v), 0)))) &
                                               cancer_before_acc == "Cancer",
                                             "Lung", cancer_before_acc_type
)]

d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_colorectal_vars, function(v) f1(get(v), 1)))) &
                                               (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_colorectal_vars, function(v) f1(get(v), 0)))) &
                                               cancer_before_acc == "Cancer",
                                             "Colorectal", cancer_before_acc_type
)]

d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_gastrointestinal_vars, function(v) f1(get(v), 1)))) &
                                               (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_gastrointestinal_vars, function(v) f1(get(v), 0)))) &
                                               cancer_before_acc == "Cancer",
                                             "Gastrointestinal Tract", cancer_before_acc_type
)]

d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_breast_vars, function(v) f1(get(v), 1)))) &
                                               (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_breast_vars, function(v) f1(get(v), 0)))) &
                                               cancer_before_acc == "Cancer",
                                             "Breast", cancer_before_acc_type
)]

d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_prostate_vars, function(v) f1(get(v), 1)))) &
                                               (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_prostate_vars, function(v) f1(get(v), 0)))) &
                                               cancer_before_acc == "Cancer",
                                             "Prostate", cancer_before_acc_type
)]
d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_blood_vars, function(v) f1(get(v), 1)))) &
                                               (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_blood_vars, function(v) f1(get(v), 0)))) &
                                               cancer_before_acc == "Cancer",
                                             "Blood", cancer_before_acc_type
)]

d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_melanoma_vars, function(v) f1(get(v), 1)))) &
                                               (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_melanoma_vars, function(v) f1(get(v), 0)))) &
                                               cancer_before_acc == "Cancer",
                                             "Melanoma", cancer_before_acc_type
)]

d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_skin_vars, function(v) f1(get(v), 1)))) &
                                               (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_skin_vars, function(v) f1(get(v), 0)))) &
                                               cancer_before_acc == "Cancer",
                                             "Other Skin", cancer_before_acc_type
)]

# d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_cns_vars, function(v) f1(get(v), 1)))) &
#                                                (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_cns_vars, function(v) f1(get(v), 0)))) &
#                                                cancer_before_acc == "Cancer",
#                                              "Central Nervous System", cancer_before_acc_type
# )]

d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_endocrine_vars, function(v) f1(get(v), 1)))) &
                                               (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_endocrine_vars, function(v) f1(get(v), 0)))) &
                                               cancer_before_acc == "Cancer",
                                             "Endocrine Gland", cancer_before_acc_type
)]

d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_other_vars, function(v) f1(get(v), 1)))) &
                                               (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_other_vars, function(v) f1(get(v), 0)))) &
                                               cancer_before_acc == "Cancer",
                                             "Other Cancer", cancer_before_acc_type
)]

d_acc_icd[, cancer_before_acc_type := ifelse(icd_ii_type_n > 1 & cancer_before_acc == "Cancer", "Multiple Primary", cancer_before_acc_type)]
d_acc_icd[, cancer_before_acc_type := ifelse(cancer_before_acc == "Healthy", "Healthy", cancer_before_acc_type)]

d_acc_icd[, cancer_before_acc_type := as.factor(cancer_before_acc_type)]
d_acc_icd[, cancer_before_acc_type := relevel(cancer_before_acc_type, ref = "Healthy")]

table(d_acc_icd$cancer_before_acc_type, useNA = "always")

# check multiple primary
nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_gynaecological == 1])
nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_genitourinary == 1])
nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_headneck == 1])
nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_lung == 1])
nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_colorectal == 1])
nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_gastrointestinal == 1])
nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_breast == 1]) # 2
nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_prostate == 1])
nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_blood == 1])
nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_melanoma == 1])
nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_skin == 1]) # 1
nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_endocrine == 1])
nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_other == 1])

# exclude those had cancer after acc and other health conditions -------------------
d_cancer_acc <- d_acc_icd[!is.na(cancer_before_acc)]

table(d_cancer_acc$cancer_before_acc, useNA = "always")
table(d_cancer_acc$cancer_before_acc_type, useNA = "always")
table(d_cancer_acc$cancer_time_since_diag, useNA = "always")
table(round(d_cancer_acc$age_diff_cancer_acc), useNA = "always")

d_cancer_acc[, age_at_acc := year(acc_startdate) - year_birth]

# Composition and ilr ----------------------
# d_cancer_acc <- readRDS(paste0(outputdir, "d_cancer_acc", ".RDS"))
d_cancer_acc <- d_cancer_acc[, -grep("ilr", names(d_cancer_acc), value = TRUE), with = FALSE]

sbp <- matrix(c(
  1, -1, -1,-1,
  0, 1, -1, -1,
  0, 0, 1, -1), ncol = 4, byrow = TRUE)

clr_cancer_acc <- complr(data = d_cancer_acc,
                         transform = "ilr",
                         parts = c("sleep_comp", "mvpa_comp", "lpa_comp", "sb_comp"),
                         sbp = sbp,
                         total = 1440)

# descriptives ----------------------------
## demographics - group by cancer vs healthy
egltable(c(
  "age", "sex", "ethnicg", "bmig", "edu", "working", 
  "smoking", "alcohol", "deprivation"
  ),
  strict = FALSE, g = "cancer_before_acc", data = d_cancer_acc)

## main variables - group by cancer vs healthy
egltable(c(
  "sleep_comp", "mvpa_comp", "lpa_comp", "sb_comp", #acomp and 0 imputed
  "sleep", "mvpa", "lpa", "sb" # raw
), strict = FALSE, g = "cancer_before_acc", data = d_cancer_acc)

## main variables - group by cancer types
egltable(c(
  "sleep_comp", "mvpa_comp", "lpa_comp", "sb_comp", #acomp and 0 imputed
  "sleep", "mvpa", "lpa", "sb" # raw
), strict = FALSE, g = "cancer_before_acc_type", data = d_cancer_acc)

egltable(c("icd_ii_subtype"), strict = FALSE, data = d_cancer_acc)

