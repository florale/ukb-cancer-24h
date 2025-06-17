source("ukb-cancer-24h-utils.R")
source(paste0(redir, "ukb_utils.R"))

# source(paste0(redir, "data/data_acc_icd.R"))
# # d_acc_icd <- readRDS(paste0(inputdir, "d_acc_icd", ".RDS"))
# 
# # var names --------------
# icd_not_cancer_any_vars <- c(
#   "icd_i_any", "icd_ix_any", "icd_v_any",
#   "icd_vii_any", "icd_viii_any", "icd_i_iii_any",
#   "icd_xi_any", "icd_xiv_any", "icd_iv_any", "icd_xiii_any", "icd_vi_any",
#   "icd_x_any"
# )
# icd_ii_fo_vars <- grep("icd_ii_.*_fo", names(d_acc_icd), value = TRUE)
# icd_ii_fl_vars <- grep("icd_ii_.*_lo", names(d_acc_icd), value = TRUE)
# 
# icd_ii_gynaecological_vars <- c("icd_ii_c53", "icd_ii_c56", "icd_ii_c54_c55")
# icd_ii_genitourinary_vars <- c("icd_ii_c64", "icd_ii_c65_c66","icd_ii_c67", "icd_ii_c68", "icd_ii_c62")
# icd_ii_headneck_vars <- c("icd_ii_c01_c14", "icd_ii_c31_c33")
# icd_ii_lung_vars <- c("icd_ii_c34")
# icd_ii_colorectal_vars <- c("icd_ii_c18_c21")
# icd_ii_gastrointestinal_vars <- c("icd_ii_c15", "icd_ii_c16", "icd_ii_c17", "icd_ii_c22", "icd_ii_c23_c24", "icd_ii_c25", "icd_ii_c26")
# icd_ii_breast_vars <- c("icd_ii_c50")
# icd_ii_prostate_vars <- c("icd_ii_c61")
# icd_ii_blood_vars <- c("icd_ii_lymphoid")
# icd_ii_melanoma_vars <- c("icd_ii_c43")
# icd_ii_skin_vars <- c("icd_ii_c44")
# icd_ii_endocrine_vars <- c("icd_ii_c73_c75")
# icd_ii_other_vars <- c("icd_ii_c30_c39", "icd_ii_c40_c41", "icd_ii_various", "icd_ii_c45", "icd_ii_c69_c72")
# icd_ii_subtype_vars <- c(
#   "icd_ii_c53", "icd_ii_c56", "icd_ii_c62", "icd_ii_c01_c14", "icd_ii_c15",
#   "icd_ii_c16", "icd_ii_c17", "icd_ii_c18_c21", "icd_ii_c22", "icd_ii_c23_c24",
#   "icd_ii_c25", "icd_ii_c26", "icd_ii_c31_c33", "icd_ii_c34", "icd_ii_c30_c39",
#   "icd_ii_c40_c41", "icd_ii_c43", "icd_ii_c44", "icd_ii_c45", "icd_ii_various", "icd_ii_c50",
#   "icd_ii_c54_c55", "icd_ii_c61", "icd_ii_c64", "icd_ii_c65_c66", "icd_ii_c67",
#   "icd_ii_c68", "icd_ii_c69_c72", "icd_ii_c73_c75", "icd_ii_lymphoid"
# )
# icd_ii_sub_vars <- c(
#   icd_ii_c53_vars, icd_ii_c56_vars, icd_ii_c62_vars, icd_ii_c01_c14_vars, icd_ii_c15_vars,
#   icd_ii_c16_vars, icd_ii_c17_vars, icd_ii_c18_c21_vars, icd_ii_c22_vars, icd_ii_c23_c24_vars,
#   icd_ii_c25_vars, icd_ii_c26_vars, icd_ii_c31_c33_vars, icd_ii_c34_vars, icd_ii_c30_c39_vars,
#   icd_ii_c40_c41_vars, icd_ii_c43_vars, icd_ii_c44_vars, icd_ii_c45_vars, icd_ii_various_vars, icd_ii_c50_vars,
#   icd_ii_c54_c55_vars, icd_ii_c61_vars, icd_ii_c64_vars, icd_ii_c65_c66_vars, icd_ii_c67_vars,
#   icd_ii_c68_vars, icd_ii_c69_c72_vars, icd_ii_c73_c75_vars, icd_ii_lymphoid_vars
# )
# icd_ii_sub_fo_vars <- paste0(icd_ii_sub_vars, "_fo")
# icd_ii_sub_diff_acc_vars <- paste0(icd_ii_sub_vars, "_diff_acc")
# 
# icd_any_vars <- c(
#   "icd_i_any", "icd_ii_any", "icd_ix_any", "icd_v_any",
#   "icd_vii_any", "icd_viii_any", "icd_i_iii_any",
#   "icd_xi_any", "icd_xiv_any", "icd_iv_any", "icd_xiii_any", "icd_vi_any",
#   "icd_x_any"
# )
# #
# # cancer data management  ---------------
# # notes for Flora: 13 cancer composite related to PA
# # (bladder, breast, colon, endometrial, oesophageal adenocarcinoma, gastric cardia, head and neck, kidney, liver, lung, myeloid leukaemia, myeloma, and rectum)
# d_acc_icd[, icd_not_cancer := ifelse((Reduce(`|`, lapply(icd_not_cancer_any_vars, function(v) f1(get(v), 1)))),
#                                      1, 0)]
# 
# d_acc_icd[, icd_any := ifelse((Reduce(`|`, lapply(icd_any_vars, function(v) f1(get(v), 1)))),
#                               1, 0)]
# 
# table(d_acc_icd$icd_not_cancer, useNA = "always")
# table(d_acc_icd$icd_ii_any, useNA = "always")
# table(d_acc_icd$icd_ii_group, useNA = "always")
# table(d_acc_icd$icd_any, useNA = "always")
# 
# # code cancer --------
# # time since first other diagnoses
# # d_acc_icd[, age_diff_fo_other_cond_acc := year(acc_startdate) - year(icd_not_ii_sub_fo)]
# d_acc_icd[, age_diff_fo_other_cond_acc := (acc_startdate - icd_not_ii_sub_fo)/365.25]
# table(round(d_acc_icd$age_diff_fo_other_cond_acc), useNA = "always")
# 
# # time since most recent other diagnoses
# d_acc_icd[, age_diff_lo_other_cond_acc := (acc_startdate - icd_not_ii_sub_lo)/365.25]
# table(round(d_acc_icd$age_diff_lo_other_cond_acc), useNA = "always")
# 
# # number of cancer diags after acc
# nrow(d_acc_icd[((acc_startdate - icd_ii_sub_fo)/365.25) <= 0])
# 
# # number of cancer diags within 1y
# nrow(d_acc_icd[((acc_startdate - icd_ii_sub_fo)/365.25) %gele% c(-1, 0)])
# 
# nrow(d_acc_icd[((acc_startdate - icd_ii_sub_fo)/365.25) < -1])
# 
# # number of other diags 1y after acc
# nrow(d_acc_icd[((acc_startdate - icd_not_ii_sub_fo)/365.25) < -1])
# 
# # number of other diags before acc
# nrow(d_acc_icd[((acc_startdate - icd_not_ii_sub_fo)/365.25) >= 0])
# 
# # number of other diags within 1y
# nrow(d_acc_icd[((acc_startdate - icd_not_ii_sub_fo)/365.25) %gele% c(-1, 0)])
# 
# # number of other diags up to 1y after acc
# nrow(d_acc_icd[((acc_startdate - icd_not_ii_sub_fo)/365.25) >= -1])
# 
# # time since any diagnoses fo
# d_acc_icd[, time_diff_acc_icd_any := (acc_startdate - icd_sub_fo)/365.25]
# d_acc_icd[time_diff_acc_icd_any %gele% c(-1, 0), time_diff_acc_icd_any := NA]
# d_acc_icd[, time_diff_acc_icd_any := ifelse(icd_any == 0 | time_diff_acc_icd_any < - 1, 0, time_diff_acc_icd_any)]
# 
# d_acc_icd[age_diff_fo_other_cond_acc %gele% c(-1, 0), age_diff_fo_other_cond_acc := NA]
# 
# # diff time between acc and cancer subtypes
# d_acc_icd[, (icd_ii_sub_diff_acc_vars) := lapply(.SD, function (x) {
#   (acc_startdate - x)/365.25
# }), .SDcols = icd_ii_sub_fo_vars]
# 
# d_acc_icd[, icd_ii_before_acc := ifelse(Reduce(`|`, lapply(icd_ii_sub_diff_acc_vars, function(v)
#   f4(get(v), 0))),
#   1, 0)]
# table(d_acc_icd$icd_ii_before_acc, useNA = "always")
# 
# # censor 1 year follow up
# d_acc_icd[, (icd_ii_sub_diff_acc_vars) := lapply(.SD, function (x) {
#   ifelse(x %gele% c(-1, 0), NA, x)
# }), .SDcols = icd_ii_sub_diff_acc_vars]
# 
# # code later positive so they wont be coded as lo (min) later, 999 is healthy
# d_acc_icd[, (icd_ii_sub_diff_acc_vars) := lapply(.SD, function (x) {
#   ifelse(x < 0, 999, x)
# }), .SDcols = icd_ii_sub_diff_acc_vars]
# 
# # main cancer predictor variables ----------------
# # factor cancer vs healthy vs other conds
# d_acc_icd[, cancer_other_before_acc := NA]
# d_acc_icd[, cancer_other_before_acc := ifelse(icd_ii_before_acc == 1, "Cancer", cancer_other_before_acc)]
# d_acc_icd[, cancer_other_before_acc := ifelse(icd_ii_before_acc == 0 & age_diff_fo_other_cond_acc > - 1, "Others", cancer_other_before_acc)]
# d_acc_icd[, cancer_other_before_acc := ifelse(icd_ii_before_acc == 0 & time_diff_acc_icd_any <= 0, "Healthy", cancer_other_before_acc)]
# d_acc_icd[, cancer_other_before_acc := as.factor(cancer_other_before_acc)]
# table(d_acc_icd$cancer_other_before_acc, useNA = "always")
# 
# # time since most **recent** cancer diagnosis at acc
# ## closet to acc by taking min abs difference between time since diag and acc
# d_acc_icd[, icd_ii_time_since_lo := do.call(pmin, c(.SD, list(na.rm = TRUE))), .SDcols = icd_ii_sub_diff_acc_vars]
# # d_acc_icd[, icd_ii_time_since_lo := ifelse(icd_ii_before_acc == 0, 0, icd_ii_time_since_lo)]
# 
# # set time since diag = 0 for other conds and healthy
# d_acc_icd[, icd_ii_time_since_lo := ifelse(cancer_other_before_acc %in% c("Healthy", "Others"), 0, icd_ii_time_since_lo)]
# d_acc_icd[, icd_ii_time_since_lo := ifelse(is.na(cancer_other_before_acc) & icd_ii_time_since_lo == 999, NA, icd_ii_time_since_lo)]
# 
# table(round(d_acc_icd$icd_ii_time_since_lo), useNA = "always")
# 
# # time since diag
# d_acc_icd[, cancer_time_since_diag_other := cut(as.numeric(icd_ii_time_since_lo),
#                                                 breaks = c(-Inf, 0, 1, 5, Inf),
#                                                 labels = c("Healthy",
#                                                            "Less than 1 year since diagnosis",
#                                                            "1-5 years since diagnosis",
#                                                            "More than 5 years since diagnosis"
#                                                 ))]
# 
# d_acc_icd[, cancer_time_since_diag_other := as.character(cancer_time_since_diag_other)]
# d_acc_icd[, cancer_time_since_diag_other := ifelse(cancer_other_before_acc == "Healthy" | cancer_time_since_diag_other == "Healthy", "Healthy", cancer_time_since_diag_other)]
# d_acc_icd[, cancer_time_since_diag_other := ifelse(cancer_other_before_acc == "Others", "Others", cancer_time_since_diag_other)]
# 
# d_acc_icd[, cancer_time_since_diag_other := as.factor(cancer_time_since_diag_other)]
# d_acc_icd[, cancer_time_since_diag_other := relevel(cancer_time_since_diag_other, ref = "Healthy")]
# 
# table(d_acc_icd$cancer_time_since_diag_other, useNA = "always")
# 
# # factor cancer vs healthy vs other conds
# # type
# d_acc_icd[, cancer_before_acc_type := NA]
# d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_gynaecological_vars, function(v) f1(get(v), 1)))) &
#                                                (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_gynaecological_vars, function(v) f1(get(v), 0)))) &
#                                                cancer_other_before_acc == "Cancer",
#                                              "Gynaecological", cancer_before_acc_type
# )]
# d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_genitourinary_vars, function(v) f1(get(v), 1)))) &
#                                                (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_genitourinary_vars, function(v) f1(get(v), 0)))) &
#                                                cancer_other_before_acc == "Cancer",
#                                              "Genitourinary", cancer_before_acc_type
# )]
# 
# d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_headneck_vars, function(v) f1(get(v), 1)))) &
#                                                (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_headneck_vars, function(v) f1(get(v), 0))))  &
#                                                cancer_other_before_acc == "Cancer",
#                                              "Head & Neck", cancer_before_acc_type
# )]
# 
# d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_lung_vars, function(v) f1(get(v), 1)))) &
#                                                (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_lung_vars, function(v) f1(get(v), 0)))) &
#                                                cancer_other_before_acc == "Cancer",
#                                              "Lung", cancer_before_acc_type
# )]
# 
# d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_colorectal_vars, function(v) f1(get(v), 1)))) &
#                                                (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_colorectal_vars, function(v) f1(get(v), 0)))) &
#                                                cancer_other_before_acc == "Cancer",
#                                              "Colorectal", cancer_before_acc_type
# )]
# 
# d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_gastrointestinal_vars, function(v) f1(get(v), 1)))) &
#                                                (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_gastrointestinal_vars, function(v) f1(get(v), 0)))) &
#                                                cancer_other_before_acc == "Cancer",
#                                              "Gastrointestinal Tract", cancer_before_acc_type
# )]
# 
# d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_breast_vars, function(v) f1(get(v), 1)))) &
#                                                (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_breast_vars, function(v) f1(get(v), 0)))) &
#                                                cancer_other_before_acc == "Cancer",
#                                              "Breast", cancer_before_acc_type
# )]
# 
# d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_prostate_vars, function(v) f1(get(v), 1)))) &
#                                                (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_prostate_vars, function(v) f1(get(v), 0)))) &
#                                                cancer_other_before_acc == "Cancer",
#                                              "Prostate", cancer_before_acc_type
# )]
# d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_blood_vars, function(v) f1(get(v), 1)))) &
#                                                (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_blood_vars, function(v) f1(get(v), 0)))) &
#                                                cancer_other_before_acc == "Cancer",
#                                              "Blood", cancer_before_acc_type
# )]
# 
# d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_melanoma_vars, function(v) f1(get(v), 1)))) &
#                                                (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_melanoma_vars, function(v) f1(get(v), 0)))) &
#                                                cancer_other_before_acc == "Cancer",
#                                              "Melanoma", cancer_before_acc_type
# )]
# 
# d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_skin_vars, function(v) f1(get(v), 1)))) &
#                                                (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_skin_vars, function(v) f1(get(v), 0)))) &
#                                                cancer_other_before_acc == "Cancer",
#                                              "Other Skin", cancer_before_acc_type
# )]
# 
# # d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_cns_vars, function(v) f1(get(v), 1)))) &
# #                                                (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_cns_vars, function(v) f1(get(v), 0)))) &
# #                                                cancer_other_before_acc == "Cancer",
# #                                              "Central Nervous System", cancer_before_acc_type
# # )]
# 
# d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_endocrine_vars, function(v) f1(get(v), 1)))) &
#                                                (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_endocrine_vars, function(v) f1(get(v), 0)))) &
#                                                cancer_other_before_acc == "Cancer",
#                                              "Endocrine Gland", cancer_before_acc_type
# )]
# 
# d_acc_icd[, cancer_before_acc_type := ifelse((Reduce(`|`, lapply(icd_ii_other_vars, function(v) f1(get(v), 1)))) &
#                                                (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_other_vars, function(v) f1(get(v), 0)))) &
#                                                cancer_other_before_acc == "Cancer",
#                                              "Other Cancer", cancer_before_acc_type
# )]
# 
# d_acc_icd[, cancer_before_acc_type := ifelse(icd_ii_type_n > 1 & cancer_other_before_acc == "Cancer", "Multiple Primary", cancer_before_acc_type)]
# d_acc_icd[, cancer_before_acc_type := ifelse(cancer_other_before_acc == "Healthy", "Healthy", cancer_before_acc_type)]
# 
# d_acc_icd[, cancer_before_acc_type := as.factor(cancer_before_acc_type)]
# d_acc_icd[, cancer_before_acc_type := relevel(cancer_before_acc_type, ref = "Healthy")]
# 
# table(d_acc_icd$cancer_before_acc_type, useNA = "always")
# 
# # check multiple primary
# nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_gynaecological == 1])
# nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_genitourinary == 1])
# nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_headneck == 1])
# nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_lung == 1])
# nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_colorectal == 1])
# nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_gastrointestinal == 1])
# nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_breast == 1]) # 2
# nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_prostate == 1])
# nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_blood == 1])
# nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_melanoma == 1])
# nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_skin == 1]) # 1
# nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_endocrine == 1])
# nrow(d_acc_icd[cancer_before_acc_type == 'Multiple Primary' & icd_ii_other == 1])
# 
# # factor cancer vs healthy vs other conds
# d_acc_icd[, cancer_before_acc_type_other := as.character(cancer_before_acc_type)]
# d_acc_icd[, cancer_before_acc_type_other := ifelse(is.na(cancer_before_acc_type) & cancer_other_before_acc == "Others",  "Others", cancer_before_acc_type_other)]
# d_acc_icd[, cancer_before_acc_type_other := as.factor(cancer_before_acc_type_other)]
# 
# table(d_acc_icd$cancer_before_acc_type_other, useNA = "always")
# 
# # exclude those had cancer after acc -------------------
# d_cancer_acc <- d_acc_icd[!is.na(cancer_other_before_acc)]
# 
# # age at acc
# d_cancer_acc[, age_at_acc := year(acc_startdate) - year_birth]
# 
# # save --------------
# table(d_cancer_acc$cancer_other_before_acc, useNA = "always")
# 
# table(d_cancer_acc$cancer_before_acc_type, useNA = "always")
# table(d_cancer_acc$cancer_before_acc_type_other, useNA = "always")
# 
# table(d_cancer_acc$cancer_time_since_diag, useNA = "always")
# table(d_cancer_acc$cancer_time_since_diag_other, useNA = "always")
# 
# table(round(d_cancer_acc$icd_ii_time_since_lo), useNA = "always")
# identical(nrow(d_cancer_acc[icd_ii_time_since_lo == 0]), (nrow(d_cancer_acc[cancer_before_acc_type_other %in% c("Healthy", "Others")])))
# 
# saveRDS(d_cancer_acc, paste0(inputdir, "d_cancer_acc", ".RDS"))

# Composition and ilr ----------------------
d_cancer_acc <- readRDS(paste0(inputdir, "d_cancer_acc", ".RDS"))
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


