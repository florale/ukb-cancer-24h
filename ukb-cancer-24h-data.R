d_acc_icd <- readRDS("d_acc_icd.RDS")
source("/Users/oliviadelia/Desktop/GitHub/ukbiobank/ukb_utils.R")
source("/Users/oliviadelia/Desktop/GitHub/ukb-cancer-24h/ukb-cancer-24h-data.R")

# exclude those had cancer after
d_acc_icd_ii_before_acc <- d_acc_icd[icd_ii_before_acc != 0 | is.na(icd_ii_before_acc)]

# # icd ii comorbidity with other long term conditions
# d_acc_icd_ii_before_acc[, icd_ii_comorbid := NA]
# d_acc_icd_ii_before_acc[, icd_ii_comorbid := ifelse(ltc == 0,  "Healthy", icd_ii_comorbid)]
# d_acc_icd_ii_before_acc[, icd_ii_comorbid := ifelse(icd_ii_n == 0 & ltc == 1, "No ICD II, other ICD only", icd_ii_comorbid)]
# 
# d_acc_icd_ii_before_acc[, icd_ii_comorbid := ifelse(icd_ii_n > 0 & ltc == 1 & Reduce(`&`, lapply(grep(paste0(ltc_vars[ltc_vars %nin% ltc_icd_ii_vars], ".*fo", collapse = "|"), names(d_icd_all), value = T), function(v)
#   is.na(get(v)))), 
#   "ICD II only", 
#   icd_ii_comorbid)]
# 
# d_acc_icd_ii_before_acc[, icd_ii_comorbid := ifelse(icd_ii_n > 0 & ltc == 1 & Reduce(`|`, lapply(grep(paste0(ltc_vars[ltc_vars %nin% ltc_icd_ii_vars], ".*fo", collapse = "|"), names(d_icd_all), value = T), function(v)
#   !is.na(get(v)))), 
#   "ICD II comorbidity", 
#   icd_ii_comorbid)]
# 
# table(d_acc_icd_ii_before_acc$icd_ii_comorbid, useNA = "always")
# 
# d_acc_icd_ii_before_acc[, icd_ii_comorbid := factor(icd_ii_comorbid, ordered = TRUE,
#                                                     levels = c("Healthy",
#                                                                "ICD II only",
#                                                                "ICD II comorbidity",
#                                                                "No ICD II, other ICD only"))]
# 
# egltable(c("age", "sex", "ethnicg", "white", "bmi", "edu", "working", "deprivation", 
#            "current_drinker", "never_smoked", "smoking", "alcohol",
#            
#            "icd_ii_comorbid",
#            
#            "sleep_comp", "mvpa_comp", "lpa_comp", "sb_comp", #acomp and 0 imputed
#            "sleep", "mvpa", "lpa", "sb" # raw
# ), strict = FALSE, data = d_acc_icd_ii_before_acc)
# 
# 
# clr_acc_icd_ii <- complr(data = d_acc_icd_ii_before_acc[, -colnames(ilr_acc), with = FALSE],
#                          transform = "ilr",
#                          parts = c("sleep_comp", "mvpa_comp", "lpa_comp", "sb_comp"),
#                          sbp = sbp,
#                          total = 1440,
#                          idvar = "eid",
#                          shape = "wide")
# saveRDS(clr_acc_icd_ii, paste0(inputdir, "d_acc_icd_ii_before_acc", ".RDS"))

## Cancer subgroup with ICD-10 Codes 
icd_ii_fo_vars <- grep("icd_ii_.*_fo", names(d_acc_icd), value = TRUE)
icd_ii_vars <- str_remove(icd_ii_fo_vars, "_fo")

icd_ii_c53_vars <- grep("C53", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c53 := ifelse(Reduce(`|`, lapply(icd_ii_c53_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c53, useNA = "always")

icd_ii_c56_vars <- grep("C56", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c56 := ifelse(Reduce(`|`, lapply(icd_ii_c56_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c56, useNA = "always")

icd_ii_c62_vars <- grep("C62", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c62 := ifelse(Reduce(`|`, lapply(icd_ii_c62_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c62, useNA = "always")

icd_ii_c01_c14_vars <- grep("C0[1-9]|C1[0-4]", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c01_c14 := ifelse(Reduce(`|`, lapply(icd_ii_c01_c14_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c01_c14, useNA = "always")

icd_ii_c15_vars <- grep("C15", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c15 := ifelse(Reduce(`|`, lapply(icd_ii_c15_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c15, useNA = "always")

icd_ii_c16_vars <- grep("C16", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c16 := ifelse(Reduce(`|`, lapply(icd_ii_c16_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c16, useNA = "always")

icd_ii_c17_vars <- grep("C17", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c17 := ifelse(Reduce(`|`, lapply(icd_ii_c17_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c17, useNA = "always")

icd_ii_c18_c21_vars <- grep("C1[8-9]|C21", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c18_c21 := ifelse(Reduce(`|`, lapply(icd_ii_c18_c21_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c18_c21, useNA = "always")

icd_ii_c22_vars <- grep("C22", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c22 := ifelse(Reduce(`|`, lapply(icd_ii_c22_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c22, useNA = "always")

icd_ii_c23_c24_vars <- grep("C2[3-4]", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c23_c24 := ifelse(Reduce(`|`, lapply(icd_ii_c23_c24_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c23_c24, useNA = "always")

icd_ii_c25_vars <- grep("C25", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c25 := ifelse(Reduce(`|`, lapply(icd_ii_c25_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c25, useNA = "always")

icd_ii_c26_vars <- grep("C26", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c26 := ifelse(Reduce(`|`, lapply(icd_ii_c26_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c26, useNA = "always")

icd_ii_c31_c33_vars <- grep("C3[1-3]", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c31_c33 := ifelse(Reduce(`|`, lapply(icd_ii_c31_c33_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c31_c33, useNA = "always")

icd_ii_c34_vars <- grep("C34", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c34 := ifelse(Reduce(`|`, lapply(icd_ii_c34_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c34, useNA = "always")

icd_ii_c30_c39_vars <- grep("C30|C3[7-9]", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c30_c39 := ifelse(Reduce(`|`, lapply(icd_ii_c30_c39_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c30_c39, useNA = "always")

icd_ii_c40_c41_vars <- grep("C4[0-1]", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c40_c41 := ifelse(Reduce(`|`, lapply(icd_ii_c40_c41_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c40_c41, useNA = "always")

icd_ii_c44_vars <- grep("C44", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c44 := ifelse(Reduce(`|`, lapply(icd_ii_c44_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c44, useNA = "always")

icd_ii_c44_c45_vars <- grep("C4[4-5]", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c44_c45 := ifelse(Reduce(`|`, lapply(icd_ii_c44_c45_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c44_c45, useNA = "always")

icd_ii_other_vars <- grep("C4[6-9]|C5[1-2|7-8]|C60|C63|C76|C80|C97", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_other := ifelse(Reduce(`|`, lapply(icd_ii_other_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_other, useNA = "always")

icd_ii_c50_vars <- grep("C50", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c50 := ifelse(Reduce(`|`, lapply(icd_ii_c50_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c50, useNA = "always")

icd_ii_c54_c55_vars <- grep("C5[4-5]", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c54_c55 := ifelse(Reduce(`|`, lapply(icd_ii_c54_c55_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c54_c55, useNA = "always")

icd_ii_c61_vars <- grep("C61", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c61 := ifelse(Reduce(`|`, lapply(icd_ii_c61_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c61, useNA = "always")

icd_ii_c64_vars <- grep("C64", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c64 := ifelse(Reduce(`|`, lapply(icd_ii_c64_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c64, useNA = "always")

icd_ii_c65_c66_vars <- grep("C6[5-6]", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c65_c66 := ifelse(Reduce(`|`, lapply(icd_ii_c65_c66_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c65_c66, useNA = "always")

icd_ii_c67_vars <- grep("C67", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c67 := ifelse(Reduce(`|`, lapply(icd_ii_c67_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c67, useNA = "always")

icd_ii_c68_vars <- grep("C68", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c68 := ifelse(Reduce(`|`, lapply(icd_ii_c68_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c68, useNA = "always")

icd_ii_c69_c72_vars <- grep("C69|C7[0-2]", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c69_c72 := ifelse(Reduce(`|`, lapply(icd_ii_c69_c72_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c69_c72, useNA = "always")

icd_ii_c73_c75_vars <- grep("C7[3-5]", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c73_c75 := ifelse(Reduce(`|`, lapply(icd_ii_c73_c75_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c73_c75, useNA = "always")

icd_ii_lymphoid_vars <- grep("C8[1-6]|C88|C9[0-6]", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_lymphoid := ifelse(Reduce(`|`, lapply(icd_ii_lymphoid_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_lymphoid, useNA = "always")


# Infectious diseases with ICD-10 Code 
icd_i_fo_vars <- grep("icd_i_.*_fo", names(d_acc_icd), value = TRUE)
icd_i_vars <- str_remove(icd_i_fo_vars, "_fo")

icd_i_tuberculosis_vars <- grep("A1[5-9]", icd_i_vars, perl = TRUE, value = T)
d_acc_icd[, icd_i_tuberculosis := ifelse(Reduce(`|`, lapply(icd_i_tuberculosis_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_i_tuberculosis, useNA = "always")

icd_i_hepatitis_vars <- grep("B18", icd_i_vars, perl = TRUE, value = T)
d_acc_icd[, icd_i_hepatitis := ifelse(Reduce(`|`, lapply(icd_i_hepatitis_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_i_hepatitis, useNA = "always")

icd_i_hiv_vars <- grep("B2[0-3]", icd_i_vars, perl = TRUE, value = T)
d_acc_icd[, icd_i_hiv := ifelse(Reduce(`|`, lapply(icd_i_hiv_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_i_hiv, useNA = "always")

# Diseases of the Circulatory System 
icd_ix_fo_vars <- grep("icd_ix_.*_fo", names(d_acc_icd), value = TRUE)
icd_ix_vars <- str_remove(icd_ix_fo_vars, "_fo")

icd_ix_hypertension_vars <- grep("I1[0-3]|I15", icd_ix_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ix_hypertension := ifelse(Reduce(`|`, lapply(icd_ix_hypertension_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ix_hypertension, useNA = "always")

icd_ix_valve_vars <- grep("I3[4-7]", icd_ix_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ix_valve := ifelse(Reduce(`|`, lapply(icd_ix_valve_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ix_valve, useNA = "always")

icd_ix_i20_vars <- grep("I20", icd_ix_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ix_i20 := ifelse(Reduce(`|`, lapply(icd_ix_i20_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ix_i20, useNA = "always")

icd_ix_myo_vars <- grep("I2[1-2]|I25", icd_ix_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ix_myo := ifelse(Reduce(`|`, lapply(icd_ix_myo_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ix_myo, useNA = "always")

icd_ix_i50_vars <- grep("I50", icd_ix_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ix_i50 := ifelse(Reduce(`|`, lapply(icd_ix_i50_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ix_i50, useNA = "always")

icd_ix_i26_vars <- grep("I26", icd_ix_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ix_i26 := ifelse(Reduce(`|`, lapply(icd_ix_i26_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ix_i26, useNA = "always")

icd_ix_i74_vars <- grep("I74", icd_ix_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ix_i74 := ifelse(Reduce(`|`, lapply(icd_ix_i74_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ix_i74, useNA = "always")

icd_ix_venous_vars <- grep("I8[1-2]", icd_ix_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ix_venous := ifelse(Reduce(`|`, lapply(icd_ix_venous_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ix_venous, useNA = "always")

icd_ix_i63_vars <- grep("I63", icd_ix_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ix_i63 := ifelse(Reduce(`|`, lapply(icd_ix_i63_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ix_i63, useNA = "always")

icd_ix_haemorrhagic_vars <- grep("I6[1-2]", icd_ix_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ix_haemorrhagic := ifelse(Reduce(`|`, lapply(icd_ix_haemorrhagic_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ix_haemorrhagic, useNA = "always")

icd_ix_i42_vars <- grep("I42", icd_ix_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ix_i42 := ifelse(Reduce(`|`, lapply(icd_ix_i42, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ix_i42, useNA = "always")

icd_ix_arrhythmia_vars <- grep("I4[4-5]|I47|I49", icd_ix_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ix_arrhythmia := ifelse(Reduce(`|`, lapply(icd_ix_arrhythmia_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ix_arrhythmia, useNA = "always")

icd_ix_aneurysms_vars <- grep("I7[1-2]", icd_ix_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ix_aneurysms := ifelse(Reduce(`|`, lapply(icd_ix_aneurysms_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ix_aneurysms, useNA = "always")

icd_ix_rheumatic_vars <- grep("I0[5-9]", icd_ix_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ix_rheumatic := ifelse(Reduce(`|`, lapply(icd_ix_rheumatic_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ix_rheumatic, useNA = "always")

# Mental and Behavioural Disorders 
icd_v_fo_vars <- grep("icd_v_.*_fo", names(d_acc_icd), value = TRUE)
icd_v_vars <- str_remove(icd_v_fo_vars, "_fo")

icd_v_alcohol_vars <- grep("F10", icd_v_vars, perl = TRUE, value = T)
d_acc_icd[, icd_v_alcohol := ifelse(Reduce(`|`, lapply(icd_v_alcohol_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_v_alcohol, useNA = "always")

icd_v_bipolar_vars <- grep("F31", icd_v_vars, perl = TRUE, value = T)
d_acc_icd[, icd_v_bipolar := ifelse(Reduce(`|`, lapply(icd_v_bipolar_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_v_bipolar, useNA = "always")

icd_v_depressive_vars <- grep("F3[2-3]", icd_v_vars, perl = TRUE, value = T)
d_acc_icd[, icd_v_depressive := ifelse(Reduce(`|`, lapply(icd_v_depressive_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_v_depressive, useNA = "always")

icd_v_anxiety_vars <- grep("F4[0-1]", icd_v_vars, perl = TRUE, value = T)
d_acc_icd[, icd_v_anxiety := ifelse(Reduce(`|`, lapply(icd_v_anxiety_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_v_anxiety, useNA = "always")

icd_v_stress_vars <- grep("F43", icd_v_vars, perl = TRUE, value = T)
d_acc_icd[, icd_v_stress := ifelse(Reduce(`|`, lapply(icd_v_stress_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_v_stress, useNA = "always")

icd_v_schizophrenia_vars <- grep("F20", icd_v_vars, perl = TRUE, value = T)
d_acc_icd[, icd_v_schizophrenia := ifelse(Reduce(`|`, lapply(icd_v_schizophrenia_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_v_schizophrenia, useNA = "always")

icd_v_personality_vars <- grep("F60", icd_v_vars, perl = TRUE, value = T)
d_acc_icd[, icd_v_personality := ifelse(Reduce(`|`, lapply(icd_v_personality_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_v_personality, useNA = "always")

icd_v_ocd_vars <- grep("F42", icd_v_vars, perl = TRUE, value = T)
d_acc_icd[, icd_v_ocd := ifelse(Reduce(`|`, lapply(icd_v_ocd_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_v_ocd, useNA = "always")

icd_v_ed_vars <- grep("F50", icd_v_vars, perl = TRUE, value = T)
d_acc_icd[, icd_v_ed := ifelse(Reduce(`|`, lapply(icd_v_ed_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_v_ed, useNA = "always")

# Disorders of the Eye 
icd_vii_fo_vars <- grep("icd_vii_.*_fo", names(d_acc_icd), value = TRUE)
icd_vii_vars <- str_remove(icd_vii_fo_vars, "_fo")

icd_vii_cataracts_vars <- grep("H2[5-6]", icd_vii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_vii_cataracts := ifelse(Reduce(`|`, lapply(icd_vii_cataracts_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vii_cataracts, useNA = "always")

icd_vii_retina_vars <- grep("H3[0-1]|H3[3-6]", icd_vii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_vii_retina := ifelse(Reduce(`|`, lapply(icd_vii_retina_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vii_retina, useNA = "always")

icd_vii_glaucoma_vars <- grep("H40|H42", icd_vii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_vii_glaucoma := ifelse(Reduce(`|`, lapply(icd_vii_glaucoma_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vii_glaucoma, useNA = "always")

icd_vii_muscles_vars <- grep("H49|H5[0-1]", icd_vii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_vii_muscles := ifelse(Reduce(`|`, lapply(icd_vii_muscles_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vii_muscles, useNA = "always")

icd_vii_blind_vars <- grep("H5[3-4]", icd_vii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_vii_blind := ifelse(Reduce(`|`, lapply(icd_vii_blind_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vii_blind, useNA = "always")

# Disorders of the Ear and Mastoid  
icd_viii_fo_vars <- grep("icd_viii_.*_fo", names(d_acc_icd), value = TRUE)
icd_viii_vars <- str_remove(icd_viii_fo_vars, "_fo")

icd_viii_inner_vars <- grep("H8[0-1]|H83", icd_viii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_viii_inner := ifelse(Reduce(`|`, lapply(icd_viii_inner_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_viii_inner, useNA = "always")

icd_viii_loss_vars <- grep("H9[0-1]", icd_viii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_viii_loss := ifelse(Reduce(`|`, lapply(icd_viii_loss_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_viii_loss, useNA = "always")

# Disorders of blood and the blood-forming organs and 
# certain disorders involving the immune mechanism
icd_i_fo_vars <- grep("icd_i_.*_fo", names(d_acc_icd), value = TRUE)
icd_i_vars <- str_remove(icd_i_fo_vars, "_fo")

icd_i_anaemias_vars <- grep("B5[0-3]", icd_i_vars, perl = TRUE, value = T)
d_acc_icd[, icd_i_anaemias := ifelse(Reduce(`|`, lapply(icd_i_anaemias_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_i_anaemias, useNA = "always")

icd_i_hanaemias_vars <- grep("B5[5-9]", icd_i_vars, perl = TRUE, value = T)
d_acc_icd[, icd_i_hanaemias := ifelse(Reduce(`|`, lapply(icd_i_hanaemias_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_i_hanaemias, useNA = "always")

icd_iii_fo_vars <- grep("icd_iii_.*_fo", names(d_acc_icd), value = TRUE)
icd_iii_vars <- str_remove(icd_iii_fo_vars, "_fo")

icd_iii_aanaemias_vars <- grep("D6[0-1]", icd_iii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_iii_aanaemias := ifelse(Reduce(`|`, lapply(icd_iii_aanaemias_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_iii_aanaemias, useNA = "always")

icd_iii_d63_vars <- grep("D63", icd_iii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_iii_d63 := ifelse(Reduce(`|`, lapply(icd_iii_d63_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_iii_d63, useNA = "always")

icd_iii_coagulation_vars <- grep("D6[6-8]", icd_iii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_iii_coagulation := ifelse(Reduce(`|`, lapply(icd_iii_coagulation_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_iii_coagulation, useNA = "always")

icd_iii_d69_vars <- grep("D69", icd_iii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_iii_d69 := ifelse(Reduce(`|`, lapply(icd_iii_d69_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_iii_d69, useNA = "always")

icd_iii_immunodeficiency_vars <- grep("D70|D8[0-4]", icd_iii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_iii_immunodeficiency := ifelse(Reduce(`|`, lapply(icd_iii_immunodeficiency_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_iii_immunodeficiency, useNA = "always")

icd_iii_d86_vars <- grep("D86", icd_iii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_iii_d86 := ifelse(Reduce(`|`, lapply(icd_iii_d86_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_iii_d86, useNA = "always")

# Disorders of the digestive system
icd_xi_fo_vars <- grep("icd_xi_.*_fo", names(d_acc_icd), value = TRUE)
icd_xi_vars <- str_remove(icd_xi_fo_vars, "_fo")

# Disorders of the genitourinary system
icd_xiv_fo_vars <- grep("icd_xiv_.*_fo", names(d_acc_icd), value = TRUE)
icd_xiv_vars <- str_remove(icd_xiv_fo_vars, "_fo")

icd_xiv_renal_vars <- grep("N1[8-9]", icd_xiv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiv_renal := ifelse(Reduce(`|`, lapply(icd_xiv_renal_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiv_renal, useNA = "always")

icd_xiv_urolithiasis_vars <- grep("N2[0-1]", icd_xiv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiv_urolithiasis := ifelse(Reduce(`|`, lapply(icd_xiv_urolithiasis_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiv_urolithiasis, useNA = "always")

icd_xiv_n30_vars <- grep("N30", icd_xiv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiv_n30 := ifelse(Reduce(`|`, lapply(icd_xiv_n30_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiv_n30, useNA = "always")

icd_xiv_n31_vars <- grep("N31", icd_xiv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiv_n31 := ifelse(Reduce(`|`, lapply(icd_xiv_n31_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiv_n31, useNA = "always")

icd_xiv_n40_vars <- grep("N40", icd_xiv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiv_n40 := ifelse(Reduce(`|`, lapply(icd_xiv_n40_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiv_n40, useNA = "always")

icd_xiv_n41_vars <- grep("N41", icd_xiv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiv_n41 := ifelse(Reduce(`|`, lapply(icd_xiv_n41_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiv_n41, useNA = "always")

icd_xiv_n43_vars <- grep("N43", icd_xiv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiv_n43 := ifelse(Reduce(`|`, lapply(icd_xiv_n43_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiv_n43, useNA = "always")

icd_xiv_n60_vars <- grep("N60", icd_xiv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiv_n60 := ifelse(Reduce(`|`, lapply(icd_xiv_n60_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiv_n60, useNA = "always")

icd_xiv_n61_vars <- grep("N61", icd_xiv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiv_n61 := ifelse(Reduce(`|`, lapply(icd_xiv_n61_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiv_n61, useNA = "always")

icd_xiv_breastpelvis_vars <- grep("N6[5-6]|N70|N73", icd_xiv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiv_breastpelvis := ifelse(Reduce(`|`, lapply(icd_xiv_breastpelvis_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiv_breastpelvis, useNA = "always")

icd_xiv_n70_vars <- grep("N70", icd_xiv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiv_n70 := ifelse(Reduce(`|`, lapply(icd_xiv_n70_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiv_n70, useNA = "always")

icd_xiv_n71_vars <- grep("N71", icd_xiv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiv_n71 := ifelse(Reduce(`|`, lapply(icd_xiv_n71_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiv_n71, useNA = "always")

icd_xiv_n72_vars <- grep("N72", icd_xiv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiv_n72 := ifelse(Reduce(`|`, lapply(icd_xiv_n72_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiv_n72, useNA = "always")

icd_xiv_n80_vars <- grep("N80", icd_xiv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiv_n80 := ifelse(Reduce(`|`, lapply(icd_xiv_n80_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiv_n80, useNA = "always")

icd_xiv_n81_vars <- grep("N81", icd_xiv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiv_n81 := ifelse(Reduce(`|`, lapply(icd_xiv_n81_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiv_n81, useNA = "always")

icd_xiv_n82_vars <- grep("N82", icd_xiv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiv_n82 := ifelse(Reduce(`|`, lapply(icd_xiv_n82_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiv_n82, useNA = "always")

icd_xiv_n92_vars <- grep("N92", icd_xiv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiv_n92 := ifelse(Reduce(`|`, lapply(icd_xiv_n92_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiv_n92, useNA = "always")

icd_xiv_n97_vars <- grep("N97", icd_xiv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiv_n97 := ifelse(Reduce(`|`, lapply(icd_xiv_n97_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiv_n97, useNA = "always")

icd_xiv_n46_vars <- grep("N46", icd_xiv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiv_n46 := ifelse(Reduce(`|`, lapply(icd_xiv_n46_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiv_n46, useNA = "always")
