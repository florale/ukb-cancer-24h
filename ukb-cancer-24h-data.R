d_acc_icd <- readRDS("~/Desktop/GitHub/ukb-cancer-24h/d_acc_icd.RDS")
source("/Users/oliviadelia/Desktop/GitHub/ukbiobank/ukb_utils.R")
source("/Users/oliviadelia/Desktop/GitHub/ukb-cancer-24h/ukb-cancer-24h-data.R")



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


# 2 Infectious diseases I ------------------------
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

icd_i_subtype_vars <- c(
  "icd_i_tuberculosis",
  "icd_i_hepatitis",
  "icd_i_hiv"
)

d_acc_icd[, icd_i_any := ifelse(Reduce(`|`, lapply(icd_i_subtype_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_i_any, useNA = "always")

# 3 Diseases of the Circulatory System Chapter IX  -----------------------------
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
d_acc_icd[, icd_ix_i42 := ifelse(Reduce(`|`, lapply(icd_ix_i42_vars, function(v)
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

icd_ix_subtype_vars <- c(
  "icd_ix_hypertension", "icd_ix_valve", "icd_ix_i20", "icd_ix_myo", "icd_ix_i50", 
  "icd_ix_i26", "icd_ix_i74", "icd_ix_venous", "icd_ix_i63", "icd_ix_haemorrhagic", 
  "icd_ix_i42", "icd_ix_arrhythmia",  "icd_ix_aneurysms", "icd_ix_rheumatic"
)

d_acc_icd[, icd_ix_any := ifelse(Reduce(`|`, lapply(icd_ix_subtype_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ix_any, useNA = "always")


# 4 Mental and Behavioral Disorders V ----------------
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

icd_v_subtype_vars <- c(
  "icd_v_alcohol", "icd_v_bipolar", "icd_v_depressive", "icd_v_anxiety", 
  "icd_v_stress", "icd_v_schizophrenia", "icd_v_personality", "icd_v_ocd", "icd_v_ed"
)

d_acc_icd[, icd_v_any := ifelse(Reduce(`|`, lapply(icd_v_subtype_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_v_any, useNA = "always")

# 5 Disorders of the Eye VII ----------
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

icd_vii_subtype_vars <- c(
  "icd_vii_cataracts", "icd_vii_retina", 
  "icd_vii_glaucoma", "icd_vii_muscles", "icd_vii_blind"
)

d_acc_icd[, icd_vii_any := ifelse(Reduce(`|`, lapply(icd_vii_subtype_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vii_any, useNA = "always")




# 6 Disorders of the Ear and Mastoid VIII ---------------
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

icd_viii_subtype_vars <- c(
  "icd_viii_inner",
  "icd_viii_loss"
)

d_acc_icd[, icd_viii_any := ifelse(Reduce(`|`, lapply(icd_viii_subtype_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_viii_any, useNA = "always")


# 7 Disorders of blood and the blood-forming organs and certain disorders involving the immune mechanism I & III---------
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

icd_i_iii_subtype_vars <- c(
  "icd_i_anaemias", "icd_i_hanaemias", 
  "icd_iii_aanaemias", "icd_iii_d63", "icd_iii_coagulation", 
  "icd_iii_d69", "icd_iii_immunodeficiency", "icd_iii_d86"
)

d_acc_icd[, icd_i_iii_any := ifelse(Reduce(`|`, lapply(icd_i_iii_subtype_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_i_iii_any, useNA = "always")


# 8 Disorders of the digestive system XI ----------
icd_xi_fo_vars <- grep("icd_xi_.*_fo", names(d_acc_icd), value = TRUE)
icd_xi_vars <- str_remove(icd_xi_fo_vars, "_fo")

icd_xi_k21_vars <- grep("K21", icd_xi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xi_k21 := ifelse(Reduce(`|`, lapply(icd_xi_k21_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xi_k21, useNA = "always")

icd_xi_k22_vars <- grep("K22", icd_xi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xi_k22 := ifelse(Reduce(`|`, lapply(icd_xi_k22_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xi_k22, useNA = "always")

icd_xi_ulcer_vars <- grep("K2[5-8]", icd_xi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xi_ulcer := ifelse(Reduce(`|`, lapply(icd_xi_ulcer_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xi_ulcer, useNA = "always")

icd_xi_k29_vars <- grep("K29", icd_xi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xi_k29 := ifelse(Reduce(`|`, lapply(icd_xi_k29_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xi_k29, useNA = "always")

icd_xi_intestinal_vars <- grep("K5[8-9]|K30", icd_xi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xi_intestinal := ifelse(Reduce(`|`, lapply(icd_xi_intestinal_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xi_intestinal, useNA = "always")

icd_xi_appendicitis_vars <- grep("K3[5-7]", icd_xi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xi_appendicitis := ifelse(Reduce(`|`, lapply(icd_xi_appendicitis_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xi_appendicitis, useNA = "always")

icd_xi_hernia_vars <- grep("K4[0-6]", icd_xi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xi_hernia := ifelse(Reduce(`|`, lapply(icd_xi_hernia_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xi_hernia, useNA = "always")

icd_xi_k50_vars <- grep("K50", icd_xi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xi_k50 := ifelse(Reduce(`|`, lapply(icd_xi_k50_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xi_k50, useNA = "always")

icd_xi_k51_vars <- grep("K51", icd_xi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xi_k51 := ifelse(Reduce(`|`, lapply(icd_xi_k51_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xi_k51, useNA = "always")

icd_xi_k55_vars <- grep("K55", icd_xi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xi_k55 := ifelse(Reduce(`|`, lapply(icd_xi_k55_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xi_k55, useNA = "always")

icd_xi_k57_vars <- grep("K57", icd_xi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xi_k57 := ifelse(Reduce(`|`, lapply(icd_xi_k57_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xi_k57, useNA = "always")

icd_xi_k76_vars <- grep("K76", icd_xi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xi_k76 := ifelse(Reduce(`|`, lapply(icd_xi_k76_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xi_k76, useNA = "always")

icd_xi_k81_vars <- grep("K81", icd_xi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xi_k81 := ifelse(Reduce(`|`, lapply(icd_xi_k81_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xi_k81, useNA = "always")

icd_xi_k86_vars <- grep("K86", icd_xi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xi_k86 := ifelse(Reduce(`|`, lapply(icd_xi_k86_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xi_k86, useNA = "always")

icd_xi_k90_vars <- grep("K90", icd_xi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xi_k90 := ifelse(Reduce(`|`, lapply(icd_xi_k90_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xi_k90, useNA = "always")

icd_xi_subtype_vars <- c(
  "icd_xi_k21", "icd_xi_k22", "icd_xi_ulcer", "icd_xi_k29", "icd_xi_intestinal", 
  "icd_xi_appendicitis", "icd_xi_hernia", "icd_xi_k50", "icd_xi_k51", "icd_xi_k55", 
  "icd_xi_k57", "icd_xi_k76", "icd_xi_k81", "icd_xi_k86", "icd_xi_k90"
)

d_acc_icd[, icd_xi_any := ifelse(Reduce(`|`, lapply(icd_xi_subtype_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xi_any, useNA = "always")



# 9 Disorders of the genitourinary system XIV -----------
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

icd_xiv_subtype_vars <- c(
  "icd_xiv_renal", "icd_xiv_urolithiasis", "icd_xiv_n30", "icd_xiv_n31", "icd_xiv_n40", 
  "icd_xiv_n41", "icd_xiv_n43", "icd_xiv_n60", "icd_xiv_n61", "icd_xiv_breastpelvis", 
  "icd_xiv_n70", "icd_xiv_n71", "icd_xiv_n72", "icd_xiv_n80", "icd_xiv_n81", "icd_xiv_n82", 
  "icd_xiv_n92", "icd_xiv_n97", "icd_xiv_n46"
)

d_acc_icd[, icd_xiv_any := ifelse(Reduce(`|`, lapply(icd_xiv_subtype_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiv_any, useNA = "always")



# 10 Diseases of the musculoskeletal system and connective tissues, endocrine, nutritional, and metabolic diseases IV & XIII----------
icd_iv_fo_vars <- grep("icd_iv_.*_fo", names(d_acc_icd), value = TRUE)
icd_iv_vars <- str_remove(icd_iv_fo_vars, "_fo")

icd_iv_e03_vars <- grep("E01|E03", icd_iv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_iv_e03 := ifelse(Reduce(`|`, lapply(icd_iv_e03_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_iv_e03, useNA = "always")

icd_iv_e04_vars <- grep("E04", icd_iv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_iv_e04 := ifelse(Reduce(`|`, lapply(icd_iv_e04_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_iv_e04, useNA = "always")

icd_iv_e05_vars <- grep("E05", icd_iv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_iv_e05 := ifelse(Reduce(`|`, lapply(icd_iv_e05_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_iv_e05, useNA = "always")

icd_iv_diabetes_vars <- grep("E1[0-1]|E1[3-4]", icd_iv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_iv_diabetes := ifelse(Reduce(`|`, lapply(icd_iv_e05_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_iv_diabetes, useNA = "always")

icd_iv_e21_vars <- grep("E21", icd_iv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_iv_e21 := ifelse(Reduce(`|`, lapply(icd_iv_e21_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_iv_e21, useNA = "always")

icd_iv_e22_vars <- grep("E22", icd_iv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_iv_e22 := ifelse(Reduce(`|`, lapply(icd_iv_e22_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_iv_e22, useNA = "always")

icd_iv_e23_vars <- grep("E23", icd_iv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_iv_e23 := ifelse(Reduce(`|`, lapply(icd_iv_e23_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_iv_e23, useNA = "always")

icd_iv_e24_vars <- grep("E24", icd_iv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_iv_e24 := ifelse(Reduce(`|`, lapply(icd_iv_e24_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_iv_e24, useNA = "always")

icd_iv_e27_vars <- grep("E27", icd_iv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_iv_e27 := ifelse(Reduce(`|`, lapply(icd_iv_e27_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_iv_e27, useNA = "always")

icd_iv_e78_vars <- grep("E78", icd_iv_vars, perl = TRUE, value = T)
d_acc_icd[, icd_iv_e78 := ifelse(Reduce(`|`, lapply(icd_iv_e78_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_iv_e78, useNA = "always")

icd_xiii_fo_vars <- grep("icd_xiii_.*_fo", names(d_acc_icd), value = TRUE)
icd_xiii_vars <- str_remove(icd_xiii_fo_vars, "_fo")

icd_xiii_arthritis_vars <- grep("M0[5-6]", icd_xiii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiii_arthritis := ifelse(Reduce(`|`, lapply(icd_xiii_arthritis_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiii_arthritis, useNA = "always")

icd_xiii_arthopathies_vars <- grep("M1[0-1]", icd_xiii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiii_arthopathies := ifelse(Reduce(`|`, lapply(icd_xiii_arthopathies_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiii_arthopathies, useNA = "always")

icd_xiii_arthroses_vars <- grep("M1[5-9]", icd_xiii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiii_arthroses := ifelse(Reduce(`|`, lapply(icd_xiii_arthroses_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiii_arthroses, useNA = "always")

icd_xiii_m23_vars <- grep("M23", icd_xiii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiii_m23 := ifelse(Reduce(`|`, lapply(icd_xiii_m23_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiii_m23, useNA = "always")

icd_xiii_m54_vars <- grep("M54", icd_xiii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiii_m54 := ifelse(Reduce(`|`, lapply(icd_xiii_m54_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiii_m54, useNA = "always")

icd_xiii_m75_vars <- grep("M75", icd_xiii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiii_m75 := ifelse(Reduce(`|`, lapply(icd_xiii_m75_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiii_m75, useNA = "always")

icd_xiii_osteo_vars <- grep("M8[0-1]", icd_xiii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_xiii_osteo := ifelse(Reduce(`|`, lapply(icd_xiii_osteo_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_xiii_osteo, useNA = "always")

icd_iv_xiii_subtype_vars <- c(
  "icd_iv_e03", "icd_iv_e04", "icd_iv_e05", "icd_iv_diabetes", "icd_iv_e21", 
  "icd_iv_e22", "icd_iv_e23", "icd_iv_e24", "icd_iv_e27", "icd_iv_e78", 
  "icd_xiii_arthritis", "icd_xiii_arthopathies", "icd_xiii_arthroses", "icd_xiii_m23", 
  "icd_xiii_m54", "icd_xiii_m75", "icd_xiii_osteo"
)

d_acc_icd[, icd_iv_xiii_any := ifelse(Reduce(`|`, lapply(icd_iv_xiii_subtype_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_iv_xiii_any, useNA = "always")

# 11 Diseases of the nervous system VI ----------
icd_vi_fo_vars <- grep("icd_vi_.*_fo", names(d_acc_icd), value = TRUE)
icd_vi_vars <- str_remove(icd_vi_fo_vars, "_fo")

icd_vi_n06_vars <- grep("G06", icd_vi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_vi_n06 := ifelse(Reduce(`|`, lapply(icd_vi_n06_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vi_n06, useNA = "always")

icd_vi_atrophies_vars <- grep("G1[1-4]", icd_vi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_vi_atrophies := ifelse(Reduce(`|`, lapply(icd_vi_atrophies_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vi_atrophies, useNA = "always")

icd_vi_movement_vars <- grep("G2[0-1]|G23|G25", icd_vi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_vi_movement := ifelse(Reduce(`|`, lapply(icd_vi_movement_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vi_movement, useNA = "always")

icd_vi_cns_vars <- grep("G3[0-1]", icd_vi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_vi_cns := ifelse(Reduce(`|`, lapply(icd_vi_cns_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vi_cns, useNA = "always")

icd_vi_g35_vars <- grep("G35", icd_vi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_vi_g35 := ifelse(Reduce(`|`, lapply(icd_vi_g35_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vi_g35, useNA = "always")

icd_vi_other_vars <- grep("G3[6-7]", icd_vi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_vi_other := ifelse(Reduce(`|`, lapply(icd_vi_other_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vi_other, useNA = "always")

icd_vi_g40_vars <- grep("G40", icd_vi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_vi_g40 := ifelse(Reduce(`|`, lapply(icd_vi_g40_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vi_g40, useNA = "always")

icd_vi_g43_vars <- grep("G43", icd_vi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_vi_g43 := ifelse(Reduce(`|`, lapply(icd_vi_g43_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vi_g43, useNA = "always")

icd_vi_g45_vars <- grep("G45", icd_vi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_vi_g45 := ifelse(Reduce(`|`, lapply(icd_vi_g45_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vi_g45, useNA = "always")

icd_vi_g47_vars <- grep("G47", icd_vi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_vi_g47 := ifelse(Reduce(`|`, lapply(icd_vi_g47_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vi_g47, useNA = "always")

icd_vi_mononeuropathy_vars <- grep("G54|G57", icd_vi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_vi_mononeuropathy := ifelse(Reduce(`|`, lapply(icd_vi_mononeuropathy_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vi_mononeuropathy, useNA = "always")

icd_vi_polyneuropathy_vars <- grep("G6[0-3]", icd_vi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_vi_polyneuropathy := ifelse(Reduce(`|`, lapply(icd_vi_polyneuropathy_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vi_polyneuropathy, useNA = "always")

icd_vi_myoneural_vars <- grep("G70|G73", icd_vi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_vi_myoneural := ifelse(Reduce(`|`, lapply(icd_vi_myoneural_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vi_myoneural, useNA = "always")

icd_vi_g71_vars <- grep("G71", icd_vi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_vi_g71 := ifelse(Reduce(`|`, lapply(icd_vi_g71_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vi_g71, useNA = "always")

icd_vi_paralytic_vars <- grep("G8[0-3]", icd_vi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_vi_paralytic := ifelse(Reduce(`|`, lapply(icd_vi_paralytic_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vi_paralytic, useNA = "always")

icd_vi_g91_vars <- grep("G91", icd_vi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_vi_g91 := ifelse(Reduce(`|`, lapply(icd_vi_g91_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vi_g91, useNA = "always")

icd_vi_g95_vars <- grep("G95", icd_vi_vars, perl = TRUE, value = T)
d_acc_icd[, icd_vi_g95 := ifelse(Reduce(`|`, lapply(icd_vi_g95_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vi_g95, useNA = "always")

icd_vi_subtype_vars <- c(
  "icd_vi_n06", "icd_vi_atrophies", "icd_vi_movement", "icd_vi_cns", "icd_vi_g35",
  "icd_vi_other", "icd_vi_g40", "icd_vi_g43", "icd_vi_g45", "icd_vi_g47", 
  "icd_vi_mononeuropathy", "icd_vi_polyneuropathy", "icd_vi_myoneural", "icd_vi_g71",
  "icd_vi_paralytic", "icd_vi_g91", "icd_vi_g95"
)

d_acc_icd[, icd_vi_any := ifelse(Reduce(`|`, lapply(icd_vi_subtype_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_vi_any, useNA = "always")


# 12 Diseases of the respiratory system X ----------
icd_x_fo_vars <- grep("icd_x_.*_fo", names(d_acc_icd), value = TRUE)
icd_x_vars <- str_remove(icd_x_fo_vars, "_fo")

icd_x_j30_vars <- grep("J30", icd_x_vars, perl = TRUE, value = T)
d_acc_icd[, icd_x_j30 := ifelse(Reduce(`|`, lapply(icd_x_j30_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_x_j30, useNA = "always")

icd_x_j31_vars <- grep("J31", icd_x_vars, perl = TRUE, value = T)
d_acc_icd[, icd_x_j31 := ifelse(Reduce(`|`, lapply(icd_x_j31_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_x_j31, useNA = "always")

icd_x_j32_vars <- grep("J32", icd_x_vars, perl = TRUE, value = T)
d_acc_icd[, icd_x_j32 := ifelse(Reduce(`|`, lapply(icd_x_j32_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_x_j32, useNA = "always")

icd_x_j37_vars <- grep("J37", icd_x_vars, perl = TRUE, value = T)
d_acc_icd[, icd_x_j37 := ifelse(Reduce(`|`, lapply(icd_x_j37_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_x_j37, useNA = "always")

icd_x_bronchitis_vars <- grep("J4[1-4]", icd_x_vars, perl = TRUE, value = T)
d_acc_icd[, icd_x_bronchitis := ifelse(Reduce(`|`, lapply(icd_x_bronchitis_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_x_bronchitis, useNA = "always")

icd_x_asthma_vars <- grep("J4[5-6]", icd_x_vars, perl = TRUE, value = T)
d_acc_icd[, icd_x_asthma := ifelse(Reduce(`|`, lapply(icd_x_asthma_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_x_asthma, useNA = "always")

icd_x_j47_vars <- grep("J47", icd_x_vars, perl = TRUE, value = T)
d_acc_icd[, icd_x_j47 := ifelse(Reduce(`|`, lapply(icd_x_j47_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_x_j47, useNA = "always")

icd_x_lung_vars <- grep("J6[0-4]|J6[6-8]", icd_x_vars, perl = TRUE, value = T)
d_acc_icd[, icd_x_lung := ifelse(Reduce(`|`, lapply(icd_x_lung_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_x_lung, useNA = "always")

icd_x_j70_j84_vars <- grep("J70|J84", icd_x_vars, perl = TRUE, value = T)
d_acc_icd[, icd_x_j70_j84 := ifelse(Reduce(`|`, lapply(icd_x_j70_j84_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_x_j70_j84, useNA = "always")

icd_x_j85_j86_vars <- grep("J85|J86", icd_x_vars, perl = TRUE, value = T)
d_acc_icd[, icd_x_j85_j86 := ifelse(Reduce(`|`, lapply(icd_x_j85_j86_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_x_j85_j86, useNA = "always")

icd_x_subtype_vars <- c(
  "icd_x_j30", "icd_x_j31", "icd_x_j32", "icd_x_j37", "icd_x_bronchitis", 
  "icd_x_asthma", "icd_x_j47", "icd_x_lung", "icd_x_j70_j84", "icd_x_j85_j86"
)

d_acc_icd[, icd_x_any := ifelse(Reduce(`|`, lapply(icd_x_subtype_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_x_any, useNA = "always")




# 1 Cancer Chapter 2 ---------------------
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

icd_ii_c43_vars <- grep("C43", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c43 := ifelse(Reduce(`|`, lapply(icd_ii_c43_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c43, useNA = "always")

icd_ii_c44_vars <- grep("C44", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c44 := ifelse(Reduce(`|`, lapply(icd_ii_c44_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c44, useNA = "always")

icd_ii_c45_vars <- grep("C45", icd_ii_vars, perl = TRUE, value = T)
d_acc_icd[, icd_ii_c45 := ifelse(Reduce(`|`, lapply(icd_ii_c45_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c45, useNA = "always")

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

icd_ii_subtype_vars <- c(
  "icd_ii_c53", "icd_ii_c56", "icd_ii_c62", "icd_ii_c01_c14", "icd_ii_c15", 
  "icd_ii_c16", "icd_ii_c17", "icd_ii_c18_c21", "icd_ii_c22", "icd_ii_c23_c24", 
  "icd_ii_c25", "icd_ii_c26", "icd_ii_c31_c33", "icd_ii_c34", "icd_ii_c30_c39", 
  "icd_ii_c40_c41", "icd_ii_c44", "icd_ii_c45", "icd_ii_other", "icd_ii_c50", 
  "icd_ii_c54_c55", "icd_ii_c61", "icd_ii_c64", "icd_ii_c65_c66", "icd_ii_c67", 
  "icd_ii_c68", "icd_ii_c69_c72", "icd_ii_c73_c75", "icd_ii_lymphoid" 
)

d_acc_icd[, icd_ii_any := ifelse(Reduce(`|`, lapply(icd_ii_subtype_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_any, useNA = "always")

d_acc_icd[, icd_ii_subtype_n :=  rowSums(d_acc_icd[, icd_ii_subtype_vars, with = FALSE], na.rm = FALSE)]
table(d_acc_icd$icd_ii_subtype_n, useNA = "always")

d_acc_icd[, icd_ii_subtype := NA]

# 1 Gynaecological Cancer
icd_ii_gynaecological_vars <- c("icd_ii_c53", "icd_ii_c56", "icd_ii_c54_c55")
d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c53 == 1, "Cervix", icd_ii_subtype)]
d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c56 == 1, "Ovary", icd_ii_subtype)]
d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c54_c55 == 1, "Uterus", icd_ii_subtype)]

d_acc_icd[, icd_ii_gynaecological := ifelse(Reduce(`|`, lapply(icd_ii_gynaecological_vars, function(v) f1(get(v), 1))),
                                            1, 0)]

# 2 Genitourinary Cancer 
icd_ii_genitourinary_vars <- c("icd_ii_c64", "icd_ii_c65_c66","icd_ii_c67", "icd_ii_c68", "icd_ii_c62")

d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c64 == 1, "Kidney", icd_ii_subtype)]
d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c65_c66 == 1, "Renal Pelvis & Ureter", icd_ii_subtype)]
d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c67 == 1, "Bladder", icd_ii_subtype)]
d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c68 == 1, "Unspecified Urinary Organs", icd_ii_subtype)]
d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c62 == 1, "Testis", icd_ii_subtype)]

d_acc_icd[, icd_ii_genitourinary := ifelse(Reduce(`|`, lapply(icd_ii_genitourinary_vars, function(v) f1(get(v), 1))),
                                           1, 0)]

# 3 Head and Neck 
icd_ii_headneck_vars <- c("icd_ii_c01_c14", "icd_ii_c31_c33")

d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c01_c14 == 1, "Lip, Oral, Cavity & Pharnyx", icd_ii_subtype)]
d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c31_c33 == 1, "Sinuses, Larynx & Trachea ", icd_ii_subtype)]

d_acc_icd[, icd_ii_headneck := ifelse(Reduce(`|`, lapply(icd_ii_headneck_vars, function(v) f1(get(v), 1))),
                                      1, 0)]
# 4 Lung 
icd_ii_bonchuslung_vars <- c("icd_ii_c34")

d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c34 == 1, "Bronchus & Lung", icd_ii_subtype)]

d_acc_icd[, icd_ii_bonchuslung := ifelse(Reduce(`|`, lapply(icd_ii_bonchuslung_vars, function(v) f1(get(v), 1))),
                                         1, 0)]

# 5 Colorectal  
icd_ii_colorectal_vars <- c("icd_ii_c18_c21")

d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c18_c21 == 1, "Colorectual", icd_ii_subtype)]

d_acc_icd[, icd_ii_colorectal := ifelse(Reduce(`|`, lapply(icd_ii_colorectal_vars, function(v) f1(get(v), 1))),
                                        1, 0)]

# 6 GI 
icd_ii_gastrointestinal_vars <- c("icd_ii_c15", "icd_ii_c16", "icd_ii_c17", "icd_ii_c22", "icd_ii_c23_c24", "icd_ii_c25", "icd_ii_c26")

d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c15 == 1, "Oesophagus", icd_ii_subtype)]
d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c16 == 1, "Stomach", icd_ii_subtype)]
d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c17 == 1, "Small Intestine", icd_ii_subtype)]
d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c22 == 1, "Liver & Intraheptic Bile Ducts", icd_ii_subtype)]
d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c23_c24 == 1, "Gallbladder & Biliary Tract", icd_ii_subtype)]
d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c25 == 1, "Pancreas", icd_ii_subtype)]
d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c26 == 1, "Ill-defined Digestive Organs", icd_ii_subtype)]

d_acc_icd[, icd_ii_gastrointestinal := ifelse(Reduce(`|`, lapply(icd_ii_gastrointestinal_vars, function(v) f1(get(v), 1))),
                                              1, 0)]

# 7 Breast 
icd_ii_breast_vars <- c("icd_ii_c50")

d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c50 == 1, "Breast", icd_ii_subtype)]

d_acc_icd[, icd_ii_breast := ifelse(Reduce(`|`, lapply(icd_ii_breast_vars, function(v) f1(get(v), 1))),
                                    1, 0)]

# 8 Prostate 
icd_ii_prostate_vars <- c("icd_ii_c61")

d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c61 == 1, "Prostate", icd_ii_subtype)]

d_acc_icd[, icd_ii_prostate := ifelse(Reduce(`|`, lapply(icd_ii_prostate_vars, function(v) f1(get(v), 1))),
                                      1, 0)]

# 9 Blood 
icd_ii_blood_vars <- c("icd_ii_lymphoid")

d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_lymphoid == 1, "Lymphoid, Haematopoietic & Related Tissue", icd_ii_subtype)]

d_acc_icd[, icd_ii_blood := ifelse(Reduce(`|`, lapply(icd_ii_blood_vars, function(v) f1(get(v), 1))),
                                   1, 0)]

# 10 Melanoma
icd_ii_melanoma_vars <- c("icd_ii_c43")

d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c43 == 1, "Melanoma", icd_ii_subtype)]

d_acc_icd[, icd_ii_melanoma := ifelse(Reduce(`|`, lapply(icd_ii_melanoma_vars, function(v) f1(get(v), 1))),
                                      1, 0)]

# 11 Skin
icd_ii_skin_vars <- c("icd_ii_c44")

d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c44 == 1, "Other Skin", icd_ii_subtype)]

d_acc_icd[, icd_ii_skin := ifelse(Reduce(`|`, lapply(icd_ii_skin_vars, function(v) f1(get(v), 1))),
                                  1, 0)]


# 12 CNS 
icd_ii_cns_vars <- c("icd_ii_c69_c72")

d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c69_c72 == 1, "Central Nervous System", icd_ii_subtype)]

d_acc_icd[, icd_ii_cns := ifelse(Reduce(`|`, lapply(icd_ii_cns_vars, function(v) f1(get(v), 1))),
                                 1, 0)]

# 13 Endocrine Gland
icd_ii_endocrine_vars <- c("icd_ii_c73_c75")

d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c73_c75 == 1, "Endocrine Gland", icd_ii_subtype)]

d_acc_icd[, icd_ii_endocrine := ifelse(Reduce(`|`, lapply(icd_ii_endocrine_vars, function(v) f1(get(v), 1))),
                                       1, 0)]

# 14 Other 
icd_ii_other_vars <- c("icd_ii_c30_c39", "icd_ii_c40_c41", "icd_ii_other", "icd_ii_c45")

d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c30_c39 == 1, "Ill-defined Sites", icd_ii_subtype)]
d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c40_c41 == 1, "Bone & Articular Cartilage", icd_ii_subtype)]
d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_other == 1, "Other Neoplasms", icd_ii_subtype)]
d_acc_icd[, icd_ii_subtype := ifelse(icd_ii_subtype_n == 1 & icd_ii_c45 == 1, "Mesothelioma", icd_ii_subtype)]

d_acc_icd[, icd_ii_other := ifelse(Reduce(`|`, lapply(icd_ii_other_vars, function(v) f1(get(v), 1))),
                                   1, 0)]

table(d_acc_icd$icd_ii_subtype, useNA = "always")

# number of cancer types
icd_ii_type_vars <- c(
  "icd_ii_gynaecological", "icd_ii_genitourinary", "icd_ii_headneck",
  "icd_ii_bonchuslung", "icd_ii_colorectal", "icd_ii_gastrointestinal",
  "icd_ii_breast", "icd_ii_prostate", "icd_ii_blood", "icd_ii_melanoma",
  "icd_ii_skin", "icd_ii_cns", "icd_ii_endocrine", "icd_ii_other"
)

d_acc_icd[, icd_ii_type_n :=  rowSums(d_acc_icd[, icd_ii_type_vars, with = FALSE], na.rm = FALSE)]
table(d_acc_icd$icd_ii_type_n, useNA = "always")

d_acc_icd[, icd_ii_any := NA]
d_acc_icd[, icd_ii_any := ifelse(icd_ii_type_n == 0, "No", icd_ii_any)]
d_acc_icd[, icd_ii_any := ifelse(icd_ii_type_n == 1, "One Primary", icd_ii_any)]
d_acc_icd[, icd_ii_any := ifelse(icd_ii_type_n > 1, "Multiple Primary", icd_ii_any)]

table(d_acc_icd$icd_ii_any, useNA = "always")


d_acc_icd[, icd_ii_type := NA]
d_acc_icd[, icd_ii_type := ifelse((Reduce(`|`, lapply(icd_ii_gynaecological_vars, function(v) f1(get(v), 1)))) &
                                    (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_gynaecological_vars, function(v) f1(get(v), 0)))),
                                  "Gynaecological", icd_ii_type
)]
d_acc_icd[, icd_ii_type := ifelse((Reduce(`|`, lapply(icd_ii_genitourinary_vars, function(v) f1(get(v), 1)))) &
                                    (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_genitourinary_vars, function(v) f1(get(v), 0)))),
                                  "Genitourinary", icd_ii_type
)]

d_acc_icd[, icd_ii_type := ifelse((Reduce(`|`, lapply(icd_ii_headneck_vars, function(v) f1(get(v), 1)))) &
                                    (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_headneck_vars, function(v) f1(get(v), 0)))),
                                  "Head & Neck", icd_ii_type
)]

d_acc_icd[, icd_ii_type := ifelse((Reduce(`|`, lapply(icd_ii_bonchuslung_vars, function(v) f1(get(v), 1)))) &
                                    (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_bonchuslung_vars, function(v) f1(get(v), 0)))),
                                  "Bronchus & Lung", icd_ii_type
)]

d_acc_icd[, icd_ii_type := ifelse((Reduce(`|`, lapply(icd_ii_colorectal_vars, function(v) f1(get(v), 1)))) &
                                    (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_colorectal_vars, function(v) f1(get(v), 0)))),
                                  "Colorectal", icd_ii_type
)]

d_acc_icd[, icd_ii_type := ifelse((Reduce(`|`, lapply(icd_ii_gastrointestinal_vars, function(v) f1(get(v), 1)))) &
                                    (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_gastrointestinal_vars, function(v) f1(get(v), 0)))),
                                  "Gastrointestinal Tract", icd_ii_type
)]

d_acc_icd[, icd_ii_type := ifelse((Reduce(`|`, lapply(icd_ii_breast_vars, function(v) f1(get(v), 1)))) &
                                    (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_breast_vars, function(v) f1(get(v), 0)))),
                                  "Breast", icd_ii_type
)]

d_acc_icd[, icd_ii_type := ifelse((Reduce(`|`, lapply(icd_ii_prostate_vars, function(v) f1(get(v), 1)))) &
                                    (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_prostate_vars, function(v) f1(get(v), 0)))),
                                  "Prostate", icd_ii_type
)]
d_acc_icd[, icd_ii_type := ifelse((Reduce(`|`, lapply(icd_ii_blood_vars, function(v) f1(get(v), 1)))) &
                                    (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_blood_vars, function(v) f1(get(v), 0)))),
                                  "Blood", icd_ii_type
)]

d_acc_icd[, icd_ii_type := ifelse((Reduce(`|`, lapply(icd_ii_melanoma_vars, function(v) f1(get(v), 1)))) &
                                    (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_melanoma_vars, function(v) f1(get(v), 0)))),
                                  "Melanoma", icd_ii_type
)]

d_acc_icd[, icd_ii_type := ifelse((Reduce(`|`, lapply(icd_ii_skin_vars, function(v) f1(get(v), 1)))) &
                                    (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_skin_vars, function(v) f1(get(v), 0)))),
                                  "Other Skin", icd_ii_type
)]
d_acc_icd[, icd_ii_type := ifelse((Reduce(`|`, lapply(icd_ii_cns_vars, function(v) f1(get(v), 1)))) &
                                    (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_cns_vars, function(v) f1(get(v), 0)))),
                                  "Central Nervous System", icd_ii_type
)]

d_acc_icd[, icd_ii_type := ifelse((Reduce(`|`, lapply(icd_ii_endocrine_vars, function(v) f1(get(v), 1)))) &
                                    (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_endocrine_vars, function(v) f1(get(v), 0)))),
                                  "Endocrine Gland", icd_ii_type
)]

d_acc_icd[, icd_ii_type := ifelse((Reduce(`|`, lapply(icd_ii_other_vars, function(v) f1(get(v), 1)))) &
                                    (Reduce(`&`, lapply(icd_ii_subtype_vars %snin% icd_ii_other_vars, function(v) f1(get(v), 0)))),
                                  "Other Cancer", icd_ii_type
)]

# d_acc_icd[, icd_ii_type := ifelse(icd_ii_type_n == 0, "No Cancer", icd_ii_type)]

## Any other conditions at all? ---------------

icd_not_cancer_any_vars <- c(
  "icd_i_any", "icd_ix_any", "icd_v_any",
  "icd_vii_any", "icd_viii_any", "icd_i_iii_any",
  "icd_xi_any", "icd_xiv_any", "icd_iv_xiii_any", "icd_vi_any",
  "icd_x_any"
)

d_acc_icd[, icd_not_cancer := ifelse((Reduce(`|`, lapply(icd_not_cancer_any_vars, function(v) f1(get(v), 1)))),
                                     1, 0)]
table(d_acc_icd$icd_not_cancer, useNA = "always")

d_acc_icd[, icd_ii_comorbid := NA]
d_acc_icd[, icd_ii_comorbid := ifelse(icd_ii_any %in% c("One Primary", "Multiple Primary") & icd_not_cancer == 0, "only cancer", icd_ii_comorbid)]
d_acc_icd[, icd_ii_comorbid := ifelse(icd_ii_any %in% c("One Primary", "Multiple Primary") & icd_not_cancer == 1, "cancer comorbid", icd_ii_comorbid)]
d_acc_icd[, icd_ii_comorbid := ifelse(icd_ii_any == "No" & icd_not_cancer == 1, "only other conditions", icd_ii_comorbid)]
d_acc_icd[, icd_ii_comorbid := ifelse(icd_ii_any == "No" & icd_not_cancer == 0, "healthy", icd_ii_comorbid)]

table(d_acc_icd$icd_ii_comorbid, useNA = "always")

d_acc_icd[, icd_ii_vs_others := NA]
d_acc_icd[, icd_ii_vs_others := ifelse(icd_ii_any %in% c("One Primary", "Multiple Primary"), "cancer", icd_ii_vs_others)]
d_acc_icd[, icd_ii_vs_others := ifelse(icd_ii_any == "No" & icd_not_cancer == 1, "only other conditions", icd_ii_vs_others)]
d_acc_icd[, icd_ii_vs_others := ifelse(icd_ii_any == "No" & icd_not_cancer == 0, "healthy", icd_ii_vs_others)]

table(d_acc_icd$icd_ii_vs_others, useNA = "always")

d_acc_icd[, icd_ii_type := ifelse(icd_ii_comorbid == "healthy", "Healthy", icd_ii_type)]
d_acc_icd[, icd_ii_type := ifelse(icd_ii_type_n > 1, "Multiple Primary", icd_ii_type)]

table(d_acc_icd$icd_ii_type, useNA = "always")

# excluding other conds
d_acc_icd_ii <- d_acc_icd[!is.na(icd_ii_type)]

# make factors
d_acc_icd_ii[, icd_ii_type := as.factor(icd_ii_type)]
d_acc_icd_ii[, icd_ii_type := relevel(icd_ii_type, ref = "Healthy")]

table(d_acc_icd_ii$icd_ii_type, useNA = "always")

# first occurrence of cancer diagnosis if any
d_acc_icd_ii[, icd_ii_fo := do.call(pmin, c(.SD, list(na.rm = TRUE))), .SDcols = icd_ii_fo_vars]

# last occurrence of cancer diagnosis if any
d_acc_icd_ii[, icd_ii_lo := do.call(pmax, c(.SD, list(na.rm = TRUE))), .SDcols = icd_ii_fo_vars]

d_acc_icd[, icd_ii_before_acc := ifelse(icd_ii_lo <= acc_startdate, 1, icd_ii_before_acc)]
# d_acc_icd[, icd_ii_before_acc := ifelse(icd_ii_fo > acc_startdate, 0, icd_ii_before_acc)]

# exclude those had cancer after
# d_acc_icd_ii_before_acc <- d_acc_icd_ii[icd_ii_fo < acc_startdate & !is.na(icd_ii_type)]
# table(d_acc_icd_ii_before_acc$icd_ii_type, useNA = "always")

# descriptives
egltable(c("icd_ii_type",
           "sleep_comp", "mvpa_comp", "lpa_comp", "sb_comp", #acomp and 0 imputed
           "sleep", "mvpa", "lpa", "sb" # raw
), strict = FALSE, data = d_acc_icd_ii)


egltable(c(
           "sleep_comp", "mvpa_comp", "lpa_comp", "sb_comp", #acomp and 0 imputed
           "sleep", "mvpa", "lpa", "sb" # raw
), strict = FALSE, g = "icd_ii_subtype", data = d_acc_icd_ii)


psych::describe(d_acc_icd_ii[icd_ii_subtype == "Bladder"]$sleep)
