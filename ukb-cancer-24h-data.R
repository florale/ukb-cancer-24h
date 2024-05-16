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

icd_ii_fo_vars <- grep("icd_ii_.*_fo", names(d_acc_icd), value = TRUE)
icd_ii_vars <- str_remove(icd_ii_fo_vars, "_fo")

## code into main categories --------
icd_ii_c53_vars <- grep("C53", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c56_vars <- grep("C56", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c62_vars <- grep("C62", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c01_c14_vars <- grep("C0[1-9]|C1[0-4]", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c15_vars <- grep("C15", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c16_vars <- grep("C16", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c17_vars <- grep("C17", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c18_c21_vars <- grep("C1[8-9]|C21", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c22_vars <- grep("C22", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c23_c24_vars <- grep("C2[3-4]", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c25_vars <- grep("C25", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c31_c33_vars <- grep("C3[1-3]", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c34_vars <- grep("C34", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c30_c39_vars <- grep("C30|C3[7-9]", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c40_c41_vars <- grep("C4[0-1]", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c43_vars <- grep("C43", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c44_vars <- grep("C44", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c44_c45_vars <- grep("C4[4-5]", icd_ii_vars, perl = TRUE, value = T)
icd_ii_other_vars <- grep("C4[6-9]|C5[1-2|7-8]|C60|C63|C76|C80|C97", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c50_vars <- grep("C50", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c54_c55_vars <- grep("C5[4-5]", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c61_vars <- grep("C61", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c64_vars <- grep("C64", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c65_c66_vars <- grep("C6[5-6]", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c67_vars <- grep("C67", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c68_vars <- grep("C68", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c69_c72_vars <- grep("C69|C7[0-2]", icd_ii_vars, perl = TRUE, value = T)
icd_ii_c73_c75_vars <- grep("C7[3-5]", icd_ii_vars, perl = TRUE, value = T)
icd_ii_lymphoid_vars <- grep("C8[1-6]|C88|C9[0-6]", icd_ii_vars, perl = TRUE, value = T)

d_acc_icd[ , icd_ii_c01_c14 := ifelse(Reduce(`|`, lapply(icd_ii_c01_c14_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
table(d_acc_icd$icd_ii_c01_c14, useNA = "always")

d_acc_icd[ , icd_ii_c01_c14 := ifelse(Reduce(`|`, lapply(icd_ii_c01_c14_vars, function(v)
  f1(get(v), 1))),
  1, 0)]
