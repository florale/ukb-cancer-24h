source("ukb-cancer-24h-utils.R")
source(paste0(redir, "data/data_acc_icd_ii.R"))

# exclude those had cancer after
d_acc_icd_ii_before_acc <- d_acc_icd[icd_ii_before_acc != 0 | is.na(icd_ii_before_acc)]

# icd ii comorbidity with other long term conditions
d_acc_icd_ii_before_acc[, icd_ii_comorbid := NA]
d_acc_icd_ii_before_acc[, icd_ii_comorbid := ifelse(ltc == 0,  "Healthy", icd_ii_comorbid)]
d_acc_icd_ii_before_acc[, icd_ii_comorbid := ifelse(icd_ii_n == 0 & ltc == 1, "No ICD II, other ICD only", icd_ii_comorbid)]

d_acc_icd_ii_before_acc[, icd_ii_comorbid := ifelse(icd_ii_n > 0 & ltc == 1 & Reduce(`&`, lapply(grep(paste0(ltc_vars[ltc_vars %nin% ltc_icd_ii_vars], ".*fo", collapse = "|"), names(d_icd_all), value = T), function(v)
  is.na(get(v)))), 
  "ICD II only", 
  icd_ii_comorbid)]

d_acc_icd_ii_before_acc[, icd_ii_comorbid := ifelse(icd_ii_n > 0 & ltc == 1 & Reduce(`|`, lapply(grep(paste0(ltc_vars[ltc_vars %nin% ltc_icd_ii_vars], ".*fo", collapse = "|"), names(d_icd_all), value = T), function(v)
  !is.na(get(v)))), 
  "ICD II comorbidity", 
  icd_ii_comorbid)]

table(d_acc_icd_ii_before_acc$icd_ii_comorbid, useNA = "always")

d_acc_icd_ii_before_acc[, icd_ii_comorbid := factor(icd_ii_comorbid, ordered = TRUE,
                                                    levels = c("Healthy",
                                                               "ICD II only",
                                                               "ICD II comorbidity",
                                                               "No ICD II, other ICD only"))]

egltable(c("age", "sex", "ethnicg", "white", "bmi", "edu", "working", "deprivation", 
           "current_drinker", "never_smoked", "smoking", "alcohol",
           
           "icd_ii_comorbid",
           
           "sleep_comp", "mvpa_comp", "lpa_comp", "sb_comp", #acomp and 0 imputed
           "sleep", "mvpa", "lpa", "sb" # raw
), strict = FALSE, data = d_acc_icd_ii_before_acc)


clr_acc_icd_ii <- complr(data = d_acc_icd_ii_before_acc[, -colnames(ilr_acc), with = FALSE],
                         transform = "ilr",
                         parts = c("sleep_comp", "mvpa_comp", "lpa_comp", "sb_comp"),
                         sbp = sbp,
                         total = 1440,
                         idvar = "eid",
                         shape = "wide")
saveRDS(clr_acc_icd_ii, paste0(inputdir, "d_acc_icd_ii_before_acc", ".RDS"))