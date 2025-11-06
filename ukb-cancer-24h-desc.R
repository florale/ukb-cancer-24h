source("ukb-cancer-24h-utils.R")
source(paste0(redir, "ukb_utils.R"))
source("ukb-cancer-24h-data.R")

# d_cancer_acc_desc <- d_cancer_acc[, c("eid", "age", "age_at_acc", "sex", "white", "edu", "working", 
#                                     "smoking", "alcohol", "deprivation", "deprivationg", 
#                                     "cancer_other_before_acc",
#                                     "cancer_time_since_diag_other",
#                                     "cancer_before_acc_type_other", 
#                                     "icd_ii_time_since_lo",
#                                     "icd_not_cancer",
#                                     "sleep_comp", "mvpa_comp", "lpa_comp", "sb_comp",
#                                     "sleep", "mvpa", "lpa", "sb")]
d_cancer_acc_desc <- d_cancer_acc[complete.cases(d_cancer_acc[, c("eid", "age", "age_at_acc", "sex", "white", "edu", "working", 
                                                                       "smoking", "alcohol", "deprivation", "deprivationg", 
                                                                       "cancer_other_before_acc",
                                                                       "cancer_time_since_diag_other",
                                                                       "cancer_before_acc_type_other", 
                                                                       "icd_ii_time_since_lo",
                                                                       "icd_not_cancer",
                                                                       "sleep_comp", "mvpa_comp", "lpa_comp", "sb_comp",
                                                                       "sleep", "mvpa", "lpa", "sb")])]
nrow(d_cancer_acc_desc)

# exclude bc incidences up to 1y followup
nrow(d_acc_icd) - nrow(d_cancer_acc) #94471 - 93490 #981

# exclude due to missing covariates
nrow(d_cancer_acc) - nrow(d_cancer_acc_desc) #93490 - 91352 #2138

# descriptives ----------------------------
## demographics - group by cancer vs healthy
egltable(c(
  "age", "age_at_acc", "sex", "white", "bmig", "edu", "working", 
  "smoking", "alcohol", "deprivation", "deprivationg", 
  "cancer_other_before_acc",
  "cancer_time_since_diag_other",
  "cancer_before_acc_type_other", 
  "icd_ii_time_since_lo",
  "sleep", "mvpa", "lpa", "sb"
),
strict = FALSE, data = d_cancer_acc_desc)

egltable(c(
  "age", "age_at_acc", "sex", "white", "bmig", "edu", "working", 
  "smoking", "alcohol", "deprivation", "deprivationg", 
  "cancer_other_before_acc",
  "icd_ii_time_since_lo"
),
strict = FALSE, g = "cancer_other_before_acc", data = d_cancer_acc_desc)

## main variables - group by cancer vs healthy
egltable(c(
  "sleep_comp", "mvpa_comp", "lpa_comp", "sb_comp", #acomp and 0 imputed
  "sleep", "mvpa", "lpa", "sb" # raw
), strict = FALSE, g = "cancer_other_before_acc", data = d_cancer_acc_desc)

# ## main variables - group by cancer types
# egltable(c(
#   "sleep_comp", "mvpa_comp", "lpa_comp", "sb_comp", #acomp and 0 imputed
#   "sleep", "mvpa", "lpa", "sb" # raw
# ), strict = FALSE, g = "cancer_before_acc_type_other", data = d_cancer_acc_desc)

# egltable(c("cancer_before_acc_type_other"), strict = FALSE, data = d_cancer_acc_desc)

egltable(c(
  "cancer_time_since_diag_other",
  "icd_ii_time_since_lo", "cancer_other_before_acc", "cancer_before_acc_type_other",
  "sleep", "mvpa", "lpa", "sb" # raw
),
strict = FALSE, data = d_cancer_acc_desc[cancer_other_before_acc %nin% c("Healthy", "Others")])

#  median time since cancer diag by types 
psych::describeBy(d_cancer_acc_desc$icd_ii_time_since_lo, group = d_cancer_acc_desc$cancer_before_acc_type_other)

# number of other conditions
d_acc_icd_health <- readRDS(paste0(inputdir, "d_acc_icd_health", ".RDS"))

icd_at_acc <- grep("_at_acc$", names(d_acc_icd_health), value = T)
icd_at_acc <- icd_at_acc[icd_at_acc %nin% c("one_icd_type_at_acc", "other_conds_at_acc", "icd_at_acc", "age_at_acc", "icd_not_ii_at_acc")]
icd_at_acc_not_ii_vars <- icd_at_acc[icd_at_acc %nin% c("icd_ii_at_acc")]

d_acc_icd_health[, icd_at_acc_n :=  rowSums(d_acc_icd_health[, icd_at_acc, with = FALSE] != "Healthy", na.rm = TRUE)]
table(d_acc_icd_health$icd_at_acc_n, useNA = "always")

d_acc_icd_health[, icd_not_cancer_n := rowSums(d_acc_icd_health[, icd_at_acc_not_ii_vars, with = FALSE] != "Healthy", na.rm = TRUE)]
table(d_acc_icd_health$icd_not_cancer_n, useNA = "always")

d_acc_icd_health[, icd_ii_at_acc_other := NA]
d_acc_icd_health[, icd_ii_at_acc_other := ifelse(icd_ii_at_acc == "Healthy", "Healthy", icd_ii_at_acc_other)]
d_acc_icd_health[, icd_ii_at_acc_other := ifelse(icd_ii_at_acc == "ICD II", "Cancer", icd_ii_at_acc_other)]
d_acc_icd_health[, icd_ii_at_acc_other := ifelse(is.na(icd_ii_at_acc) & is.na(icd_ii_at_acc_other), "Others", icd_ii_at_acc_other)]
table(d_acc_icd_health$icd_ii_at_acc_other, useNA = "always")

d_cancer_acc_desc <- merge(d_cancer_acc_desc, d_acc_icd_health[, .(eid, icd_ii_at_acc_other, icd_not_cancer_n)], by = "eid", all.x = TRUE)
psych::describeBy(d_cancer_acc_desc$icd_not_cancer_n, group = d_cancer_acc_desc$icd_ii_at_acc_other)
# psych::describeBy(d_acc_icd_health$icd_not_cancer_n, group = d_acc_icd_health$icd_ii_at_acc_other)

# other cond
d_cancer_acc_desc[, other_conds_at_acc := ifelse(cancer_other_before_acc == "Healthy", 0, icd_not_cancer)]
table(d_cancer_acc_desc$other_conds_at_acc, useNA = "always")

75859 - 66403

# dag -------------------
coords <- tibble::tribble(
  ~name,           ~x,  ~y,
  "cancer",         1,   0,
  "behaviours",     2,   0,
  "demographics",   1,   1,
  "lifestyle",      1.25,   1,
  "bmi",            1.5,   0.25
)

dag <-  ggdag::dagify(behaviours ~ cancer,
                      behaviours ~ demographics + lifestyle + bmi,
                      cancer ~ bmi + demographics + lifestyle + bmi,
                      # lifestyle ~ demographics,
                      bmi ~ cancer + demographics + lifestyle,
                      exposure = "cancer",
                      outcome = "behaviours",
                      latent = c("bmi"),
                      
                      labels = c("cancer" = "Cancer",
                                 "behaviours" = "24h behaviours",
                                 "demographics" = "Demographic factors\n(age, sex, deprivation,\neducation, employment)",
                                 "lifestyle" = "Lifestyle factors\n(smoking, alcohol)",
                                 # "education" = "Education",
                                 # "employment" = "Employment",
                                 "bmi" = "BMI"),
                      
                      # labels = c("cancer" = "Exposure",
                      #            "behaviours" = "Outcome",
                      #            "demographics" = "Confounders",
                      #            "lifestyle" = "Confounders",
                      #            # "education" = "Education",
                      #            # "employment" = "Employment",
                      #            "bmi" = "Mechanism"),
                      coords = coords
)

plot_dag <- ggdag::ggdag_parents(dag, "cancer",
                                 text = FALSE
)
plot_dag[["data"]]$parent <- NA
plot_dag[["data"]]$parent <- ifelse(plot_dag[["data"]]$name == "behaviours", "Outcome", plot_dag[["data"]]$parent)
plot_dag[["data"]]$parent <- ifelse(plot_dag[["data"]]$name == "cancer", "Exposure", plot_dag[["data"]]$parent)
plot_dag[["data"]]$parent <- ifelse(plot_dag[["data"]]$name == "bmi", "Mechanism (Unadjusted)", plot_dag[["data"]]$parent)
plot_dag[["data"]]$parent <- ifelse(plot_dag[["data"]]$name == "demographics", "Confounder (Adjusted)", plot_dag[["data"]]$parent)
plot_dag[["data"]]$parent <- ifelse(plot_dag[["data"]]$name == "lifestyle", "Confounder (Adjusted)", plot_dag[["data"]]$parent)
plot_dag[["data"]]$parent <- factor(plot_dag[["data"]]$parent, ordered = TRUE,
                                    levels = c(
                                      "Exposure",
                                      "Outcome",
                                      "Confounder (Adjusted)",
                                      "Mechanism (Unadjusted)"
                                    ))

(plot_dag <- plot_dag +
  geom_dag_point(aes(color = parent)) +  # Adjust node colors
  geom_dag_edges(
    arrow_directed = grid::arrow(length = grid::unit(2, "pt"), type = "closed"),
    arrow_bidirected = grid::arrow(length = grid::unit(2, "pt"), ends = "both", type = "closed")) +
  scale_color_manual(values = c(
    "Exposure" = "#978787", # DCD5CE
    "Outcome" = "#9DB3A8",
    "Confounder (Adjusted)" = "#C99696", #EAD3BF
    "Mechanism (Unadjusted)" = "#EFE3E0"
  ), drop = TRUE, name = NULL)+
  geom_label(aes(label = label),
             color = "black",
             vjust = 0,
             nudge_x = 0, nudge_y = 0.075,
             family = "Arial Narrow", size = 4,
             show.legend = NA) +

  theme_void() +
  theme(
    legend.position     = "bottom"
  ))

saveRDS(plot_dag, paste0(outputdir, "plot_dag", ".RDS"))

grDevices::png(
  file = paste0(outputdir, "dag", ".png"),
  width = 8500,
  height = 7000,
  res = 800
)
plot_dag
dev.off()

# raw compostion by groups -------------------
## sex
print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(mvpa, probs = 0.50),
  quant25 = quantile(mvpa, probs = 0.25),
  quant75 = quantile(mvpa, probs = 0.75)
), by = sex], digits = 2)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(lpa, probs = 0.50),
  quant25 = quantile(lpa, probs = 0.25),
  quant75 = quantile(lpa, probs = 0.75)
), by = sex], digits = 3)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(sb, probs = 0.50),
  quant25 = quantile(sb, probs = 0.25),
  quant75 = quantile(sb, probs = 0.75)
), by = sex], digits = 3)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(sleep, probs = 0.50),
  quant25 = quantile(sleep, probs = 0.25),
  quant75 = quantile(sleep, probs = 0.75)
), by = sex], digits = 3)

## white
print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(mvpa, probs = 0.50),
  quant25 = quantile(mvpa, probs = 0.25),
  quant75 = quantile(mvpa, probs = 0.75)
), by = white], digits = 2)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(lpa, probs = 0.50),
  quant25 = quantile(lpa, probs = 0.25),
  quant75 = quantile(lpa, probs = 0.75)
), by = white], digits = 3)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(sb, probs = 0.50),
  quant25 = quantile(sb, probs = 0.25),
  quant75 = quantile(sb, probs = 0.75)
), by = white], digits = 3)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(sleep, probs = 0.50),
  quant25 = quantile(sleep, probs = 0.25),
  quant75 = quantile(sleep, probs = 0.75)
), by = white], digits = 3)

## bmi
print(d_cancer_acc[!is.na(bmig) & cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(mvpa, probs = 0.50),
  quant25 = quantile(mvpa, probs = 0.25),
  quant75 = quantile(mvpa, probs = 0.75)
), by = bmig], digits = 2)

print(d_cancer_acc[!is.na(bmig) & cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(lpa, probs = 0.50),
  quant25 = quantile(lpa, probs = 0.25),
  quant75 = quantile(lpa, probs = 0.75)
), by = bmig], digits = 3)

print(d_cancer_acc[!is.na(bmig) & cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(sb, probs = 0.50),
  quant25 = quantile(sb, probs = 0.25),
  quant75 = quantile(sb, probs = 0.75)
), by = bmig], digits = 3)

print(d_cancer_acc[!is.na(bmig) & cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(sleep, probs = 0.50),
  quant25 = quantile(sleep, probs = 0.25),
  quant75 = quantile(sleep, probs = 0.75)
), by = bmig], digits = 3)

## edu
print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(mvpa, probs = 0.50),
  quant25 = quantile(mvpa, probs = 0.25),
  quant75 = quantile(mvpa, probs = 0.75)
), by = edu], digits = 2)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(lpa, probs = 0.50),
  quant25 = quantile(lpa, probs = 0.25),
  quant75 = quantile(lpa, probs = 0.75)
), by = edu], digits = 3)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(sb, probs = 0.50),
  quant25 = quantile(sb, probs = 0.25),
  quant75 = quantile(sb, probs = 0.75)
), by = edu], digits = 3)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(sleep, probs = 0.50),
  quant25 = quantile(sleep, probs = 0.25),
  quant75 = quantile(sleep, probs = 0.75)
), by = edu], digits = 3)

## working
print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(mvpa, probs = 0.50),
  quant25 = quantile(mvpa, probs = 0.25),
  quant75 = quantile(mvpa, probs = 0.75)
), by = working], digits = 2)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(lpa, probs = 0.50),
  quant25 = quantile(lpa, probs = 0.25),
  quant75 = quantile(lpa, probs = 0.75)
), by = working], digits = 3)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(sb, probs = 0.50),
  quant25 = quantile(sb, probs = 0.25),
  quant75 = quantile(sb, probs = 0.75)
), by = working], digits = 3)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(sleep, probs = 0.50),
  quant25 = quantile(sleep, probs = 0.25),
  quant75 = quantile(sleep, probs = 0.75)
), by = working], digits = 3)

## deprivation
print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(mvpa, probs = 0.50),
  quant25 = quantile(mvpa, probs = 0.25),
  quant75 = quantile(mvpa, probs = 0.75)
), by = deprivationg], digits = 2)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(lpa, probs = 0.50),
  quant25 = quantile(lpa, probs = 0.25),
  quant75 = quantile(lpa, probs = 0.75)
), by = deprivationg], digits = 3)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(sb, probs = 0.50),
  quant25 = quantile(sb, probs = 0.25),
  quant75 = quantile(sb, probs = 0.75)
), by = deprivationg], digits = 3)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(sleep, probs = 0.50),
  quant25 = quantile(sleep, probs = 0.25),
  quant75 = quantile(sleep, probs = 0.75)
), by = deprivationg], digits = 3)

## smoking
print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(mvpa, probs = 0.50),
  quant25 = quantile(mvpa, probs = 0.25),
  quant75 = quantile(mvpa, probs = 0.75)
), by = smoking], digits = 1)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(lpa, probs = 0.50),
  quant25 = quantile(lpa, probs = 0.25),
  quant75 = quantile(lpa, probs = 0.75)
), by = smoking], digits = 3)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(sb, probs = 0.50),
  quant25 = quantile(sb, probs = 0.25),
  quant75 = quantile(sb, probs = 0.75)
), by = smoking], digits = 3)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(sleep, probs = 0.50),
  quant25 = quantile(sleep, probs = 0.25),
  quant75 = quantile(sleep, probs = 0.75)
), by = smoking], digits = 3)

## age group
print(d_cancer_acc_desc[age_at_acc >= 65 & cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(mvpa, probs = 0.50),
  quant25 = quantile(mvpa, probs = 0.25),
  quant75 = quantile(mvpa, probs = 0.75)
)], digits = 1)

print(d_cancer_acc_desc[age_at_acc >= 65 & cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(lpa, probs = 0.50),
  quant25 = quantile(lpa, probs = 0.25),
  quant75 = quantile(lpa, probs = 0.75)
)], digits = 3)

print(d_cancer_acc_desc[age_at_acc >= 65 & cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(sb, probs = 0.50),
  quant25 = quantile(sb, probs = 0.25),
  quant75 = quantile(sb, probs = 0.75)
)], digits = 3)

print(d_cancer_acc_desc[age_at_acc >= 65 & cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(sleep, probs = 0.50),
  quant25 = quantile(sleep, probs = 0.25),
  quant75 = quantile(sleep, probs = 0.75)
)], digits = 3)

print(d_cancer_acc_desc[age_at_acc < 65 & cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(mvpa, probs = 0.50),
  quant25 = quantile(mvpa, probs = 0.25),
  quant75 = quantile(mvpa, probs = 0.75)
)], digits = 1)

print(d_cancer_acc_desc[age_at_acc < 65 & cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(lpa, probs = 0.50),
  quant25 = quantile(lpa, probs = 0.25),
  quant75 = quantile(lpa, probs = 0.75)
)], digits = 3)

print(d_cancer_acc_desc[age_at_acc < 65 & cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(sb, probs = 0.50),
  quant25 = quantile(sb, probs = 0.25),
  quant75 = quantile(sb, probs = 0.75)
)], digits = 3)

print(d_cancer_acc_desc[age_at_acc < 65 & cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(sleep, probs = 0.50),
  quant25 = quantile(sleep, probs = 0.25),
  quant75 = quantile(sleep, probs = 0.75)
)], digits = 3)

## time since diag
print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(mvpa, probs = 0.50),
  quant25 = quantile(mvpa, probs = 0.25),
  quant75 = quantile(mvpa, probs = 0.75)
), by = cancer_time_since_diag_other], digits = 1)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(lpa, probs = 0.50),
  quant25 = quantile(lpa, probs = 0.25),
  quant75 = quantile(lpa, probs = 0.75)
), by = cancer_time_since_diag_other], digits = 3)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(sb, probs = 0.50),
  quant25 = quantile(sb, probs = 0.25),
  quant75 = quantile(sb, probs = 0.75)
), by = cancer_time_since_diag_other], digits = 3)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(sleep, probs = 0.50),
  quant25 = quantile(sleep, probs = 0.25),
  quant75 = quantile(sleep, probs = 0.75)
), by = cancer_time_since_diag_other], digits = 3)

## type
print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(mvpa, probs = 0.50),
  quant25 = quantile(mvpa, probs = 0.25),
  quant75 = quantile(mvpa, probs = 0.75)
), by = cancer_before_acc_type_other], digits = 1)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(lpa, probs = 0.50),
  quant25 = quantile(lpa, probs = 0.25),
  quant75 = quantile(lpa, probs = 0.75)
), by = cancer_before_acc_type_other], digits = 3)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(sb, probs = 0.50),
  quant25 = quantile(sb, probs = 0.25),
  quant75 = quantile(sb, probs = 0.75)
), by = cancer_before_acc_type_other], digits = 3)

print(d_cancer_acc_desc[cancer_other_before_acc == "Cancer", .(
  quant50 = quantile(sleep, probs = 0.50),
  quant25 = quantile(sleep, probs = 0.25),
  quant75 = quantile(sleep, probs = 0.75)
), by = cancer_before_acc_type_other], digits = 3)

# nobs
egltable(c(
  "sex", "white", "edu", "working", 
  "smoking", "alcohol", "deprivationg"
),
strict = FALSE, data = d_cancer_acc_desc[cancer_other_before_acc == "Cancer"])

egltable(c("bmig"
),
strict = FALSE, data = d_cancer_acc[cancer_other_before_acc == "Cancer"])

nrow(d_cancer_acc_desc[age_at_acc < 65 & cancer_other_before_acc == "Cancer"])
nrow(d_cancer_acc_desc[age_at_acc < 65 & cancer_other_before_acc == "Cancer"])/10152*100
nrow(d_cancer_acc_desc[age_at_acc >= 65 & cancer_other_before_acc == "Cancer"])
nrow(d_cancer_acc_desc[age_at_acc >= 65 & cancer_other_before_acc == "Cancer"])/10152*100
