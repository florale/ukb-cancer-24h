source("ukb-cancer-24h-data.R")

fit_quantile_cancer <- brmcoda(clr_cancer_acc,
                               bf(mvbind(ilr1, ilr2, ilr3) ~ cancer_before_acc,
                                  quantile = 0.5),
                               # save_pars = save_pars(all = TRUE),
                               warmup = 500, chains = 4, cores = 4, backend = "cmdstanr",
                               family = asym_laplace()
)
saveRDS(fit_quantile_cancer, paste0(outputdir, "fit_quantile_cancer", ".RDS"))
