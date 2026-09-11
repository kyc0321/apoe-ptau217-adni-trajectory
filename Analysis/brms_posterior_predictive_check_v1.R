suppressPackageStartupMessages({library(brms)})
OUT <- "Analysis/outputs"
models <- c("ADAS_pT217","MMSE_pT217","Hippo_pT217","ADAS_GFAP","MMSE_GFAP","Hippo_GFAP")
for (m in models) {
  fit <- readRDS(file.path(OUT, paste0("brmsfit_", m, ".rds")))
  yrep_bm <- posterior_predict(fit, resp = "bm", ndraws = 500)
  y_bm <- get_y(fit, resp = "bm")
  yrep_oc <- posterior_predict(fit, resp = "outc", ndraws = 500)
  y_oc <- get_y(fit, resp = "outc")
  bp_mean_bm <- mean(apply(yrep_bm, 1, mean) >= mean(y_bm, na.rm=TRUE))
  bp_sd_bm   <- mean(apply(yrep_bm, 1, sd)   >= sd(y_bm, na.rm=TRUE))
  bp_mean_oc <- mean(apply(yrep_oc, 1, mean) >= mean(y_oc, na.rm=TRUE))
  bp_sd_oc   <- mean(apply(yrep_oc, 1, sd)   >= sd(y_oc, na.rm=TRUE))
  cat(sprintf("%-14s bm: Bayes-p(mean)=%.2f Bayes-p(sd)=%.2f | outc: Bayes-p(mean)=%.2f Bayes-p(sd)=%.2f\n",
              m, bp_mean_bm, bp_sd_bm, bp_mean_oc, bp_sd_oc))
  rm(fit, yrep_bm, yrep_oc); gc(verbose=FALSE)
}
