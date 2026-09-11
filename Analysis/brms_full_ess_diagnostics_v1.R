suppressPackageStartupMessages({library(brms); library(posterior)})
OUT <- "Analysis/outputs"
models <- c("ADAS_pT217","MMSE_pT217","Hippo_pT217","ADAS_GFAP","MMSE_GFAP","Hippo_GFAP")
for (m in models) {
  fit <- readRDS(file.path(OUT, paste0("brmsfit_", m, ".rds")))
  ess <- summarise_draws(as_draws_array(fit), "ess_bulk", "ess_tail")
  cat(sprintf("%-14s min ess_bulk (all params) = %.0f | min ess_tail (all params) = %.0f\n",
              m, min(ess$ess_bulk, na.rm=TRUE), min(ess$ess_tail, na.rm=TRUE)))
  rm(fit); gc(verbose=FALSE)
}
