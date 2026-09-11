suppressPackageStartupMessages({library(brms)})
OUT <- "Analysis/outputs"
models <- c("ADAS_pT217","MMSE_pT217","Hippo_pT217","ADAS_GFAP","MMSE_GFAP","Hippo_GFAP")

cat("=== Priors (from first model, brms defaults are formula/data-dependent but consistent across these 6) ===\n")
fit1 <- readRDS(file.path(OUT, paste0("brmsfit_", models[1], ".rds")))
print(prior_summary(fit1))

cat("\n=== Per-model diagnostics ===\n")
for (m in models) {
  fp <- file.path(OUT, paste0("brmsfit_", m, ".rds"))
  if (!file.exists(fp)) { cat(m, ": file not found\n"); next }
  fit <- readRDS(fp)
  rhats <- brms::rhat(fit)
  ess_b <- brms::neff_ratio(fit) * (ndraws(fit))  # approx bulk ESS via neff_ratio*ndraws is rough; use posterior pkg instead below
  np <- nuts_params(fit)
  n_div <- sum(np$Parameter == "divergent__" & np$Value == 1)
  n_total <- sum(np$Parameter == "divergent__")
  smry <- as.data.frame(summary(fit)$fixed)
  ess_bulk_range <- range(smry$Bulk_ESS, na.rm=TRUE)
  ess_tail_range <- range(smry$Tail_ESS, na.rm=TRUE)
  rhat_range <- range(smry$Rhat, na.rm=TRUE)
  cat(sprintf("%-14s Rhat range [%.3f, %.3f] | Bulk_ESS range [%.0f, %.0f] | Tail_ESS range [%.0f, %.0f] | divergences %d/%d\n",
              m, rhat_range[1], rhat_range[2], ess_bulk_range[1], ess_bulk_range[2],
              ess_tail_range[1], ess_tail_range[2], n_div, n_total))
  rm(fit); gc(verbose=FALSE)
}
