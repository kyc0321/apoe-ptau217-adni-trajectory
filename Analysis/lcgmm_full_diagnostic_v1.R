# LCGMM full diagnostic — entropy, mean posterior class probability, bootstrap LRT
# Author: Kuan-Yu Chu | 2026-05-26
# Issue #8: original table reports BIC only; reviewer-grade reporting needs
#   mean posterior class probability per class AND parametric bootstrap LRT.

suppressPackageStartupMessages({ library(dplyr); library(lcmm) })
set.seed(20260525)

OUT <- "/Users/kuan-yuchu/Desktop/LOAD/Analysis/outputs"
biom_long <- utils::read.csv(file.path(OUT, "analytic_biom_long.csv"), stringsAsFactors = FALSE)

# Compute entropy and BLRT manually since lcmm doesn't expose them directly
calc_entropy <- function(pp) {
  # pp = matrix N x K of posterior class probabilities
  K <- ncol(pp); N <- nrow(pp)
  H <- -sum(pp * log(pp + 1e-12)) / N    # per-subject avg entropy
  Hmax <- log(K)
  if (K == 1) return(1)
  1 - H / Hmax
}

fit_lcgmm <- function(yvar, K_max = 4) {
  dat <- biom_long %>% filter(!is.na(.data[[yvar]])) %>% as.data.frame()
  res <- list()
  cat("\n=== LCGMM diagnostic:", yvar, "===\n")
  # k = 1 baseline
  m1 <- hlme(as.formula(paste(yvar, "~ time_yrs")),
             random = ~time_yrs, subject = "RID", data = dat, ng = 1)
  res[[1]] <- list(k = 1, m = m1, ll = m1$loglik, AIC = m1$AIC, BIC = m1$BIC,
                   entropy = 1, min_pct = 100, mean_pp = 1)
  cat("  k=1: loglik =", round(m1$loglik, 1), "  BIC =", round(m1$BIC, 1), "\n")

  for (k in 2:K_max) {
    mk <- tryCatch(hlme(as.formula(paste(yvar, "~ time_yrs")),
                       mixture = ~time_yrs, random = ~time_yrs,
                       subject = "RID", data = dat, ng = k, B = m1),
                   error = function(e) NULL)
    if (is.null(mk)) { cat("  k=", k, " failed\n"); next }
    pp <- as.matrix(mk$pprob[, -c(1,2)])
    cls <- mk$pprob$class
    cls_pct <- 100 * as.numeric(table(cls)) / length(cls)
    # Mean posterior class probability per class (subjects assigned to that class)
    mean_pp_per_class <- sapply(seq_len(k), function(c) {
      idx <- which(cls == c)
      if (length(idx) == 0) NA else mean(pp[idx, c])
    })
    entropy <- calc_entropy(pp)
    res[[k]] <- list(k = k, m = mk, ll = mk$loglik, AIC = mk$AIC, BIC = mk$BIC,
                     entropy = entropy, min_pct = min(cls_pct),
                     mean_pp = mean(mean_pp_per_class, na.rm = TRUE),
                     mean_pp_per_class = mean_pp_per_class)
    cat("  k=", k, ": BIC =", round(mk$BIC, 1),
        " entropy =", round(entropy, 3),
        " min class % =", round(min(cls_pct), 1),
        " mean post prob =", round(mean(mean_pp_per_class, na.rm=TRUE), 3), "\n")
  }
  res
}

cat("Fitting LCGMM for log_gfap and log_ptau...\n")
res_g <- fit_lcgmm("log_gfap")
res_p <- fit_lcgmm("log_ptau")

# Build summary table
mk_row <- function(bm, r) {
  data.frame(biomarker = bm, k = r$k,
             loglik = r$ll, AIC = r$AIC, BIC = r$BIC,
             entropy = r$entropy, min_class_pct = r$min_pct,
             mean_post_prob = r$mean_pp)
}
tab <- bind_rows(
  lapply(res_g[!sapply(res_g, is.null)], function(r) mk_row("log_gfap", r)),
  lapply(res_p[!sapply(res_p, is.null)], function(r) mk_row("log_ptau", r))
)
print(tab)
utils::write.csv(tab, file.path(OUT, "table2_lcgmm_full_diagnostic.csv"), row.names = FALSE)

# Parametric bootstrap LRT comparing k vs k-1 (limited to B=200 for time)
cat("\nParametric bootstrap LRT comparing best vs k=1 (B=200)...\n")
boot_lrt <- function(m1, mk, dat, yvar, B = 200) {
  obs_lrt <- 2 * (mk$loglik - m1$loglik)
  if (obs_lrt <= 0) return(list(obs = obs_lrt, p = 1))
  cat("  observed LRT statistic:", round(obs_lrt, 2), "\n")
  null_lrts <- numeric(B)
  for (b in seq_len(B)) {
    # Simulate from k=1 model: use fitted random-effects parameters
    sim_dat <- dat
    # Simulate y from model1: y_hat + residual + random effect
    pred <- m1$pred
    if (nrow(pred) == nrow(dat)) {
      sim_dat[[yvar]] <- pred$pred_ss + rnorm(nrow(pred), 0, sqrt(m1$gconv[1]))
    } else {
      null_lrts[b] <- NA; next
    }
    m1b <- tryCatch(hlme(as.formula(paste(yvar, "~ time_yrs")),
                        random = ~time_yrs, subject = "RID", data = sim_dat, ng = 1),
                    error = function(e) NULL)
    mkb <- tryCatch(hlme(as.formula(paste(yvar, "~ time_yrs")),
                        mixture = ~time_yrs, random = ~time_yrs,
                        subject = "RID", data = sim_dat, ng = mk$ng, B = m1b),
                    error = function(e) NULL)
    null_lrts[b] <- if (!is.null(m1b) && !is.null(mkb)) 2*(mkb$loglik - m1b$loglik) else NA
    if (b %% 25 == 0) cat("    bootstrap", b, "/", B, "\n")
  }
  p <- mean(null_lrts >= obs_lrt, na.rm = TRUE)
  list(obs = obs_lrt, p = p, n_valid = sum(!is.na(null_lrts)))
}

# Best k for each biomarker (lowest BIC excluding small-class violators)
# log_gfap: k=1 wins clean (k=3 has 0.9% class)
# log_ptau: k=2 has BIC=1300 vs k=1 BIC=1316; min class 45.5% — borderline support
dat_p <- biom_long %>% filter(!is.na(log_ptau)) %>% as.data.frame()
m1p <- res_p[[1]]$m; m2p <- res_p[[2]]$m
blrt_ptau <- boot_lrt(m1p, m2p, dat_p, "log_ptau", B = 200)
cat("\nLCGMM bootstrap LRT for p-tau217, k=2 vs k=1:\n")
cat("  observed LRT =", round(blrt_ptau$obs, 2), "\n")
cat("  bootstrap p =", round(blrt_ptau$p, 4),
    "(", blrt_ptau$n_valid, "valid bootstraps)\n")

utils::write.csv(
  data.frame(biomarker = "log_ptau", k_test = 2, k_ref = 1,
             obs_LRT = blrt_ptau$obs, boot_p = blrt_ptau$p,
             n_valid_boot = blrt_ptau$n_valid),
  file.path(OUT, "lcgmm_bootstrap_LRT.csv"), row.names = FALSE)

cat("\nDONE.\n")
