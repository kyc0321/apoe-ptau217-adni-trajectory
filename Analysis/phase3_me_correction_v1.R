# Phase III — Measurement-error correction via Multiple Imputation of BLUP slopes
# Author: Kuan-Yu Chu | 2026-05-25
# Strategy: For each RID, draw 100 plausible slopes from N(BLUP_mean, BLUP_SE),
#           refit the outcome LME each draw, pool via Rubin's rules.
# This replaces the JMbayes2 step (JMbayes2 is for longitudinal+survival, not pure longitudinal-longitudinal).

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(lme4); library(lmerTest)
})
set.seed(20260525)

DATA <- "/Users/kuan-yuchu/Desktop/LOAD/ADNI_csv"
OUT  <- "/Users/kuan-yuchu/Desktop/LOAD/Analysis/outputs"

# ---- Load assembled long datasets from master pipeline ----
biom_long <- utils::read.csv(file.path(OUT, "analytic_biom_long.csv"),
                             stringsAsFactors = FALSE)
biom_long$dx3 <- factor(biom_long$dx3, levels = c("CN","MCI","AD"))

# ---- Refit Phase II to obtain BLUP MEAN and conditional SE ----
phase2 <- function(yvar) {
  d <- biom_long %>% filter(!is.na(.data[[yvar]]), !is.na(BMI_z))
  fml <- as.formula(paste0(yvar,
    " ~ time_yrs*APOE4 + age_z + sex + edu_z + dx3 + BMI_z + HTN_base + DM_base + (1 + time_yrs | RID)"))
  m <- lmerTest::lmer(fml, data = d, REML = TRUE,
                      control = lmerControl(optimizer="bobyqa"))
  re <- ranef(m, condVar = TRUE)$RID
  cv <- attr(re, "postVar")  # array: 2 x 2 x n_RID
  fe_time <- fixef(m)["time_yrs"]
  data.frame(
    RID = as.integer(rownames(re)),
    slope_mean = re[, "time_yrs"] + fe_time,
    slope_se   = sqrt(cv[2, 2, ])
  )
}
cat("Refitting Phase II for GFAP...\n")
blup_gfap <- phase2("log_gfap"); names(blup_gfap)[2:3] <- c("slope_gfap","se_gfap")
cat("Refitting Phase II for p-tau217...\n")
blup_ptau <- phase2("log_ptau"); names(blup_ptau)[2:3] <- c("slope_ptau","se_ptau")
blups <- full_join(blup_gfap, blup_ptau, by = "RID")
utils::write.csv(blups, file.path(OUT, "blup_slopes_with_SE.csv"), row.names = FALSE)
cat("Saved BLUP with conditional SE\n")
cat("GFAP slope SE  median (IQR):",
    sprintf("%.4f (%.4f-%.4f)", median(blups$se_gfap, na.rm=TRUE),
            quantile(blups$se_gfap,.25,na.rm=TRUE),
            quantile(blups$se_gfap,.75,na.rm=TRUE)), "\n")
cat("pT217 slope SE median (IQR):",
    sprintf("%.4f (%.4f-%.4f)", median(blups$se_ptau, na.rm=TRUE),
            quantile(blups$se_ptau,.25,na.rm=TRUE),
            quantile(blups$se_ptau,.75,na.rm=TRUE)), "\n")

# ---- Build outcome long tables (re-read & rebuild) ----
attr_tbl <- utils::read.csv(file.path(OUT, "analytic_baseline_attr.csv"),
                            stringsAsFactors = FALSE)
attr_tbl$dx3 <- factor(attr_tbl$dx3, levels = c("CN","MCI","AD"))
zsc <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)

adas <- utils::read.csv(file.path(DATA, "ADAS_24May2026.csv"),
                        stringsAsFactors = FALSE, na.strings = c("","NA"))
mmse <- utils::read.csv(file.path(DATA, "MMSE_24May2026.csv"),
                        stringsAsFactors = FALSE, na.strings = c("","NA"))
fs   <- utils::read.csv(file.path(DATA, "ADSP_PHC_T1_FS_24May2026.csv"),
                        stringsAsFactors = FALSE, na.strings = c("","NA"))
load("/tmp/adnimerge_unpack/ADNIMERGE/data/adnimerge.rdata")

bl_dates <- biom_long %>% group_by(RID) %>%
  summarise(bl_date = as.Date(min(EXAMDATE)), .groups = "drop")

build_out <- function(df, ycol, datecol) {
  out <- df %>% mutate(RID = as.integer(RID),
                       EXAMDATE = as.Date(.data[[datecol]]),
                       y = suppressWarnings(as.numeric(.data[[ycol]]))) %>%
    filter(!is.na(EXAMDATE), !is.na(y), RID %in% attr_tbl$RID) %>%
    left_join(bl_dates, by = "RID") %>%
    mutate(time_yrs = as.numeric(EXAMDATE - bl_date)/365.25) %>%
    filter(!is.na(time_yrs)) %>%
    left_join(attr_tbl, by = "RID")
  out$age_z <- zsc(out$age_bl); out$edu_z <- zsc(out$edu)
  out
}
adas_long <- build_out(adas, "TOTAL13",  "VISDATE")
mmse_long <- build_out(mmse, "MMSCORE",  "VISDATE")
fs_long   <- fs %>% mutate(hippo_ratio = (Left.Hippocampus_combat + Right.Hippocampus_combat)/EstimatedTotalIntraCranialVol_combat) %>%
  build_out("hippo_ratio", "PHC_SCANDATE")

# ---- Multiple-imputation measurement-error correction ----
M <- 100   # imputations
cat("\nRunning MI measurement-error correction (M =", M, ")...\n")

me_correct <- function(outc_df, ycol, blup_mean_col, blup_se_col, M = 100) {
  blup_sub <- blups %>% filter(!is.na(.data[[blup_mean_col]]),
                               !is.na(.data[[blup_se_col]]))
  # FIX: standardise using FIXED constants from the BLUP MEAN distribution,
  # not the per-imputation distribution, so MI variance reflects pure measurement
  # error rather than per-draw rescaling noise.
  scale_mean <- mean(blup_sub[[blup_mean_col]], na.rm = TRUE)
  scale_sd   <- sd(  blup_sub[[blup_mean_col]], na.rm = TRUE)
  betas <- numeric(M); ses <- numeric(M)
  for (m in seq_len(M)) {
    s <- blup_sub
    s$slope_imp <- rnorm(nrow(s), s[[blup_mean_col]], s[[blup_se_col]])
    s$slope_z   <- (s$slope_imp - scale_mean) / scale_sd
    d <- outc_df %>% inner_join(s[, c("RID","slope_z")], by = "RID")
    fml <- as.formula(paste0(ycol,
      " ~ time_yrs*slope_z + APOE4 + age_z + sex + edu_z + dx3 + (1 + time_yrs | RID)"))
    fit <- tryCatch(lmerTest::lmer(fml, data = d, REML = TRUE,
                                   control = lmerControl(optimizer="bobyqa")),
                    error = function(e) NULL)
    if (is.null(fit)) { betas[m] <- NA; ses[m] <- NA; next }
    cf <- summary(fit)$coefficients
    ix <- grep("time_yrs:slope_z", rownames(cf))
    if (length(ix) == 0) { betas[m] <- NA; ses[m] <- NA; next }
    betas[m] <- cf[ix, "Estimate"]; ses[m] <- cf[ix, "Std. Error"]
  }
  # Rubin's rules
  betas <- betas[!is.na(betas)]; ses <- ses[!is.na(ses)]
  if (length(betas) < 2) return(NULL)
  b_pool <- mean(betas)
  vw <- mean(ses^2)               # within-imputation variance
  vb <- var(betas)                # between-imputation variance
  vt <- vw + (1 + 1/length(betas)) * vb
  se_pool <- sqrt(vt)
  # Approximate degrees of freedom (Barnard & Rubin 1999, simplified)
  riv <- (1 + 1/length(betas))*vb / vw
  df  <- (length(betas) - 1) * (1 + 1/riv)^2
  t_stat <- b_pool / se_pool
  p_val  <- 2 * pt(-abs(t_stat), df = max(df, 5))
  data.frame(outcome = ycol, predictor = blup_mean_col,
             beta_per_SD = b_pool, SE = se_pool,
             CI_lo = b_pool - 1.96*se_pool, CI_hi = b_pool + 1.96*se_pool,
             p_value = p_val, M_used = length(betas),
             riv = riv, df = df)
}

results <- bind_rows(
  me_correct(adas_long, "y", "slope_gfap", "se_gfap", M),
  me_correct(adas_long, "y", "slope_ptau", "se_ptau", M),
  me_correct(mmse_long, "y", "slope_gfap", "se_gfap", M),
  me_correct(mmse_long, "y", "slope_ptau", "se_ptau", M),
  me_correct(fs_long,   "y", "slope_gfap", "se_gfap", M),
  me_correct(fs_long,   "y", "slope_ptau", "se_ptau", M)
)
results$outcome <- rep(c("adas13","adas13","mmse","mmse","hippo_ratio","hippo_ratio"),
                      length.out = nrow(results))
results$q_BH   <- p.adjust(results$p_value, method = "BH")
results$q_Bonf <- pmin(results$p_value * 6, 1)

utils::write.csv(results, file.path(OUT, "table4_ME_corrected_primary.csv"), row.names = FALSE)
cat("\n=== ME-corrected Phase III results (primary) ===\n")
print(results)
cat("\nSaved table4_ME_corrected_primary.csv\n")
