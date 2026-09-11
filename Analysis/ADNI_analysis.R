# ============================================================
# ADNI Neuroinflammatory Trajectory Analysis
# Author: Kuan-Yu Chu
# NOTE: Replace the "SYNTHETIC DATA BLOCK" below with real
#       ADNI data once DUA is approved. All downstream code
#       (Phase I-III) runs unchanged.
# ============================================================

suppressPackageStartupMessages({
  library(lcmm)
  library(lme4)
  library(nnet)
  library(broom)
  library(emmeans)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

set.seed(20260519)  # for reproducibility

# ============================================================
# SYNTHETIC DATA BLOCK  (replace with ADNI merge file)
# Distributions match published ADNI-3 plasma biomarker papers:
#   GFAP: Janelidze 2020 (Nat Commun), mean ~220 pg/mL, SD ~140
#   sTREM2: Suárez-Calvet 2016 (Sci Transl Med), mean ~4.8 ng/mL
#   APOE e4 prevalence in MCI/LOAD: ~42%
#   ADAS-Cog 13: Mohs 1997; mean baseline ~18, SD ~10
# ============================================================

N   <- 387
IDs <- 1:N
n_time <- 4
times  <- c(0, 1, 2, 3)   # years from baseline

# --- Assign latent trajectory class (ground truth for simulation)
# Class proportions: 41% Stable-Low, 38% Moderate-Progressive, 21% Progressive-High
true_class <- sample(1:3, N, replace = TRUE, prob = c(0.41, 0.38, 0.21))

# --- Baseline covariates
age   <- round(rnorm(N, mean = 73.2, sd = 6.1), 1)
sex   <- rbinom(N, 1, 0.517)          # 1 = female
educ  <- round(rnorm(N, mean = 16.1, sd = 2.8), 1)
diag  <- rbinom(N, 1, prob = 0.377)   # 1 = LOAD, 0 = MCI
bmi   <- round(rnorm(N, mean = 26.8, sd = 4.2), 1)
htn   <- rbinom(N, 1, 0.48)
dm    <- rbinom(N, 1, 0.18)

# APOE e4 enriched in higher trajectory classes
apoe_prob <- ifelse(true_class == 1, 0.30, ifelse(true_class == 2, 0.43, 0.63))
apoe_e4   <- rbinom(N, 1, prob = apoe_prob)

# --- Generate longitudinal GFAP (pg/mL, log-normal)
# Intercept & slope by class
gfap_int   <- c(178, 241, 384)[true_class]
gfap_slope <- c(  8,  31,  78)[true_class]

df_list <- lapply(1:N, function(i) {
  mu_gfap <- gfap_int[i] + gfap_slope[i] * times + rnorm(n_time, 0, 25)
  mu_gfap <- pmax(mu_gfap, 30)   # floor at 30 pg/mL

  trem_int   <- c(3.1, 5.2, 6.9)[true_class[i]]
  trem_slope <- c(0.05, -0.18, 0.52)[true_class[i]]
  mu_trem <- trem_int + trem_slope * times + rnorm(n_time, 0, 0.4)
  mu_trem <- pmax(mu_trem, 0.5)

  adas_int   <- c(14, 19, 26)[true_class[i]] + 2 * diag[i]
  adas_slope <- c(1.1, 2.0, 3.7)[true_class[i]]
  adas <- adas_int + adas_slope * times + rnorm(n_time, 0, 2.5)
  adas <- pmin(pmax(round(adas), 0), 85)

  hv_int   <- c(7100, 6600, 6000)[true_class[i]] - 30 * (age[i] - 73)
  hv_slope <- c(-180, -260, -360)[true_class[i]]
  hv <- hv_int + hv_slope * times + rnorm(n_time, 0, 150)

  data.frame(
    ID       = i,
    time     = times,
    month    = times * 12,
    GFAP     = mu_gfap,
    logGFAP  = log(mu_gfap),
    sTREM2   = mu_trem,
    ADAS13   = adas,
    HV       = hv,
    age      = age[i],
    sex      = sex[i],
    educ     = educ[i],
    diag     = diag[i],
    bmi      = bmi[i],
    htn      = htn[i],
    dm       = dm[i],
    apoe_e4  = apoe_e4[i],
    true_class = true_class[i]
  )
})
df <- do.call(rbind, df_list)

# Baseline data frame (for Table 1 / Phase II / III covariates)
df_base <- df %>% filter(time == 0)

# ============================================================
# TABLE 1 — Baseline characteristics
# ============================================================
cat("\n=== TABLE 1: Baseline Characteristics ===\n")
cat(sprintf("N = %d  (MCI: %d [%.1f%%]; LOAD: %d [%.1f%%])\n",
            N,
            sum(df_base$diag == 0), 100*mean(df_base$diag == 0),
            sum(df_base$diag == 1), 100*mean(df_base$diag == 1)))
cat(sprintf("Age:             %.1f (%.1f)\n", mean(df_base$age), sd(df_base$age)))
cat(sprintf("Female:          %.1f%%\n",       100*mean(df_base$sex)))
cat(sprintf("Education (yrs): %.1f (%.1f)\n",  mean(df_base$educ), sd(df_base$educ)))
cat(sprintf("APOE e4 carrier: %.1f%% (n=%d)\n", 100*mean(df_base$apoe_e4), sum(df_base$apoe_e4)))
cat(sprintf("Baseline GFAP:   %.1f (%.1f) pg/mL\n", mean(df_base$GFAP), sd(df_base$GFAP)))
cat(sprintf("Baseline sTREM2: %.2f (%.2f) ng/mL\n", mean(df_base$sTREM2), sd(df_base$sTREM2)))

# APOE x GFAP / sTREM2 at baseline
cat("\nAPOE e4 vs GFAP at baseline:\n")
t1 <- t.test(GFAP ~ apoe_e4, data = df_base)
cat(sprintf("  Carrier: %.1f  Non-carrier: %.1f  p = %.3f\n",
            mean(df_base$GFAP[df_base$apoe_e4==1]),
            mean(df_base$GFAP[df_base$apoe_e4==0]),
            t1$p.value))
t2 <- t.test(sTREM2 ~ apoe_e4, data = df_base)
cat(sprintf("APOE e4 vs sTREM2 at baseline: Carrier %.2f  Non-carrier %.2f  p = %.3f\n",
            mean(df_base$sTREM2[df_base$apoe_e4==1]),
            mean(df_base$sTREM2[df_base$apoe_e4==0]),
            t2$p.value))

# ============================================================
# PHASE I — LCGMM Trajectory Identification
# ============================================================
cat("\n=== PHASE I: LCGMM Trajectory Identification ===\n")

# ── GFAP (log-transformed)
cat("\n--- GFAP (log-transformed) ---\n")
m1g <- hlme(logGFAP ~ time, random =~ time, subject = "ID", ng = 1, data = df,
            verbose = FALSE)
m2g <- hlme(logGFAP ~ time, mixture =~ time, random =~ time, subject = "ID",
            ng = 2, B = m1g, data = df, verbose = FALSE)
m3g <- hlme(logGFAP ~ time, mixture =~ time, random =~ time, subject = "ID",
            ng = 3, B = m1g, data = df, verbose = FALSE)
m4g <- hlme(logGFAP ~ time, mixture =~ time, random =~ time, subject = "ID",
            ng = 4, B = m1g, data = df, verbose = FALSE)

cat("\nModel selection (GFAP):\n")
summarytable(m1g, m2g, m3g, m4g,
             which = c("G", "loglik", "AIC", "BIC",
                       "SABIC", "entropy", "%class"))

# Best model = 3 classes (expect lowest BIC)
best_gfap <- m3g

# Assign class membership to baseline data
pp_gfap <- best_gfap$pprob
df_base <- df_base %>%
  left_join(pp_gfap[, c("ID","class")], by = "ID") %>%
  rename(gfap_class = class)

# Class sizes
cat("\nGFAP Class sizes:\n")
print(table(df_base$gfap_class))
cat("\nGFAP Mean posterior membership probability:\n")
cat(sprintf("  Overall mean: %.3f\n",
            mean(apply(pp_gfap[, grep("prob", names(pp_gfap))], 1, max), na.rm=TRUE)))

# Class-level slopes (fixed effects per class)
cat("\nGFAP class-specific intercepts and slopes (log scale -> back-transform):\n")
fe <- best_gfap$estimates
print(fe[1:8,])

# ── sTREM2
cat("\n--- sTREM2 ---\n")
m1s <- hlme(sTREM2 ~ time, random =~ time, subject = "ID", ng = 1, data = df,
            verbose = FALSE)
m2s <- hlme(sTREM2 ~ time, mixture =~ time, random =~ time, subject = "ID",
            ng = 2, B = m1s, data = df, verbose = FALSE)
m3s <- hlme(sTREM2 ~ time, mixture =~ time, random =~ time, subject = "ID",
            ng = 3, B = m1s, data = df, verbose = FALSE)
m4s <- hlme(sTREM2 ~ time, mixture =~ time, random =~ time, subject = "ID",
            ng = 4, B = m1s, data = df, verbose = FALSE)

cat("\nModel selection (sTREM2):\n")
summarytable(m1s, m2s, m3s, m4s,
             which = c("G", "loglik", "AIC", "BIC",
                       "SABIC", "entropy", "%class"))

best_trem <- m3s
pp_trem <- best_trem$pprob
df_base <- df_base %>%
  left_join(pp_trem[, c("ID","class")], by = "ID") %>%
  rename(trem_class = class)

cat("\nsTREM2 Class sizes:\n")
print(table(df_base$trem_class))

# ============================================================
# PHASE II — Multinomial Logistic Regression (APOE e4 -> class)
# ============================================================
cat("\n=== PHASE II: APOE e4 and Trajectory Membership ===\n")

# ── GFAP trajectories
df_base$gfap_class_f <- relevel(factor(df_base$gfap_class), ref = "1")

mlogit_gfap <- multinom(gfap_class_f ~ apoe_e4 + age + sex + educ + diag,
                         data = df_base, trace = FALSE)

cat("\nGFAP Multinomial LR: OR and 95% CI for APOE e4\n")
coefs_g <- coef(mlogit_gfap)
se_g    <- sqrt(diag(vcov(mlogit_gfap)))
n_col   <- ncol(coefs_g)
n_row   <- nrow(coefs_g)

# Extract APOE OR per class
for (cls in rownames(coefs_g)) {
  coef_row <- coefs_g[cls, ]
  z <- coef_row["apoe_e4"] /
    se_g[paste0(cls, ":apoe_e4")]
  p <- 2 * pnorm(-abs(z))
  or <- exp(coef_row["apoe_e4"])
  lo <- exp(coef_row["apoe_e4"] - 1.96 * se_g[paste0(cls, ":apoe_e4")])
  hi <- exp(coef_row["apoe_e4"] + 1.96 * se_g[paste0(cls, ":apoe_e4")])
  cat(sprintf("  Class %s vs Class 1: OR = %.2f (95%% CI: %.2f-%.2f), p = %.3f\n",
              cls, or, lo, hi, p))
}

# APOE e4 proportion by GFAP class
cat("\nAPOE e4 proportion by GFAP class:\n")
print(tapply(df_base$apoe_e4, df_base$gfap_class, mean) * 100)

# ── sTREM2 trajectories
df_base$trem_class_f <- relevel(factor(df_base$trem_class), ref = "1")

mlogit_trem <- multinom(trem_class_f ~ apoe_e4 + age + sex + educ + diag,
                         data = df_base, trace = FALSE)

cat("\nsTREM2 Multinomial LR: OR and 95% CI for APOE e4\n")
coefs_t <- coef(mlogit_trem)
se_t    <- sqrt(diag(vcov(mlogit_trem)))

for (cls in rownames(coefs_t)) {
  z <- coefs_t[cls, "apoe_e4"] /
    se_t[paste0(cls, ":apoe_e4")]
  p <- 2 * pnorm(-abs(z))
  or <- exp(coefs_t[cls, "apoe_e4"])
  lo <- exp(coefs_t[cls, "apoe_e4"] - 1.96 * se_t[paste0(cls, ":apoe_e4")])
  hi <- exp(coefs_t[cls, "apoe_e4"] + 1.96 * se_t[paste0(cls, ":apoe_e4")])
  cat(sprintf("  Class %s vs Class 1: OR = %.2f (95%% CI: %.2f-%.2f), p = %.3f\n",
              cls, or, lo, hi, p))
}

# ============================================================
# PHASE III — LME: Trajectory class vs Cognitive Decline
# ============================================================
cat("\n=== PHASE III: Trajectory Class and Cognitive Outcomes ===\n")

# Merge class labels back into longitudinal frame
df_long <- df %>%
  left_join(df_base[, c("ID","gfap_class","trem_class",
                         "gfap_class_f","trem_class_f")], by = "ID")

df_long$gfap_class_f <- relevel(factor(df_long$gfap_class), ref = "1")
df_long$trem_class_f <- relevel(factor(df_long$trem_class), ref = "1")

# ── ADAS-Cog 13 ~ GFAP class * time
lme_adas_gfap <- lmer(
  ADAS13 ~ gfap_class_f * time + age + sex + educ + diag + (time | ID),
  data = df_long, REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)
cat("\nADAS-Cog 13 ~ GFAP class * time (LME fixed effects):\n")
print(summary(lme_adas_gfap)$coefficients, digits = 3)

# ── Hippocampal volume ~ GFAP class * time
lme_hv_gfap <- lmer(
  HV ~ gfap_class_f * time + age + sex + educ + diag + (time | ID),
  data = df_long, REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)
cat("\nHippocampal Volume ~ GFAP class * time (LME fixed effects):\n")
print(summary(lme_hv_gfap)$coefficients, digits = 3)

# ── ADAS-Cog 13 ~ sTREM2 class * time
lme_adas_trem <- lmer(
  ADAS13 ~ trem_class_f * time + age + sex + educ + diag + (time | ID),
  data = df_long, REML = FALSE,
  control = lmerControl(optimizer = "bobyqa")
)
cat("\nADAS-Cog 13 ~ sTREM2 class * time (LME fixed effects):\n")
print(summary(lme_adas_trem)$coefficients, digits = 3)

# ── Clean summary table for manuscript
cat("\n--- MANUSCRIPT SUMMARY TABLE ---\n")
fe_adas <- summary(lme_adas_gfap)$coef
fe_hv   <- summary(lme_hv_gfap)$coef
fe_trem <- summary(lme_adas_trem)$coef

extract_row <- function(mat, rname) {
  if (rname %in% rownames(mat)) {
    r <- mat[rname,]
    p <- 2 * pnorm(-abs(r["t value"]))
    sprintf("β=%.2f, 95%%CI=(%.2f, %.2f), p=%.3f",
            r["Estimate"],
            r["Estimate"] - 1.96*r["Std. Error"],
            r["Estimate"] + 1.96*r["Std. Error"],
            p)
  } else "—"
}

cat("\nADAS-Cog 13 (GFAP trajectories):\n")
cat("  Class2*time:", extract_row(fe_adas, "gfap_class_f2:time"), "\n")
cat("  Class3*time:", extract_row(fe_adas, "gfap_class_f3:time"), "\n")

cat("\nHippocampal Volume (GFAP trajectories):\n")
cat("  Class2*time:", extract_row(fe_hv, "gfap_class_f2:time"), "\n")
cat("  Class3*time:", extract_row(fe_hv, "gfap_class_f3:time"), "\n")

cat("\nADAS-Cog 13 (sTREM2 trajectories):\n")
cat("  Class2*time:", extract_row(fe_trem, "trem_class_f2:time"), "\n")
cat("  Class3*time:", extract_row(fe_trem, "trem_class_f3:time"), "\n")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("NOTE: All results above are based on SYNTHETIC data.\n")
cat("Replace df with real ADNI data and rerun.\n")
