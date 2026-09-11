# Diagnostic battery — issues #4, #6, #9, #10 from manuscript critical review
# Author: Kuan-Yu Chu | 2026-05-26
# (#4) Per-dx hippocampal volume slope to verify biological plausibility
# (#6) Count negative time_yrs in outcome tables (baseline anchoring sanity)
# (#9) VIF for Phase II fixed effects (collinearity)
# (#10) BLUP normality QQ-plot + Shapiro test

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr)
  library(lme4); library(lmerTest); library(nlme); library(car); library(ggplot2)
})
set.seed(20260525)

DATA <- "/Users/kuan-yuchu/Desktop/LOAD/ADNI_csv"
OUT  <- "/Users/kuan-yuchu/Desktop/LOAD/Analysis/outputs"
FIG  <- "/Users/kuan-yuchu/Desktop/LOAD/figures/v4"

biom_long <- utils::read.csv(file.path(OUT, "analytic_biom_long.csv"), stringsAsFactors = FALSE)
attr_tbl  <- utils::read.csv(file.path(OUT, "analytic_baseline_attr.csv"), stringsAsFactors = FALSE)
biom_long$dx3 <- factor(biom_long$dx3, levels = c("CN","MCI","AD"))
attr_tbl$dx3  <- factor(attr_tbl$dx3,  levels = c("CN","MCI","AD"))

zsc <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

# Rebuild outcome long tables
adas <- utils::read.csv(file.path(DATA, "ADAS_24May2026.csv"), stringsAsFactors = FALSE, na.strings = c("","NA"))
mmse <- utils::read.csv(file.path(DATA, "MMSE_24May2026.csv"), stringsAsFactors = FALSE, na.strings = c("","NA"))
fs   <- utils::read.csv(file.path(DATA, "ADSP_PHC_T1_FS_24May2026.csv"), stringsAsFactors = FALSE, na.strings = c("","NA"))

bl_dates <- biom_long %>% group_by(RID) %>%
  summarise(bl_date = as.Date(min(EXAMDATE)), .groups = "drop")

build_outc <- function(df, ycol, datecol) {
  df %>% mutate(RID = as.integer(RID), EXAMDATE = as.Date(.data[[datecol]]),
                y = suppressWarnings(as.numeric(.data[[ycol]]))) %>%
    filter(!is.na(EXAMDATE), !is.na(y), RID %in% attr_tbl$RID) %>%
    left_join(bl_dates, by = "RID") %>%
    mutate(time_yrs = as.numeric(EXAMDATE - bl_date)/365.25) %>%
    filter(!is.na(time_yrs)) %>%
    left_join(attr_tbl, by = "RID")
}
adas_long <- build_outc(adas, "TOTAL13",  "VISDATE")
mmse_long <- build_outc(mmse, "MMSCORE",  "VISDATE")
fs_long   <- fs %>% mutate(hippo_ratio = (Left.Hippocampus_combat + Right.Hippocampus_combat) /
                             EstimatedTotalIntraCranialVol_combat,
                           hippo_mm3   = Left.Hippocampus_combat + Right.Hippocampus_combat) %>%
  build_outc("hippo_mm3", "PHC_SCANDATE")

# ============================================================
# (#6) NEGATIVE time_yrs CHECK
# ============================================================
cat("============================================================\n")
cat("#6  NEGATIVE time_yrs CHECK\n")
cat("============================================================\n")
for (nm in c("adas_long","mmse_long","fs_long")) {
  d <- get(nm)
  n_neg <- sum(d$time_yrs < 0)
  pct_neg <- 100 * n_neg / nrow(d)
  rid_neg <- length(unique(d$RID[d$time_yrs < 0]))
  cat(sprintf("%-10s : %d/%d rows (%.1f%%) < 0; %d unique RIDs with any negative\n",
              nm, n_neg, nrow(d), pct_neg, rid_neg))
  cat(sprintf("              earliest time = %.2f yrs; latest = %.2f yrs\n",
              min(d$time_yrs), max(d$time_yrs)))
}
cat("\nInterpretation: negative time_yrs reflects outcome assessments BEFORE the first\n")
cat("biomarker measurement (ADNI cognitive/imaging often pre-date plasma assays from ADNI3/4).\n")
cat("These pre-baseline observations are still informative for outcome trajectory estimation in LME.\n\n")

# ============================================================
# (#4) PER-DX HIPPOCAMPAL VOLUME SLOPE (biological plausibility check)
# ============================================================
cat("============================================================\n")
cat("#4  PER-DX HIPPOCAMPAL VOLUME SLOPE (mm³/year)\n")
cat("============================================================\n")
# Use only post-baseline (time_yrs >= 0) for a clean within-cohort slope
fs_pos <- fs_long %>% filter(time_yrs >= 0)
fit_per_dx <- function(d, dx_lab) {
  d_dx <- d %>% filter(dx3 == dx_lab)
  if (nrow(d_dx) < 20) return(data.frame(dx = dx_lab, n_obs = nrow(d_dx),
                                          n_id = length(unique(d_dx$RID)),
                                          slope_mm3_yr = NA, SE = NA, p = NA))
  m <- tryCatch(lmer(hippo_mm3 ~ time_yrs + (1 + time_yrs | RID), data = d_dx, REML = TRUE,
                     control = lmerControl(optimizer="bobyqa")),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  cf <- summary(m)$coefficients
  data.frame(dx = dx_lab,
             n_obs = nrow(d_dx), n_id = length(unique(d_dx$RID)),
             slope_mm3_yr = cf["time_yrs","Estimate"],
             SE = cf["time_yrs","Std. Error"],
             p = cf["time_yrs","Pr(>|t|)"])
}
hippo_rates <- bind_rows(lapply(levels(fs_pos$dx3), function(x) fit_per_dx(fs_pos, x)))
cat("Per-dx hippocampal volume slope (mm³/year, post-baseline only):\n")
print(hippo_rates)
cat("\nBiological expectation (literature):\n")
cat("  CN  : −50 to −150 mm³/yr  (normal aging atrophy)\n")
cat("  MCI : −200 to −400 mm³/yr (accelerated)\n")
cat("  AD  : −300 to −600 mm³/yr (further accelerated)\n\n")
utils::write.csv(hippo_rates, file.path(OUT, "hippo_slope_per_dx.csv"), row.names = FALSE)

# ============================================================
# (#9) VIF FOR PHASE II FIXED EFFECTS
# ============================================================
cat("============================================================\n")
cat("#9  VIF for Phase II fixed effects (collinearity check)\n")
cat("============================================================\n")
# Use a non-mixed model with same fixed-effect structure to compute VIF
d_p2 <- biom_long %>% filter(!is.na(log_ptau), !is.na(BMI_z)) %>%
  mutate(age_z = zsc(age_bl), edu_z = zsc(edu))
m_lm <- lm(log_ptau ~ time_yrs + APOE4 + age_z + sex + edu_z + dx3 +
                       BMI_z + HTN_base + DM_base, data = d_p2)
vifs <- car::vif(m_lm)
cat("VIFs (fixed-effect-only OLS proxy for Phase II model):\n")
print(round(vifs, 2))
cat(sprintf("Max VIF: %.2f — %s\n",
            max(vifs[,1] %||% vifs),
            ifelse(max(vifs[,1] %||% vifs) < 5, "PASS (<5)", "WARNING (≥5)")))
utils::write.csv(data.frame(term = rownames(vifs), VIF = vifs[,1] %||% vifs),
                 file.path(OUT, "vif_phase2.csv"), row.names = FALSE)

# ============================================================
# (#10) BLUP NORMALITY (QQ-plot + Shapiro)
# ============================================================
cat("\n============================================================\n")
cat("#10  BLUP NORMALITY CHECK\n")
cat("============================================================\n")
blups <- utils::read.csv(file.path(OUT, "blup_slopes_with_SE.csv"), stringsAsFactors = FALSE)
for (col in c("slope_gfap","slope_ptau")) {
  x <- blups[[col]]; x <- x[!is.na(x)]
  sh <- shapiro.test(x)
  cat(sprintf("%-12s  n = %d  Shapiro W = %.4f  p = %.4g  (%s)\n",
              col, length(x), sh$statistic, sh$p.value,
              ifelse(sh$p.value < 0.05, "rejects normality", "consistent with normality")))
}
cat("\nNote: brms primary analysis does NOT assume marginal BLUP normality;\n")
cat("HMC samples directly from the joint posterior. Frequentist BLUP-based\n")
cat("sensitivities (S-G, S-H) inherit normality assumptions.\n\n")

# QQ-plot
qq_df <- bind_rows(
  data.frame(biomarker = "GFAP slope BLUP",    value = blups$slope_gfap),
  data.frame(biomarker = "p-tau217 slope BLUP", value = blups$slope_ptau)
) %>% filter(!is.na(value))

p_qq <- ggplot(qq_df, aes(sample = value)) +
  stat_qq(alpha = 0.6) + stat_qq_line(colour = "red") +
  facet_wrap(~ biomarker, scales = "free") +
  labs(title = "Supplementary Figure. QQ-plot of per-participant BLUP slopes",
       x = "Theoretical quantiles", y = "Sample quantiles") +
  theme_minimal(base_size = 11)
ggsave(file.path(FIG, "SuppFig_BLUP_QQ.png"), p_qq, dpi = 300, width = 8, height = 4, bg = "white")
ggsave(file.path(FIG, "SuppFig_BLUP_QQ.pdf"), p_qq, width = 8, height = 4)
cat("Saved: SuppFig_BLUP_QQ.png/pdf\n\n")

cat("============================================================\n")
cat(" DIAGNOSTIC BATTERY COMPLETE  ", as.character(Sys.time()), "\n")
cat("============================================================\n")
