# Phase III — Multivariate Bayesian Joint Longitudinal Model (brms, Path A)
# Author: Kuan-Yu Chu | 2026-05-25
# Approach: bivariate response with shared random effects across responses
#   bf(log_biom ~ time*APOE4 + covs + (1 + time | p | RID))
#   bf(outcome  ~ time       + covs + (1 + time | p | RID))
#   The cross-response random-slope correlation = joint-model association.
# 6 fits: GFAP and pT217 paired with ADAS-Cog 13, MMSE, hippo ratio.
# Hippo ratio scaled x1000 for numerical stability.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(brms)
})
options(brms.backend = "rstan", brms.file_refit = "never")
set.seed(20260525)

DATA <- "/Users/kuan-yuchu/Desktop/LOAD/ADNI_csv"
OUT  <- "/Users/kuan-yuchu/Desktop/LOAD/Analysis/outputs"

# ---- Load assembled long datasets ----
biom_long <- utils::read.csv(file.path(OUT, "analytic_biom_long.csv"),
                             stringsAsFactors = FALSE)
attr_tbl  <- utils::read.csv(file.path(OUT, "analytic_baseline_attr.csv"),
                             stringsAsFactors = FALSE)
biom_long$dx3 <- factor(biom_long$dx3, levels = c("CN","MCI","AD"))
attr_tbl$dx3  <- factor(attr_tbl$dx3,  levels = c("CN","MCI","AD"))

zsc <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
attr_tbl <- attr_tbl %>%
  mutate(age_z = zsc(age_bl), edu_z = zsc(edu))

# Build outcome long
adas <- utils::read.csv(file.path(DATA, "ADAS_24May2026.csv"),
                        stringsAsFactors = FALSE, na.strings = c("","NA"))
mmse <- utils::read.csv(file.path(DATA, "MMSE_24May2026.csv"),
                        stringsAsFactors = FALSE, na.strings = c("","NA"))
fs   <- utils::read.csv(file.path(DATA, "ADSP_PHC_T1_FS_24May2026.csv"),
                        stringsAsFactors = FALSE, na.strings = c("","NA"))

bl_dates <- biom_long %>% group_by(RID) %>%
  summarise(bl_date = as.Date(min(EXAMDATE)), .groups = "drop")

build_outc <- function(df, ycol, datecol) {
  df %>% mutate(RID = as.integer(RID),
                EXAMDATE = as.Date(.data[[datecol]]),
                y = suppressWarnings(as.numeric(.data[[ycol]]))) %>%
    filter(!is.na(EXAMDATE), !is.na(y), RID %in% attr_tbl$RID) %>%
    left_join(bl_dates, by = "RID") %>%
    mutate(time_yrs = as.numeric(EXAMDATE - bl_date)/365.25) %>%
    filter(!is.na(time_yrs))
}
adas_long <- build_outc(adas, "TOTAL13",  "VISDATE") %>% transmute(RID, time_yrs, adas = y)
mmse_long <- build_outc(mmse, "MMSCORE",  "VISDATE") %>% transmute(RID, time_yrs, mmse = y)
fs_long   <- fs %>% mutate(hippo_ratio = (Left.Hippocampus_combat + Right.Hippocampus_combat) /
                             EstimatedTotalIntraCranialVol_combat) %>%
  build_outc("hippo_ratio", "PHC_SCANDATE") %>%
  transmute(RID, time_yrs, hippo_x1k = y * 1000)

# Drop NA BMI_z rows from biom_long; keep covariates per RID
biom_red <- biom_long %>% filter(!is.na(BMI_z)) %>%
  transmute(RID, time_yrs, log_gfap, log_ptau)

# ---- Generic joint-model fitter ----
fit_brms_joint <- function(bm_col, outc_df, outc_col, lab) {
  cat("\n=== Fitting:", lab, "===\n  t0:", as.character(Sys.time()), "\n")
  # Each RID's biomarker rows
  bm_df <- biom_red %>% select(RID, time_yrs, bm = all_of(bm_col)) %>%
    filter(!is.na(bm))
  # Each RID's outcome rows
  outc_df <- outc_df %>% select(RID, time_yrs, outc = all_of(outc_col)) %>%
    filter(!is.na(outc))
  # Intersection
  rid_keep <- intersect(unique(bm_df$RID), unique(outc_df$RID))
  bm_df   <- bm_df  [bm_df$RID   %in% rid_keep, ]
  outc_df <- outc_df[outc_df$RID %in% rid_keep, ]
  # Stack with subset indicators so brms can use each response on its own rows.
  # Use 0 as placeholder (will not be used because subset filter excludes).
  bm_rows <- bm_df   %>% mutate(outc = 0, has_bm = TRUE,  has_oc = FALSE)
  oc_rows <- outc_df %>% mutate(bm   = 0, has_bm = FALSE, has_oc = TRUE)
  d <- bind_rows(bm_rows, oc_rows) %>%
    left_join(attr_tbl %>% select(RID, age_z, sex, edu_z, dx3, APOE4),
              by = "RID")
  d$RID <- factor(d$RID)
  cat("  bm rows:", nrow(bm_df), " outc rows:", nrow(outc_df),
      " unique RIDs:", length(rid_keep), "\n")

  # Bivariate formula; subset() restricts each response to rows with real data
  f_bm <- bf(bm   | subset(has_bm) ~ time_yrs + APOE4 + age_z + sex + dx3 +
               (1 + time_yrs | p | RID))
  f_oc <- bf(outc | subset(has_oc) ~ time_yrs + APOE4 + age_z + sex + edu_z + dx3 +
               (1 + time_yrs | p | RID))

  fit <- tryCatch(
    brm(f_bm + f_oc + set_rescor(FALSE),
        data = d, family = gaussian(),
        chains = 4, iter = 2000, warmup = 800, cores = 4,
        seed = 20260525, refresh = 0,
        control = list(adapt_delta = 0.95, max_treedepth = 12)),
    error = function(e) { message("  brms err: ", e$message); NULL })

  if (is.null(fit)) return(list(label = lab, ok = FALSE, mod = NULL))

  # Extract cross-response random-slope correlation: cor_RID__bm_time_yrs__outc_time_yrs
  ps <- as_draws_df(fit)
  # Find correlation parameter (named like "cor_RID__bm_time_yrs__outc_time_yrs")
  cor_nm <- grep("^cor_RID__bm_time_yrs__outc_time_yrs$", names(ps), value = TRUE)
  if (length(cor_nm) == 0) {
    cor_nm <- grep("^cor_RID.*time_yrs.*time_yrs$", names(ps), value = TRUE)[1]
  }
  if (length(cor_nm) == 0 || is.na(cor_nm)) {
    cat("  WARNING: cross-slope correlation parameter not found in posterior\n")
    return(list(label = lab, ok = TRUE, mod = fit,
                summary = data.frame(label = lab, note = "cor param not found",
                                     mean = NA, lo = NA, hi = NA,
                                     p_dir = NA, rhat = NA)))
  }
  draws <- ps[[cor_nm]]
  summ <- data.frame(
    label = lab,
    parameter = cor_nm,
    mean  = mean(draws),
    median = median(draws),
    lo    = quantile(draws, 0.025),
    hi    = quantile(draws, 0.975),
    p_dir = mean(draws < 0),      # P(rho < 0) for negative associations
    rhat  = max(brms::rhat(fit)[grep("cor_RID", names(brms::rhat(fit)))]),
    n_id  = length(rid_keep),
    n_bm  = nrow(bm_df),
    n_oc  = nrow(outc_df)
  )
  cat("  rho =", round(summ$mean, 3),
      " [", round(summ$lo, 3), ",", round(summ$hi, 3), "]\n",
      "  P(rho<0) =", round(summ$p_dir, 3), "\n")
  list(label = lab, ok = TRUE, mod = fit, summary = summ)
}

# ---- Run 6 joint fits ----
pairs <- list(
  list(bm = "log_gfap", outc = adas_long, oc = "adas",       lab = "ADAS_GFAP"),
  list(bm = "log_ptau", outc = adas_long, oc = "adas",       lab = "ADAS_pT217"),
  list(bm = "log_gfap", outc = mmse_long, oc = "mmse",       lab = "MMSE_GFAP"),
  list(bm = "log_ptau", outc = mmse_long, oc = "mmse",       lab = "MMSE_pT217"),
  list(bm = "log_gfap", outc = fs_long,   oc = "hippo_x1k",  lab = "Hippo_GFAP"),
  list(bm = "log_ptau", outc = fs_long,   oc = "hippo_x1k",  lab = "Hippo_pT217")
)

all_summ <- list()
for (p in pairs) {
  res <- fit_brms_joint(p$bm, p$outc, p$oc, p$lab)
  if (!is.null(res$mod)) saveRDS(res$mod, file.path(OUT, paste0("brmsfit_", p$lab, ".rds")))
  if (!is.null(res$summary)) all_summ[[p$lab]] <- res$summary
}

if (length(all_summ) > 0) {
  out_tbl <- bind_rows(all_summ)
  utils::write.csv(out_tbl, file.path(OUT, "table4_brms_joint_primary.csv"), row.names = FALSE)
  cat("\n=== brms joint-model results ===\n")
  print(out_tbl)
  cat("\nSaved table4_brms_joint_primary.csv\n")
} else {
  cat("\nNo successful fits.\n")
}
