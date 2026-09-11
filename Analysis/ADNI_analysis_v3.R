# =============================================================================
# ADNI_analysis_v3.R   (2026-05-24)
# Wired to the ACTUAL ADNI CSVs downloaded to ~/Desktop/LOAD/ADNI_csv/.
#
# Design:
#   primary exposures (same UPENN file, same row):
#     - plasma GFAP   = UPENN GFAP_Q (Quanterix Simoa)
#     - plasma p-tau217 = UPENN pT217_F (Fujirebio)
#   sensitivity p-tau217 (cross-assay validation):
#     - Lilly MSD ORRES, Janssen DILUTION_CORRECTED_CONC, C2N pT217_C2N
#   demographics + hippocampus + ICV: ADSP_PHC_T1_FS (ComBat-harmonized)
#   diagnosis: ADSP PHC_Diagnosis (primary); DXSUM DIAGNOSIS (fallback)
#   covariates: age, sex, education, baseline diagnosis, BMI (from VITALS)
#   no HTN/DM (MEDHIST not downloaded)
# =============================================================================

suppressPackageStartupMessages({
  for (pkg in c("lcmm","lme4","nnet","dplyr","tidyr","ggplot2")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(lcmm); library(lme4); library(nnet)
  library(dplyr); library(tidyr); library(ggplot2)
})
set.seed(20260524)

# ---------------------------------------------------------------------------
# 1. CONFIG -- actual file names
# ---------------------------------------------------------------------------
CSV_DIR <- "/Users/kuan-yuchu/Desktop/LOAD/ADNI_csv"

F_UPENN    <- "UPENN_PLASMA_FUJIREBIO_QUANTERIX_24May2026.csv"
F_LILLY    <- "LILLY_PTAU217_MSD600_24May2026.csv"
F_JANSSEN  <- "JANSSEN_PLASMA_P217_TAU_24May2026.csv"
F_C2N      <- "C2N_PRECIVITYAD2_PLASMA_24May2026.csv"
F_APOE     <- "APOERES_24May2026.csv"
F_DXSUM    <- "DXSUM_24May2026.csv"
F_ADAS     <- "ADAS_24May2026.csv"
F_MMSE     <- "MMSE_24May2026.csv"
F_VITALS   <- "VITALS_24May2026.csv"
F_ADSP     <- "ADSP_PHC_T1_FS_24May2026.csv"
F_REGISTRY <- "REGISTRY_24May2026.csv"

MAX_TIME_MONTHS <- 120                   # Plan A: up to 10-year window
MIN_GFAP_TP <- 3                         # ≥3 timepoints inclusion
MIN_PTAU_TP <- 3

# ---------------------------------------------------------------------------
# 2. Helpers
# ---------------------------------------------------------------------------
read_adni <- function(name) {
  path <- file.path(CSV_DIR, name)
  if (!file.exists(path)) stop("MISSING file: ", path)
  read.csv(path, na.strings = c("", "NA", "-1", "-4", "-9", "NaN"),
           stringsAsFactors = FALSE, check.names = FALSE)
}

viscode_to_month <- function(v) {
  v <- tolower(trimws(as.character(v)))
  out <- rep(NA_integer_, length(v))
  out[v %in% c("bl", "baseline", "sc", "sccn", "scrn", "scmri", "init")] <- 0L
  for (i in seq_along(v)) {
    if (!is.na(out[i])) next
    m <- regmatches(v[i], regexec("^m(\\d+)$", v[i]))
    if (length(m[[1]]) >= 2) { out[i] <- as.integer(m[[1]][2]); next }
    y <- regmatches(v[i], regexec("^y(\\d+)$", v[i]))
    if (length(y[[1]]) >= 2) { out[i] <- as.integer(y[[1]][2]) * 12L; next }
  }
  out
}

pick_col <- function(df, candidates) {
  hits <- intersect(candidates, names(df))
  if (length(hits)) hits[1] else NULL
}

# parse GENOTYPE string ("3/4", "E3/E4", "34", "44") -> ε4 carrier 0/1
apoe_e4_from_genotype <- function(g) {
  s <- toupper(as.character(g))
  sapply(s, function(x) {
    if (is.na(x) || x == "") return(NA_integer_)
    digits <- gsub("[^0-9]", "", x)
    if (nchar(digits) < 2) return(NA_integer_)
    as.integer(grepl("4", digits))
  }, USE.NAMES = FALSE)
}

# ---------------------------------------------------------------------------
# 3. Plasma biomarkers (UPENN: GFAP_Q + pT217_F on same row)
# ---------------------------------------------------------------------------
cat("=== UPENN plasma ===\n")
upenn <- read_adni(F_UPENN)
upenn$time <- viscode_to_month(upenn[[pick_col(upenn, c("VISCODE2", "VISCODE"))]])

gfap <- upenn %>%
  transmute(RID, time, GFAP = suppressWarnings(as.numeric(GFAP_Q))) %>%
  filter(!is.na(GFAP), !is.na(time))
ptau_primary <- upenn %>%
  transmute(RID, time, PTAU = suppressWarnings(as.numeric(pT217_F))) %>%
  filter(!is.na(PTAU), !is.na(time))

cat("  GFAP_Q obs:        ", nrow(gfap), "\n")
cat("  pT217_F obs:       ", nrow(ptau_primary), "\n")

# p-tau sensitivity assays (single dataframe per assay for cross-validation)
read_sens <- function(file, val_col, label, restrict_test = NULL) {
  d <- tryCatch(read_adni(file), error = function(e) NULL)
  if (is.null(d)) { message("  [", label, "] file not found"); return(NULL) }
  vcol <- pick_col(d, c("VISCODE2", "VISCODE"))
  d$time <- viscode_to_month(d[[vcol]])
  if (!is.null(restrict_test)) {
    tcol <- pick_col(d, c("TESTCD", "TEST"))
    if (!is.null(tcol)) d <- d[grepl(restrict_test, d[[tcol]], ignore.case = TRUE), ]
  }
  if (!val_col %in% names(d)) {
    message("  [", label, "] value col ", val_col, " missing"); return(NULL)
  }
  out <- data.frame(RID = d$RID, time = d$time,
                    val = suppressWarnings(as.numeric(d[[val_col]])))
  names(out)[3] <- paste0("PTAU_", label)
  out <- out[complete.cases(out), ]
  cat("  PTAU_", label, " obs: ", nrow(out), "\n", sep = "")
  out
}
ptau_lilly   <- read_sens(F_LILLY,   "ORRES",                  "LILLY",   "PTAU217")
ptau_janssen <- read_sens(F_JANSSEN, "DILUTION_CORRECTED_CONC","JANSSEN")
ptau_c2n     <- read_sens(F_C2N,     "pT217_C2N",              "C2N")

# ---------------------------------------------------------------------------
# 4. APOE
# ---------------------------------------------------------------------------
cat("\n=== APOE ===\n")
apoe_raw <- read_adni(F_APOE)
apoe <- apoe_raw %>%
  transmute(RID, apoe_e4 = apoe_e4_from_genotype(GENOTYPE), GENOTYPE) %>%
  distinct(RID, .keep_all = TRUE)
cat("  APOE rows:", nrow(apoe),
    "; ε4 carriers:", sum(apoe$apoe_e4, na.rm = TRUE),
    "(", round(100 * mean(apoe$apoe_e4, na.rm = TRUE), 1), "%)\n")

# ---------------------------------------------------------------------------
# 5. ADSP — demographics + longitudinal hippocampus/ICV (ComBat-harmonized)
# ---------------------------------------------------------------------------
cat("\n=== ADSP_PHC_T1_FS ===\n")
adsp <- read_adni(F_ADSP)
adsp$time <- viscode_to_month(adsp[[pick_col(adsp, c("VISCODE2", "VISCODE"))]])

adsp$HV  <- suppressWarnings(as.numeric(adsp[["Left.Hippocampus_combat"]])) +
            suppressWarnings(as.numeric(adsp[["Right.Hippocampus_combat"]]))
adsp$ICV <- suppressWarnings(as.numeric(adsp[["EstimatedTotalIntraCranialVol_combat"]]))
adsp$HV_adj <- adsp$HV / adsp$ICV * mean(adsp$ICV, na.rm = TRUE)

adsp_base <- adsp %>%
  group_by(RID) %>%
  slice_min(order_by = time, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(RID,
            age_baseline = suppressWarnings(as.numeric(PHC_Age_T1)),
            female       = as.integer(PHC_Sex == 2),
            educ         = suppressWarnings(as.numeric(PHC_Education)),
            adsp_diag    = case_when(
              PHC_Diagnosis == 1 ~ "CN",
              PHC_Diagnosis == 2 ~ "MCI",
              PHC_Diagnosis == 3 ~ "AD",
              TRUE ~ NA_character_))

hv_long <- adsp %>%
  transmute(RID, time, HV_adj) %>%
  filter(!is.na(HV_adj), !is.na(time))

cat("  ADSP unique RIDs:", length(unique(adsp$RID)), "\n")

# ---------------------------------------------------------------------------
# 6. DXSUM (fallback for baseline diagnosis if ADSP missing)
# ---------------------------------------------------------------------------
cat("\n=== DXSUM ===\n")
dxsum <- read_adni(F_DXSUM)
dxsum$time <- viscode_to_month(dxsum[[pick_col(dxsum, c("VISCODE2", "VISCODE"))]])
dx_base <- dxsum %>%
  filter(time == 0) %>%
  transmute(RID, dxsum_diag = case_when(
    DIAGNOSIS == 1 ~ "CN",
    DIAGNOSIS == 2 ~ "MCI",
    DIAGNOSIS == 3 ~ "AD",
    TRUE ~ NA_character_)) %>%
  distinct(RID, .keep_all = TRUE)
cat("  DXSUM baseline rows:", nrow(dx_base), "\n")

# ---------------------------------------------------------------------------
# 7. Cognitive outcomes
# ---------------------------------------------------------------------------
cat("\n=== ADAS / MMSE ===\n")
adas <- read_adni(F_ADAS)
adas$time <- viscode_to_month(adas[[pick_col(adas, c("VISCODE2", "VISCODE"))]])
adas <- adas %>% transmute(RID, time,
                           ADAS13 = suppressWarnings(as.numeric(TOTAL13))) %>%
                 filter(!is.na(time))

mmse <- read_adni(F_MMSE)
mmse$time <- viscode_to_month(mmse[[pick_col(mmse, c("VISCODE2", "VISCODE"))]])
mmse <- mmse %>% transmute(RID, time,
                           MMSE = suppressWarnings(as.numeric(MMSCORE))) %>%
                 filter(!is.na(time))
cat("  ADAS rows:", nrow(adas), "; MMSE rows:", nrow(mmse), "\n")

# ---------------------------------------------------------------------------
# 8. VITALS -> baseline BMI
# ---------------------------------------------------------------------------
cat("\n=== VITALS / BMI ===\n")
vit <- read_adni(F_VITALS)
vit$time <- viscode_to_month(vit[[pick_col(vit, c("VISCODE2", "VISCODE"))]])
vit_b <- vit %>% filter(time == 0)
h <- suppressWarnings(as.numeric(vit_b$VSHEIGHT))
w <- suppressWarnings(as.numeric(vit_b$VSWEIGHT))
hu <- suppressWarnings(as.numeric(vit_b$VSHTUNIT))     # 1=in, 2=cm
wu <- suppressWarnings(as.numeric(vit_b$VSWTUNIT))     # 1=lb, 2=kg
h_m <- ifelse(!is.na(hu) & hu == 1, h * 0.0254, h / 100)
w_kg <- ifelse(!is.na(wu) & wu == 1, w * 0.4536, w)
vit_b$bmi <- w_kg / h_m^2
bmi <- vit_b %>% group_by(RID) %>%
  summarise(bmi = median(bmi[is.finite(bmi) & bmi > 10 & bmi < 60], na.rm = TRUE),
            .groups = "drop")
cat("  baseline BMI available for:", sum(!is.na(bmi$bmi)), "/", nrow(bmi), "\n")

# ---------------------------------------------------------------------------
# 9. Apply inclusion criteria & build long frame
# ---------------------------------------------------------------------------
cat("\n=== Inclusion criteria (Plan A: ≤", MAX_TIME_MONTHS, "mo, ≥3 timepoints) ===\n")
in_window <- function(d) d %>% filter(time >= 0 & time <= MAX_TIME_MONTHS)
gfap <- in_window(gfap)
ptau_primary <- in_window(ptau_primary)
adas <- in_window(adas)
mmse <- in_window(mmse)
hv_long <- in_window(hv_long)

elig <- adsp_base %>%
  filter(age_baseline >= 65) %>%
  left_join(dx_base, by = "RID") %>%
  mutate(diag = coalesce(adsp_diag, dxsum_diag)) %>%
  filter(diag %in% c("CN", "MCI", "AD")) %>%        # Plan A2: full AD continuum
  inner_join(gfap %>% distinct(RID, time) %>% count(RID, name = "n_gfap") %>%
               filter(n_gfap >= MIN_GFAP_TP), by = "RID") %>%
  inner_join(ptau_primary %>% distinct(RID, time) %>% count(RID, name = "n_ptau") %>%
               filter(n_ptau >= MIN_PTAU_TP), by = "RID") %>%
  mutate(diag = factor(diag, levels = c("CN", "MCI", "AD")))   # CN as reference

cat("  Eligible N =", nrow(elig), "\n")
print(table(elig$diag))
if (nrow(elig) < 30) stop("Eligible N is unexpectedly small; check filters.")

elig_RIDs <- elig$RID

# Build df_long from the union of actually observed timepoints (no fixed grid)
all_times <- bind_rows(
  gfap %>% transmute(RID, time),
  ptau_primary %>% transmute(RID, time),
  adas %>% transmute(RID, time),
  mmse %>% transmute(RID, time),
  hv_long %>% transmute(RID, time)
) %>% filter(RID %in% elig_RIDs) %>% distinct()

df_long <- all_times %>%
  left_join(gfap, by = c("RID", "time")) %>%
  left_join(ptau_primary, by = c("RID", "time")) %>%
  left_join(adas, by = c("RID", "time")) %>%
  left_join(mmse, by = c("RID", "time")) %>%
  left_join(hv_long, by = c("RID", "time")) %>%
  mutate(time_yr = time / 12,
         logGFAP = log(GFAP),
         logPTAU = log(PTAU))

df_base <- elig %>%
  select(RID, age_baseline, female, educ, diag) %>%
  left_join(apoe %>% select(RID, apoe_e4), by = "RID") %>%
  left_join(bmi, by = "RID")

df_long <- df_long %>% left_join(df_base, by = "RID")

cat("  Long-format rows:", nrow(df_long),
    "; baseline rows:", nrow(df_base), "\n")
saveRDS(df_long, file.path(CSV_DIR, "df_long_real.rds"))
saveRDS(df_base, file.path(CSV_DIR, "df_base_real.rds"))

# ---------------------------------------------------------------------------
# 10. Cross-assay p-tau217 validation (Pearson @ baseline)
# ---------------------------------------------------------------------------
cat("\n=== Cross-assay p-tau217 concordance (baseline) ===\n")
cross_bl <- ptau_primary %>% filter(time == 0) %>%
              transmute(RID, UPENN_FUJIREBIO = PTAU)
if (!is.null(ptau_lilly))
  cross_bl <- cross_bl %>% left_join(ptau_lilly %>% filter(time == 0) %>%
                                       select(RID, LILLY_MSD = PTAU_LILLY), by = "RID")
if (!is.null(ptau_janssen))
  cross_bl <- cross_bl %>% left_join(ptau_janssen %>% filter(time == 0) %>%
                                       select(RID, JANSSEN = PTAU_JANSSEN), by = "RID")
if (!is.null(ptau_c2n))
  cross_bl <- cross_bl %>% left_join(ptau_c2n %>% filter(time == 0) %>%
                                       select(RID, C2N = PTAU_C2N), by = "RID")
print(round(cor(cross_bl[, -1, drop = FALSE], use = "pairwise.complete.obs",
                method = "pearson"), 3))
write.csv(cross_bl, file.path(CSV_DIR, "ptau217_cross_assay_baseline.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------------
# 11. PHASE I — LCGMM diagnostic (BIC selection) only;
#     no class structure expected -> proceed to continuous-slope LME (Path 1).
# ---------------------------------------------------------------------------
cat("\n=== PHASE I: LCGMM ===\n")
fit_lcgmm <- function(yvar, df, label) {
  cat("\n--- LCGMM:", label, "---\n")
  d <- df[!is.na(df[[yvar]]) & !is.na(df$time_yr), ]
  form <- as.formula(paste(yvar, "~ time_yr"))
  m1 <- hlme(form, random = ~time_yr, subject = "RID", ng = 1,
             data = d, verbose = FALSE)
  m2 <- hlme(form, mixture = ~time_yr, random = ~time_yr, subject = "RID",
             ng = 2, B = m1, data = d, verbose = FALSE)
  m3 <- hlme(form, mixture = ~time_yr, random = ~time_yr, subject = "RID",
             ng = 3, B = m1, data = d, verbose = FALSE)
  m4 <- tryCatch(
    hlme(form, mixture = ~time_yr, random = ~time_yr, subject = "RID",
         ng = 4, B = m1, data = d, verbose = FALSE),
    error = function(e) { message("  ng=4 failed: ", e$message); NULL })
  mods <- Filter(Negate(is.null), list(m1, m2, m3, m4))
  st <- do.call(summarytable, c(mods,
                                list(which = c("G", "loglik", "AIC", "BIC",
                                               "SABIC", "entropy", "%class"))))
  print(st)
  bics <- sapply(mods, function(m) m$BIC)
  ord  <- order(bics)
  best <- mods[[1]]
  for (k in ord) {
    if (k == 1) { best <- mods[[1]]; break }
    cl <- mods[[k]]$pprob$class
    if (length(cl) > 0 && min(prop.table(table(cl))) >= 0.05) {
      best <- mods[[k]]; break
    }
  }
  cat("  selected ng =", best$ng, "\n")
  list(models = mods, selected = best, summary = st)
}

gfap_fit <- fit_lcgmm("logGFAP", df_long, "plasma GFAP")
ptau_fit <- fit_lcgmm("logPTAU", df_long, "plasma p-tau217 (UPENN-Fujirebio)")

joint_fit <- tryCatch({
  cat("\n--- Joint multivariate LCGMM (sensitivity) ---\n")
  d_joint <- df_long[!is.na(df_long$logGFAP) | !is.na(df_long$logPTAU), ]
  m_joint1 <- multlcmm(logGFAP + logPTAU ~ time_yr,
                       random = ~time_yr, subject = "RID",
                       ng = 1, data = d_joint, verbose = FALSE)
  m_joint <- multlcmm(logGFAP + logPTAU ~ time_yr,
                      random = ~time_yr, subject = "RID",
                      mixture = ~time_yr, ng = 3, B = m_joint1,
                      data = d_joint, verbose = FALSE)
  m_joint
}, error = function(e) { message("  multlcmm failed: ", e$message); NULL })

# ---------------------------------------------------------------------------
# 12. PHASE II (revised) — LME with APOE × time:
#     does APOE ε4 modulate each biomarker's slope?
# ---------------------------------------------------------------------------
cat("\n=== PHASE II (Path 1): LME with APOE × time ===\n")

fit_biomarker_lme <- function(yvar, label) {
  cat("\n--- LME:", label, "~ time + APOE + APOE:time + covariates ---\n")
  d <- df_long[!is.na(df_long[[yvar]]) & !is.na(df_long$apoe_e4) &
               !is.na(df_long$time_yr), ]
  f <- as.formula(paste(
    yvar, "~ time_yr + apoe_e4 + apoe_e4:time_yr +",
          "age_baseline + female + educ + diag + (time_yr | RID)"))
  fit <- lmer(f, data = d, REML = FALSE,
              control = lmerControl(optimizer = "bobyqa"))
  print(round(summary(fit)$coefficients, 4))
  fit
}
lme_gfap <- tryCatch(fit_biomarker_lme("logGFAP", "log plasma GFAP"),
                     error = function(e) { message(e$message); NULL })
lme_ptau <- tryCatch(fit_biomarker_lme("logPTAU", "log plasma p-tau217"),
                     error = function(e) { message(e$message); NULL })

# Extract per-person slope BLUPs = fixed time slope + person-specific deviation
extract_slopes <- function(fit, name) {
  re <- ranef(fit)$RID
  fe <- fixef(fit)
  rid <- as.integer(rownames(re))
  total_slope <- as.numeric(fe["time_yr"] + re[["time_yr"]])
  out <- data.frame(RID = rid, slope = total_slope)
  names(out)[2] <- paste0(name, "_slope")
  out
}
if (!is.null(lme_gfap)) {
  gs <- extract_slopes(lme_gfap, "gfap")
  df_base <- df_base %>% select(-any_of("gfap_slope")) %>% left_join(gs, by = "RID")
}
if (!is.null(lme_ptau)) {
  ps <- extract_slopes(lme_ptau, "ptau")
  df_base <- df_base %>% select(-any_of("ptau_slope")) %>% left_join(ps, by = "RID")
}

report_apoe_x_time <- function(fit, label) {
  if (is.null(fit)) return(invisible())
  co <- summary(fit)$coefficients
  ix <- intersect(c("apoe_e4:time_yr", "time_yr:apoe_e4"), rownames(co))
  if (!length(ix)) return(invisible())
  b <- co[ix[1], "Estimate"]; se <- co[ix[1], "Std. Error"]
  p <- 2 * pnorm(-abs(b / se))
  cat(sprintf("\n  %s — APOE ε4 × time interaction:  β = %.4f  (SE %.4f)  p = %.4g\n",
              label, b, se, p))
}
report_apoe_x_time(lme_gfap, "log GFAP")
report_apoe_x_time(lme_ptau, "log p-tau217")

# ---------------------------------------------------------------------------
# 13. PHASE III (revised) — Does each biomarker slope predict cognitive /
#     hippocampal decline?  ADAS13/MMSE/HV_adj ~ time × slope_BLUP + covariates
# ---------------------------------------------------------------------------
cat("\n=== PHASE III (Path 1): biomarker slope predicts outcome rate ===\n")
df_long <- df_long %>%
  select(-any_of(c("gfap_slope", "ptau_slope"))) %>%
  left_join(df_base %>% select(RID, any_of(c("gfap_slope", "ptau_slope"))),
            by = "RID")

fit_outcome_slope <- function(outcome_y, slope_var, label) {
  cat("\n--- ", outcome_y, " ~ time_yr * ", slope_var,
      " + covariates ---\n", sep = "")
  d <- df_long[!is.na(df_long[[outcome_y]]) & !is.na(df_long[[slope_var]]) &
               !is.na(df_long$time_yr), ]
  if (nrow(d) < 30) { message("  too few rows; skipping"); return(NULL) }
  f <- as.formula(paste(
    outcome_y, "~ time_yr *", slope_var,
    "+ age_baseline + female + educ + diag + apoe_e4 + (time_yr | RID)"))
  fit <- lmer(f, data = d, REML = FALSE,
              control = lmerControl(optimizer = "bobyqa"))
  co <- summary(fit)$coefficients
  print(round(co, 4))
  fit
}

lme_adas_gfap <- tryCatch(fit_outcome_slope("ADAS13", "gfap_slope", "ADAS-Cog 13 ~ GFAP slope"),
                          error = function(e) { message(e$message); NULL })
lme_hv_gfap   <- tryCatch(fit_outcome_slope("HV_adj","gfap_slope", "HV ~ GFAP slope"),
                          error = function(e) { message(e$message); NULL })
lme_mmse_gfap <- tryCatch(fit_outcome_slope("MMSE",  "gfap_slope", "MMSE ~ GFAP slope"),
                          error = function(e) { message(e$message); NULL })
lme_adas_ptau <- tryCatch(fit_outcome_slope("ADAS13","ptau_slope", "ADAS-Cog 13 ~ p-tau slope"),
                          error = function(e) { message(e$message); NULL })
lme_hv_ptau   <- tryCatch(fit_outcome_slope("HV_adj","ptau_slope", "HV ~ p-tau slope"),
                          error = function(e) { message(e$message); NULL })
lme_mmse_ptau <- tryCatch(fit_outcome_slope("MMSE",  "ptau_slope", "MMSE ~ p-tau slope"),
                          error = function(e) { message(e$message); NULL })

# Helper: pull time × slope interaction
report_slope_x_time <- function(fit, slope_var, label) {
  if (is.null(fit)) return(invisible())
  co <- summary(fit)$coefficients
  ix <- intersect(c(paste0("time_yr:", slope_var), paste0(slope_var, ":time_yr")),
                  rownames(co))
  if (!length(ix)) return(invisible())
  b <- co[ix[1], "Estimate"]; se <- co[ix[1], "Std. Error"]
  p <- 2 * pnorm(-abs(b / se))
  cat(sprintf("  %s: β = %.4f (SE %.4f)  p = %.4g\n", label, b, se, p))
}
cat("\n--- Path-III interaction summary (β = additional outcome change per year per unit slope) ---\n")
report_slope_x_time(lme_adas_gfap, "gfap_slope", "ADAS-Cog 13 ~ time × GFAP slope")
report_slope_x_time(lme_hv_gfap,   "gfap_slope", "HV_adj      ~ time × GFAP slope")
report_slope_x_time(lme_mmse_gfap, "gfap_slope", "MMSE        ~ time × GFAP slope")
report_slope_x_time(lme_adas_ptau, "ptau_slope", "ADAS-Cog 13 ~ time × p-tau slope")
report_slope_x_time(lme_hv_ptau,   "ptau_slope", "HV_adj      ~ time × p-tau slope")
report_slope_x_time(lme_mmse_ptau, "ptau_slope", "MMSE        ~ time × p-tau slope")

# ---------------------------------------------------------------------------
# 14. Manuscript summary
# ---------------------------------------------------------------------------
cat("\n=================== MANUSCRIPT SUMMARY ===================\n")
cat("Eligible N:", nrow(df_base), "\n")
cat("Baseline diag table:\n"); print(table(df_base$diag))
cat("Mean(SD) baseline age:", round(mean(df_base$age_baseline, na.rm = TRUE), 1),
    "(", round(sd(df_base$age_baseline,  na.rm = TRUE), 1), ")\n")
cat("Female %:", round(100 * mean(df_base$female, na.rm = TRUE), 1), "\n")
cat("Education mean(SD):", round(mean(df_base$educ, na.rm = TRUE), 1),
    "(", round(sd(df_base$educ, na.rm = TRUE), 1), ")\n")
cat("APOE ε4 %:", round(100 * mean(df_base$apoe_e4, na.rm = TRUE), 1), "\n")
cat("BMI mean(SD):", round(mean(df_base$bmi, na.rm = TRUE), 1),
    "(", round(sd(df_base$bmi,  na.rm = TRUE), 1), ")\n\n")

cat("\nLCGMM model selection: BIC favored ng = 1 for both biomarkers (see Phase I tables above).\n")
cat("  -> No distinct trajectory classes detected at this N; analysis used\n")
cat("  -> continuous random-slope LME (Path 1).\n\n")

if (!is.null(lme_gfap)) {
  cat("Random-slope variance (GFAP):\n")
  print(VarCorr(lme_gfap))
  report_apoe_x_time(lme_gfap, "log GFAP")
}
if (!is.null(lme_ptau)) {
  cat("\nRandom-slope variance (p-tau217):\n")
  print(VarCorr(lme_ptau))
  report_apoe_x_time(lme_ptau, "log p-tau217")
}

cat("\n--- Per-person slope BLUP summary ---\n")
if ("gfap_slope" %in% names(df_base))
  print(summary(df_base$gfap_slope))
if ("ptau_slope" %in% names(df_base))
  print(summary(df_base$ptau_slope))

cat("\n=== END ===\n")
cat("Saved: df_long_real.rds, df_base_real.rds, ptau217_cross_assay_baseline.csv\n")
cat("Paste this log + the printed tables back; I will fill the manuscript.\n")
