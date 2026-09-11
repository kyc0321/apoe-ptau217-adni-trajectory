# ========================================================================
# ADNI Plasma GFAP / p-tau217 — Master Re-run Pipeline (v1)
# ------------------------------------------------------------------------
# Author : Kuan-Yu Chu
# Date   : 2026-05-25
# SAP    : /Mauuscripts/SAP_v1.md (locked 2026-02-12; amended 2026-05-25)
# Data   : ADNI 2026-05-24 release
# Inclusion: age >=65, dx in {CN,MCI,AD}, APOE available,
#            >=2 plasma GFAP AND >=2 plasma p-tau217 (UPENN-Fujirebio)
# Phases:
#   I   LCGMM (exploratory only)
#   II  LME (primary): APOE x time on each biomarker
#   III JMbayes2 joint models (primary): biomarker slope -> outcomes
#   IV  Cross-assay p-tau217 concordance
# ========================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(ggplot2)
  library(data.table); library(lme4); library(lmerTest); library(nlme)
  library(lcmm); library(JMbayes2)
})
set.seed(20260525)

DATA <- "/Users/kuan-yuchu/Desktop/LOAD/ADNI_csv"
OUT  <- "/Users/kuan-yuchu/Desktop/LOAD/Analysis/outputs"
LOG  <- "/Users/kuan-yuchu/Desktop/LOAD/Analysis/logs"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(LOG, showWarnings = FALSE, recursive = TRUE)

rc <- function(f) utils::read.csv(file.path(DATA, f),
                                  stringsAsFactors = FALSE, na.strings = c("", "NA"))
wc <- function(x, f) utils::write.csv(x, file.path(OUT, f), row.names = FALSE)

cat("==========================================\n",
    " MASTER PIPELINE v1  ", as.character(Sys.time()), "\n",
    "==========================================\n", sep = "")

# ------------------------------------------------------------------------
# 1. LOAD DATA
# ------------------------------------------------------------------------
upenn  <- rc("UPENN_PLASMA_FUJIREBIO_QUANTERIX_24May2026.csv")
apoe   <- rc("APOERES_24May2026.csv")
vitals <- rc("VITALS_24May2026.csv")
adas   <- rc("ADAS_24May2026.csv")
mmse   <- rc("MMSE_24May2026.csv")
fs     <- rc("ADSP_PHC_T1_FS_24May2026.csv")
labd   <- rc("LABDATA_25May2026.csv")
load("/tmp/adnimerge_unpack/ADNIMERGE/data/adnimerge.rdata")
comorb <- utils::read.csv(file.path(OUT, "comorbidities.csv"),
                          stringsAsFactors = FALSE)

# Optional cross-assay
lilly  <- tryCatch(rc("LILLY_PTAU217_MSD600_24May2026.csv"), error = function(e) NULL)
c2n    <- tryCatch(rc("C2N_PRECIVITYAD2_PLASMA_24May2026.csv"), error = function(e) NULL)
jans   <- tryCatch(rc("JANSSEN_PLASMA_P217_TAU_24May2026.csv"), error = function(e) NULL)

# ------------------------------------------------------------------------
# 2. PLASMA BIOMARKERS — long form
# ------------------------------------------------------------------------
upenn <- upenn %>%
  mutate(RID = as.integer(RID),
         EXAMDATE = as.Date(EXAMDATE),
         pT217_F = suppressWarnings(as.numeric(pT217_F)),
         GFAP_Q  = suppressWarnings(as.numeric(GFAP_Q)),
         pT217_F = ifelse(pT217_F <= 0, NA, pT217_F),
         GFAP_Q  = ifelse(GFAP_Q  <= 0, NA, GFAP_Q),
         log_ptau = log(pT217_F),
         log_gfap = log(GFAP_Q))

bl_dates <- upenn %>% group_by(RID) %>%
  summarise(bl_date = min(EXAMDATE[!is.na(pT217_F) | !is.na(GFAP_Q)], na.rm = TRUE),
            .groups = "drop")
upenn <- upenn %>% left_join(bl_dates, by = "RID") %>%
  mutate(time_yrs = as.numeric(EXAMDATE - bl_date) / 365.25)

counts <- upenn %>% group_by(RID) %>%
  summarise(n_gfap = sum(!is.na(GFAP_Q)),
            n_ptau = sum(!is.na(pT217_F)), .groups = "drop")

# ------------------------------------------------------------------------
# 3. APOE + baseline demographics
# ------------------------------------------------------------------------
apoe <- apoe %>% mutate(RID = as.integer(RID)) %>%
  filter(!is.na(GENOTYPE), GENOTYPE != "") %>%
  distinct(RID, .keep_all = TRUE) %>%
  mutate(APOE4 = as.integer(str_count(GENOTYPE, "4") >= 1)) %>%
  select(RID, GENOTYPE, APOE4)

adni_bl <- adnimerge %>% filter(VISCODE == "bl") %>%
  transmute(RID = as.integer(RID),
            age_bl = as.numeric(AGE),
            sex = PTGENDER,
            edu = as.numeric(PTEDUCAT),
            dx_bl = DX.bl) %>%
  distinct(RID, .keep_all = TRUE) %>%
  mutate(dx3 = case_when(
    dx_bl %in% c("CN","SMC") ~ "CN",
    dx_bl %in% c("EMCI","LMCI","MCI") ~ "MCI",
    dx_bl %in% c("AD","Dementia") ~ "AD",
    TRUE ~ NA_character_))

# ------------------------------------------------------------------------
# 4. TIME-VARYING BMI (data.table rolling join within 365 days)
# ------------------------------------------------------------------------
vit <- vitals %>% mutate(
  RID = as.integer(RID),
  EXAMDATE = as.Date(VISDATE),
  VSWEIGHT = suppressWarnings(as.numeric(VSWEIGHT)),
  VSHEIGHT = suppressWarnings(as.numeric(VSHEIGHT)),
  VSWTUNIT = as.integer(VSWTUNIT),
  VSHTUNIT = as.integer(VSHTUNIT),
  wt_kg = ifelse(!is.na(VSWTUNIT) & VSWTUNIT == 1, VSWEIGHT*0.4536, VSWEIGHT),
  ht_m  = ifelse(!is.na(VSHTUNIT) & VSHTUNIT == 1, VSHEIGHT*0.0254, VSHEIGHT/100),
  BMI   = wt_kg/(ht_m^2)
) %>% filter(BMI > 10, BMI < 70, !is.na(EXAMDATE)) %>%
  select(RID, EXAMDATE, BMI)

dt_vit <- as.data.table(vit); setkey(dt_vit, RID, EXAMDATE)

attach_tv <- function(biom_df, src_dt, val_col, win_days = 365) {
  # For each biom row, find closest src row within +/- win_days (either direction)
  dt_b <- as.data.table(biom_df)[, .(RID, EXAMDATE, .row = .I)]
  setkey(dt_b, RID, EXAMDATE)
  fwd <- src_dt[dt_b, on = c("RID","EXAMDATE"), roll = win_days,  rollends = c(TRUE,TRUE)]
  bwd <- src_dt[dt_b, on = c("RID","EXAMDATE"), roll = -win_days, rollends = c(TRUE,TRUE)]
  # Pick whichever exists (preferring smaller absolute date diff)
  out <- biom_df
  vfwd <- fwd[[val_col]]
  vbwd <- bwd[[val_col]]
  picked <- ifelse(is.na(vfwd), vbwd,
            ifelse(is.na(vbwd), vfwd,
                   ifelse(abs(as.numeric(fwd$EXAMDATE - dt_b$EXAMDATE)) <=
                          abs(as.numeric(bwd$EXAMDATE - dt_b$EXAMDATE)), vfwd, vbwd)))
  out[[paste0(val_col, "_tv")]] <- picked
  out
}

biom_long <- upenn %>% select(RID, EXAMDATE, time_yrs, log_gfap, log_ptau, pT217_F, GFAP_Q)
# Baseline BMI per RID (closest VITALS within +/-2y of biomarker baseline)
bmi_bl <- biom_long %>% distinct(RID) %>%
  left_join(bl_dates, by = "RID") %>%
  rowwise() %>%
  mutate(BMI_bl = {
    cand <- vit[vit$RID == RID & abs(as.numeric(vit$EXAMDATE - bl_date)) <= 730, ]
    if (nrow(cand) == 0) NA_real_ else cand$BMI[which.min(abs(as.numeric(cand$EXAMDATE - bl_date)))]
  }) %>% ungroup() %>% select(RID, BMI_bl)
biom_long <- biom_long %>% left_join(bmi_bl, by = "RID")
biom_long$BMI_tv <- biom_long$BMI_bl   # alias for downstream code expecting BMI_tv
cat("RIDs with baseline BMI:", sum(!is.na(bmi_bl$BMI_bl)), "/", nrow(bmi_bl), "\n")

# ------------------------------------------------------------------------
# 5. eGFR — CKD-EPI 2021 race-free, time-varying
# ------------------------------------------------------------------------
crea <- labd %>% transmute(
  RID = as.integer(RID),
  EXAMDATE = as.Date(EXAMDATE),
  scr = suppressWarnings(as.numeric(RCT392))
) %>% filter(!is.na(EXAMDATE), !is.na(scr), scr > 0, scr < 15)

cat("Creatinine rows:", nrow(crea), "; unique RIDs:", length(unique(crea$RID)), "\n")

# eGFR strategy: BASELINE-ONLY per RID (most-recent creatinine within +/-2 years of bl_date).
# Time-varying eGFR is impractical: ADNI creatinine only at ~8% of biomarker visits.
# This follows standard practice in AD plasma-biomarker trajectory papers.
crea_bl <- biom_long %>% distinct(RID) %>%
  left_join(bl_dates, by = "RID") %>%
  rowwise() %>%
  mutate(
    scr_bl = {
      cand <- crea[crea$RID == RID & abs(as.numeric(crea$EXAMDATE - bl_date)) <= 730, ]
      if (nrow(cand) == 0) NA_real_ else cand$scr[which.min(abs(as.numeric(cand$EXAMDATE - bl_date)))]
    }
  ) %>% ungroup() %>% select(RID, scr_bl)

biom_long <- biom_long %>%
  left_join(crea_bl, by = "RID") %>%
  left_join(adni_bl %>% select(RID, age_bl, sex), by = "RID") %>%
  mutate(
    is_female = as.integer(sex == "Female"),
    k = ifelse(is_female == 1, 0.7, 0.9),
    a = ifelse(is_female == 1, -0.241, -0.302),
    eGFR_bl = 142 *
      pmin(scr_bl/k, 1)^a *
      pmax(scr_bl/k, 1)^(-1.200) *
      (0.9938^age_bl) *
      ifelse(is_female == 1, 1.012, 1)
  )
n_rid_egfr <- biom_long %>% distinct(RID, eGFR_bl) %>% filter(!is.na(eGFR_bl)) %>% nrow()
cat("RIDs with baseline eGFR (within 2y of biomarker bl):", n_rid_egfr, "/",
    length(unique(biom_long$RID)), "\n")

# ------------------------------------------------------------------------
# 6. APPLY INCLUSION (SAP §4 revised: >=2/>=2)
# ------------------------------------------------------------------------
elig <- counts %>%
  filter(n_gfap >= 2, n_ptau >= 2) %>%
  inner_join(apoe, by = "RID") %>%
  inner_join(adni_bl, by = "RID") %>%
  filter(age_bl >= 65, !is.na(dx3))

analytic_rid <- elig$RID
cat("\n=== INCLUSION RESULT ===\n")
cat("Analytic N =", length(analytic_rid), "\n")
print(table(elig$dx3, useNA = "ifany"))

# Restrict biom_long to analytic
biom_long <- biom_long %>% filter(RID %in% analytic_rid)

# Build per-participant attributes
attr_tbl <- elig %>%
  left_join(comorb, by = "RID") %>%
  mutate(HTN_base    = ifelse(is.na(HTN_base), 0L, HTN_base),
         DM_base     = ifelse(is.na(DM_base),  0L, DM_base),
         stroke_ever = ifelse(is.na(stroke_ever), 0L, stroke_ever)) %>%
  select(RID, age_bl, sex, edu, dx3, APOE4, HTN_base, DM_base, stroke_ever, n_gfap, n_ptau)

# Merge attr to long
biom_long <- biom_long %>% select(-age_bl, -sex) %>%
  left_join(attr_tbl, by = "RID")

# ------------------------------------------------------------------------
# 7. TABLE 1 — Baseline characteristics
# ------------------------------------------------------------------------
bl_biom <- biom_long %>% group_by(RID) %>%
  filter(time_yrs == min(time_yrs)) %>% ungroup() %>%
  select(RID, GFAP_bl = GFAP_Q, pT217_bl = pT217_F, BMI_bl = BMI_tv, eGFR_bl = eGFR_bl)

t1_dat <- attr_tbl %>% left_join(bl_biom, by = "RID")

table1 <- tibble::tribble(
  ~Characteristic, ~Value,
  "N", as.character(nrow(t1_dat)),
  "Age, mean (SD)", sprintf("%.1f (%.1f)", mean(t1_dat$age_bl), sd(t1_dat$age_bl)),
  "Female, n (%)", sprintf("%d (%.1f%%)", sum(t1_dat$sex=="Female"), 100*mean(t1_dat$sex=="Female")),
  "Education, mean (SD)", sprintf("%.1f (%.1f)", mean(t1_dat$edu, na.rm=TRUE), sd(t1_dat$edu, na.rm=TRUE)),
  "Diagnosis CN, n (%)", sprintf("%d (%.1f%%)", sum(t1_dat$dx3=="CN"), 100*mean(t1_dat$dx3=="CN")),
  "Diagnosis MCI, n (%)", sprintf("%d (%.1f%%)", sum(t1_dat$dx3=="MCI"), 100*mean(t1_dat$dx3=="MCI")),
  "Diagnosis AD, n (%)", sprintf("%d (%.1f%%)", sum(t1_dat$dx3=="AD"), 100*mean(t1_dat$dx3=="AD")),
  "APOE4+, n (%)", sprintf("%d (%.1f%%)", sum(t1_dat$APOE4==1), 100*mean(t1_dat$APOE4==1)),
  "HTN, n (%)", sprintf("%d (%.1f%%)", sum(t1_dat$HTN_base==1), 100*mean(t1_dat$HTN_base==1)),
  "DM, n (%)", sprintf("%d (%.1f%%)", sum(t1_dat$DM_base==1), 100*mean(t1_dat$DM_base==1)),
  "Stroke ever, n (%)", sprintf("%d (%.1f%%)", sum(t1_dat$stroke_ever==1), 100*mean(t1_dat$stroke_ever==1)),
  "Baseline BMI, mean (SD)", sprintf("%.1f (%.1f)", mean(t1_dat$BMI_bl, na.rm=TRUE), sd(t1_dat$BMI_bl, na.rm=TRUE)),
  "Baseline eGFR, mean (SD)", sprintf("%.1f (%.1f)", mean(t1_dat$eGFR_bl, na.rm=TRUE), sd(t1_dat$eGFR_bl, na.rm=TRUE)),
  "Baseline GFAP pg/mL, median (IQR)", sprintf("%.1f (%.1f-%.1f)", median(t1_dat$GFAP_bl, na.rm=TRUE),
            quantile(t1_dat$GFAP_bl,.25,na.rm=TRUE), quantile(t1_dat$GFAP_bl,.75,na.rm=TRUE)),
  "Baseline p-tau217 pg/mL, median (IQR)", sprintf("%.3f (%.3f-%.3f)", median(t1_dat$pT217_bl, na.rm=TRUE),
            quantile(t1_dat$pT217_bl,.25,na.rm=TRUE), quantile(t1_dat$pT217_bl,.75,na.rm=TRUE)),
  "GFAP measurements, median (IQR)", sprintf("%.0f (%.0f-%.0f)", median(t1_dat$n_gfap),
            quantile(t1_dat$n_gfap,.25), quantile(t1_dat$n_gfap,.75)),
  "p-tau217 measurements, median (IQR)", sprintf("%.0f (%.0f-%.0f)", median(t1_dat$n_ptau),
            quantile(t1_dat$n_ptau,.25), quantile(t1_dat$n_ptau,.75))
)
wc(table1, "table1_baseline.csv")
cat("Saved Table 1\n")

# ------------------------------------------------------------------------
# 8. PHASE I — LCGMM (exploratory only)
# ------------------------------------------------------------------------
cat("\n=== PHASE I: LCGMM (exploratory) ===\n")
lcgmm_results <- list()
for (bm in c("log_gfap","log_ptau")) {
  dat_lc <- biom_long %>% filter(!is.na(.data[[bm]])) %>% as.data.frame()
  res <- tibble::tibble(biomarker=character(), k=integer(),
                        loglik=double(), AIC=double(), BIC=double(),
                        entropy=double(), min_class_pct=double())
  m0 <- tryCatch(
    hlme(as.formula(paste(bm, "~ time_yrs")), random = ~time_yrs,
         subject = "RID", data = dat_lc, ng = 1),
    error = function(e) NULL)
  if (!is.null(m0)) {
    res <- res %>% add_row(biomarker = bm, k = 1L,
                           loglik = as.numeric(m0$loglik),
                           AIC = m0$AIC, BIC = m0$BIC,
                           entropy = 1, min_class_pct = 100)
  }
  for (k in 2:4) {
    mk <- tryCatch(
      hlme(as.formula(paste(bm, "~ time_yrs")),
           mixture = ~time_yrs, random = ~time_yrs,
           subject = "RID", data = dat_lc, ng = k,
           B = m0),
      error = function(e) NULL)
    if (!is.null(mk)) {
      pp <- mk$pprob
      tab_pct <- 100*as.numeric(table(pp$class)) / nrow(pp)
      res <- res %>% add_row(biomarker = bm, k = k,
                             loglik = as.numeric(mk$loglik),
                             AIC = mk$AIC, BIC = mk$BIC,
                             entropy = NA_real_,
                             min_class_pct = min(tab_pct))
    }
  }
  lcgmm_results[[bm]] <- res
}
table2 <- bind_rows(lcgmm_results)
wc(table2, "table2_lcgmm.csv")
cat("Saved Table 2 (LCGMM exploratory)\n")
print(table2)

# ------------------------------------------------------------------------
# 9. PHASE II — Primary LME with Tier 1/2 covariates
# ------------------------------------------------------------------------
cat("\n=== PHASE II: LME (primary) ===\n")
# Standardise covariates within analytic cohort
zsc <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
biom_long <- biom_long %>% mutate(
  age_z   = zsc(age_bl),
  edu_z   = zsc(edu),
  BMI_z   = zsc(BMI_tv),
  eGFR_z  = zsc(eGFR_bl),
  dx3     = factor(dx3, levels = c("CN","MCI","AD"))
)

# Two biomarkers x two outcomes-of-interest at this phase
phase2_fit <- function(yvar, include_egfr = FALSE) {
  rhs <- "time_yrs*APOE4 + age_z + sex + edu_z + dx3 + BMI_z + HTN_base + DM_base"
  if (include_egfr) rhs <- paste(rhs, "+ eGFR_z")
  fml <- as.formula(paste0(yvar, " ~ ", rhs, " + (1 + time_yrs | RID)"))
  d <- biom_long %>% filter(!is.na(.data[[yvar]]), !is.na(BMI_z))
  if (include_egfr) d <- d %>% filter(!is.na(eGFR_z))
  m <- lmerTest::lmer(fml, data = d, REML = TRUE,
                      control = lmerControl(optimizer="bobyqa"))
  list(model = m, n_obs = nrow(d), n_id = length(unique(d$RID)))
}
phase2_extract <- function(fit) {
  cf <- as.data.frame(summary(fit$model)$coefficients)
  cf$term <- rownames(cf)
  ci <- tryCatch(confint(fit$model, parm = "beta_", method = "Wald"),
                 error = function(e) NULL)
  if (!is.null(ci)) {
    ci <- as.data.frame(ci); ci$term <- rownames(ci)
    cf <- cf %>% left_join(ci, by = "term")
  }
  cf$n_obs <- fit$n_obs; cf$n_id <- fit$n_id
  cf$singular <- isSingular(fit$model)
  cf
}

phase2_gfap <- phase2_fit("log_gfap", include_egfr = FALSE)
phase2_ptau <- phase2_fit("log_ptau", include_egfr = FALSE)
table3 <- bind_rows(
  phase2_extract(phase2_gfap) %>% mutate(biomarker = "log_gfap", model = "primary_noEGFR"),
  phase2_extract(phase2_ptau) %>% mutate(biomarker = "log_ptau", model = "primary_noEGFR")
)
# eGFR sensitivity in subset that has eGFR
phase2_gfap_e <- tryCatch(phase2_fit("log_gfap", include_egfr = TRUE), error = function(e) NULL)
phase2_ptau_e <- tryCatch(phase2_fit("log_ptau", include_egfr = TRUE), error = function(e) NULL)
if (!is.null(phase2_gfap_e))
  table3 <- bind_rows(table3, phase2_extract(phase2_gfap_e) %>%
                        mutate(biomarker = "log_gfap", model = "sensitivity_withEGFR"))
if (!is.null(phase2_ptau_e))
  table3 <- bind_rows(table3, phase2_extract(phase2_ptau_e) %>%
                        mutate(biomarker = "log_ptau", model = "sensitivity_withEGFR"))
wc(table3, "table3_lme_phase2.csv")
cat("Saved Table 3 (Phase II LME, primary + eGFR sensitivity)\n")
key_rows <- table3[grepl("time_yrs:APOE4|^APOE4|^time_yrs$", table3$term), ]
print(key_rows[, c("biomarker","model","term","Estimate","Std. Error","Pr(>|t|)","singular","n_obs","n_id")])

# Extract per-participant slope BLUPs for Phase III sensitivity
slope_blup <- function(fit, bm_name) {
  re <- ranef(fit$model)$RID
  fe_time <- fixef(fit$model)["time_yrs"]
  data.frame(RID = as.integer(rownames(re)),
             slope = re[, "time_yrs"] + fe_time)
}
blup_gfap <- slope_blup(phase2_gfap, "GFAP")
blup_ptau <- slope_blup(phase2_ptau, "pT217")
names(blup_gfap)[2] <- "slope_gfap"
names(blup_ptau)[2] <- "slope_ptau"
blups <- full_join(blup_gfap, blup_ptau, by = "RID")
wc(blups, "blup_slopes.csv")

# ------------------------------------------------------------------------
# 10. OUTCOME long tables
# ------------------------------------------------------------------------
cat("\n=== Building outcome long tables ===\n")
# ADAS-Cog 13
adas_long <- adas %>% transmute(
  RID = as.integer(RID),
  EXAMDATE = as.Date(VISDATE),
  adas13 = suppressWarnings(as.numeric(TOTAL13))
) %>% filter(!is.na(EXAMDATE), !is.na(adas13), RID %in% analytic_rid) %>%
  left_join(bl_dates, by = "RID") %>%
  mutate(time_yrs = as.numeric(EXAMDATE - bl_date)/365.25) %>%
  filter(!is.na(time_yrs)) %>% left_join(attr_tbl, by = "RID")

mmse_long <- mmse %>% transmute(
  RID = as.integer(RID),
  EXAMDATE = as.Date(VISDATE),
  mmse = suppressWarnings(as.numeric(MMSCORE))
) %>% filter(!is.na(EXAMDATE), !is.na(mmse), RID %in% analytic_rid) %>%
  left_join(bl_dates, by = "RID") %>%
  mutate(time_yrs = as.numeric(EXAMDATE - bl_date)/365.25) %>%
  filter(!is.na(time_yrs)) %>% left_join(attr_tbl, by = "RID")

fs_long <- fs %>% transmute(
  RID = as.integer(RID),
  EXAMDATE = as.Date(PHC_SCANDATE),
  hippo = Left.Hippocampus_combat + Right.Hippocampus_combat,
  ICV   = EstimatedTotalIntraCranialVol_combat,
  hippo_ratio = hippo/ICV
) %>% filter(!is.na(hippo_ratio), RID %in% analytic_rid) %>%
  left_join(bl_dates, by = "RID") %>%
  mutate(time_yrs = as.numeric(EXAMDATE - bl_date)/365.25) %>%
  filter(!is.na(time_yrs)) %>% left_join(attr_tbl, by = "RID")

cat("ADAS rows:", nrow(adas_long), " RIDs:", length(unique(adas_long$RID)), "\n")
cat("MMSE rows:", nrow(mmse_long), " RIDs:", length(unique(mmse_long$RID)), "\n")
cat("Hippo rows:", nrow(fs_long),  " RIDs:", length(unique(fs_long$RID)),  "\n")

# z-standardise covars in outcome tables
for (d_name in c("adas_long","mmse_long","fs_long")) {
  d <- get(d_name)
  d$age_z <- zsc(d$age_bl); d$edu_z <- zsc(d$edu)
  d$dx3   <- factor(d$dx3, levels = c("CN","MCI","AD"))
  assign(d_name, d)
}

# ------------------------------------------------------------------------
# 11. PHASE III — Two-stage BLUP-as-predictor (sensitivity S-G); JMbayes2 follows
# ------------------------------------------------------------------------
cat("\n=== PHASE III (two-stage sensitivity first) ===\n")
two_stage <- function(outc_df, ycol, slope_col) {
  d <- outc_df %>% left_join(blups, by = "RID") %>%
    filter(!is.na(.data[[slope_col]]))
  d$slope_z <- zsc(d[[slope_col]])
  fml <- as.formula(paste0(
    ycol, " ~ time_yrs*slope_z + APOE4 + age_z + sex + edu_z + dx3 + (1 + time_yrs | RID)"))
  m <- tryCatch(lmerTest::lmer(fml, data = d, REML = TRUE,
                               control = lmerControl(optimizer = "bobyqa")),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  cf <- summary(m)$coefficients
  ix <- grep("time_yrs:slope_z", rownames(cf))
  if (length(ix) == 0) return(NULL)
  est <- cf[ix, "Estimate"]; se <- cf[ix, "Std. Error"]
  data.frame(outcome = ycol, predictor = slope_col,
             beta_per_SD = est, SE = se,
             CI_lo = est - 1.96*se, CI_hi = est + 1.96*se,
             p_value = cf[ix, "Pr(>|t|)"],
             n_obs = nrow(d), n_id = length(unique(d$RID)))
}

table4_ts <- bind_rows(
  two_stage(adas_long, "adas13", "slope_gfap"),
  two_stage(adas_long, "adas13", "slope_ptau"),
  two_stage(mmse_long, "mmse",   "slope_gfap"),
  two_stage(mmse_long, "mmse",   "slope_ptau"),
  two_stage(fs_long,   "hippo_ratio", "slope_gfap"),
  two_stage(fs_long,   "hippo_ratio", "slope_ptau")
)
table4_ts$q_BH <- p.adjust(table4_ts$p_value, method = "BH")
table4_ts$q_Bonf <- pmin(table4_ts$p_value * 6, 1)
wc(table4_ts, "table4_twostage_sensitivity.csv")
cat("Saved Table 4 (two-stage sensitivity)\n")
print(table4_ts)

# ------------------------------------------------------------------------
# 12. PHASE III — JMbayes2 joint models (primary)
# ------------------------------------------------------------------------
cat("\n=== PHASE III: JMbayes2 joint models (primary) ===\n")
fit_jm <- function(outc_df, ycol, bm_long_df, bm_col) {
  # Biomarker submodel data
  bd <- bm_long_df %>% filter(!is.na(.data[[bm_col]]),
                              !is.na(BMI_z)) %>%
    as.data.frame()
  bd$bm <- bd[[bm_col]]
  # Outcome submodel data
  od <- outc_df %>% as.data.frame()
  od$y <- od[[ycol]]
  od <- od[!is.na(od$y), ]
  # Restrict both to intersection of RIDs and require >=2 obs/RID in each
  rid_bd <- as.integer(names(which(table(bd$RID) >= 2)))
  rid_od <- as.integer(names(which(table(od$RID) >= 2)))
  rid_keep <- intersect(rid_bd, rid_od)
  bd <- bd[bd$RID %in% rid_keep, ]
  od <- od[od$RID %in% rid_keep, ]
  cat("    JM data: bm n_obs=", nrow(bd), " n_id=", length(unique(bd$RID)),
      "; outc n_obs=", nrow(od), " n_id=", length(unique(od$RID)),
      " (shared RIDs=", length(rid_keep), ")\n", sep = "")
  if (length(rid_keep) < 30) return(NULL)
  bm_lme <- tryCatch(
    nlme::lme(bm ~ time_yrs + APOE4 + age_z + sex + dx3 + BMI_z,
              random = ~ time_yrs | RID, data = bd,
              control = nlme::lmeControl(opt = "optim", returnObject = TRUE)),
    error = function(e) { message("bm_lme err: ", e$message); NULL })
  if (is.null(bm_lme)) return(NULL)
  outc_lme <- tryCatch(
    nlme::lme(y ~ time_yrs + APOE4 + age_z + sex + edu_z + dx3,
              random = ~ time_yrs | RID, data = od,
              control = nlme::lmeControl(opt = "optim", returnObject = TRUE)),
    error = function(e) { message("outc_lme err: ", e$message); NULL })
  if (is.null(outc_lme)) return(NULL)
  jm_fit <- tryCatch(
    JMbayes2::jm(outc_lme, bm_lme, time_var = "time_yrs",
                 functional_forms = list("bm" = ~ slope(bm)),
                 n_iter = 4000L, n_burnin = 1500L, n_chains = 2L,
                 cores = 2L),
    error = function(e) { message("jm err (slope): ", e$message); NULL })
  if (is.null(jm_fit)) {
    # fallback: default association (current value of bm)
    jm_fit <- tryCatch(
      JMbayes2::jm(outc_lme, bm_lme, time_var = "time_yrs",
                   n_iter = 4000L, n_burnin = 1500L, n_chains = 2L,
                   cores = 2L),
      error = function(e) { message("jm err (default): ", e$message); NULL })
  }
  jm_fit
}

extract_jm <- function(jm_fit, outcome, predictor) {
  if (is.null(jm_fit)) {
    return(data.frame(outcome=outcome, predictor=predictor,
                      mean=NA, lo=NA, hi=NA, P=NA, rhat=NA, note="failed"))
  }
  s <- summary(jm_fit)
  # Pull alpha (association parameter) coefficients
  al <- s$Outcome[grepl("slope\\(bm\\)", rownames(s$Outcome)), , drop = FALSE]
  if (nrow(al) == 0) {
    return(data.frame(outcome=outcome, predictor=predictor,
                      mean=NA, lo=NA, hi=NA, P=NA, rhat=NA, note="no_alpha"))
  }
  data.frame(outcome=outcome, predictor=predictor,
             mean = al[1,"Mean"],
             lo   = al[1,"2.5%"],
             hi   = al[1,"97.5%"],
             P    = al[1,"P"],
             rhat = ifelse("Rhat" %in% colnames(al), al[1,"Rhat"], NA),
             note = "ok")
}

jm_results <- list()
pairs <- list(
  list(out=adas_long, yc="adas13",       bmcol="log_gfap", lab="ADAS-Cog13 ~ GFAP"),
  list(out=adas_long, yc="adas13",       bmcol="log_ptau", lab="ADAS-Cog13 ~ pT217"),
  list(out=mmse_long, yc="mmse",         bmcol="log_gfap", lab="MMSE ~ GFAP"),
  list(out=mmse_long, yc="mmse",         bmcol="log_ptau", lab="MMSE ~ pT217"),
  list(out=fs_long,   yc="hippo_ratio",  bmcol="log_gfap", lab="Hippo ~ GFAP"),
  list(out=fs_long,   yc="hippo_ratio",  bmcol="log_ptau", lab="Hippo ~ pT217")
)
for (p in pairs) {
  cat("\n  fitting:", p$lab, "...\n")
  t0 <- Sys.time()
  jm <- fit_jm(p$out, p$yc, biom_long, p$bmcol)
  cat("    elapsed:", round(as.numeric(difftime(Sys.time(), t0, units="mins")),2), "min\n")
  jm_results[[p$lab]] <- extract_jm(jm, p$yc, p$bmcol)
  saveRDS(jm, file.path(OUT, paste0("jm_", gsub("[~ ]","_", p$lab), ".rds")))
}
table4_jm <- bind_rows(jm_results)
if (nrow(table4_jm) > 0 && "P" %in% colnames(table4_jm)) {
  table4_jm$q_BH   <- p.adjust(table4_jm$P, method = "BH")
  table4_jm$q_Bonf <- pmin(table4_jm$P * 6, 1)
}
wc(table4_jm, "table4_jointmodel_primary.csv")
cat("\nSaved Table 4 (joint model PRIMARY)\n")
print(table4_jm)

# ------------------------------------------------------------------------
# 13. PHASE IV — Cross-assay p-tau217 (descriptive)
# ------------------------------------------------------------------------
cat("\n=== PHASE IV: cross-assay p-tau217 ===\n")
ptau_bl <- biom_long %>% group_by(RID) %>%
  filter(time_yrs == min(time_yrs)) %>% ungroup() %>%
  transmute(RID, EXAMDATE, pT217_F = pT217_F)

corx <- function(d, col, label) {
  if (is.null(d) || is.null(col) || !(col %in% names(d)))
    return(data.frame(assay=label, n=0, pearson=NA, spearman=NA))
  d$val <- suppressWarnings(as.numeric(d[[col]]))
  d <- d[!is.na(d$val) & d$val > 0 & as.integer(d$RID) %in% analytic_rid, ]
  if (nrow(d) == 0) return(data.frame(assay=label, n=0, pearson=NA, spearman=NA))
  d$EXAMDATE <- as.Date(d$EXAMDATE)
  d_bl <- d %>% mutate(RID = as.integer(RID)) %>%
    group_by(RID) %>% filter(EXAMDATE == min(EXAMDATE)) %>% ungroup() %>%
    select(RID, val)
  merged <- inner_join(ptau_bl, d_bl, by = "RID")
  if (nrow(merged) < 5)
    return(data.frame(assay=label, n=nrow(merged), pearson=NA, spearman=NA))
  data.frame(assay = label, n = nrow(merged),
             pearson  = cor(merged$pT217_F, merged$val, use="complete.obs"),
             spearman = cor(merged$pT217_F, merged$val, use="complete.obs", method="spearman"))
}

# Lilly: long-format keyed by TESTCD; filter to PTAU217 and use ORRES
lilly_w <- if (!is.null(lilly)) {
  lilly %>% filter(grepl("PTAU217", TESTCD, ignore.case = TRUE) |
                   grepl("p-?tau-?217", TEST, ignore.case = TRUE)) %>%
    transmute(RID = as.integer(RID), EXAMDATE = as.Date(EXAMDATE),
              lilly_val = suppressWarnings(as.numeric(ORRES)))
} else NULL

xa_results <- bind_rows(
  corx(lilly_w, "lilly_val",  "Lilly MSD600"),
  corx(c2n,     "pT217_C2N",  "C2N PrecivityAD2"),
  corx(jans,    "DILUTION_CORRECTED_CONC", "Janssen")
)
wc(xa_results, "table_cross_assay.csv")
cat("Cross-assay results:\n"); print(xa_results)

# ------------------------------------------------------------------------
# 14. Save analytic long dataset for downstream figures
# ------------------------------------------------------------------------
wc(biom_long, "analytic_biom_long.csv")
wc(attr_tbl,  "analytic_baseline_attr.csv")

cat("\n==========================================\n",
    " PIPELINE COMPLETE  ", as.character(Sys.time()), "\n",
    "==========================================\n", sep = "")
