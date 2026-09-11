# ADNI Plasma GFAP / p-tau217 — Joint longitudinal model pipeline (v1)
# Author: Kuan-Yu Chu | 2026-05-25
# SAP: SAP_v1.md (locked 2026-02-12; amended 2026-05-25)
# Inclusion: age >=65, dx in {CN,MCI,AD}, APOE available, >=2 GFAP AND >=2 pT217 (revised)
# Primary framework: LME (Phase II), JMbayes2 joint model (Phase III)
# LCGMM: exploratory diagnostic only (Phase I)

# --------- 0. Environment ---------
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(ggplot2)
  library(lme4); library(nlme); library(lcmm)
  # library(JMbayes2)   # uncomment after install
  # library(lmerTest)   # uncomment after install
  # library(brms)       # uncomment after install
})
set.seed(20260525)

DATA <- "/Users/kuan-yuchu/Desktop/LOAD/ADNI_csv"
OUT  <- "/Users/kuan-yuchu/Desktop/LOAD/Analysis/outputs"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

read_csv  <- function(f) utils::read.csv(f, stringsAsFactors = FALSE, na.strings = c("", "NA"))
write_csv <- function(x, f) utils::write.csv(x, f, row.names = FALSE)

# --------- 1. Load data ---------
upenn <- read_csv(file.path(DATA, "UPENN_PLASMA_FUJIREBIO_QUANTERIX_24May2026.csv"))
apoe  <- read_csv(file.path(DATA, "APOERES_24May2026.csv"))
vitals<- read_csv(file.path(DATA, "VITALS_24May2026.csv"))
adas  <- read_csv(file.path(DATA, "ADAS_24May2026.csv"))
mmse  <- read_csv(file.path(DATA, "MMSE_24May2026.csv"))
fs    <- read_csv(file.path(DATA, "ADSP_PHC_T1_FS_24May2026.csv"))
load("/tmp/adnimerge_unpack/ADNIMERGE/data/adnimerge.rdata")

# --------- 2. Build biomarker long table ---------
upenn <- upenn %>%
  mutate(
    RID = as.integer(RID),
    EXAMDATE = as.Date(EXAMDATE),
    pT217_F = suppressWarnings(as.numeric(pT217_F)),
    GFAP_Q  = suppressWarnings(as.numeric(GFAP_Q)),
    pT217_F = ifelse(pT217_F <= 0, NA, pT217_F),
    GFAP_Q  = ifelse(GFAP_Q  <= 0, NA, GFAP_Q),
    log_ptau  = log(pT217_F),
    log_gfap  = log(GFAP_Q)
  )

# Baseline per RID (earliest EXAMDATE)
bl_dates <- upenn %>% group_by(RID) %>%
  summarise(bl_date = min(EXAMDATE[!is.na(pT217_F) | !is.na(GFAP_Q)], na.rm = TRUE),
            .groups = "drop")

upenn <- upenn %>% left_join(bl_dates, by = "RID") %>%
  mutate(time_yrs = as.numeric(EXAMDATE - bl_date) / 365.25)

# Per-RID measurement counts
counts <- upenn %>% group_by(RID) %>%
  summarise(n_gfap = sum(!is.na(GFAP_Q)),
            n_ptau = sum(!is.na(pT217_F)), .groups = "drop")

# --------- 3. APOE + baseline demographics ---------
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

# --------- 4. Time-varying BMI ---------
vitals <- vitals %>% mutate(
  RID = as.integer(RID),
  EXAMDATE = as.Date(EXAMDATE),
  VSWEIGHT = suppressWarnings(as.numeric(VSWEIGHT)),
  VSHEIGHT = suppressWarnings(as.numeric(VSHEIGHT)),
  VSWTUNIT = as.integer(VSWTUNIT),  # 1=lbs, 2=kg in ADNI
  VSHTUNIT = as.integer(VSHTUNIT)   # 1=in, 2=cm
)
vitals <- vitals %>%
  mutate(
    wt_kg = ifelse(!is.na(VSWTUNIT) & VSWTUNIT == 1, VSWEIGHT * 0.4536, VSWEIGHT),
    ht_m  = ifelse(!is.na(VSHTUNIT) & VSHTUNIT == 1, VSHEIGHT * 0.0254, VSHEIGHT / 100),
    BMI   = wt_kg / (ht_m^2)
  ) %>%
  filter(BMI > 10, BMI < 70) %>%   # sanity
  select(RID, EXAMDATE, BMI)

# LOCF within 12 months: for each biomarker visit, take most recent BMI within +/-365 days
attach_bmi <- function(biom_long, vitals_long) {
  out <- biom_long
  out$BMI_tv <- NA_real_
  for (i in seq_len(nrow(out))) {
    cand <- vitals_long[vitals_long$RID == out$RID[i] &
                          abs(as.numeric(vitals_long$EXAMDATE - out$EXAMDATE[i])) <= 365, ]
    if (nrow(cand) > 0) {
      diffs <- abs(as.numeric(cand$EXAMDATE - out$EXAMDATE[i]))
      out$BMI_tv[i] <- cand$BMI[which.min(diffs)]
    }
  }
  out
}
# NOTE: this loop is O(N*V); for N=325 with ~3 visits = ~1000 iters * vitals ~17k rows -> ~17M comparisons.
# Acceptable on M-series Mac (~15s). If slow, switch to data.table rolling join.

# --------- 5. eGFR (PENDING LAB DATA DOWNLOAD) ---------
# Action: download "Laboratory Data — Clinical Chemistry" from ida.loni.usc.edu
# Expected columns: RID, EXAMDATE, RCT11 or CREATININE (serum creatinine mg/dL)
# CKD-EPI 2021 race-free formula:
#   eGFR = 142 * min(Scr/k,1)^a * max(Scr/k,1)^-1.200 * 0.9938^age * (1.012 if female)
#   k = 0.7 (F), 0.9 (M); a = -0.241 (F), -0.302 (M)
# If lab data unavailable, set eGFR_tv = NA and omit from primary models per SAP §10.4.
egfr <- NULL  # placeholder

# --------- 6. Stroke / HTN / DM (PENDING MEDICAL HISTORY UNPACK) ---------
# Source: Medical_History.zip in ADNI_csv/. Unpack and parse MEDHIST table
# Stroke: MHCAT == "Vascular" with MHDESC matching /stroke|CVA|cerebrovascular/i
# HTN/DM: similar matching
# Placeholder:
htn_dm_stroke <- data.frame(RID = integer(), HTN_base = integer(),
                            DM_base = integer(), stroke_ever = integer())

# --------- 7. Hippocampal volume (ratio) ---------
fs <- fs %>% mutate(
  RID = as.integer(RID),
  EXAMDATE = as.Date(EXAMDATE),
  hippo = Left.Hippocampus_combat + Right.Hippocampus_combat,
  ICV   = EstimatedTotalIntraCranialVol_combat,
  hippo_ratio = hippo / ICV
) %>% filter(!is.na(hippo_ratio))

# --------- 8. ADAS-Cog 13 / MMSE long ---------
# Inspect headers separately; placeholder column names — verify:
# ADAS: TOTAL13 or TOTSCORE; MMSE: MMSCORE
# (Add after running head() on ADAS, MMSE CSVs.)

# --------- 9. Apply inclusion (SAP §4 revised) ---------
analytic_rid <- counts %>%
  filter(n_gfap >= 2, n_ptau >= 2) %>%
  inner_join(apoe, by = "RID") %>%
  inner_join(adni_bl, by = "RID") %>%
  filter(age_bl >= 65, !is.na(dx3)) %>%
  pull(RID)

cat("Analytic N =", length(analytic_rid), "(expected ~325)\n")

# Build long dataset
biom_long <- upenn %>%
  filter(RID %in% analytic_rid) %>%
  select(RID, EXAMDATE, time_yrs, log_gfap, log_ptau)

# --------- 10. Phase I — LCGMM (EXPLORATORY ONLY) ---------
# (placeholder — see SAP §6.2 for criteria)

# --------- 11. Phase II — primary LME (per SAP §6.3) ---------
# Once BMI_tv + HTN/DM + (eGFR_tv) merged, fit:
#   m_gfap <- lmer(log_gfap ~ time_yrs*APOE4 + age_z + sex + edu_z + dx3
#                  + BMI_tv_z + HTN_base + DM_base + (1 + time_yrs | RID),
#                  data = dat, REML = TRUE)
#   confint(m_gfap, method = "profile")
#   isSingular(m_gfap)

# --------- 12. Phase III — JMbayes2 joint model (per SAP §6.4) ---------
# requires JMbayes2 install:
#   install.packages("JMbayes2")
# Skeleton:
#   library(JMbayes2)
#   biom_lme  <- lme(log_ptau ~ time_yrs*APOE4 + age_z + sex + dx3 + BMI_tv_z,
#                    random = ~ time_yrs | RID, data = dat, control = lmeControl(opt = "optim"))
#   outc_lme  <- lme(adas13 ~ time_yrs + APOE4 + age_z + sex + dx3,
#                    random = ~ time_yrs | RID, data = outc_dat)
#   jm_fit    <- jm(outc_lme, biom_lme, time_var = "time_yrs",
#                   functional_forms = ~ slope(log_ptau):time_yrs,
#                   n_iter = 6000, n_burnin = 2000, n_chains = 4)
#   summary(jm_fit)

cat("Skeleton complete. Required next steps:\n")
cat("  - install JMbayes2, brms, lmerTest, data.table\n")
cat("  - download ADNI Lab Chemistry (creatinine) for eGFR\n")
cat("  - unpack Medical_History.zip for HTN/DM/stroke\n")
cat("  - confirm ADAS/MMSE outcome column names\n")
