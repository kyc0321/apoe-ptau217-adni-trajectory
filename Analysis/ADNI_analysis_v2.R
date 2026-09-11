# =============================================================================
# ADNI_analysis_v2.R
# Option C: dual-biomarker blood-based trajectory pipeline
#   exposure: plasma GFAP  (astrocytic reactivity)
#   exposure: plasma p-tau (preferred p-tau217; falls back to p-tau181)
#   genetic modifier: APOE e4
#   outcomes: ADAS-Cog 13, MMSE, hippocampal volume
# Real-data ingestion (replaces the synthetic block in v1).
# =============================================================================

suppressPackageStartupMessages({
  for (pkg in c("lcmm","lme4","nnet","dplyr","tidyr","ggplot2"))
    if (!requireNamespace(pkg, quietly=TRUE))
      install.packages(pkg, repos="https://cloud.r-project.org")
  library(lcmm); library(lme4); library(nnet)
  library(dplyr); library(tidyr); library(ggplot2)
})
set.seed(20260523)

# -----------------------------------------------------------------------------
# 1. CONFIGURATION -- ** EDIT FILE NAMES TO MATCH YOUR LONI DOWNLOADS **
# -----------------------------------------------------------------------------
CSV_DIR <- "~/Desktop/LOAD/ADNI_csv"

FILES <- list(
  gfap      = "PLASMA_GFAP.csv",       # plasma GFAP (Quanterix Simoa)
  ptau217   = "PLASMA_PTAU217.csv",    # plasma p-tau217 (Lilly MSD) -- preferred
  ptau181   = "PLASMA_PTAU181.csv",    # fallback if p-tau217 not released for cohort
  apoe      = "APOERES.csv",
  demog     = "PTDEMOG.csv",
  diag      = "DXSUM_PDXCONV_ADNIALL.csv",
  vitals    = "VITALS.csv",
  medhist   = "MEDHIST.csv",
  adas      = "ADAS_ADNIGO23.csv",
  mmse      = "MMSE.csv",
  freesurf  = "UCSFFSX51_ADNI.csv",
  registry  = "REGISTRY.csv"
)

PTAU_PREFERENCE <- c("ptau217","ptau181")   # use first one that's available

# -----------------------------------------------------------------------------
# 2. HELPERS
# -----------------------------------------------------------------------------
read_adni <- function(name) {
  path <- file.path(CSV_DIR, FILES[[name]])
  if (!file.exists(path)) { message("  [MISSING] ", path); return(NULL) }
  read.csv(path, na.strings = c("", "NA", "-1", "-4", "-9", "NaN"),
           stringsAsFactors = FALSE, check.names = FALSE)
}

# Visit code -> months since baseline
viscode_to_month <- function(v) {
  v <- tolower(trimws(as.character(v)))
  out <- suppressWarnings(as.integer(NA))
  out <- ifelse(v %in% c("bl","baseline","scrn","sc","sccn"), 0L, NA_integer_)
  m <- regmatches(v, regexec("^m(\\d+)$", v))
  for (i in seq_along(v)) {
    if (is.na(out[i]) && length(m[[i]]) >= 2)
      out[i] <- as.integer(m[[i]][2])
  }
  out
}

pick_col <- function(df, candidates) {
  hits <- candidates[candidates %in% names(df)]
  if (length(hits) == 0) return(NULL)
  hits[1]
}

# -----------------------------------------------------------------------------
# 3. INGEST EACH TABLE -> CANONICAL FRAME (RID, time, var)
# -----------------------------------------------------------------------------
cat("=== Reading ADNI CSVs ===\n")

## 3a. Plasma GFAP
gfap_raw <- read_adni("gfap")
stopifnot(!is.null(gfap_raw))
gfap_visit <- pick_col(gfap_raw, c("VISCODE2","VISCODE","VISIT"))
gfap_val   <- pick_col(gfap_raw, c("PLASMA_GFAP","GFAP","GFAP_PLASMA","PLASMAGFAP"))
stopifnot(!is.null(gfap_visit), !is.null(gfap_val))
gfap <- data.frame(
  RID  = gfap_raw$RID,
  time = viscode_to_month(gfap_raw[[gfap_visit]]),
  GFAP = suppressWarnings(as.numeric(gfap_raw[[gfap_val]]))
) %>% filter(!is.na(GFAP), !is.na(time))

## 3b. Plasma p-tau (try p-tau217 first, fall back to p-tau181)
ptau_used <- NULL
ptau <- NULL
for (which in PTAU_PREFERENCE) {
  raw <- read_adni(which)
  if (is.null(raw)) next
  pv <- pick_col(raw, c("VISCODE2","VISCODE","VISIT"))
  cand <- if (which == "ptau217")
    c("PLASMA_PTAU217","PTAU217","PTAU_217","P_TAU217","PTau217")
  else
    c("PLASMA_PTAU181","PTAU181","PTAU_181","P_TAU181","PTau181")
  pcol <- pick_col(raw, cand)
  if (is.null(pv) || is.null(pcol)) {
    message("  [", which, "] visit or value column not found; skipping")
    next
  }
  ptau <- data.frame(
    RID  = raw$RID,
    time = viscode_to_month(raw[[pv]]),
    PTAU = suppressWarnings(as.numeric(raw[[pcol]]))
  ) %>% filter(!is.na(PTAU), !is.na(time))
  ptau_used <- which
  cat("  using plasma", ptau_used, "(n_obs =", nrow(ptau), ")\n")
  break
}
stopifnot(!is.null(ptau))

## 3c. APOE
apoe_raw <- read_adni("apoe")
ap1 <- pick_col(apoe_raw, c("APGEN1","APOE1","APOEGEN1"))
ap2 <- pick_col(apoe_raw, c("APGEN2","APOE2","APOEGEN2"))
apoe <- data.frame(
  RID = apoe_raw$RID,
  apoe_e4 = as.integer(apoe_raw[[ap1]] == 4 | apoe_raw[[ap2]] == 4)
)

## 3d. Demographics (sex, education, DOB)
dem_raw <- read_adni("demog")
sex_col <- pick_col(dem_raw, c("PTGENDER","SEX"))
edu_col <- pick_col(dem_raw, c("PTEDUCAT","EDUC","EDUCATION"))
dob_col <- pick_col(dem_raw, c("PTDOB","DOB"))
demog <- data.frame(
  RID  = dem_raw$RID,
  female = as.integer(dem_raw[[sex_col]] %in% c(2, "2", "F", "Female", "female")),
  educ   = suppressWarnings(as.numeric(dem_raw[[edu_col]])),
  DOB    = dem_raw[[dob_col]]   # may need parsing later
)

## 3e. Clinical diagnosis (longitudinal -> we take baseline)
dx_raw <- read_adni("diag")
dx_vis <- pick_col(dx_raw, c("VISCODE2","VISCODE","VISIT"))
dx_val <- pick_col(dx_raw, c("DIAGNOSIS","DXCURREN","DXCHANGE","DX_bl"))
dx_long <- data.frame(
  RID = dx_raw$RID,
  time = viscode_to_month(dx_raw[[dx_vis]]),
  dx_code = dx_raw[[dx_val]]
) %>% filter(!is.na(time))
# Map: 1=CN, 2=MCI, 3=AD (ADNI convention varies; treat 2 and 3 as MCI/AD)
dx_base <- dx_long %>% filter(time == 0) %>%
  mutate(diag = case_when(
    dx_code %in% c(2,"2","MCI","LMCI","EMCI") ~ "MCI",
    dx_code %in% c(3,"3","AD","Dementia","LOAD") ~ "AD",
    TRUE ~ NA_character_)) %>%
  select(RID, diag)

## 3f. ADAS-Cog 13
adas_raw <- read_adni("adas")
av <- pick_col(adas_raw, c("VISCODE2","VISCODE","VISIT"))
a13 <- pick_col(adas_raw, c("TOTAL13","TOTSCORE","ADAS13","ADASTOTAL13"))
adas <- data.frame(
  RID = adas_raw$RID,
  time = viscode_to_month(adas_raw[[av]]),
  ADAS13 = suppressWarnings(as.numeric(adas_raw[[a13]]))
) %>% filter(!is.na(time))

## 3g. MMSE
mm_raw <- read_adni("mmse")
mv <- pick_col(mm_raw, c("VISCODE2","VISCODE","VISIT"))
ms <- pick_col(mm_raw, c("MMSCORE","MMSE","MMSE_TOTAL"))
mmse <- data.frame(
  RID = mm_raw$RID,
  time = viscode_to_month(mm_raw[[mv]]),
  MMSE = suppressWarnings(as.numeric(mm_raw[[ms]]))
) %>% filter(!is.na(time))

## 3h. Vitals (BMI)
vit_raw <- read_adni("vitals")
hv_col <- pick_col(vit_raw, c("VSHEIGHT","HEIGHT"))
wv_col <- pick_col(vit_raw, c("VSWEIGHT","WEIGHT"))
hu_col <- pick_col(vit_raw, c("VSHTUNIT"))         # 1=inches/2=cm in ADNI
wu_col <- pick_col(vit_raw, c("VSWTUNIT"))         # 1=lbs/2=kg
vit_base <- vit_raw %>%
  filter(viscode_to_month(if (!is.null(pick_col(vit_raw,c("VISCODE2","VISCODE"))))
                          vit_raw[[pick_col(vit_raw,c("VISCODE2","VISCODE"))]] else "bl") == 0)
calc_bmi <- function(h, w, hu, wu) {
  h_m <- ifelse(!is.na(hu) & hu == 1, h*0.0254, h/100)   # inches->m or cm->m
  w_kg <- ifelse(!is.na(wu) & wu == 1, w*0.4536, w)       # lbs->kg
  w_kg / h_m^2
}
vit <- data.frame(
  RID = vit_base$RID,
  bmi = with(vit_base,
             calc_bmi(suppressWarnings(as.numeric(vit_base[[hv_col]])),
                      suppressWarnings(as.numeric(vit_base[[wv_col]])),
                      if (!is.null(hu_col)) vit_base[[hu_col]] else NA,
                      if (!is.null(wu_col)) vit_base[[wu_col]] else NA))
) %>% group_by(RID) %>% summarise(bmi = median(bmi, na.rm=TRUE), .groups="drop")

## 3i. Medical history (HTN, DM) -- optional, may be free-text
mh_raw <- read_adni("medhist")
htn <- dm <- data.frame(RID = integer(0), v = integer(0))
if (!is.null(mh_raw)) {
  mh_term <- pick_col(mh_raw, c("MHTERM","MHCM","MHTERMC","MED_HISTORY"))
  if (!is.null(mh_term)) {
    mh_raw$txt <- tolower(as.character(mh_raw[[mh_term]]))
    htn <- mh_raw %>% group_by(RID) %>%
      summarise(htn = as.integer(any(grepl("hypertens|high blood pressure|htn", txt))), .groups="drop")
    dm <- mh_raw %>% group_by(RID) %>%
      summarise(dm  = as.integer(any(grepl("diabet|dm type|t2dm", txt))), .groups="drop")
  }
}

## 3j. FreeSurfer hippocampal volume + ICV
fs_raw <- read_adni("freesurf")
fv <- pick_col(fs_raw, c("VISCODE2","VISCODE","VISIT"))
# ADNI UCSFFSX/FSL hippocampal subfield columns: ST29SV (left) + ST88SV (right) total hippocampal volume; ICV: ST10CV
lh <- pick_col(fs_raw, c("ST29SV","LEFT_HIPPOCAMPUS","HIPPO_L"))
rh <- pick_col(fs_raw, c("ST88SV","RIGHT_HIPPOCAMPUS","HIPPO_R"))
icv<- pick_col(fs_raw, c("ST10CV","ICV","INTRACRANIALVOLUME"))
fs <- data.frame(
  RID = fs_raw$RID,
  time = viscode_to_month(fs_raw[[fv]]),
  HV  = suppressWarnings(as.numeric(fs_raw[[lh]])) +
        suppressWarnings(as.numeric(fs_raw[[rh]])),
  ICV = suppressWarnings(as.numeric(fs_raw[[icv]]))
) %>% filter(!is.na(time))
fs$HV_adj <- fs$HV / fs$ICV * mean(fs$ICV, na.rm=TRUE)  # ICV-normalized

## 3k. Registry (for baseline age)
reg_raw <- read_adni("registry")
reg <- NULL
if (!is.null(reg_raw)) {
  rv <- pick_col(reg_raw, c("VISCODE2","VISCODE"))
  rd <- pick_col(reg_raw, c("EXAMDATE","USERDATE"))
  reg <- data.frame(
    RID = reg_raw$RID,
    time = viscode_to_month(reg_raw[[rv]]),
    examdate = reg_raw[[rd]]
  ) %>% filter(!is.na(time))
}

# Baseline age = examdate(bl) - DOB.  DOB in ADNI PTDEMOG is mm/yyyy (no day);
# we assign day 15 for parsing.
parse_dob <- function(x) {
  s <- as.character(x); s[s==""] <- NA
  # try "mm/yyyy" or "yyyy-mm-dd" etc.
  out <- suppressWarnings(as.Date(paste0("15/", s), format="%d/%m/%Y"))
  if (all(is.na(out))) out <- suppressWarnings(as.Date(s))
  out
}
demog$DOB_d <- parse_dob(demog$DOB)
bl_dates <- reg %>% filter(time == 0) %>% select(RID, examdate)
bl_dates$examdate <- suppressWarnings(as.Date(bl_dates$examdate))
age_bl <- demog %>% left_join(bl_dates, by="RID") %>%
  mutate(age_baseline = as.numeric(difftime(examdate, DOB_d, units="days"))/365.25) %>%
  select(RID, age_baseline)

# -----------------------------------------------------------------------------
# 4. APPLY INCLUSION CRITERIA AND BUILD LONG FRAME
# -----------------------------------------------------------------------------
cat("\n=== Applying inclusion criteria ===\n")

# Restrict to <=36 months window for biomarkers / outcomes
within_36 <- function(d) d %>% filter(time <= 36)

gfap <- within_36(gfap)
ptau <- within_36(ptau)
adas <- within_36(adas)
mmse <- within_36(mmse)
fs   <- within_36(fs)

# Eligibility: age>=65, baseline MCI/AD, >=3 GFAP timepoints, >=3 p-tau timepoints
elig_RIDs <- age_bl %>% filter(age_baseline >= 65) %>%
  inner_join(dx_base, by="RID") %>%
  filter(diag %in% c("MCI","AD")) %>%
  inner_join(gfap %>% count(RID, name="n_gfap") %>% filter(n_gfap >= 3), by="RID") %>%
  inner_join(ptau %>% count(RID, name="n_ptau") %>% filter(n_ptau >= 3), by="RID") %>%
  pull(RID)

cat("Eligible N =", length(elig_RIDs), "\n")
stopifnot(length(elig_RIDs) >= 50)   # sanity check

# Build longitudinal frame
df_long <- expand.grid(
  RID = elig_RIDs,
  time = c(0,12,24,36)
) %>%
  left_join(gfap, by = c("RID","time")) %>%
  left_join(ptau, by = c("RID","time")) %>%
  left_join(adas, by = c("RID","time")) %>%
  left_join(mmse, by = c("RID","time")) %>%
  left_join(fs %>% select(RID, time, HV_adj), by = c("RID","time")) %>%
  mutate(time_yr = time/12,
         logGFAP = log(GFAP),
         logPTAU = log(PTAU))

# Baseline covariates
df_base <- data.frame(RID = elig_RIDs) %>%
  left_join(age_bl, by="RID") %>%
  left_join(dx_base, by="RID") %>%
  left_join(apoe, by="RID") %>%
  left_join(demog %>% select(RID, female, educ), by="RID") %>%
  left_join(vit, by="RID") %>%
  left_join(htn, by="RID") %>%
  left_join(dm,  by="RID")

df_long <- df_long %>% left_join(df_base, by="RID")

cat("Long-format rows:", nrow(df_long), "; baseline rows:", nrow(df_base), "\n")
saveRDS(df_long, file.path(CSV_DIR, "df_long_real.rds"))
saveRDS(df_base, file.path(CSV_DIR, "df_base_real.rds"))

# -----------------------------------------------------------------------------
# 5. PHASE I -- LCGMM (separate for GFAP and p-tau; multlcmm joint as sensitivity)
# -----------------------------------------------------------------------------
cat("\n=== PHASE I: LCGMM ===\n")
df_lc <- df_long %>% filter(!is.na(logGFAP) | !is.na(logPTAU))

fit_lcgmm <- function(yvar, df, label) {
  cat("\n--- LCGMM for", label, "---\n")
  d <- df %>% filter(!is.na(.data[[yvar]]))
  form <- as.formula(paste(yvar, "~ time_yr"))
  m1 <- hlme(form, random=~time_yr, subject="RID", ng=1, data=d, verbose=FALSE)
  m2 <- hlme(form, mixture=~time_yr, random=~time_yr, subject="RID", ng=2, B=m1, data=d, verbose=FALSE)
  m3 <- hlme(form, mixture=~time_yr, random=~time_yr, subject="RID", ng=3, B=m1, data=d, verbose=FALSE)
  m4 <- hlme(form, mixture=~time_yr, random=~time_yr, subject="RID", ng=4, B=m1, data=d, verbose=FALSE)
  st <- summarytable(m1,m2,m3,m4,
                     which=c("G","loglik","AIC","BIC","SABIC","entropy","%class"))
  print(st)
  # Select smallest BIC meeting min-5% class size
  bics <- c(m1$BIC,m2$BIC,m3$BIC,m4$BIC)
  mods <- list(m1,m2,m3,m4)
  best <- NULL
  for (k in order(bics)) {
    pp <- mods[[k]]$pprob
    if (k == 1) { best <- mods[[k]]; break }
    cl <- pp$class
    if (min(prop.table(table(cl))) >= 0.05) { best <- mods[[k]]; break }
  }
  cat("  selected ng =", best$ng, "\n")
  list(models = mods, selected = best, summary = st)
}

gfap_fit <- fit_lcgmm("logGFAP", df_lc, "plasma GFAP")
ptau_fit <- fit_lcgmm("logPTAU", df_lc, paste("plasma", toupper(ptau_used)))

# Joint multivariate LCGMM (sensitivity) -- may not converge with many classes
joint_fit <- tryCatch({
  cat("\n--- Joint multivariate LCGMM (sensitivity) ---\n")
  multlcmm(logGFAP + logPTAU ~ time_yr,
           random = ~time_yr, subject = "RID",
           mixture = ~time_yr, ng = 3,
           data = df_lc, verbose = FALSE)
}, error = function(e) { message("  joint multlcmm failed: ", e$message); NULL })

# Attach class memberships to baseline frame
df_base$gfap_class <- gfap_fit$selected$pprob$class[match(df_base$RID, gfap_fit$selected$pprob$RID)]
df_base$ptau_class <- ptau_fit$selected$pprob$class[match(df_base$RID, ptau_fit$selected$pprob$RID)]
df_base$gfap_class_f <- factor(df_base$gfap_class); df_base$gfap_class_f <- relevel(df_base$gfap_class_f, ref="1")
df_base$ptau_class_f <- factor(df_base$ptau_class); df_base$ptau_class_f <- relevel(df_base$ptau_class_f, ref="1")

# -----------------------------------------------------------------------------
# 6. PHASE II -- APOE e4 -> trajectory class (multinomial logit)
# -----------------------------------------------------------------------------
cat("\n=== PHASE II: APOE e4 vs trajectory class ===\n")
run_multinom <- function(class_var) {
  fit <- multinom(as.formula(paste(class_var,
                  "~ apoe_e4 + age_baseline + female + educ + diag")),
                  data = df_base, trace = FALSE)
  co <- coef(fit); se <- sqrt(diag(vcov(fit)))
  out <- data.frame()
  for (cls in rownames(co)) {
    b <- co[cls, "apoe_e4"]; s <- se[paste0(cls, ":apoe_e4")]
    p <- 2*pnorm(-abs(b/s))
    out <- rbind(out, data.frame(class=cls, OR=exp(b),
                                 lo=exp(b-1.96*s), hi=exp(b+1.96*s), p=p))
  }
  out
}
gfap_mlogit <- run_multinom("gfap_class_f");  cat("\nGFAP:\n"); print(gfap_mlogit)
ptau_mlogit <- run_multinom("ptau_class_f");  cat("\np-tau:\n"); print(ptau_mlogit)

# -----------------------------------------------------------------------------
# 7. PHASE III -- LME: class * time on ADAS-Cog 13, MMSE, HV
# -----------------------------------------------------------------------------
cat("\n=== PHASE III: trajectory class vs longitudinal outcomes ===\n")
df_long <- df_long %>%
  left_join(df_base %>% select(RID, gfap_class_f, ptau_class_f), by="RID")

run_lme <- function(yvar, class_var, label) {
  cat("\n", yvar, "~", class_var, "* time_yr\n")
  d <- df_long %>% filter(!is.na(.data[[yvar]]), !is.na(.data[[class_var]]))
  f <- as.formula(paste(yvar, "~", class_var,
                        "* time_yr + age_baseline + female + educ + diag + (time_yr | RID)"))
  fit <- lmer(f, data=d, REML=FALSE,
              control=lmerControl(optimizer="bobyqa"))
  co <- summary(fit)$coefficients
  print(round(co, 3))
  fit
}

lme_adas_gfap <- run_lme("ADAS13", "gfap_class_f", "ADAS13 ~ GFAP class")
lme_hv_gfap   <- run_lme("HV_adj","gfap_class_f", "HV ~ GFAP class")
lme_mmse_gfap <- run_lme("MMSE",  "gfap_class_f", "MMSE ~ GFAP class")
lme_adas_ptau <- run_lme("ADAS13","ptau_class_f", "ADAS13 ~ p-tau class")
lme_hv_ptau   <- run_lme("HV_adj","ptau_class_f", "HV ~ p-tau class")
lme_mmse_ptau <- run_lme("MMSE",  "ptau_class_f", "MMSE ~ p-tau class")

# -----------------------------------------------------------------------------
# 8. SUMMARY for manuscript
# -----------------------------------------------------------------------------
cat("\n=================== MANUSCRIPT SUMMARY ===================\n")
cat("Plasma p-tau marker used:", ptau_used, "\n")
cat("Eligible N:", length(elig_RIDs), "\n")
cat("Baseline n by diagnosis:\n"); print(table(df_base$diag))
cat("APOE e4 carrier %:", round(100*mean(df_base$apoe_e4, na.rm=TRUE),1), "\n")
cat("\nGFAP LCGMM selected ng =", gfap_fit$selected$ng, "\n")
cat("p-tau LCGMM selected ng =", ptau_fit$selected$ng, "\n")
cat("GFAP class sizes:\n"); print(table(df_base$gfap_class))
cat("p-tau class sizes:\n"); print(table(df_base$ptau_class))
if (!is.null(joint_fit)) {
  cat("\nJoint multlcmm converged. ng =", joint_fit$ng, "\n")
  cat("Concordance (cross-tab joint vs GFAP class):\n")
  jc <- joint_fit$pprob$class[match(df_base$RID, joint_fit$pprob$RID)]
  print(table(GFAP_only = df_base$gfap_class, Joint = jc))
}
cat("\n=== END ===\nPaste this console output back; I will fill the manuscript.\n")
