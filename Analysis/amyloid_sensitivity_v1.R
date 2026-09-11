# Amyloid-stratified sensitivity analysis (Path A response to reviewer issue #3)
# Author: Kuan-Yu Chu | 2026-05-25
# A+ definition: baseline plasma AB42/AB40 ratio < 0.085 (Fujirebio Lumipulse)
# Literature support: Janelidze 2021, Brum 2024; conservative cutoff yields ~50% A+ here
# Rerun Phase II LME (both biomarkers) + Phase III brms joint (p-tau217 × 3 outcomes only)
# GFAP-paired joint models in A+ subset skipped — full-cohort GFAP was null; A+ subset would not change conclusion

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr)
  library(lme4); library(lmerTest); library(nlme); library(brms)
})
options(brms.backend = "rstan")
set.seed(20260525)

DATA <- "/Users/kuan-yuchu/Desktop/LOAD/ADNI_csv"
OUT  <- "/Users/kuan-yuchu/Desktop/LOAD/Analysis/outputs"

# ---- Load assembled long datasets + add A+/A- flag ----
biom_long <- utils::read.csv(file.path(OUT, "analytic_biom_long.csv"), stringsAsFactors = FALSE)
attr_tbl  <- utils::read.csv(file.path(OUT, "analytic_baseline_attr.csv"), stringsAsFactors = FALSE)
biom_long$dx3 <- factor(biom_long$dx3, levels = c("CN","MCI","AD"))
attr_tbl$dx3  <- factor(attr_tbl$dx3,  levels = c("CN","MCI","AD"))

# Baseline AB42/40 per RID from UPENN file
upenn <- utils::read.csv(file.path(DATA, "UPENN_PLASMA_FUJIREBIO_QUANTERIX_24May2026.csv"),
                         stringsAsFactors = FALSE)
upenn <- upenn %>% mutate(
  RID = as.integer(RID),
  EXAMDATE = as.Date(EXAMDATE),
  AB42_AB40_F = suppressWarnings(as.numeric(AB42_AB40_F)),
  AB42_AB40_F = ifelse(AB42_AB40_F <= 0 | AB42_AB40_F > 1, NA, AB42_AB40_F)
)
ab_bl <- upenn %>% filter(!is.na(AB42_AB40_F)) %>%
  group_by(RID) %>% slice_min(EXAMDATE, with_ties = FALSE) %>% ungroup() %>%
  transmute(RID, ab42_40_bl = AB42_AB40_F,
            A_pos = as.integer(AB42_AB40_F < 0.085))

attr_tbl_a <- attr_tbl %>% left_join(ab_bl, by = "RID")
cat("=== Amyloid status in analytic cohort (N =", nrow(attr_tbl_a), ") ===\n")
cat("Baseline AB42/40 available for:", sum(!is.na(attr_tbl_a$ab42_40_bl)), "/", nrow(attr_tbl_a), "\n")
cat("A+ (AB42/40 < 0.085) at baseline:",
    sum(attr_tbl_a$A_pos == 1, na.rm = TRUE), " (",
    round(100 * mean(attr_tbl_a$A_pos == 1, na.rm = TRUE), 1), "%)\n", sep = "")
cat("By diagnosis:\n")
print(table(attr_tbl_a$dx3, ifelse(attr_tbl_a$A_pos == 1, "A+", "A-"), useNA = "ifany"))

# Add A_pos to biom_long
biom_long_a <- biom_long %>%
  left_join(attr_tbl_a %>% select(RID, A_pos, ab42_40_bl), by = "RID")

# Restrict to A+ subset
biom_apos <- biom_long_a %>% filter(A_pos == 1)
attr_apos <- attr_tbl_a %>% filter(A_pos == 1)
cat("\nA+ subset for Phase II/III: n_id =", length(unique(biom_apos$RID)),
    "; n_obs (any biom) =", nrow(biom_apos), "\n\n")

# ---- Phase II LME on A+ subset ----
phase2_apos <- function(yvar) {
  d <- biom_apos %>% filter(!is.na(.data[[yvar]]), !is.na(BMI_z))
  fml <- as.formula(paste0(yvar,
    " ~ time_yrs*APOE4 + age_z + sex + edu_z + dx3 + BMI_z + HTN_base + DM_base + (1 + time_yrs | RID)"))
  m <- lmerTest::lmer(fml, data = d, REML = TRUE,
                      control = lmerControl(optimizer="bobyqa"))
  cf <- as.data.frame(summary(m)$coefficients)
  cf$term <- rownames(cf); cf$n_obs <- nrow(d); cf$n_id <- length(unique(d$RID))
  ci <- tryCatch(confint(m, parm = "beta_", method = "Wald"), error = function(e) NULL)
  if (!is.null(ci)) { ci <- as.data.frame(ci); ci$term <- rownames(ci); cf <- left_join(cf, ci, by="term") }
  cf
}
cat("=== Phase II A+ — log_gfap ===\n")
p2g <- phase2_apos("log_gfap")
print(p2g[grepl("time_yrs:APOE4|^APOE4$|^time_yrs$", p2g$term),
          c("term","Estimate","Std. Error","Pr(>|t|)","n_obs","n_id")])
cat("\n=== Phase II A+ — log_ptau ===\n")
p2p <- phase2_apos("log_ptau")
print(p2p[grepl("time_yrs:APOE4|^APOE4$|^time_yrs$", p2p$term),
          c("term","Estimate","Std. Error","Pr(>|t|)","n_obs","n_id")])
phase2_apos_out <- bind_rows(p2g %>% mutate(biomarker = "log_gfap"),
                              p2p %>% mutate(biomarker = "log_ptau"))
utils::write.csv(phase2_apos_out, file.path(OUT, "phase2_amyloid_pos.csv"), row.names = FALSE)

# ---- Phase III brms joint on A+ subset (p-tau217 × 3 outcomes only) ----
# Build outcome long tables restricted to A+ RIDs
adas <- utils::read.csv(file.path(DATA, "ADAS_24May2026.csv"), stringsAsFactors = FALSE, na.strings = c("","NA"))
mmse <- utils::read.csv(file.path(DATA, "MMSE_24May2026.csv"), stringsAsFactors = FALSE, na.strings = c("","NA"))
fs   <- utils::read.csv(file.path(DATA, "ADSP_PHC_T1_FS_24May2026.csv"), stringsAsFactors = FALSE, na.strings = c("","NA"))

zsc <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
apos_rid <- attr_apos$RID
bl_dates <- biom_apos %>% group_by(RID) %>%
  summarise(bl_date = as.Date(min(EXAMDATE)), .groups = "drop")

build_outc <- function(df, ycol, datecol) {
  df %>% mutate(RID = as.integer(RID), EXAMDATE = as.Date(.data[[datecol]]),
                y = suppressWarnings(as.numeric(.data[[ycol]]))) %>%
    filter(!is.na(EXAMDATE), !is.na(y), RID %in% apos_rid) %>%
    left_join(bl_dates, by = "RID") %>%
    mutate(time_yrs = as.numeric(EXAMDATE - bl_date)/365.25) %>%
    filter(!is.na(time_yrs))
}
adas_long <- build_outc(adas, "TOTAL13",  "VISDATE")  %>% transmute(RID, time_yrs, outc = y)
mmse_long <- build_outc(mmse, "MMSCORE",  "VISDATE")  %>% transmute(RID, time_yrs, outc = y)
fs_long   <- fs %>% mutate(hippo_ratio = (Left.Hippocampus_combat + Right.Hippocampus_combat) /
                             EstimatedTotalIntraCranialVol_combat) %>%
  build_outc("hippo_ratio", "PHC_SCANDATE") %>% transmute(RID, time_yrs, outc = y * 1000)

attr_apos_z <- attr_apos %>% mutate(age_z = zsc(age_bl), edu_z = zsc(edu))
biom_red <- biom_apos %>% filter(!is.na(BMI_z)) %>%
  transmute(RID, time_yrs, bm = log_ptau) %>% filter(!is.na(bm))

fit_one <- function(outc_df, lab) {
  cat("\n=== Fitting A+:", lab, "===  ", as.character(Sys.time()), "\n")
  rid_keep <- intersect(unique(biom_red$RID), unique(outc_df$RID))
  bm_df <- biom_red[biom_red$RID %in% rid_keep, ]
  oc_df <- outc_df[outc_df$RID %in% rid_keep, ]
  bm_rows <- bm_df %>% mutate(outc = 0, has_bm = TRUE, has_oc = FALSE)
  oc_rows <- oc_df %>% mutate(bm = 0,  has_bm = FALSE, has_oc = TRUE)
  d <- bind_rows(bm_rows, oc_rows) %>%
    left_join(attr_apos_z %>% select(RID, age_z, sex, edu_z, dx3, APOE4), by = "RID")
  d$RID <- factor(d$RID)
  cat("  n_id=", length(rid_keep), " | bm rows=", nrow(bm_df), " | outc rows=", nrow(oc_df), "\n", sep="")
  f_bm <- bf(bm   | subset(has_bm) ~ time_yrs + APOE4 + age_z + sex + dx3 + (1 + time_yrs | p | RID))
  f_oc <- bf(outc | subset(has_oc) ~ time_yrs + APOE4 + age_z + sex + edu_z + dx3 + (1 + time_yrs | p | RID))
  fit <- tryCatch(
    brm(f_bm + f_oc + set_rescor(FALSE), data = d, family = gaussian(),
        chains = 4, iter = 4000, warmup = 1500, cores = 4,
        seed = 20260525, refresh = 0,
        control = list(adapt_delta = 0.97, max_treedepth = 13)),
    error = function(e) { message(" brms err: ", e$message); NULL })
  if (is.null(fit)) return(NULL)
  ps <- as_draws_df(fit)
  cor_nm <- grep("^cor_RID__bm_time_yrs__outc_time_yrs$", names(ps), value = TRUE)
  draws <- ps[[cor_nm]]
  summ <- data.frame(label = lab, parameter = cor_nm,
                     mean = mean(draws), median = median(draws),
                     lo = quantile(draws, .025), hi = quantile(draws, .975),
                     p_dir = mean(draws < 0),
                     rhat = max(brms::rhat(fit)[grep("cor_RID", names(brms::rhat(fit)))]),
                     n_id = length(rid_keep), n_bm = nrow(bm_df), n_oc = nrow(oc_df))
  cat("  rho =", round(summ$mean, 3), " [", round(summ$lo, 3), ", ", round(summ$hi, 3),
      "]  Rhat=", round(summ$rhat, 3), "\n", sep = "")
  saveRDS(fit, file.path(OUT, paste0("brmsfit_Apos_", lab, ".rds")))
  summ
}

results <- bind_rows(
  fit_one(adas_long, "ADAS_pT217_Apos"),
  fit_one(mmse_long, "MMSE_pT217_Apos"),
  fit_one(fs_long,   "Hippo_pT217_Apos")
)
utils::write.csv(results, file.path(OUT, "table4_brms_amyloid_pos.csv"), row.names = FALSE)
cat("\n=== A+ subset brms joint results ===\n"); print(results)
