# =============================================================================
# mediation_analysis_v1.R
# Formal causal mediation analysis: APOE e4 -> amyloid positivity (A+/A-) ->
# plasma p-tau217 slope.
#
# Motivated by Brain Communications Reviewer 2 (BRAINCOM-2026-676, decision
# 27-Jul-2026): "mediation is inferred from attenuation within a subgroup
# instead of a formal mediation analysis." The manuscript's own Limitations
# section (Sec 4.6) had already flagged this as a gap and proposed exactly
# this fix ("a counterfactual-based decomposition ... to quantify the
# proportion mediated"); this script implements it using the amyloid
# ascertainment already available (plasma AB42/40), i.e. without new PET
# data (Tier 1 scope; PET-based re-ascertainment is a separate, larger
# undertaking noted as a limitation below).
#
# Design: mediation::mediate() (Imai, Keele & Tingley 2010) with a binary
# mediator (A+/A-) and continuous outcome (per-participant BLUP slope of
# log plasma p-tau217, from the primary LME). Nonparametric bootstrap
# (BCa), sims = 5000.
#
# Caveat carried into the write-up: the outcome is a BLUP (best linear
# unbiased predictor) slope, which is subject to shrinkage/measurement
# error -- the manuscript's own Sec 3.5c already flags BLUP non-normality
# for this reason. This is a standard, disclosed limitation of any
# mediation analysis that must use a per-subject longitudinal summary
# rather than a raw repeated-measures outcome, and is noted explicitly in
# the results reported here.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(mediation)
})

DATA <- "/Users/kuan-yuchu/Desktop/LOAD-BC/ADNI_csv"
OUT  <- "/Users/kuan-yuchu/Desktop/LOAD-BC/Analysis/outputs"
set.seed(20260728)

# ---- 1. Per-participant BLUP slope for log p-tau217 (already computed) -----
blup <- read.csv(file.path(OUT, "blup_slopes_with_SE.csv"))

# ---- 2. Baseline covariates + APOE4 (one row per RID) ----------------------
biom_long <- read.csv(file.path(OUT, "analytic_biom_long.csv"))
covars <- biom_long %>%
  distinct(RID, .keep_all = TRUE) %>%
  transmute(RID, APOE4 = as.integer(APOE4), age_z, sex = factor(sex),
            edu_z, dx3 = factor(dx3, levels = c("CN","MCI","AD")),
            BMI_z, HTN_base = as.integer(HTN_base), DM_base = as.integer(DM_base))

# ---- 3. Baseline amyloid positivity (A+/A-), same definition as the ------
#         existing amyloid-stratified sensitivity (Sec 3.5): plasma
#         AB42/40 < 0.085 (Fujirebio Lumipulse; Janelidze 2021, Brum 2024) --
upenn <- read.csv(file.path(DATA, "UPENN_PLASMA_FUJIREBIO_QUANTERIX_24May2026.csv"))
upenn <- upenn %>% mutate(
  RID = as.integer(RID),
  EXAMDATE = as.Date(EXAMDATE),
  AB42_AB40_F = suppressWarnings(as.numeric(AB42_AB40_F)),
  AB42_AB40_F = ifelse(AB42_AB40_F <= 0 | AB42_AB40_F > 1, NA, AB42_AB40_F)
)
ab_bl <- upenn %>% filter(!is.na(AB42_AB40_F)) %>%
  group_by(RID) %>% slice_min(EXAMDATE, with_ties = FALSE) %>% ungroup() %>%
  transmute(RID, A_pos = as.integer(AB42_AB40_F < 0.085))

# ---- 4. Assemble the mediation analytic dataset ----------------------------
med_dat <- blup %>%
  inner_join(covars, by = "RID") %>%
  inner_join(ab_bl, by = "RID")
need_cols <- c("slope_ptau","APOE4","A_pos","age_z","sex","edu_z","dx3",
               "BMI_z","HTN_base","DM_base")
med_dat <- med_dat[complete.cases(med_dat[, need_cols]), ]

cat("=== Mediation analytic sample ===\n")
cat("n =", nrow(med_dat), "(of", nrow(blup), "with a BLUP slope,",
    nrow(ab_bl), "with baseline AB42/40)\n")
cat("APOE4+:", sum(med_dat$APOE4), " | A+:", sum(med_dat$A_pos), "\n\n")

# ---- 5. Mediator model: A+ ~ APOE4 + covariates (logistic) ------------------
mediator_model <- glm(
  A_pos ~ APOE4 + age_z + sex + edu_z + dx3 + BMI_z + HTN_base + DM_base,
  data = med_dat, family = binomial(link = "logit")
)

# ---- 6. Outcome model: p-tau217 slope ~ APOE4 + A+ + covariates (linear) ----
outcome_model <- lm(
  slope_ptau ~ APOE4 + A_pos + age_z + sex + edu_z + dx3 + BMI_z + HTN_base + DM_base,
  data = med_dat
)

cat("=== Mediator model (A+ ~ APOE4 + covariates) ===\n")
print(summary(mediator_model)$coefficients["APOE4", ])
cat("\n=== Outcome model (p-tau217 slope ~ APOE4 + A+ + covariates) ===\n")
print(summary(outcome_model)$coefficients[c("APOE4","A_pos"), ])

# ---- 7. Formal causal mediation analysis (bootstrap) ------------------------
cat("\n=== Fitting mediation::mediate() (nonparametric bootstrap, sims=5000) ===\n")
med_out <- mediate(
  mediator_model, outcome_model,
  treat = "APOE4", mediator = "A_pos",
  boot = TRUE, boot.ci.type = "bca", sims = 5000
)

med_summary <- summary(med_out)
sink(file.path(OUT, "mediation_summary_v1.txt"))
print(med_summary)
sink()

result_row <- data.frame(
  quantity = c("ACME","ADE","Total Effect","Prop. Mediated"),
  estimate = c(med_out$d0, med_out$z0, med_out$tau.coef, med_out$n0),
  ci_lower = c(med_out$d0.ci[1], med_out$z0.ci[1], med_out$tau.ci[1], med_out$n0.ci[1]),
  ci_upper = c(med_out$d0.ci[2], med_out$z0.ci[2], med_out$tau.ci[2], med_out$n0.ci[2]),
  p_value  = c(med_out$d0.p, med_out$z0.p, med_out$tau.p, med_out$n0.p)
)
write.csv(result_row, file.path(OUT, "table_mediation_v1.csv"), row.names = FALSE)

cat("\n=== Mediation results (ACME=indirect via A+, ADE=direct, Total, Prop.Mediated) ===\n")
print(result_row)
cat("\nFull summary written to", file.path(OUT, "mediation_summary_v1.txt"), "\n")
cat("Table written to", file.path(OUT, "table_mediation_v1.csv"), "\n")
