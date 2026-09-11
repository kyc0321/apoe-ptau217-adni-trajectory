# mediation_analysis_v4_continuous_amyloid.R
# Sensitivity analysis for Reviewer #1 comment 5 (Neurological Sciences,
# decision 2026-09-10): the manuscript dichotomizes amyloid status at
# AB42/40 < 0.085, but this is "a relatively noisy mediator for a process
# that evolves continuously." Re-runs the CORRECTED (v2, dx3 included --
# the primary specification, since v3 showed the collider concern doesn't
# change the story) mediation analysis with baseline AB42/40 modeled
# CONTINUOUSLY as the mediator instead of the binary A+/A- indicator.
#
# NOTE: AB42/40 is INVERSELY related to amyloid pathology (lower ratio =
# more amyloid), so the mediator-model and outcome-model coefficient SIGNS
# on AB42_AB40_F are expected to flip relative to the A_pos (binary)
# versions -- this is expected, not an error. Proportion mediated is
# compared in magnitude, not raw sign.

suppressPackageStartupMessages({
  library(dplyr); library(mediation)
})

BASE <- "/Users/kuan-yuchu/Desktop/#2 LOAD-BC-NS"
DATA <- file.path(BASE, "ADNI_csv")
OUT  <- file.path(BASE, "Analysis/outputs")
set.seed(20260728)

blup <- read.csv(file.path(OUT, "blup_slopes_with_SE_v2_fixed.csv"))
biom_long <- read.csv(file.path(OUT, "analytic_biom_long.csv"))
covars <- biom_long %>%
  distinct(RID, .keep_all = TRUE) %>%
  transmute(RID, APOE4 = as.integer(APOE4), age_z, sex = factor(sex),
            edu_z, dx3 = factor(dx3, levels = c("CN","MCI","AD")),
            BMI_z, HTN_base = as.integer(HTN_base), DM_base = as.integer(DM_base))

upenn <- read.csv(file.path(DATA, "UPENN_PLASMA_FUJIREBIO_QUANTERIX_24May2026.csv"))
upenn <- upenn %>% mutate(
  RID = as.integer(RID), EXAMDATE = as.Date(EXAMDATE),
  AB42_AB40_F = suppressWarnings(as.numeric(AB42_AB40_F)),
  AB42_AB40_F = ifelse(AB42_AB40_F <= 0 | AB42_AB40_F > 1, NA, AB42_AB40_F)
)
ab_bl <- upenn %>% filter(!is.na(AB42_AB40_F)) %>%
  group_by(RID) %>% slice_min(EXAMDATE, with_ties = FALSE) %>% ungroup() %>%
  transmute(RID, AB4240_cont = AB42_AB40_F)   # keep continuous, no dichotomization

med_dat <- blup %>% inner_join(covars, by = "RID") %>% inner_join(ab_bl, by = "RID")
need_cols <- c("slope_ptau","APOE4","AB4240_cont","age_z","sex","edu_z","dx3","BMI_z","HTN_base","DM_base")
med_dat <- med_dat[complete.cases(med_dat[, need_cols]), ]
cat("=== Continuous-mediator sample n =", nrow(med_dat), "===\n")
cat("AB42/40 summary:\n"); print(summary(med_dat$AB4240_cont))
cat("\n")

mediator_model_cont <- lm(
  AB4240_cont ~ APOE4 + age_z + sex + edu_z + dx3 + BMI_z + HTN_base + DM_base,
  data = med_dat
)
outcome_model_cont <- lm(
  slope_ptau ~ APOE4 + AB4240_cont + age_z + sex + edu_z + dx3 + BMI_z + HTN_base + DM_base,
  data = med_dat
)

cat("=== Mediator model (AB42/40 ~ APOE4 + covariates) ===\n")
print(summary(mediator_model_cont)$coefficients["APOE4", ])
cat("\n=== Outcome model (p-tau217 slope ~ APOE4 + AB42/40 + covariates) ===\n")
print(summary(outcome_model_cont)$coefficients[c("APOE4","AB4240_cont"), ])

cat("\n=== Fitting mediation::mediate(), continuous mediator, sims=5000 ===\n")
med_out_cont <- mediate(mediator_model_cont, outcome_model_cont, treat = "APOE4", mediator = "AB4240_cont",
                         boot = TRUE, boot.ci.type = "bca", sims = 5000)

res_cont <- data.frame(
  quantity = c("ACME","ADE","Total Effect","Prop. Mediated"),
  estimate = c(med_out_cont$d0, med_out_cont$z0, med_out_cont$tau.coef, med_out_cont$n0),
  ci_lower = c(med_out_cont$d0.ci[1], med_out_cont$z0.ci[1], med_out_cont$tau.ci[1], med_out_cont$n0.ci[1]),
  ci_upper = c(med_out_cont$d0.ci[2], med_out_cont$z0.ci[2], med_out_cont$tau.ci[2], med_out_cont$n0.ci[2]),
  p_value  = c(med_out_cont$d0.p, med_out_cont$z0.p, med_out_cont$tau.p, med_out_cont$n0.p)
)
cat("\n=== Mediation results, CONTINUOUS AB42/40 mediator ===\n")
print(res_cont)
write.csv(res_cont, file.path(OUT, "table_mediation_v4_continuous.csv"), row.names = FALSE)
cat("\nSaved table_mediation_v4_continuous.csv\n")

cat("\n=== FOR COMPARISON: v2 (binary A+/A-, WITH dx3, primary spec) ===\n")
print(read.csv(file.path(OUT, "table_mediation_v2_fixed.csv")))
