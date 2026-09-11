# mediation_analysis_v3_no_dx3.R
# Sensitivity analysis for Reviewer #1 comment 4 (Neurological Sciences,
# decision 2026-09-10): baseline diagnosis (dx3: CN/MCI/AD) may lie
# downstream of APOE4/amyloid/tau pathology, making it a potential collider
# rather than a legitimate confounder when included in the mediator and
# outcome models. Re-runs the CORRECTED (v2) mediation analysis with dx3
# dropped from both models, to check whether the 6.4% proportion-mediated
# result depends on conditioning on a variable that may itself be a
# downstream consequence of the exposure/mediator.
#
# NOTE: models are written out explicitly (not built via paste+as.formula)
# because mediation::mediate() re-evaluates/updates the model formulas
# internally during bootstrap resampling, and a dynamically-pasted formula
# object loses access to variables (like a local `rhs` string) that lived
# only in the constructing function's environment -- caught as a real
# "object 'rhs' not found" error on first attempt.

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
  transmute(RID, A_pos = as.integer(AB42_AB40_F < 0.085))

med_dat <- blup %>% inner_join(covars, by = "RID") %>% inner_join(ab_bl, by = "RID")
need_cols <- c("slope_ptau","APOE4","A_pos","age_z","sex","edu_z","BMI_z","HTN_base","DM_base")
med_dat <- med_dat[complete.cases(med_dat[, need_cols]), ]
cat("=== Sensitivity sample n =", nrow(med_dat), "===\n\n")

summarize_med <- function(med_out, label) {
  cat("===", label, "===\n")
  res <- data.frame(
    quantity = c("ACME","ADE","Total Effect","Prop. Mediated"),
    estimate = c(med_out$d0, med_out$z0, med_out$tau.coef, med_out$n0),
    ci_lower = c(med_out$d0.ci[1], med_out$z0.ci[1], med_out$tau.ci[1], med_out$n0.ci[1]),
    ci_upper = c(med_out$d0.ci[2], med_out$z0.ci[2], med_out$tau.ci[2], med_out$n0.ci[2]),
    p_value  = c(med_out$d0.p, med_out$z0.p, med_out$tau.p, med_out$n0.p)
  )
  print(res)
  cat("\n")
  res
}

# ---- WITHOUT dx3 (sensitivity: dx3 excluded as a possible collider) --------
mediator_model_nodx3 <- glm(
  A_pos ~ APOE4 + age_z + sex + edu_z + BMI_z + HTN_base + DM_base,
  data = med_dat, family = binomial(link = "logit")
)
outcome_model_nodx3 <- lm(
  slope_ptau ~ APOE4 + A_pos + age_z + sex + edu_z + BMI_z + HTN_base + DM_base,
  data = med_dat
)
med_out_nodx3 <- mediate(mediator_model_nodx3, outcome_model_nodx3, treat = "APOE4", mediator = "A_pos",
                          boot = TRUE, boot.ci.type = "bca", sims = 5000)
res_no_dx3 <- summarize_med(med_out_nodx3, "WITHOUT dx3 (collider-concern sensitivity)")

write.csv(res_no_dx3, file.path(OUT, "table_mediation_v3_no_dx3.csv"), row.names = FALSE)
cat("Saved table_mediation_v3_no_dx3.csv\n\n")

cat("=== FOR COMPARISON: v2 (WITH dx3) result, re-printed from saved file ===\n")
print(read.csv(file.path(OUT, "table_mediation_v2_fixed.csv")))
