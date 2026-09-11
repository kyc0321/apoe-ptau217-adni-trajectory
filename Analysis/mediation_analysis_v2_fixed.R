# mediation_analysis_v2_fixed.R
# Re-run of mediation_analysis_v1.R using the corrected BLUP slopes
# (blup_slopes_with_SE_v2_fixed.csv, from phase3_me_correction_v2_fixed.R),
# which restore the APOE4-associated portion of the slope that v1's BLUP
# reconstruction had silently omitted (see that script's header for the bug).
# Design otherwise IDENTICAL to v1: mediation::mediate(), binary mediator
# (A+/A- via plasma AB42/40 < 0.085), nonparametric BCa bootstrap, sims=5000.
# Paths also fixed to the current "#2 LOAD-BC-NS" folder location (v1 used a
# stale pre-rename path).

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
  RID = as.integer(RID),
  EXAMDATE = as.Date(EXAMDATE),
  AB42_AB40_F = suppressWarnings(as.numeric(AB42_AB40_F)),
  AB42_AB40_F = ifelse(AB42_AB40_F <= 0 | AB42_AB40_F > 1, NA, AB42_AB40_F)
)
ab_bl <- upenn %>% filter(!is.na(AB42_AB40_F)) %>%
  group_by(RID) %>% slice_min(EXAMDATE, with_ties = FALSE) %>% ungroup() %>%
  transmute(RID, A_pos = as.integer(AB42_AB40_F < 0.085))

med_dat <- blup %>%
  inner_join(covars, by = "RID") %>%
  inner_join(ab_bl, by = "RID")
need_cols <- c("slope_ptau","APOE4","A_pos","age_z","sex","edu_z","dx3",
               "BMI_z","HTN_base","DM_base")
med_dat <- med_dat[complete.cases(med_dat[, need_cols]), ]

cat("=== Mediation analytic sample (v2, fixed BLUP) ===\n")
cat("n =", nrow(med_dat), "(of", nrow(blup), "with a BLUP slope,",
    nrow(ab_bl), "with baseline AB42/40)\n")
cat("APOE4+:", sum(med_dat$APOE4), " | A+:", sum(med_dat$A_pos), "\n\n")

mediator_model <- glm(
  A_pos ~ APOE4 + age_z + sex + edu_z + dx3 + BMI_z + HTN_base + DM_base,
  data = med_dat, family = binomial(link = "logit")
)
outcome_model <- lm(
  slope_ptau ~ APOE4 + A_pos + age_z + sex + edu_z + dx3 + BMI_z + HTN_base + DM_base,
  data = med_dat
)

cat("=== Mediator model (A+ ~ APOE4 + covariates) ===\n")
print(summary(mediator_model)$coefficients["APOE4", ])
cat("\n=== Outcome model (p-tau217 slope ~ APOE4 + A+ + covariates) ===\n")
print(summary(outcome_model)$coefficients[c("APOE4","A_pos"), ])

cat("\n=== Fitting mediation::mediate() (nonparametric bootstrap, sims=5000) ===\n")
med_out <- mediate(
  mediator_model, outcome_model,
  treat = "APOE4", mediator = "A_pos",
  boot = TRUE, boot.ci.type = "bca", sims = 5000
)
med_summary <- summary(med_out)
sink(file.path(OUT, "mediation_summary_v2_fixed.txt"))
print(med_summary)
sink()

result_row <- data.frame(
  quantity = c("ACME","ADE","Total Effect","Prop. Mediated"),
  estimate = c(med_out$d0, med_out$z0, med_out$tau.coef, med_out$n0),
  ci_lower = c(med_out$d0.ci[1], med_out$z0.ci[1], med_out$tau.ci[1], med_out$n0.ci[1]),
  ci_upper = c(med_out$d0.ci[2], med_out$z0.ci[2], med_out$tau.ci[2], med_out$n0.ci[2]),
  p_value  = c(med_out$d0.p, med_out$z0.p, med_out$tau.p, med_out$n0.p)
)
write.csv(result_row, file.path(OUT, "table_mediation_v2_fixed.csv"), row.names = FALSE)

cat("\n=== Mediation results, v2 FIXED (ACME=indirect via A+, ADE=direct, Total, Prop.Mediated) ===\n")
print(result_row)

cat("\n\n=== FOR COMPARISON: original v1 (buggy) results, re-printed from saved file ===\n")
print(read.csv(file.path(OUT, "table_mediation_v1.csv")))
