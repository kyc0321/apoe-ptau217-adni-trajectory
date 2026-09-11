suppressPackageStartupMessages({library(dplyr); library(lme4); library(lmerTest)})
BASE <- "/Users/kuan-yuchu/Desktop/#2 LOAD-BC-NS"
DATA <- file.path(BASE, "ADNI_csv")
OUT  <- file.path(BASE, "Analysis/outputs")

biom_long <- read.csv(file.path(OUT, "analytic_biom_long.csv"))
biom_long$dx3 <- factor(biom_long$dx3, levels = c("CN","MCI","AD"))

upenn <- read.csv(file.path(DATA, "UPENN_PLASMA_FUJIREBIO_QUANTERIX_24May2026.csv"))
upenn <- upenn %>% mutate(RID = as.integer(RID), EXAMDATE = as.Date(EXAMDATE),
  AB42_AB40_F = suppressWarnings(as.numeric(AB42_AB40_F)),
  AB42_AB40_F = ifelse(AB42_AB40_F <= 0 | AB42_AB40_F > 1, NA, AB42_AB40_F))
ab_bl <- upenn %>% filter(!is.na(AB42_AB40_F)) %>% group_by(RID) %>%
  slice_min(EXAMDATE, with_ties = FALSE) %>% ungroup() %>%
  transmute(RID, A_pos = as.integer(AB42_AB40_F < 0.085))

biom_long_a <- biom_long %>% left_join(ab_bl, by = "RID")

cat("=== 1. Reproduce original A+ subgroup Phase II refit (log_ptau) ===\n")
d_apos <- biom_long_a %>% filter(A_pos == 1, !is.na(log_ptau), !is.na(BMI_z))
cat("n_obs =", nrow(d_apos), " n_id =", length(unique(d_apos$RID)), "\n")
m_apos <- lmerTest::lmer(log_ptau ~ time_yrs*APOE4 + age_z + sex + edu_z + dx3 + BMI_z + HTN_base + DM_base +
                          (1 + time_yrs | RID), data = d_apos, REML = TRUE, control = lmerControl(optimizer="bobyqa"))
print(summary(m_apos)$coefficients["time_yrs:APOE4", ])

cat("\n=== 2. Full-cohort 3-way interaction test: time x APOE4 x A_pos (does amyloid MODIFY the APOE-time slope?) ===\n")
d_full <- biom_long_a %>% filter(!is.na(log_ptau), !is.na(BMI_z), !is.na(A_pos))
cat("n_obs =", nrow(d_full), " n_id =", length(unique(d_full$RID)), "\n")
m_full3way <- lmerTest::lmer(log_ptau ~ time_yrs*APOE4*A_pos + age_z + sex + edu_z + dx3 + BMI_z + HTN_base + DM_base +
                              (1 + time_yrs | RID), data = d_full, REML = TRUE, control = lmerControl(optimizer="bobyqa"))
cat("\n-- key terms --\n")
co <- summary(m_full3way)$coefficients
print(co[grepl("time_yrs:APOE4|time_yrs:A_pos|time_yrs:APOE4:A_pos", rownames(co)), ])

cat("\n=== 3. For comparison: A- (amyloid-negative) subgroup refit ===\n")
d_aneg <- biom_long_a %>% filter(A_pos == 0, !is.na(log_ptau), !is.na(BMI_z))
cat("n_obs =", nrow(d_aneg), " n_id =", length(unique(d_aneg$RID)), "\n")
m_aneg <- lmerTest::lmer(log_ptau ~ time_yrs*APOE4 + age_z + sex + edu_z + dx3 + BMI_z + HTN_base + DM_base +
                          (1 + time_yrs | RID), data = d_aneg, REML = TRUE, control = lmerControl(optimizer="bobyqa"))
print(summary(m_aneg)$coefficients["time_yrs:APOE4", ])
