suppressPackageStartupMessages({library(dplyr); library(lme4); library(lmerTest)})
BASE <- "/Users/kuan-yuchu/Desktop/#2 LOAD-BC-NS"
DATA <- file.path(BASE, "ADNI_csv")
OUT  <- file.path(BASE, "Analysis/outputs")

apoe_raw <- read.csv(file.path(DATA, "APOERES_24May2026.csv"), stringsAsFactors = FALSE)
apoe <- apoe_raw %>%
  filter(GENOTYPE != "", !is.na(GENOTYPE)) %>%
  mutate(RID = as.integer(RID)) %>%
  distinct(RID, .keep_all = TRUE) %>%   # one genotype per person (should be constant across visits)
  mutate(
    a1 = as.integer(substr(GENOTYPE, 1, 1)),
    a2 = as.integer(substr(GENOTYPE, 3, 3)),
    apoe4_dose = (a1 == 4) + (a2 == 4),
    apoe4_homo = as.integer(a1 == 4 & a2 == 4),
    apoe2_carrier = as.integer(a1 == 2 | a2 == 2)
  ) %>%
  select(RID, GENOTYPE, a1, a2, apoe4_dose, apoe4_homo, apoe2_carrier)

biom_long <- read.csv(file.path(OUT, "analytic_biom_long.csv"))
biom_long$dx3 <- factor(biom_long$dx3, levels = c("CN","MCI","AD"))
d <- biom_long %>% filter(!is.na(log_ptau), !is.na(BMI_z)) %>%
  left_join(apoe, by = "RID")

cat("=== Genotype merge check ===\n")
cat("n_obs with genotype available:", sum(!is.na(d$GENOTYPE)), "/", nrow(d), "\n")
cat("n_obs where binary APOE4 flag disagrees with genotype-derived dose>0:",
    sum((d$APOE4==1) != (d$apoe4_dose>0), na.rm=TRUE), "\n\n")

cat("=== Genotype distribution, unique participants in analytic sample ===\n")
gt_id <- d %>% distinct(RID, .keep_all = TRUE)
print(table(gt_id$GENOTYPE, useNA = "ifany"))
cat("\napoe4_dose distribution (0/1/2 alleles):\n")
print(table(gt_id$apoe4_dose, useNA = "ifany"))
cat("\nepsilon4 homozygotes:", sum(gt_id$apoe4_homo, na.rm=TRUE), "\n")
cat("epsilon2 carriers:", sum(gt_id$apoe2_carrier, na.rm=TRUE), "\n\n")

cat("=== Model A: APOE4 dose (0/1/2) instead of binary carrier status ===\n")
d_dose <- d %>% filter(!is.na(apoe4_dose))
m_dose <- lmerTest::lmer(log_ptau ~ time_yrs*apoe4_dose + age_z + sex + edu_z + dx3 + BMI_z + HTN_base + DM_base +
                          (1 + time_yrs | RID), data = d_dose, REML = TRUE, control = lmerControl(optimizer="bobyqa"))
print(summary(m_dose)$coefficients["time_yrs:apoe4_dose", ])
cat("n_obs =", nrow(d_dose), " n_id =", length(unique(d_dose$RID)), "\n\n")

cat("=== Model B: primary model EXCLUDING e4 homozygotes (heterozygote vs non-carrier only) ===\n")
d_nohomo <- d %>% filter(apoe4_homo == 0 | is.na(apoe4_homo))
m_nohomo <- lmerTest::lmer(log_ptau ~ time_yrs*APOE4 + age_z + sex + edu_z + dx3 + BMI_z + HTN_base + DM_base +
                            (1 + time_yrs | RID), data = d_nohomo, REML = TRUE, control = lmerControl(optimizer="bobyqa"))
print(summary(m_nohomo)$coefficients["time_yrs:APOE4", ])
cat("n_obs =", nrow(d_nohomo), " n_id =", length(unique(d_nohomo$RID)),
    " (excluded", sum(gt_id$apoe4_homo==1, na.rm=TRUE), "homozygotes)\n\n")

cat("=== Model C: heterozygote-only vs homozygote-only dose-response check (3-level factor) ===\n")
d_dose$dose_f <- factor(d_dose$apoe4_dose, levels=c(0,1,2), labels=c("0 alleles","1 allele (het)","2 alleles (homo)"))
m_dosef <- lmerTest::lmer(log_ptau ~ time_yrs*dose_f + age_z + sex + edu_z + dx3 + BMI_z + HTN_base + DM_base +
                           (1 + time_yrs | RID), data = d_dose, REML = TRUE, control = lmerControl(optimizer="bobyqa"))
co <- summary(m_dosef)$coefficients
print(co[grepl("time_yrs:dose_f", rownames(co)), c("Estimate","Std. Error","Pr(>|t|)")])
