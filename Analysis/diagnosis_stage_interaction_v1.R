suppressPackageStartupMessages({library(dplyr); library(lme4); library(lmerTest)})
BASE <- "/Users/kuan-yuchu/Desktop/#2 LOAD-BC-NS"
OUT  <- file.path(BASE, "Analysis/outputs")
biom_long <- read.csv(file.path(OUT, "analytic_biom_long.csv"))
biom_long$dx3 <- factor(biom_long$dx3, levels = c("CN","MCI","AD"))
d <- biom_long %>% filter(!is.na(log_ptau), !is.na(BMI_z))

cat("=== Model D: additive stage-dependent slopes, time x APOE4 + time x dx3, NO 3-way ===\n")
m_additive <- lmerTest::lmer(log_ptau ~ time_yrs*APOE4 + time_yrs*dx3 + age_z + sex + edu_z + BMI_z + HTN_base + DM_base +
                              (1 + time_yrs | RID), data = d, REML = FALSE, control = lmerControl(optimizer="bobyqa"))
co <- summary(m_additive)$coefficients
cat("time_yrs:APOE4 with stage-dependent slopes included:\n")
print(co[grepl("^time_yrs:APOE4$", rownames(co)), c("Estimate","Std. Error","Pr(>|t|)")])

cat("\n=== LRT: does the 3-way term add anything BEYOND additive stage-dependent slopes? ===\n")
m_3way_ml <- lmerTest::lmer(log_ptau ~ time_yrs*APOE4*dx3 + age_z + sex + edu_z + BMI_z + HTN_base + DM_base +
                             (1 + time_yrs | RID), data = d, REML = FALSE, control = lmerControl(optimizer="bobyqa"))
print(anova(m_additive, m_3way_ml))
