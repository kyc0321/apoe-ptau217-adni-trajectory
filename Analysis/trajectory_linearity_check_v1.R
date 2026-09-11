suppressPackageStartupMessages({library(dplyr); library(lme4); library(lmerTest)})
BASE <- "/Users/kuan-yuchu/Desktop/#2 LOAD-BC-NS"
OUT  <- file.path(BASE, "Analysis/outputs")

biom_long <- read.csv(file.path(OUT, "analytic_biom_long.csv"))
biom_long$dx3 <- factor(biom_long$dx3, levels = c("CN","MCI","AD"))

d <- biom_long %>% filter(!is.na(log_ptau), !is.na(BMI_z))
cat("n_obs =", nrow(d), " n_id =", length(unique(d$RID)), "\n")
obs_per_id <- d %>% count(RID)
cat("Measurements per participant: "); print(table(obs_per_id$n))

cat("\n=== Primary (linear time) model, for reference ===\n")
m_lin <- lmerTest::lmer(log_ptau ~ time_yrs*APOE4 + age_z + sex + edu_z + dx3 + BMI_z + HTN_base + DM_base +
                         (1 + time_yrs | RID), data = d, REML = TRUE, control = lmerControl(optimizer="bobyqa"))
print(summary(m_lin)$coefficients["time_yrs:APOE4", ])

cat("\n=== Sensitivity: quadratic time added (fixed effect only -- random quadratic slope not identifiable with median 2-3 obs/person) ===\n")
d$time_yrs2 <- d$time_yrs^2
m_quad <- lmerTest::lmer(log_ptau ~ time_yrs*APOE4 + time_yrs2*APOE4 + age_z + sex + edu_z + dx3 + BMI_z + HTN_base + DM_base +
                          (1 + time_yrs | RID), data = d, REML = TRUE, control = lmerControl(optimizer="bobyqa"))
co <- summary(m_quad)$coefficients
print(co[grepl("time_yrs|time_yrs2", rownames(co)), c("Estimate","Std. Error","Pr(>|t|)")])

cat("\n=== Model comparison (does adding quadratic time improve fit? refit ML for LRT) ===\n")
m_lin_ml  <- lmerTest::lmer(log_ptau ~ time_yrs*APOE4 + age_z + sex + edu_z + dx3 + BMI_z + HTN_base + DM_base +
                            (1 + time_yrs | RID), data = d, REML = FALSE, control = lmerControl(optimizer="bobyqa"))
m_quad_ml <- lmerTest::lmer(log_ptau ~ time_yrs*APOE4 + time_yrs2*APOE4 + age_z + sex + edu_z + dx3 + BMI_z + HTN_base + DM_base +
                            (1 + time_yrs | RID), data = d, REML = FALSE, control = lmerControl(optimizer="bobyqa"))
print(anova(m_lin_ml, m_quad_ml))

cat("\n=== Restricted cubic spline sensitivity (3 knots on time_yrs, fixed-effect only) ===\n")
if (!requireNamespace("splines", quietly=TRUE)) install.packages("splines")
library(splines)
m_spline <- lmerTest::lmer(log_ptau ~ ns(time_yrs, df=2)*APOE4 + age_z + sex + edu_z + dx3 + BMI_z + HTN_base + DM_base +
                            (1 + time_yrs | RID), data = d, REML = TRUE, control = lmerControl(optimizer="bobyqa"))
co_s <- summary(m_spline)$coefficients
print(co_s[grepl("ns\\(time_yrs", rownames(co_s)), c("Estimate","Std. Error","Pr(>|t|)")])
