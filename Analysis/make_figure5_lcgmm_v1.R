# Figure 5 — LCGMM 2-class p-tau217 trajectory visualisation
# Author: Kuan-Yu Chu | 2026-05-26

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(lcmm)
  library(patchwork); library(showtext)
})
showtext_auto()
set.seed(20260525)

OUT <- "/Users/kuan-yuchu/Desktop/LOAD/Analysis/outputs"
FIG <- "/Users/kuan-yuchu/Desktop/LOAD/figures/v4"
DATA <- "/Users/kuan-yuchu/Desktop/LOAD/ADNI_csv"

biom_long <- utils::read.csv(file.path(OUT, "analytic_biom_long.csv"), stringsAsFactors = FALSE)
attr_tbl  <- utils::read.csv(file.path(OUT, "analytic_baseline_attr.csv"), stringsAsFactors = FALSE)
biom_long$dx3 <- factor(biom_long$dx3, levels = c("CN","MCI","AD"))
attr_tbl$dx3  <- factor(attr_tbl$dx3,  levels = c("CN","MCI","AD"))

# Refit LCGMM k=2 on log_ptau
dat <- biom_long %>% filter(!is.na(log_ptau)) %>% as.data.frame()
m1 <- hlme(log_ptau ~ time_yrs, random = ~time_yrs, subject = "RID", data = dat, ng = 1)
m2 <- hlme(log_ptau ~ time_yrs, mixture = ~time_yrs, random = ~time_yrs,
           subject = "RID", data = dat, ng = 2, B = m1)
cat("Class assignments:\n"); print(table(m2$pprob$class))

# Build per-RID class table
classes <- m2$pprob %>% transmute(RID = RID, class = factor(class, levels = c(1,2)),
                                  prob_c1 = `prob1`, prob_c2 = `prob2`)

# Order classes by mean slope so colour mapping is consistent (Slow vs Fast)
slopes_per_class <- biom_long %>% inner_join(classes, by = "RID") %>%
  group_by(class, RID) %>%
  summarise(n = n(), slope = if (n >= 2) coef(lm(log_ptau ~ time_yrs))[2] else NA_real_,
            .groups = "drop") %>%
  filter(!is.na(slope)) %>% group_by(class) %>%
  summarise(mean_slope = mean(slope), .groups = "drop")
cat("Mean per-class empirical slope:\n"); print(slopes_per_class)
# Decide which class is "Fast"
fast_class <- slopes_per_class$class[which.max(slopes_per_class$mean_slope)]
classes$class_lab <- factor(ifelse(as.character(classes$class) == as.character(fast_class),
                                   "Steep-slope class", "Gentle-slope class"),
                            levels = c("Gentle-slope class", "Steep-slope class"))

# Merge class label into long data
dat_p <- biom_long %>% inner_join(classes %>% select(RID, class_lab), by = "RID") %>%
  inner_join(attr_tbl %>% select(RID, dx3, APOE4), by = "RID") %>%
  rename(dx3 = dx3.x)

# Marginal LME fit per class for population-average line
library(lme4)
mk_fit <- function(d_one) {
  fit <- lmer(log_ptau ~ time_yrs + (1 + time_yrs | RID), data = d_one, REML = TRUE,
              control = lmerControl(optimizer = "bobyqa"))
  predict(fit, newdata = data.frame(time_yrs = seq(0, 10, 0.2)), re.form = NA)
}
new_df <- expand.grid(time_yrs = seq(0, 10, 0.2),
                      class_lab = levels(dat_p$class_lab))
for (lv in levels(dat_p$class_lab)) {
  new_df$log_ptau[new_df$class_lab == lv] <-
    mk_fit(dat_p %>% filter(class_lab == lv, !is.na(log_ptau)) %>% as.data.frame())
}

# ---- PANEL A: spaghetti plot coloured by class with marginal fits ----
panelA <- ggplot(dat_p %>% filter(!is.na(log_ptau)),
                 aes(x = time_yrs, y = log_ptau, group = RID, color = class_lab)) +
  geom_line(alpha = 0.18) +
  geom_line(data = new_df, aes(group = class_lab), size = 1.5, inherit.aes = TRUE) +
  scale_color_manual(values = c("Gentle-slope class" = "#2980b9",
                                "Steep-slope class"  = "#c0392b")) +
  labs(x = "Time from plasma baseline (years)",
       y = "log plasma p-tau217 (pg/mL)", color = NULL,
       title = "A. LCGMM 2-class p-tau217 trajectories (full diagnostic supported)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

# ---- PANEL B: per-class composition by dx, APOE, Amyloid ----
# Add amyloid status (re-derive from UPENN AB42/40)
upenn <- utils::read.csv(file.path(DATA, "UPENN_PLASMA_FUJIREBIO_QUANTERIX_24May2026.csv"),
                         stringsAsFactors = FALSE)
ab_bl <- upenn %>% mutate(RID = as.integer(RID), EXAMDATE = as.Date(EXAMDATE),
                          AB42_AB40_F = suppressWarnings(as.numeric(AB42_AB40_F))) %>%
  filter(!is.na(AB42_AB40_F), AB42_AB40_F > 0, AB42_AB40_F < 1) %>%
  group_by(RID) %>% slice_min(EXAMDATE, with_ties = FALSE) %>% ungroup() %>%
  transmute(RID, A_pos = factor(ifelse(AB42_AB40_F < 0.085, "A+", "A−"),
                                levels = c("A−","A+")))

attr_class <- attr_tbl %>%
  inner_join(classes %>% select(RID, class_lab), by = "RID") %>%
  left_join(ab_bl, by = "RID")

comp <- bind_rows(
  attr_class %>% mutate(Var = "Diagnosis", Lvl = as.character(dx3)),
  attr_class %>% mutate(Var = "APOE ε4",    Lvl = ifelse(APOE4 == 1, "ε4+", "ε4−")),
  attr_class %>% filter(!is.na(A_pos)) %>%
    mutate(Var = "Amyloid",   Lvl = as.character(A_pos))
) %>% count(class_lab, Var, Lvl) %>%
  group_by(class_lab, Var) %>% mutate(pct = 100*n/sum(n)) %>% ungroup()

comp$Var <- factor(comp$Var, levels = c("Diagnosis","APOE ε4","Amyloid"))

panelB <- ggplot(comp, aes(x = Lvl, y = pct, fill = class_lab)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  facet_wrap(~ Var, scales = "free_x") +
  scale_fill_manual(values = c("Gentle-slope class" = "#2980b9",
                               "Steep-slope class"  = "#c0392b")) +
  labs(x = NULL, y = "% within class", fill = NULL,
       title = "B. Class composition: diagnosis, APOE ε4 carriage, and amyloid status (plasma Aβ42/40 < 0.085)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))

fig5 <- panelA / panelB
ggsave(file.path(FIG, "Figure5_LCGMM_pT217_classes.png"), fig5, dpi = 300,
       width = 9, height = 7.5, bg = "white")
ggsave(file.path(FIG, "Figure5_LCGMM_pT217_classes.pdf"), fig5,
       width = 9, height = 7.5)
cat("Saved: Figure5_LCGMM_pT217_classes.png/pdf\n")

# Cross-tabs for caption
cat("\n=== Cross-tabulation: LCGMM class × clinical/biological strata ===\n")
cat("\nDiagnosis:\n")
print(round(prop.table(table(attr_class$class_lab, attr_class$dx3), margin = 1) * 100, 1))
cat("\nAPOE ε4 carrier %:\n")
print(round(prop.table(table(attr_class$class_lab, attr_class$APOE4), margin = 1) * 100, 1))
cat("\nA+ %:\n")
print(round(prop.table(table(attr_class$class_lab, attr_class$A_pos), margin = 1) * 100, 1))

# Save class assignments for downstream use
utils::write.csv(classes, file.path(OUT, "lcgmm_pT217_class_assignments.csv"), row.names = FALSE)
cat("Saved: lcgmm_pT217_class_assignments.csv\n")
