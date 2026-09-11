# Generate Figures 1-4 for manuscript (v4)
# Author: Kuan-Yu Chu | 2026-05-25
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
  library(showtext)
})
# Use showtext to embed Unicode glyphs (Greek ρ etc.) into PDF/PNG cleanly
showtext_auto()
font_add_google <- tryCatch(showtext::font_add_google, error = function(e) NULL)
OUT  <- "/Users/kuan-yuchu/Desktop/LOAD/Analysis/outputs"
DATA <- "/Users/kuan-yuchu/Desktop/LOAD/ADNI_csv"
FIG  <- "/Users/kuan-yuchu/Desktop/LOAD/figures/v4"
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

biom_long <- utils::read.csv(file.path(OUT, "analytic_biom_long.csv"), stringsAsFactors = FALSE)
attr_tbl  <- utils::read.csv(file.path(OUT, "analytic_baseline_attr.csv"), stringsAsFactors = FALSE)
tab4_jm   <- utils::read.csv(file.path(OUT, "table4_brms_joint_primary.csv"), stringsAsFactors = FALSE)
tab4_refit<- if (file.exists(file.path(OUT, "table4_brms_GFAP_refit.csv")))
  utils::read.csv(file.path(OUT, "table4_brms_GFAP_refit.csv"), stringsAsFactors = FALSE) else NULL
xa        <- utils::read.csv(file.path(OUT, "table_cross_assay.csv"), stringsAsFactors = FALSE)

# Merge refit GFAP rows into tab4_jm
if (!is.null(tab4_refit)) {
  tab4_jm <- tab4_jm[!grepl("GFAP$", tab4_jm$label), ]
  tab4_jm <- rbind(tab4_jm, tab4_refit)
}

# ---- Save helper (showtext_auto enables Unicode/Greek across PNG and PDF) ----
sv <- function(p, name, w = 6, h = 4.5) {
  ggsave(file.path(FIG, paste0(name, ".png")), p, dpi = 300, width = w, height = h, bg = "white")
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h)
}

# ============================================================
# Figure 1 — STROBE participant flow
# ============================================================
# Static rectangles via ggplot
flow_df <- data.frame(
  step = factor(c("ADNI plasma-biomarker registry\n(UPENN-Fujirebio Lumipulse, all)",
                  "Age ≥ 65 at baseline",
                  "Baseline dx ∈ {CN, MCI, AD}",
                  "APOE genotype available",
                  "≥ 2 GFAP and ≥ 2 p-tau217 measurements\n(SAP §4 revised threshold)",
                  "Analytic cohort\nN = 325 (CN 181 / MCI 106 / AD 38)"),
                 levels = rev(c("ADNI plasma-biomarker registry\n(UPENN-Fujirebio Lumipulse, all)",
                                "Age ≥ 65 at baseline",
                                "Baseline dx ∈ {CN, MCI, AD}",
                                "APOE genotype available",
                                "≥ 2 GFAP and ≥ 2 p-tau217 measurements\n(SAP §4 revised threshold)",
                                "Analytic cohort\nN = 325 (CN 181 / MCI 106 / AD 38)"))),
  n = c(1593, 910, 886, 886, 325, 325),
  is_final = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)
)
flow_df$y <- as.numeric(flow_df$step)
fig1 <- ggplot(flow_df, aes(x = 1, y = y)) +
  geom_tile(aes(fill = is_final), width = 0.9, height = 0.7, color = "black") +
  geom_text(aes(label = paste0(step, "\nN = ", n)), size = 3.3, lineheight = 1) +
  scale_fill_manual(values = c(`FALSE` = "grey95", `TRUE` = "#dceaf7")) +
  scale_x_continuous(limits = c(0.4, 1.6), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.3, 6.7), expand = c(0, 0)) +
  geom_segment(data = data.frame(y_from = 5:2, y_to = (5:2) - 0.65),
               aes(x = 1, xend = 1, y = y_from - 0.35, yend = y_to),
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
               inherit.aes = FALSE) +
  theme_void() + theme(legend.position = "none",
                       plot.title = element_text(face = "bold", hjust = 0.5)) +
  labs(title = "Figure 1. STROBE participant flow")
sv(fig1, "Figure1_STROBE_flow", w = 6, h = 7)

# ============================================================
# Figure 2 — Spaghetti plots, GFAP & p-tau217 over time, by APOE4
# ============================================================
biom_long$APOE4_lab <- factor(ifelse(biom_long$APOE4 == 1, "APOE ε4 carrier", "Non-carrier"),
                              levels = c("Non-carrier", "APOE ε4 carrier"))

# Marginal predictions: simple LME (no covariates) for plotting
library(lme4)
mk_fit <- function(y) lmer(as.formula(paste0(y, " ~ time_yrs*APOE4 + (1 + time_yrs | RID)")),
                           data = biom_long %>% filter(!is.na(.data[[y]])), REML = TRUE,
                           control = lmerControl(optimizer = "bobyqa"))
m_g <- mk_fit("log_gfap"); m_p <- mk_fit("log_ptau")
new_df <- expand.grid(time_yrs = seq(0, 10, 0.2),
                      APOE4 = c(0, 1))
new_df$pred_gfap <- predict(m_g, newdata = new_df, re.form = NA)
new_df$pred_ptau <- predict(m_p, newdata = new_df, re.form = NA)
new_df$APOE4_lab <- factor(ifelse(new_df$APOE4 == 1, "APOE ε4 carrier", "Non-carrier"),
                           levels = c("Non-carrier", "APOE ε4 carrier"))

plot_spag <- function(yvar, ylab, fitvar) {
  ggplot(biom_long %>% filter(!is.na(.data[[yvar]])),
         aes(x = time_yrs, y = .data[[yvar]], group = RID, color = APOE4_lab)) +
    geom_line(alpha = 0.15) +
    geom_line(data = new_df, aes(x = time_yrs, y = .data[[fitvar]],
                                 color = APOE4_lab, group = APOE4_lab),
              size = 1.4, inherit.aes = FALSE) +
    scale_color_manual(values = c("Non-carrier" = "#3a6ea5", "APOE ε4 carrier" = "#c0392b")) +
    labs(x = "Time from baseline (years)", y = ylab, color = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank())
}
g_panel <- plot_spag("log_gfap", "log plasma GFAP (pg/mL)", "pred_gfap") +
  ggtitle("A. Plasma GFAP")
p_panel <- plot_spag("log_ptau", "log plasma p-tau217 (pg/mL)", "pred_ptau") +
  ggtitle("B. Plasma p-tau217")

library(patchwork)
fig2 <- (g_panel | p_panel) +
  plot_annotation(title = "Figure 2. Individual longitudinal trajectories by APOE ε4 status",
                  subtitle = "Light lines = individual RIDs; heavy lines = APOE-stratified marginal LME fits.")
sv(fig2, "Figure2_spaghetti", w = 10, h = 4.5)

# ============================================================
# Figure 3 — Forest plot of brms joint-model ρ values
#   Panel A: Full cohort (GFAP + p-tau217 × 3 outcomes)
#   Panel B: p-tau217 sensitivities (Full / Amyloid+ / Post-baseline only) × 3 outcomes
# ============================================================
tab4_jm$Outcome <- factor(sub("_.*", "", tab4_jm$label),
                         levels = c("ADAS", "MMSE", "Hippo"),
                         labels = c("ADAS-Cog 13", "MMSE", "Hippocampal ratio"))
tab4_jm$Biomarker <- ifelse(grepl("GFAP", tab4_jm$label), "Plasma GFAP", "Plasma p-tau217")
tab4_jm$Biomarker <- factor(tab4_jm$Biomarker, levels = c("Plasma GFAP", "Plasma p-tau217"))

panelA <- ggplot(tab4_jm, aes(x = mean, y = Outcome, color = Biomarker)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, size = 0.7,
                 position = position_dodge(width = 0.5)) +
  geom_point(size = 3.2, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("Plasma GFAP" = "#7f8c8d", "Plasma p-tau217" = "#c0392b")) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(x = expression(paste("Cross-response random-slope correlation ", rho, " (95% CrI)")),
       y = NULL, color = NULL,
       title = "A. Full analytic cohort (N = 282–296 per pair)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

# Panel B: p-tau217 across sensitivity specifications
build_sens <- function() {
  full <- tab4_jm %>% filter(Biomarker == "Plasma p-tau217") %>%
    transmute(Outcome, Spec = "Full cohort",
              mean, lo, hi, n_id)
  apos_path <- file.path(OUT, "table4_brms_amyloid_pos.csv")
  pbl_path  <- file.path(OUT, "table4_brms_postbaseline.csv")
  rows <- list(full)
  if (file.exists(apos_path)) {
    a <- utils::read.csv(apos_path, stringsAsFactors = FALSE)
    a$Outcome <- factor(c("ADAS-Cog 13","MMSE","Hippocampal ratio"),
                       levels = c("ADAS-Cog 13", "MMSE", "Hippocampal ratio"))
    rows[[length(rows)+1]] <- a %>% transmute(Outcome, Spec = "Amyloid-positive subset",
                                              mean, lo, hi, n_id)
  }
  if (file.exists(pbl_path)) {
    p <- utils::read.csv(pbl_path, stringsAsFactors = FALSE)
    p$Outcome <- factor(c("ADAS-Cog 13","MMSE","Hippocampal ratio"),
                       levels = c("ADAS-Cog 13", "MMSE", "Hippocampal ratio"))
    rows[[length(rows)+1]] <- p %>% transmute(Outcome, Spec = "Post-baseline only",
                                              mean, lo, hi, n_id)
  }
  bind_rows(rows)
}
sens_df <- build_sens()
sens_df$Spec <- factor(sens_df$Spec,
                       levels = c("Full cohort", "Amyloid-positive subset", "Post-baseline only"))

panelB <- ggplot(sens_df, aes(x = mean, y = Outcome, color = Spec)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.20, size = 0.7,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 3.2, position = position_dodge(width = 0.6)) +
  scale_color_manual(values = c("Full cohort" = "#c0392b",
                                "Amyloid-positive subset" = "#2980b9",
                                "Post-baseline only" = "#27ae60")) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  labs(x = expression(paste("Cross-response random-slope correlation ", rho, " (95% CrI)")),
       y = NULL, color = NULL,
       title = "B. Plasma p-tau217: sensitivity to specification") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

library(patchwork)
fig3 <- (panelA / panelB) + plot_annotation(
  title = "Figure 3. Joint-model biomarker–outcome co-trajectory correlations",
  subtitle = "ρ from brms multivariate joint LME. Panel A: full-cohort primary. Panel B: p-tau217 sensitivities (Amyloid+ stratification, Post-baseline only).")
sv(fig3, "Figure3_brms_forest", w = 8.5, h = 7)

# ============================================================
# Figure 4 — Cross-assay p-tau217 scatter at baseline
# ============================================================
# Reload assay data + baseline UPENN ptau
analytic_rid <- attr_tbl$RID
ptau_bl <- biom_long %>% group_by(RID) %>%
  filter(time_yrs == min(time_yrs, na.rm = TRUE)) %>% ungroup() %>%
  transmute(RID = as.integer(RID), pT217_F = exp(log_ptau))

c2n  <- utils::read.csv(file.path(DATA, "C2N_PRECIVITYAD2_PLASMA_24May2026.csv"),
                       stringsAsFactors = FALSE, na.strings = c("","NA"))
lilly <- utils::read.csv(file.path(DATA, "LILLY_PTAU217_MSD600_24May2026.csv"),
                         stringsAsFactors = FALSE, na.strings = c("","NA"))
jans  <- utils::read.csv(file.path(DATA, "JANSSEN_PLASMA_P217_TAU_24May2026.csv"),
                         stringsAsFactors = FALSE, na.strings = c("","NA"))

prep_assay <- function(d, col, lab) {
  d %>% mutate(RID = as.integer(RID), val = suppressWarnings(as.numeric(.data[[col]]))) %>%
    filter(!is.na(val), val > 0, RID %in% analytic_rid) %>%
    group_by(RID) %>% slice_min(as.Date(EXAMDATE), with_ties = FALSE) %>% ungroup() %>%
    transmute(RID, val, assay = lab)
}
c2n_d <- prep_assay(c2n, "pT217_C2N", "C2N PrecivityAD2 (n = 208)")
ll_d  <- prep_assay(lilly %>% filter(TESTCD == "PTAU217"), "ORRES", "Lilly MSD600 (n = 23)")
jn_d  <- prep_assay(jans, "DILUTION_CORRECTED_CONC", "Janssen (n = 14)")

xa_df <- bind_rows(c2n_d, ll_d, jn_d) %>% inner_join(ptau_bl, by = "RID")

fig4 <- ggplot(xa_df, aes(x = pT217_F, y = val)) +
  geom_point(alpha = 0.6, color = "#34495e") +
  geom_smooth(method = "lm", se = TRUE, color = "#c0392b", fill = "#fadbd8", size = 0.7) +
  facet_wrap(~ assay, scales = "free_y") +
  labs(x = "UPENN-Fujirebio Lumipulse plasma p-tau217 (pg/mL), baseline",
       y = "Other-assay value",
       title = "Figure 4. Cross-assay concordance for plasma p-tau217 at baseline",
       subtitle = "Restricted to analytic cohort (N = 325).") +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))
sv(fig4, "Figure4_crossassay", w = 11, h = 4)

cat("Figures saved to:", FIG, "\n")
list.files(FIG)
