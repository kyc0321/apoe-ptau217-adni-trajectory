# Figures for LOAD manuscript v4
suppressPackageStartupMessages({
  library(lme4); library(dplyr); library(ggplot2); library(tidyr)
})

CSV_DIR <- "/Users/kuan-yuchu/Desktop/LOAD/ADNI_csv"
FIG_DIR <- "/Users/kuan-yuchu/Desktop/LOAD/figures"
dir.create(FIG_DIR, showWarnings = FALSE)

df_long <- readRDS(file.path(CSV_DIR, "df_long_real.rds"))
df_base <- readRDS(file.path(CSV_DIR, "df_base_real.rds"))
df_long$diag <- factor(df_long$diag, levels = c("CN","MCI","AD"))
df_base$diag <- factor(df_base$diag, levels = c("CN","MCI","AD"))

# Re-fit biomarker LMEs to extract BLUPs
fit_b <- function(y) {
  d <- df_long[!is.na(df_long[[y]]) & !is.na(df_long$apoe_e4), ]
  f <- as.formula(paste(y, "~ time_yr + apoe_e4 + apoe_e4:time_yr +",
                        "age_baseline + female + educ + diag + (time_yr | RID)"))
  lmer(f, data = d, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
}
m_gfap <- fit_b("logGFAP"); m_ptau <- fit_b("logPTAU")
slope <- function(fit, name) {
  re <- ranef(fit)$RID; fe <- fixef(fit)
  data.frame(RID = as.integer(rownames(re)),
             slope = as.numeric(fe["time_yr"] + re[["time_yr"]]))
}
gs <- slope(m_gfap); ps <- slope(m_ptau)
df_base <- df_base %>% select(-any_of(c("gfap_slope","ptau_slope"))) %>%
  left_join(gs %>% rename(gfap_slope = slope), by = "RID") %>%
  left_join(ps %>% rename(ptau_slope = slope), by = "RID")

# Outcome rate BLUPs (random slope from unconditional LME)
out_slope <- function(y) {
  d <- df_long[!is.na(df_long[[y]]), ]
  f <- as.formula(paste(y, "~ time_yr + age_baseline + female + educ + diag + (time_yr | RID)"))
  fit <- tryCatch(lmer(f, data = d, REML = FALSE,
                       control = lmerControl(optimizer = "bobyqa")),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  re <- ranef(fit)$RID; fe <- fixef(fit)
  if (!"time_yr" %in% colnames(re)) {
    return(data.frame(RID = as.integer(rownames(re)),
                      rate = rep(as.numeric(fe["time_yr"]), nrow(re))))
  }
  data.frame(RID = as.integer(rownames(re)),
             rate = as.numeric(fe["time_yr"] + re[["time_yr"]]))
}
adas_r <- out_slope("ADAS13"); mmse_r <- out_slope("MMSE"); hv_r <- out_slope("HV_adj")
df_base <- df_base %>%
  left_join(adas_r %>% rename(adas_rate = rate), by = "RID") %>%
  left_join(mmse_r %>% rename(mmse_rate = rate), by = "RID") %>%
  left_join(hv_r   %>% rename(hv_rate   = rate), by = "RID")

theme_pub <- theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom",
        panel.grid.minor = element_blank())

save_tiff <- function(p, fn, w, h) {
  path <- file.path(FIG_DIR, fn)
  ggsave(path, p, width = w, height = h, dpi = 300,
         device = "tiff", compression = "lzw")
  message("Saved: ", path)
}

# ---- FIGURE 2: spaghetti GFAP + p-tau by APOE ----
df_long$apoe_lab <- factor(ifelse(df_long$apoe_e4 == 1, "APOE ε4 carrier",
                                  "Non-carrier"),
                           levels = c("Non-carrier", "APOE ε4 carrier"))
mk_spag <- function(yvar, lab, ylim_pad = NULL) {
  d <- df_long[!is.na(df_long[[yvar]]) & !is.na(df_long$apoe_lab), ]
  yvals <- d[[yvar]]
  ggplot(d, aes(time_yr, .data[[yvar]], group = RID, color = apoe_lab)) +
    geom_line(alpha = 0.30, linewidth = 0.35) +
    geom_smooth(aes(group = apoe_lab), method = "lm", se = TRUE,
                linewidth = 1.0, alpha = 0.20) +
    scale_color_manual(values = c("Non-carrier" = "#2C7BB6",
                                  "APOE ε4 carrier" = "#D7191C")) +
    labs(x = "Years from baseline", y = lab,
         color = NULL) +
    theme_pub
}
p2a <- mk_spag("logGFAP", "log plasma GFAP (pg/mL)") +
  ggtitle("A. Plasma GFAP")
p2b <- mk_spag("logPTAU", "log plasma p-tau217 (pg/mL)") +
  ggtitle("B. Plasma p-tau217 (UPENN-Fujirebio)")
fig2 <- gridExtra::grid.arrange(p2a, p2b, ncol = 2)
save_tiff(fig2, "Figure2_spaghetti_GFAP_pTau.tiff", 11, 4.5)

# ---- FIGURE 3: slope BLUPs vs outcome rates (6 panels) ----
panel <- function(xv, yv, xlab, ylab, title) {
  d <- df_base[!is.na(df_base[[xv]]) & !is.na(df_base[[yv]]), ]
  d$apoe_lab <- factor(ifelse(d$apoe_e4 == 1, "ε4 carrier", "Non-carrier"),
                       levels = c("Non-carrier", "ε4 carrier"))
  r <- suppressWarnings(cor(d[[xv]], d[[yv]], use = "complete.obs"))
  ggplot(d, aes(.data[[xv]], .data[[yv]], color = apoe_lab)) +
    geom_smooth(method = "lm", se = TRUE, color = "black",
                fill = "gray80", linewidth = 0.7) +
    geom_point(alpha = 0.7, size = 1.8) +
    scale_color_manual(values = c("Non-carrier" = "#2C7BB6",
                                  "ε4 carrier" = "#D7191C")) +
    labs(x = xlab, y = ylab, title = title,
         subtitle = sprintf("Pearson r = %.2f", r),
         color = NULL) +
    theme_pub + theme(legend.position = "none")
}
p3a <- panel("gfap_slope","adas_rate","GFAP slope (log-units/yr)",
             "ADAS-Cog 13 rate (pts/yr)", "A. GFAP slope vs ADAS-Cog 13")
p3b <- panel("gfap_slope","mmse_rate","GFAP slope (log-units/yr)",
             "MMSE rate (pts/yr)",         "B. GFAP slope vs MMSE")
p3c <- panel("gfap_slope","hv_rate","GFAP slope (log-units/yr)",
             "Hippocampal volume rate (mm³/yr)", "C. GFAP slope vs HV")
p3d <- panel("ptau_slope","adas_rate","p-tau217 slope (log-units/yr)",
             "ADAS-Cog 13 rate (pts/yr)", "D. p-tau217 slope vs ADAS-Cog 13")
p3e <- panel("ptau_slope","mmse_rate","p-tau217 slope (log-units/yr)",
             "MMSE rate (pts/yr)",         "E. p-tau217 slope vs MMSE")
p3f <- panel("ptau_slope","hv_rate","p-tau217 slope (log-units/yr)",
             "Hippocampal volume rate (mm³/yr)", "F. p-tau217 slope vs HV")
fig3 <- gridExtra::grid.arrange(p3a, p3b, p3c, p3d, p3e, p3f, ncol = 3)
save_tiff(fig3, "Figure3_slope_BLUPs_vs_outcomes.tiff", 13, 8)

# ---- FIGURE 4: cross-assay p-tau217 concordance at baseline ----
cb <- read.csv(file.path(CSV_DIR, "ptau217_cross_assay_baseline.csv"),
               check.names = FALSE)
# Drop ADNI sentinel values (e.g. -4.0 = below LOQ / not collected). p-tau is positive.
for (col in c("UPENN_FUJIREBIO","LILLY_MSD","JANSSEN","C2N"))
  if (col %in% names(cb)) cb[[col]][cb[[col]] <= 0] <- NA
mk_cross <- function(x, y, xl, yl, title) {
  d <- cb[!is.na(cb[[x]]) & !is.na(cb[[y]]), ]
  r <- cor(d[[x]], d[[y]], use = "complete.obs")
  ggplot(d, aes(.data[[x]], .data[[y]])) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
    geom_smooth(method = "lm", se = TRUE, color = "#D7191C",
                fill = "#fde0dd", linewidth = 0.7) +
    geom_point(alpha = 0.6, size = 1.8, color = "#2C7BB6") +
    labs(x = xl, y = yl, title = title,
         subtitle = sprintf("Pearson r = %.2f (n = %d)", r, nrow(d))) +
    theme_pub
}
p4a <- mk_cross("UPENN_FUJIREBIO","C2N",
                "UPENN-Fujirebio (pg/mL)", "C2N PrecivityAD2 (pg/mL)",
                "A. UPENN-Fujirebio vs C2N")
p4b <- mk_cross("UPENN_FUJIREBIO","LILLY_MSD",
                "UPENN-Fujirebio (pg/mL)", "Lilly MSD600 (pg/mL)",
                "B. UPENN-Fujirebio vs Lilly MSD")
fig4 <- gridExtra::grid.arrange(p4a, p4b, ncol = 2)
save_tiff(fig4, "Figure4_ptau217_cross_assay.tiff", 10, 4.5)

cat("\nAll figures saved to:", FIG_DIR, "\n")
