# Sensitivity: re-fit Phase II + III LME with random-INTERCEPT only
# (drops the random slope; resolves the singular boundary warning).
suppressPackageStartupMessages({ library(lme4); library(dplyr) })

CSV_DIR <- "/Users/kuan-yuchu/Desktop/LOAD/ADNI_csv"
df_long <- readRDS(file.path(CSV_DIR, "df_long_real.rds"))
df_base <- readRDS(file.path(CSV_DIR, "df_base_real.rds"))

df_long$diag <- factor(df_long$diag, levels = c("CN","MCI","AD"))
df_base$diag <- factor(df_base$diag, levels = c("CN","MCI","AD"))

fit_lme_intercept <- function(yvar, formula_rhs, data, label) {
  cat("\n---", label, "---\n")
  d <- data[complete.cases(data[, c(yvar, all.vars(as.formula(paste("~", formula_rhs))))]), ]
  f <- as.formula(paste(yvar, "~", formula_rhs, "+ (1 | RID)"))
  fit <- lmer(f, data = d, REML = FALSE,
              control = lmerControl(optimizer = "bobyqa"))
  co <- summary(fit)$coefficients
  print(round(co, 4))
  invisible(fit)
}

cat("\n=== PHASE II sensitivity (random intercept only) ===\n")
g <- fit_lme_intercept("logGFAP",
  "time_yr + apoe_e4 + apoe_e4:time_yr + age_baseline + female + educ + diag",
  df_long, "log GFAP")
p <- fit_lme_intercept("logPTAU",
  "time_yr + apoe_e4 + apoe_e4:time_yr + age_baseline + female + educ + diag",
  df_long, "log p-tau217")

# Re-fit random-slope LMEs to extract per-person slope BLUPs (used unchanged)
quick_blup <- function(yvar) {
  d <- df_long[!is.na(df_long[[yvar]]), ]
  f <- as.formula(paste(yvar,
    "~ time_yr + apoe_e4 + apoe_e4:time_yr + age_baseline + female + educ + diag + (time_yr | RID)"))
  fit <- lmer(f, data = d, REML = FALSE,
              control = lmerControl(optimizer = "bobyqa"))
  re <- ranef(fit)$RID; fe <- fixef(fit)
  data.frame(RID = as.integer(rownames(re)),
             slope = as.numeric(fe["time_yr"] + re[["time_yr"]]))
}
gs <- quick_blup("logGFAP"); ps <- quick_blup("logPTAU")
df_long <- df_long %>%
  select(-any_of(c("gfap_slope","ptau_slope"))) %>%
  left_join(gs %>% rename(gfap_slope = slope), by = "RID") %>%
  left_join(ps %>% rename(ptau_slope = slope), by = "RID")

cat("\n=== PHASE III sensitivity (random intercept only) ===\n")
for (out in c("ADAS13","MMSE","HV_adj")) {
  for (sv in c("gfap_slope","ptau_slope")) {
    fit_lme_intercept(out,
      paste("time_yr *", sv, "+ age_baseline + female + educ + diag + apoe_e4"),
      df_long, paste(out, "~ time × ", sv))
  }
}

cat("\n=== END SENSITIVITY ===\n")
