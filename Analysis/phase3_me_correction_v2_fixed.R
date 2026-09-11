# Phase III — Measurement-error correction via Multiple Imputation of BLUP slopes
# v2: FIXES a bug in phase3_me_correction_v1.R caught via Neurological Sciences
# Reviewer #1 (decision 2026-09-10), comment on the mediation analysis.
#
# THE BUG (v1, line ~32): the model formula includes a time_yrs*APOE4 fixed-effect
# interaction, but the BLUP slope was reconstructed as
#     slope_mean = ranef(m)$RID[,"time_yrs"] + fixef(m)["time_yrs"]
# i.e. random deviation + the FIXED MAIN EFFECT ONLY. ranef() returns each
# person's random deviation net of their FULL fixed-effects prediction, which
# for an APOE4 carrier includes fixef(m)["time_yrs:APOE4"] as well as the main
# effect. Omitting that interaction term when adding the fixed portion back in
# means every APOE4 carrier's reconstructed slope is missing their own
# APOE4-associated acceleration (~0.03/year in the primary model, non-trivial
# over a multi-year follow-up). This is the most likely reason the downstream
# mediation analysis (mediation_analysis_v1.R) found an ACME/total effect near
# zero despite a significant APOE4×time interaction in the primary LME: the
# outcome variable it was handed had already had most of the APOE4 signal
# silently stripped out during BLUP reconstruction, before mediation() ever
# saw it.
#
# THE FIX: add APOE4_i * fixef(time_yrs:APOE4) back in per person.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(lme4); library(lmerTest)
})
set.seed(20260525)

# NOTE: v1 used stale hardcoded paths ("/Users/kuan-yuchu/Desktop/LOAD/...")
# left over from before the project folder was renamed to "#2 LOAD-BC-NS" --
# fixed here to the current, real location.
BASE <- "/Users/kuan-yuchu/Desktop/#2 LOAD-BC-NS"
DATA <- file.path(BASE, "ADNI_csv")
OUT  <- file.path(BASE, "Analysis/outputs")

biom_long <- utils::read.csv(file.path(OUT, "analytic_biom_long.csv"), stringsAsFactors = FALSE)
biom_long$dx3 <- factor(biom_long$dx3, levels = c("CN","MCI","AD"))

phase2 <- function(yvar) {
  d <- biom_long %>% filter(!is.na(.data[[yvar]]), !is.na(BMI_z))
  fml <- as.formula(paste0(yvar,
    " ~ time_yrs*APOE4 + age_z + sex + edu_z + dx3 + BMI_z + HTN_base + DM_base + (1 + time_yrs | RID)"))
  m <- lmerTest::lmer(fml, data = d, REML = TRUE,
                      control = lmerControl(optimizer = "bobyqa"))

  re <- ranef(m, condVar = TRUE)$RID
  cv <- attr(re, "postVar")
  fe <- fixef(m)
  fe_time <- fe["time_yrs"]
  fe_int  <- fe["time_yrs:APOE4"]

  ridlist <- as.integer(rownames(re))
  # per-person APOE4 status, one row per RID (should be time-invariant, but
  # take the first observed value defensively rather than assume)
  apoe_by_rid <- d %>% distinct(RID, .keep_all = TRUE) %>% select(RID, APOE4)
  apoe_vec <- apoe_by_rid$APOE4[match(ridlist, apoe_by_rid$RID)]
  if (anyNA(apoe_vec)) stop("APOE4 status missing for some RIDs in ranef() output -- investigate before proceeding")

  data.frame(
    RID = ridlist,
    # THE FIX is this line: + apoe_vec * fe_int added back in, vs. v1 which omitted it
    slope_mean = re[, "time_yrs"] + fe_time + apoe_vec * fe_int,
    slope_se   = sqrt(cv[2, 2, ])
  )
}

cat("Refitting Phase II for GFAP...\n")
blup_gfap <- phase2("log_gfap"); names(blup_gfap)[2:3] <- c("slope_gfap","se_gfap")
cat("Refitting Phase II for p-tau217...\n")
blup_ptau <- phase2("log_ptau"); names(blup_ptau)[2:3] <- c("slope_ptau","se_ptau")
blups <- full_join(blup_gfap, blup_ptau, by = "RID")

utils::write.csv(blups, file.path(OUT, "blup_slopes_with_SE_v2_fixed.csv"), row.names = FALSE)
cat("Saved BLUP with conditional SE (v2, fixed) -> blup_slopes_with_SE_v2_fixed.csv\n")

# ---- quick before/after comparison, split by APOE4 status ----
old <- utils::read.csv(file.path(OUT, "blup_slopes_with_SE.csv"))
cmp <- blups %>% left_join(old %>% select(RID, slope_ptau_old = slope_ptau), by = "RID") %>%
  left_join(biom_long %>% distinct(RID, .keep_all = TRUE) %>% select(RID, APOE4), by = "RID")
cat("\n=== Mean p-tau217 slope, v1 (buggy) vs v2 (fixed), by APOE4 ===\n")
print(cmp %>% group_by(APOE4) %>%
        summarise(n = n(), mean_v1 = mean(slope_ptau_old, na.rm = TRUE),
                  mean_v2_fixed = mean(slope_ptau, na.rm = TRUE),
                  mean_diff = mean(slope_ptau - slope_ptau_old, na.rm = TRUE)))
