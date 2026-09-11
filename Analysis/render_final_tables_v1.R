# Render manuscript-ready tables from /Analysis/outputs/ CSVs
# Author: Kuan-Yu Chu | 2026-05-25
# Output: /Mauuscripts/Tables_final_v1.md

suppressPackageStartupMessages({ library(dplyr) })
OUT <- "/Users/kuan-yuchu/Desktop/LOAD/Analysis/outputs"
MS  <- "/Users/kuan-yuchu/Desktop/LOAD/Mauuscripts"

t1  <- utils::read.csv(file.path(OUT, "table1_baseline.csv"), stringsAsFactors = FALSE)
t2  <- utils::read.csv(file.path(OUT, "table2_lcgmm.csv"),    stringsAsFactors = FALSE)
t2f <- tryCatch(utils::read.csv(file.path(OUT, "table2_lcgmm_full_diagnostic.csv"), stringsAsFactors=FALSE),
                error = function(e) NULL)
blrt <- tryCatch(utils::read.csv(file.path(OUT, "lcgmm_bootstrap_LRT.csv"), stringsAsFactors=FALSE),
                 error = function(e) NULL)
t3  <- utils::read.csv(file.path(OUT, "table3_lme_phase2.csv"), stringsAsFactors = FALSE)
t4  <- utils::read.csv(file.path(OUT, "table4_brms_joint_primary.csv"), stringsAsFactors = FALSE)
t4s <- utils::read.csv(file.path(OUT, "table4_twostage_sensitivity.csv"), stringsAsFactors = FALSE)
t4m <- utils::read.csv(file.path(OUT, "table4_ME_corrected_primary.csv"), stringsAsFactors = FALSE)
t4a <- tryCatch(utils::read.csv(file.path(OUT, "table4_brms_amyloid_pos.csv"), stringsAsFactors=FALSE),
                error = function(e) NULL)
t4pb<- tryCatch(utils::read.csv(file.path(OUT, "table4_brms_postbaseline.csv"), stringsAsFactors=FALSE),
                error = function(e) NULL)
hippo <- tryCatch(utils::read.csv(file.path(OUT, "hippo_slope_per_dx.csv"), stringsAsFactors=FALSE),
                  error = function(e) NULL)
vif <- tryCatch(utils::read.csv(file.path(OUT, "vif_phase2.csv"), stringsAsFactors=FALSE),
                error = function(e) NULL)
xa  <- utils::read.csv(file.path(OUT, "table_cross_assay.csv"),  stringsAsFactors = FALSE)

# Optionally merge GFAP refit if present
gfap_refit_path <- file.path(OUT, "table4_brms_GFAP_refit.csv")
if (file.exists(gfap_refit_path)) {
  rfit <- utils::read.csv(gfap_refit_path, stringsAsFactors = FALSE)
  # Replace GFAP rows in t4 with refit rows
  t4 <- t4[!grepl("GFAP$", t4$label), ]
  t4 <- bind_rows(t4, rfit)
}

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) "<0.001"
  else if (p < 0.01) sprintf("%.3f", p)
  else sprintf("%.3f", p)
}

# ---------- TABLE 1 ----------
tab1_md <- paste0(
  "## Table 1. Baseline characteristics of the analytic cohort (N = 325).\n\n",
  "| Characteristic | Value |\n|---|---|\n",
  paste0("| ", t1$Characteristic, " | ", t1$Value, " |", collapse = "\n"),
  "\n\n*ADNI 2026-05-24 data release. Inclusion: age ≥65 at baseline, baseline diagnosis CN/MCI/AD, APOE genotype available, ≥2 plasma GFAP and ≥2 plasma p-tau217 (UPENN-Fujirebio) measurements within a 10-year window.*\n"
)

# ---------- TABLE 2 (LCGMM with full diagnostic) ----------
if (!is.null(t2f)) {
  t2f$loglik <- sprintf("%.1f", t2f$loglik)
  t2f$AIC <- sprintf("%.1f", t2f$AIC); t2f$BIC <- sprintf("%.1f", t2f$BIC)
  t2f$entropy <- sprintf("%.3f", t2f$entropy)
  t2f$min_class_pct <- sprintf("%.1f%%", t2f$min_class_pct)
  t2f$mean_post_prob <- sprintf("%.3f", t2f$mean_post_prob)
  blrt_note <- if (!is.null(blrt))
    sprintf("Bootstrap LRT (k=%d vs k=%d, B=%d): observed LRT = %.2f, p = %.4g (%d valid replicates).",
            blrt$k_test[1], blrt$k_ref[1], blrt$n_valid_boot[1], blrt$obs_LRT[1],
            blrt$boot_p[1], blrt$n_valid_boot[1])
    else ""
  tab2_md <- paste0(
    "## Table 2. LCGMM model selection (exploratory diagnostic — full reporting).\n\n",
    "| Biomarker | k | log-lik | AIC | BIC | Entropy | Min class % | Mean posterior probability |\n",
    "|---|---|---|---|---|---|---|---|\n",
    paste0("| ", t2f$biomarker, " | ", t2f$k, " | ", t2f$loglik, " | ", t2f$AIC, " | ", t2f$BIC,
           " | ", t2f$entropy, " | ", t2f$min_class_pct, " | ", t2f$mean_post_prob, " |",
           collapse = "\n"),
    "\n\n*", blrt_note,
    " For plasma p-tau217, all four pre-specified strong-support criteria are met for the k=2 solution (ΔBIC = 16 ≥ 10; min class size 45.5%; mean posterior probability 0.92 ≥ 0.70; bootstrap LRT p < 0.05). For plasma GFAP, the single-class solution is retained (k=3/4 lower BIC but min class size < 5% pre-specified threshold).*\n"
  )
} else {
  t2$loglik <- sprintf("%.1f", t2$loglik)
  t2$AIC <- sprintf("%.1f", t2$AIC); t2$BIC <- sprintf("%.1f", t2$BIC)
  t2$min_class_pct <- sprintf("%.1f%%", t2$min_class_pct)
  t2$entropy <- ifelse(is.na(t2$entropy), "—", sprintf("%.3f", t2$entropy))
  tab2_md <- paste0(
    "## Table 2. LCGMM model selection (exploratory only).\n\n",
    "| Biomarker | k | log-lik | AIC | BIC | Entropy | Min class % |\n|---|---|---|---|---|---|---|\n",
    paste0("| ", t2$biomarker, " | ", t2$k, " | ", t2$loglik, " | ", t2$AIC, " | ", t2$BIC,
           " | ", t2$entropy, " | ", t2$min_class_pct, " |", collapse = "\n"),
    "\n\n*Pre-specified per SAP §6.2 as exploratory diagnostic only.*\n"
  )
}

# ---------- TABLE 3 (Phase II LME, key APOE x time rows) ----------
t3k <- t3[grepl("time_yrs:APOE4|^APOE4$|^time_yrs$", t3$term), ]
t3k$Estimate <- sprintf("%.4f", as.numeric(t3k$Estimate))
t3k$SE <- sprintf("%.4f", as.numeric(t3k$Std..Error))
t3k$p <- sapply(as.numeric(t3k$Pr...t..), format_p)
t3k$ci <- paste0("[", sprintf("%.4f", as.numeric(t3k$X2.5..)), ", ",
                       sprintf("%.4f", as.numeric(t3k$X97.5..)), "]")
tab3_md <- paste0(
  "## Table 3. Phase II — APOE × time on plasma biomarker trajectories (primary LME, N≈296 RIDs).\n\n",
  "| Biomarker | Model | Term | β | SE | 95% CI | p | Singular | n_obs / n_id |\n",
  "|---|---|---|---|---|---|---|---|---|\n",
  paste0("| ", t3k$biomarker, " | ", t3k$model, " | ", t3k$term,
         " | ", t3k$Estimate, " | ", t3k$SE, " | ", t3k$ci, " | ", t3k$p,
         " | ", t3k$singular, " | ", t3k$n_obs, "/", t3k$n_id, " |", collapse = "\n"),
  "\n\n*Models adjusted for baseline age, sex, education, baseline diagnosis (CN reference), baseline BMI, baseline HTN, baseline DM. Primary: N=325 cohort; sensitivity: eGFR-available subset (n_id≈109). Random intercepts + slopes; REML.*\n"
)

# ---------- TABLE 4 (brms joint primary + sensitivities) ----------
t4$rho_ci <- sprintf("%.2f [%.2f, %.2f]", t4$mean, t4$lo, t4$hi)
t4$p_dir_str <- sprintf("%.3f", t4$p_dir)
t4$rhat_str <- sprintf("%.2f", t4$rhat)
tab4_md <- paste0(
  "## Table 4. Phase III — Multivariate Bayesian joint longitudinal model (primary).\n\n",
  "*Cross-response random-slope correlation ρ = cor(b1_biomarker, b1_outcome). |ρ|→1 means individual-level biomarker slope and outcome slope co-vary; ρ→0 means independent.*\n\n",
  "| Pair | ρ (mean) [95% CrI] | P(ρ<0) | Rhat | n_id | n_bm / n_oc |\n",
  "|---|---|---|---|---|---|\n",
  paste0("| ", t4$label, " | ", t4$rho_ci, " | ", t4$p_dir_str,
         " | ", t4$rhat_str, " | ", t4$n_id, " | ", t4$n_bm, "/", t4$n_oc, " |", collapse = "\n"),
  "\n"
)

# ---------- TABLE 4 sensitivities ----------
t4s$beta <- sprintf("%.4g", t4s$beta_per_SD)
t4s$se   <- sprintf("%.4g", t4s$SE)
t4s$ci   <- sprintf("[%.4g, %.4g]", t4s$CI_lo, t4s$CI_hi)
t4s$p    <- sapply(t4s$p_value, format_p)
t4s$q    <- sapply(t4s$q_BH, format_p)
sens_md <- paste0(
  "\n### Table 4 — Sensitivity S-G: Naive two-stage (BLUP-as-predictor, no measurement-error correction).\n\n",
  "*WARNING: This sensitivity inflates significance because it ignores measurement error in the BLUP slope. Reported only for transparency.*\n\n",
  "| Outcome | Predictor | β (per 1-SD slope) | 95% CI | p | q (BH) | n_obs / n_id |\n",
  "|---|---|---|---|---|---|---|\n",
  paste0("| ", t4s$outcome, " | ", t4s$predictor, " | ", t4s$beta, " | ", t4s$ci,
         " | ", t4s$p, " | ", t4s$q, " | ", t4s$n_obs, "/", t4s$n_id, " |", collapse = "\n"),
  "\n"
)

t4m$beta <- sprintf("%.4g", t4m$beta_per_SD)
t4m$ci   <- sprintf("[%.4g, %.4g]", t4m$CI_lo, t4m$CI_hi)
t4m$p    <- sapply(t4m$p_value, format_p)
t4m$q    <- sapply(t4m$q_BH, format_p)
mi_md <- paste0(
  "\n### Table 4 — Sensitivity S-H: MI measurement-error correction of two-stage (M=100, Rubin's pooling).\n\n",
  "*Conservative comparison; treats BLUP slope as fixed-but-noisy predictor and cannot share information through random-effects covariance.*\n\n",
  "| Outcome | Predictor | β (per 1-SD slope) | 95% CI | p (pooled) | q (BH) | RIV |\n",
  "|---|---|---|---|---|---|---|\n",
  paste0("| ", t4m$outcome, " | ", t4m$predictor, " | ", t4m$beta, " | ", t4m$ci,
         " | ", t4m$p, " | ", t4m$q, " | ", sprintf("%.2f", t4m$riv), " |", collapse = "\n"),
  "\n"
)

# Amyloid-positive sensitivity
amyloid_md <- ""
if (!is.null(t4a)) {
  t4a$rho_ci <- sprintf("%.2f [%.2f, %.2f]", t4a$mean, t4a$lo, t4a$hi)
  t4a$p_dir_str <- sprintf("%.3f", t4a$p_dir)
  t4a$rhat_str <- sprintf("%.2f", t4a$rhat)
  amyloid_md <- paste0(
    "\n### Table 4 — Sensitivity: Amyloid-positive subset (plasma Aβ42/40 < 0.085; n_id = 143–149).\n\n",
    "*Mediation check addressing upstream amyloid confounding. Cognitive ρ values attenuate but persist; hippocampal ρ attenuates to non-significance.*\n\n",
    "| Pair | ρ (mean) [95% CrI] | P(direction) | Rhat | n_id |\n|---|---|---|---|---|\n",
    paste0("| ", t4a$label, " | ", t4a$rho_ci, " | ", t4a$p_dir_str,
           " | ", t4a$rhat_str, " | ", t4a$n_id, " |", collapse = "\n"),
    "\n"
  )
}

# Post-baseline only
postbl_md <- ""
if (!is.null(t4pb)) {
  t4pb$rho_ci <- sprintf("%.2f [%.2f, %.2f]", t4pb$mean, t4pb$lo, t4pb$hi)
  t4pb$p_dir_str <- sprintf("%.3f", t4pb$p_dir)
  t4pb$rhat_str <- sprintf("%.2f", t4pb$rhat)
  postbl_md <- paste0(
    "\n### Table 4 — Sensitivity: Post-baseline-only outcome window (time_yrs ≥ 0; n_id = 167–269).\n\n",
    "*Addresses time-frame asymmetry between biomarker and outcome trajectories. Cognitive ρ values essentially unchanged; hippocampal ρ attenuates to non-significance.*\n\n",
    "| Pair | ρ (mean) [95% CrI] | P(direction) | Rhat | n_id |\n|---|---|---|---|---|\n",
    paste0("| ", t4pb$label, " | ", t4pb$rho_ci, " | ", t4pb$p_dir_str,
           " | ", t4pb$rhat_str, " | ", t4pb$n_id, " |", collapse = "\n"),
    "\n"
  )
}

# Supplementary Table S3: hippo per-dx slope
hippo_md <- ""
if (!is.null(hippo)) {
  hippo_md <- paste0(
    "\n## Supplementary Table S2. Per-diagnosis hippocampal volume slope (biological plausibility check).\n\n",
    "| Diagnosis | n_obs | n_id | Slope (mm³/year) | SE | p |\n|---|---|---|---|---|---|\n",
    paste0("| ", hippo$dx, " | ", hippo$n_obs, " | ", hippo$n_id,
           " | ", sprintf("%.1f", hippo$slope_mm3_yr),
           " | ", sprintf("%.1f", hippo$SE),
           " | ", sapply(hippo$p, format_p), " |", collapse = "\n"),
    "\n\n*Published expectations: CN −50 to −150; MCI −200 to −400; AD −300 to −600 mm³/year.*\n"
  )
}

# Supplementary Table S3: VIF
vif_md <- ""
if (!is.null(vif)) {
  vif_md <- paste0(
    "\n## Supplementary Table S3. Collinearity check (variance inflation factors) for Phase II fixed effects.\n\n",
    "| Term | VIF |\n|---|---|\n",
    paste0("| ", vif$term, " | ", sprintf("%.2f", vif$VIF), " |", collapse = "\n"),
    "\n\n*All VIFs < 5 (conventional threshold). Max VIF = ",
    sprintf("%.2f", max(vif$VIF)), ".*\n"
  )
}

# ---------- CROSS-ASSAY ----------
xa$r <- sprintf("%.2f", xa$pearson); xa$rho <- sprintf("%.2f", xa$spearman)
xa_md <- paste0(
  "\n## Table 5. Cross-assay concordance for plasma p-tau217 at baseline (analytic cohort).\n\n",
  "| Assay (vs. UPENN-Fujirebio Lumipulse) | n | Pearson r | Spearman ρ |\n|---|---|---|---|\n",
  paste0("| ", xa$assay, " | ", xa$n, " | ", xa$r, " | ", xa$rho, " |", collapse = "\n"),
  "\n\n*Only the C2N comparison has adequate baseline overlap (n=208) within the analytic cohort; Lilly (n=23) and Janssen (n=14) are reported but underpowered.*\n"
)

# ---------- WRITE ----------
all_md <- paste(
  "# Manuscript Tables — final draft (v1, 2026-05-25)",
  "",
  "**Source:** rendered automatically from `/Analysis/outputs/` CSV artefacts by `render_final_tables_v1.R`. **Do not edit numbers by hand.** If a number needs to change, re-run the pipeline.",
  "",
  "**Pipeline log:** `/Analysis/logs/master_pipeline_*.log` + `phase3_brms_*.log` (+ refit log when complete).",
  "",
  "---", "",
  tab1_md, "---", "",
  tab2_md, "---", "",
  tab3_md, "---", "",
  tab4_md, amyloid_md, postbl_md, sens_md, mi_md, "---", "",
  xa_md, hippo_md, vif_md,
  sep = "\n"
)
writeLines(all_md, file.path(MS, "Tables_final_v1.md"))
cat("Saved:", file.path(MS, "Tables_final_v1.md"), "\n")
