# R Re-run Checklist — LOAD Plasma GFAP / p-tau217 Trajectory Analysis (v1)

**Owner:** Kuan-Yu Chu
**Independent code reviewer:** Dr. Feng-Ping Chu (to be confirmed)
**Pipeline files (current):** `ADNI_analysis_v3.R`, `ADNI_sensitivity_intercept.R`
**Pipeline files (to be added):** `ADNI_jointmodel_v1.R`, `ADNI_renal_BMI_sensitivity_v1.R`, `ADNI_LCGMM_exploratory_v1.R`
**Target:** all numbers in Abstract / Tables 1–4 / Figures 2–4 must be regenerated from this re-run.

> **Reviewer workflow:** Each item below has (i) an actionable step, (ii) a sanity check, and (iii) a required artefact (CSV / plot / log) that Dr. Chu should sign off on before the manuscript is finalised. Use the checkboxes.

---

## Step 0 — Environment & reproducibility

- [ ] **0.1** Pin R version: `R 4.5.3`. Capture `sessionInfo()` to `/Analysis/logs/sessionInfo_<date>.txt`.
- [ ] **0.2** Use `renv::init()` (or `groundhog::groundhog.library(..., date = "2026-05-01")`) so package versions are frozen. Commit `renv.lock`.
- [ ] **0.3** Required packages (versions ≥): `lcmm 1.9.2`, `lme4 1.1-35`, `lmerTest 3.1-3`, `nlme 3.1-164`, `JMbayes2 0.5-0`, `simex 1.8`, `mice 3.16`, `nnet 7.3`, `ggplot2`, `dplyr`, `tidyr`, `broom.mixed`, `performance`, `DHARMa`, `gtsummary`.
- [ ] **0.4** Set a global RNG seed at top of every script: `set.seed(20260525)`.
- [ ] **0.5** All file paths read from a single `config.yml`; no hard-coded user paths.
- [ ] **0.6** `git init` (or move to existing repo); commit raw scripts before re-run starts. Reviewer reviews **commit-by-commit**, not the final blob.

## Step 1 — Data load & participant flow (drives Figure 1, Table 1)

- [ ] **1.1** Load raw ADNI CSVs from `/ADNI_csv/`: UPENN biomarker panel (GFAP, p-tau217 pT217_F + Lilly/Janssen/C2N columns), APOE genotype, ADNIMERGE, ADSP-harmonised FreeSurfer (`*_combat`), vital signs (BMI), labs (serum creatinine → eGFR via CKD-EPI 2021 race-free).
- [ ] **1.2** Apply inclusion filters **sequentially** and log N at each step (this becomes the **STROBE Figure 1 flow diagram**):
  - age ≥65 at baseline
  - baseline dx ∈ {CN, MCI, AD}
  - APOE genotype available
  - **≥2 plasma GFAP measurements** (relaxed from ≥3 per SAP §4 amendment 2026-05-25)
  - **≥2 plasma p-tau217 (UPENN-Fujirebio) measurements** (relaxed from ≥3)
  - non-missing covariates (age, sex, edu, dx, BMI). eGFR pending lab download — if available, add as final gate.
- [ ] **1.3** Sanity check: **final N should be ≈ 325** per `outputs/inclusion_sensitivity.csv`. Diagnosis breakdown approximately CN 181 / MCI 106 / AD 38; APOE4+ ≈ 46%.
- [ ] **1.4** Compute `follow_up_months` per participant; report **median (IQR)** — replace the misleading "120 months max" in Table 1.
- [ ] **1.5** Compute baseline GFAP and p-tau217 mean (SD), overall and APOE-stratified — add as new rows to Table 1.
- [ ] **1.6** **Artefact:** `/Analysis/outputs/table1_baseline.csv` + `/Analysis/outputs/figure1_strobe_flow.png`.

## Step 2 — Variable construction & QC

- [ ] **2.1** Log-transform GFAP and p-tau217: `log(value)` (natural log). Confirm no zeros / negatives; treat below-LOQ sentinels per assay documentation (drop, *not* impute as half-LOQ).
- [ ] **2.2** Hippocampal volume: `hippo = Left.Hippocampus_combat + Right.Hippocampus_combat`; ICV-normalise as **ratio** `hippo / EstimatedTotalIntraCranialVol_combat` (preferred) and as residualised version (regress hippo on ICV; keep residual). Report which is used.
  - ✅ **Unit-verification result (2026-05-25):** raw ADSP-FreeSurfer data confirmed in mm³ (median L=3,274; R=3,361 mm³; ICV=1,527,336 mm³). The v3 β=−3,439 mm³/yr "per unit log-slope" is **not a unit error** — it reflects the unstandardised slope predictor whose actual range is 0.07–0.12, giving a typical effect of ≈ −310 mm³/yr (biologically plausible). **Fix is reporting unit, not data:** all Phase III β must be re-scaled to **per 1-SD log-slope** before tabulation (Step 5.5). Retain ICV-ratio outcome as primary.
- [ ] **2.3** Compute eGFR (CKD-EPI 2021) at baseline; carry as covariate.
- [ ] **2.4** Standardise all continuous predictors (age, edu, BMI, eGFR, baseline biomarker) using `scale()` so coefficients are comparable.
- [ ] **2.5** **Artefact:** `/Analysis/outputs/qc_distributions.pdf` (histograms + QQ-plots of log-GFAP, log-p-tau217, hippo ratio, eGFR).

## Step 3 — Phase I: LCGMM (now **exploratory only**)

- [ ] **3.1** Re-label all script comments / output filenames / Table 2 caption to read **"Exploratory trajectory-class diagnostic (LCGMM)"**.
- [ ] **3.2** Fit k = 1–4 with class-specific linear slopes and subject-level random intercept+slope, on log-GFAP and log-p-tau217 separately.
- [ ] **3.3** Report BIC, AIC, log-lik, entropy, **mean posterior class probability**, **min class size %**, and **parametric bootstrap LRT** (`lcmm::permutest` or manual bootstrap, B≥500) per k.
- [ ] **3.4** Pre-specify: a class solution is "supported" only if ΔBIC ≥ 10 vs k=1 *and* min class size ≥ 5% *and* mean posterior prob ≥ 0.70 *and* bootstrap LRT p<0.05. (Document that none is expected to meet this bar at N=77.)
- [ ] **3.5** Add an explicit caveat sentence in Table 2 footnote: *"At N=77, statistical power for detecting genuine latent classes is limited; null findings should not be interpreted as evidence of homogeneity."*
- [ ] **3.6** **Artefact:** `/Analysis/outputs/table2_lcgmm.csv` + `/Analysis/outputs/lcgmm_diagnostics.pdf`.

## Step 4 — Phase II: Primary LME (APOE × time)

- [ ] **4.1** For each biomarker, fit:
  ```
  log_biom ~ time*APOE4 + age_z + sex + edu_z + dx
             + BMI_tv_z + eGFR_tv_z         # Tier 1 time-varying
             + HTN_base + DM_base            # Tier 2 baseline binary
             + (1 + time | RID)
  ```
  using `lmer()` with REML; report fixed effects with **95% CI** (profile or Wald), random-effect variances, and `isSingular()` flag. Time-varying BMI and eGFR carried as LOCF within 12 months of biomarker visit. If eGFR data unavailable at re-run time, omit `eGFR_tv_z` and note in §10 transparency.
- [ ] **4.2** Sensitivity A: drop random slope (`(1 | RID)`).
- [ ] **4.3** Sensitivity B: **Bayesian fit with weakly informative priors** (`brms::brm`, 4 chains, 4000 iters) to resolve singular boundary issues from REML — this should be the preferred sensitivity, not the random-intercept-only model.
- [ ] **4.4** Sensitivity C: drop eGFR; drop BMI; both — confirm direction of effects unchanged.
- [ ] **4.5** Residual diagnostics with `DHARMa::simulateResiduals()`; attach plot.
- [ ] **4.6** Compute per-participant EB BLUPs of random slope: `ranef(fit)$RID[, "time"] + fixef(fit)["time"]`. Save to `/Analysis/outputs/blup_slopes.csv` with columns: RID, GFAP_slope, GFAP_slope_SE, ptau_slope, ptau_slope_SE.
- [ ] **4.7** **Artefact:** `/Analysis/outputs/table3_lme_apoe.csv` + residual diagnostic PDF.

## Step 5 — Phase III: Slope-as-predictor (joint model preferred)

- [ ] **5.1** Primary: **joint longitudinal-outcome model** using `JMbayes2::jm()` — one shared-parameter model per (biomarker, outcome) pair — to propagate biomarker measurement error into the outcome regression. Six models total (2 biomarkers × 3 outcomes: ADAS-Cog 13, MMSE, hippo ratio).
- [ ] **5.2** Sensitivity 1: two-stage BLUP-as-predictor (the current v3 approach) — report side-by-side with joint-model β; expect joint-model SEs to widen.
- [ ] **5.3** Sensitivity 2: SIMEX correction on the two-stage model using `simex::simex()` with the BLUP SE as measurement-error SD.
- [ ] **5.4** All Phase-III models adjusted for age, sex, edu, dx, APOE4, BMI, eGFR.
- [ ] **5.5** Report **β per 1-SD log-slope** with 95% CI (not "per unit") — re-scale all coefficients before populating Table 4.
- [ ] **5.6** Multiple-testing: Benjamini–Hochberg FDR across the 6 Phase-III tests as primary; Bonferroni as conservative. Pre-specify in SAP.
- [ ] **5.7** **Artefact:** `/Analysis/outputs/table4_slope_outcome.csv` + `/Analysis/outputs/joint_model_diagnostics.pdf`.

## Step 6 — Cross-assay p-tau217 concordance (Figure 4)

- [ ] **6.1** Restrict to baseline visits of the **analytic N=77** (current v3 reports n=489, which mixes non-analytic participants — fix this).
- [ ] **6.2** Pearson and Spearman correlations between UPENN-Fujirebio and each of Lilly MSD600, Janssen, C2N PrecivityAD2.
- [ ] **6.3** **Bland–Altman plots** for each pair — add as supplementary figure (current Figure 4 only shows correlations).
- [ ] **6.4** If Janssen overlap remains insufficient at N=77, drop from main text and state explicitly.
- [ ] **6.5** **Artefact:** `/Analysis/outputs/figure4_crossassay.png` + `/Analysis/outputs/supp_blandaltman.pdf`.

## Step 7 — Subgroup & robustness

- [ ] **7.1** Stratify Phase II and III by baseline diagnosis (CN vs MCI; AD n=4 excluded from stratified models, noted explicitly).
- [ ] **7.2** Stratify by sex.
- [ ] **7.3** Influence diagnostics: refit primary models excluding top 5% Cook's-distance participants; report stability.
- [ ] **7.4** **Artefact:** `/Analysis/outputs/supp_subgroup_forest.pdf`.

## Step 8 — Figures (regenerate, do not reuse v3)

- [ ] **8.1** **Figure 1** — STROBE participant-flow diagram (auto-generated from Step 1.2 log; do **not** leave as "to be drawn manually").
- [ ] **8.2** **Figure 2** — Spaghetti plots A (GFAP) / B (p-tau217), x = months from baseline, y = ln(pg/mL), coloured by APOE4, with APOE-stratified marginal LME fits (`ggeffects::predict_response`). Legend must specify: units, transformation, model, n per group.
- [ ] **8.3** **Figure 3** — 2×3 panel grid of BLUP slope vs outcome rate (with shrinkage-corrected joint-model overlay); each panel labelled A–F with Pearson r and joint-model β (95% CI).
- [ ] **8.4** **Figure 4** — Cross-assay scatter (UPENN vs C2N; UPENN vs Lilly), restricted to N=77 baseline; Janssen moved to supplement.
- [ ] **8.5** Save all figures at 300 dpi, both PNG and PDF, under `/figures/v4/`.

## Step 9 — Tables: regenerate from CSV artefacts

- [ ] **9.1** Use `gtsummary` or `flextable` to render Tables 1–4 directly from the CSVs in `/Analysis/outputs/` — eliminate any hand-typed numbers in the .docx.
- [ ] **9.2** Add **95% CIs** to every β in Tables 3 and 4.
- [ ] **9.3** Add **q (BH-FDR)** column to Table 4.

## Step 10 — Reproducibility package & code review handoff

- [ ] **10.1** Push final pipeline to a private GitHub repo; tag `v1.0-submission`.
- [ ] **10.2** Generate Zenodo DOI; update Data Availability statement.
- [ ] **10.3** Send Dr. Feng-Ping Chu: (a) repo link, (b) this checklist, (c) `/Analysis/outputs/` artefacts, (d) `sessionInfo` log.
- [ ] **10.4** Reviewer (Dr. Chu) sign-off log at `/Analysis/code_review_signoff.md` with date, commit SHA reviewed, and explicit `LGTM` per checklist section.
- [ ] **10.5** Only after sign-off: copy regenerated numbers into Abstract / Results / Tables; **do not edit narrative numbers earlier**, to avoid divergence between text and outputs.

---

## Red-flag triage before any rerun

These three items override the rest of the checklist — fix before anything else:

1. **Hippocampal-volume reporting unit** (Step 2.2) — ✅ resolved 2026-05-25: data units verified mm³; v3 implausibility caused by unstandardised slope predictor, not raw data. Fix is re-scale β to per 1-SD log-slope at Step 5.5.
2. **eGFR confounding** (Steps 1.1, 2.3, 4.1) — ⏳ blocked: ADNI lab/biospecimen CSV with serum creatinine not yet downloaded. **Action required by Kuan-Yu Chu**: download from ida.loni.usc.edu (search "Laboratory Data" → "Clinical Chemistry"). Until then, primary models omit eGFR and SAP §10 transparency applies.
3. **Two-stage BLUP-as-predictor measurement error** (Step 5.1) — ⏳ in progress: JMbayes2 to be installed and joint-model script written (Step 5).
