# Plasma GFAP and p-tau217 Trajectories in Older Adults — ADNI Joint Longitudinal Analysis

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](LICENSE.md)
[![STROBE](https://img.shields.io/badge/Reporting-STROBE-green)](LOAD_NeurolSci_R1_Submission/05_STROBE_Checklist.docx)

Reproducibility repository for the manuscript:

> **APOE ε4 Accelerates Plasma p-tau217 Trajectory: A Bayesian Joint Longitudinal Analysis of ADNI**
>
> Kuan-Yu Chu, DDS, MHA, PhD (sole author, 2026). Submitted to *Neurological Sciences* (Springer); revision (R2) in response to peer review.

> **Note on DOI:** a prior version of this README referenced a Zenodo DOI (10.5281/zenodo.20393272) that was never actually minted — the repository had been initialised locally but nothing had been committed or archived. This is corrected as of the v1.1.0 release below; see [Data and Code Availability](#data-and-code-availability) for the real, resolving DOI once published.

---

## Summary

This repository contains the **complete analysis pipeline, statistical analysis plan (SAP), and audit trail** for a longitudinal observational ADNI secondary analysis (2026-05-24 data release; N = 325 eligible, N = 296 analytic after complete-case restriction on Tier 1/2 covariates; CN 181 / MCI 106 / AD 38; up to 10 years follow-up).

**Key findings (as of the R2 revision):**

- APOE ε4 carriage was associated with a steeper plasma p-tau217 slope (time×APOE4 interaction β = 0.0307 log-units/year, p = 0.002) but **not** plasma GFAP slope (p = 0.149).
- A formal counterfactual mediation analysis found this effect is **predominantly direct** (ADE = 0.029, p < 0.001), with a small (6.4%), only partially robust indirect pathway through baseline amyloid positivity (ACME p = 0.045; attenuates to marginal significance, p = 0.088, when amyloid is modeled continuously rather than dichotomized). This reverses the original manuscript's "primarily amyloid-mediated" conclusion, following a code bug found and fixed during peer review (see [Audit trail](#audit-trail-and-the-r2-mediation-bug) below).
- Six separate bivariate Bayesian joint longitudinal models (`brms`, one per biomarker-outcome pair — not a single multivariate model) estimated individual-level co-trajectory between p-tau217 slope and outcomes: ADAS-Cog 13 ρ = +0.67, MMSE ρ = −0.69, hippocampal ratio ρ = −0.39 (all excluding zero). The three GFAP-paired models did not fully converge on all parameters (max Rhat up to 1.36) and their null findings are reported as provisional.

---

## Audit trail and the R2 mediation bug

During peer review, a reviewer's concern about the mediation analysis's construction led to finding a real coding bug in `phase3_me_correction_v1.R`: the per-participant slope used as the mediation outcome was reconstructed as `random_deviation + fixed_main_effect_of_time`, omitting the `APOE4 × time` fixed-interaction term the underlying model actually estimates. This silently stripped most of the APOE4 signal out of the mediation outcome before `mediate()` ever saw it, producing a spurious near-zero total effect in the original submission.

The fix (`phase3_me_correction_v2_fixed.R`) and its full downstream consequences (corrected mediation results, a resolved apparent contradiction with an amyloid-subgroup finding, and a parallel hippocampal-ratio finding that turned out to be the same statistical artifact) are logged chronologically in `Analysis/code_self_review_log.md`, with every new script from this revision listed there.

---

## Repository structure

```
LOAD-BC-NS/
├── LOAD_NeurolSci_R1_Submission/       # Current submission package (R2, Neurological Sciences)
│   ├── 01_Title_Page.docx
│   ├── 02_Cover_Letter.docx
│   ├── 03_Main_Text.docx               # Blinded manuscript
│   ├── 04_Supplementary_Materials.docx
│   ├── 05_STROBE_Checklist.docx
│   ├── 07_Response_to_Reviewer1.docx   # Point-by-point response
│   ├── Supplementary_File_S1_SAP.docx
│   └── Figures/                        # Main + Supplementary figures (TIFF)
│
├── Analysis/
│   ├── ADNI_master_pipeline_v1.R       # Phases I-IV master script
│   ├── phase3_brms_joint_v1.R          # Primary Phase III: 6 bivariate brms joint models
│   ├── phase3_brms_gfap_refit_v1.R     # GFAP refit with longer chains (did not fully converge; see Table S7)
│   ├── phase3_me_correction_v1.R       # SUPERSEDED (has the R2 bug); kept for audit trail
│   ├── phase3_me_correction_v2_fixed.R # R2 fix: corrected BLUP reconstruction
│   ├── mediation_analysis_v1.R         # SUPERSEDED (uses buggy BLUPs)
│   ├── mediation_analysis_v2_fixed.R   # R2 primary mediation result
│   ├── mediation_analysis_v3_no_dx3.R  # R2 collider sensitivity
│   ├── mediation_analysis_v4_continuous_amyloid.R  # R2 continuous-mediator sensitivity
│   ├── amyloid_APOE_time_interaction_v1.R          # R2: formal 3-way interaction (resolves Sec 3.5 artifact)
│   ├── hippocampal_amyloid_bygroup_joint_model_v1.R # R2: by-group pooled brms model
│   ├── hippocampal_amyloid_bygroup_extract_diff_v1.R
│   ├── diagnosis_stage_interaction_v1.R            # R2: time x diagnosis, 3-way sensitivity
│   ├── apoe_allele_dose_sensitivity_v1.R           # R2: real APOE genotype dose-response
│   ├── trajectory_linearity_check_v1.R             # R2: quadratic/spline sensitivity
│   ├── brms_rhat_ess_divergence_v1.R               # R2: full convergence diagnostics
│   ├── brms_full_ess_diagnostics_v1.R
│   ├── brms_posterior_predictive_check_v1.R
│   ├── hippocampal_bygroup_full_diagnostics_v1.R
│   ├── inclusion_sensitivity_v1.R
│   ├── parse_comorbidities_v1.R
│   ├── render_final_tables_v1.R
│   ├── make_figures_v4.R
│   ├── code_self_review_log.md         # Chronological audit trail, including the R2 bug
│   └── outputs/                        # CSV/RDS artefacts (whitelisted in .gitignore where de-identified)
│
├── ADNI_csv/                            # NOT IN REPO — see "Obtaining the data" below
├── .gitignore
├── CITATION.cff
├── LICENSE.md
└── README.md                            # This file
```

---

## Obtaining the data

This repository **does not include raw ADNI participant-level data**, in accordance with the [ADNI Data Use Agreement](https://adni.loni.usc.edu/data-samples/access-data/).

To reproduce the analysis:

1. Apply for ADNI access at https://ida.loni.usc.edu/
2. Download the following files (2026-05-24 release):
   - `UPENN_PLASMA_FUJIREBIO_QUANTERIX_24May2026.csv`
   - `APOERES_24May2026.csv`
   - `ADAS_24May2026.csv` and `MMSE_24May2026.csv`
   - `VITALS_24May2026.csv`
   - `ADSP_PHC_T1_FS_24May2026.csv`
   - `LABDATA_25May2026.csv` (for serum creatinine, column `RCT392`)
   - `Medical_History.zip` (contains `RECMHIST_24May2026.csv`, `INITHEALTH_24May2026.csv`)
   - `LILLY_PTAU217_MSD600_24May2026.csv`, `C2N_PRECIVITYAD2_PLASMA_24May2026.csv`, `JANSSEN_PLASMA_P217_TAU_24May2026.csv` (cross-assay)
   - `ADNIMERGE_0.0.1.tar.gz` (R package)
3. Place these in a local folder `ADNI_csv/` at the repository root.

---

## Reproducing the analysis

### Environment

- R 4.5.3 (verified)
- Required R packages: `dplyr`, `tidyr`, `stringr`, `lme4`, `lmerTest`, `nlme`, `lcmm`, `brms` (with `rstan` backend), `posterior`, `mediation`, `splines`, `ggplot2`

### Core pipeline (R1)

```bash
Rscript Analysis/inclusion_sensitivity_v1.R
Rscript Analysis/parse_comorbidities_v1.R
Rscript Analysis/ADNI_master_pipeline_v1.R
Rscript Analysis/phase3_brms_joint_v1.R
Rscript Analysis/phase3_brms_gfap_refit_v1.R
Rscript Analysis/render_final_tables_v1.R
Rscript Analysis/make_figures_v4.R
```

### R2 revision scripts (run after the core pipeline; supersede/extend it)

```bash
Rscript Analysis/phase3_me_correction_v2_fixed.R      # corrected BLUPs (fixes the R2 bug)
Rscript Analysis/mediation_analysis_v2_fixed.R        # primary corrected mediation result
Rscript Analysis/mediation_analysis_v3_no_dx3.R       # collider sensitivity
Rscript Analysis/mediation_analysis_v4_continuous_amyloid.R
Rscript Analysis/amyloid_APOE_time_interaction_v1.R
Rscript Analysis/diagnosis_stage_interaction_v1.R
Rscript Analysis/apoe_allele_dose_sensitivity_v1.R
Rscript Analysis/trajectory_linearity_check_v1.R
Rscript Analysis/hippocampal_amyloid_bygroup_joint_model_v1.R   # slow: ~10+ min, by-group brms fit
Rscript Analysis/hippocampal_amyloid_bygroup_extract_diff_v1.R
Rscript Analysis/brms_rhat_ess_divergence_v1.R
Rscript Analysis/brms_full_ess_diagnostics_v1.R
Rscript Analysis/brms_posterior_predictive_check_v1.R
Rscript Analysis/hippocampal_bygroup_full_diagnostics_v1.R
```

All outputs are written to `Analysis/outputs/`.

---

## Statistical Analysis Plan

The pre-specified SAP (`Supplementary_File_S1_SAP.docx`) was locked on **2026-02-12**, with 10 amendments dated **2026-05-25**, all committed in writing before inspecting outcome-modelling results. The R2 revision's additional sensitivity analyses (listed above) were **not** part of the original SAP and are explicitly labeled as reviewer-motivated, post-hoc, and exploratory in the manuscript (Methods, Reproducibility subsection), consistent with STROBE's recommendation to distinguish pre-specified from post-hoc analyses.

---

## Data and Code Availability

- **Analysis code:** This GitHub repository (released under MIT; see [LICENSE.md](LICENSE.md))
- **Manuscript text, SAP, figures, and tables:** This repository (released under CC BY 4.0)
- **brms fit objects (`.rds`):** Archived on Zenodo (DOI to be added here once minted for this release — see the note at the top of this README)
- **Raw ADNI data:** Available only to qualified investigators via https://ida.loni.usc.edu/

---

## Citation

See [CITATION.cff](CITATION.cff) for machine-readable metadata.

---

## ADNI Acknowledgement

Data used in this analysis were obtained from the Alzheimer's Disease Neuroimaging Initiative (ADNI) database (adni.loni.usc.edu). ADNI investigators contributed to the design and implementation of ADNI and/or provided data but did not participate in analysis or writing of this report. A complete listing of ADNI investigators is available at https://adni.loni.usc.edu/about/governance/.

ADNI is funded by NIH Grant U01 AG024904 with contributions from many private and public partners listed at the ADNI website.

---

## Contact

Kuan-Yu Chu, DDS, MHA, PhD
Taipei, Taiwan
ORCID: [0000-0003-0574-4220](https://orcid.org/0000-0003-0574-4220)
kyc0321@gmail.com
