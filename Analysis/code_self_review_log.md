# Code Self-Review Log

**Purpose:** Chronological audit trail of every analytical decision, debug step, unit/scale verification, and commit. This log replaces external code review (per SAP §9, single-author paper).

**Project:** Plasma GFAP / p-tau217 trajectory analysis in ADNI older adults
**Author:** Kuan-Yu Chu
**SAP:** `Mauuscripts/SAP_v1.md` (locked 2026-02-12; 10 amendments dated 2026-05-25)

---

## Entry log

| Date | Commit (after first push) | Stage | Decision / verification | Outcome |
|---|---|---|---|---|
| 2026-02-12 | — | SAP lock | Three hypotheses H1/H2/H3 pre-specified; LCGMM primary; JMbayes2 primary for Phase III | SAP v1 locked |
| 2026-05-25 | TBD | Sensitivity 1 | Inclusion-criterion sensitivity (≥1/2/3/4 measurements) on 2026-05-24 ADNI release | N = 35 at ≥3; N = 325 at ≥2 → SAP §4 amendment 4: threshold relaxed to ≥2 |
| 2026-05-25 | TBD | Methods 1 | LCGMM repositioned from primary to exploratory only (SAP §6.2) | Amendment 1 |
| 2026-05-25 | TBD | Cohort assembly | VITALS date column verified as `VISDATE` (not `EXAMDATE`) | Bug fixed; pipeline rerun |
| 2026-05-25 | TBD | Unit verification | Hippocampal volume confirmed as raw mm³ (L median 3,274; R median 3,361; ICV 1,527,336); v3 implausible β = −3,439 mm³/yr identified as reporting-unit issue (unstandardised slope predictor), not data error | Slope predictors standardised to per 1-SD throughout |
| 2026-05-25 | TBD | Covariate 1 | Tier 1 (BMI, eGFR) + Tier 2 (HTN, DM) + Tier 3 stroke sensitivity exclusion adopted | SAP §5.4 amendment 3 |
| 2026-05-25 | TBD | Comorbidity regex | Stroke regex `\\bSTROKE\\b\|TIA` initially produced 12.0% false-positive rate due to "ESSENTIAL", "HOSPTIAL", "HIATIAL" matching un-bounded TIA; tightened to `\\bTIA\\b\|TRANSIENT ISCHEMIC` → 4.6% (consistent with literature) | Bug fixed in `parse_comorbidities_v1.R` |
| 2026-05-25 | TBD | Phase III primary | JMbayes2::jm() attempted on 6 (biomarker, outcome) pairs; all failed with "zero-length argument" because JMbayes2 requires a survival component our design does not have | SAP §10 amendment 6: pivot to MI-on-BLUP |
| 2026-05-25 | TBD | Phase III primary | MI-on-BLUP two-stage attempted; all 6 effects null (RIV 0.74–0.93); identified BLUP-as-fixed-but-noisy treatment as over-conservative because MI cannot share information through random-effects covariance | SAP §10 amendment 9: pivot to brms multivariate joint |
| 2026-05-25 | TBD | Phase III primary | brms multivariate joint model with `\| p \|` shared random effects across responses; ADAS_GFAP / MMSE_GFAP / Hippo_GFAP initial Rhat 1.05–1.36 (poor); refit with iter = 6000, warmup = 2000, adapt_delta = 0.99 → all Rhat ≤ 1.05 | Final Phase III primary |
| 2026-05-25 | TBD | Confounder 1 | eGFR found available for only ~50/325 RIDs within 2 years of biomarker baseline (ADNI3/4 reduced routine chemistry collection); demoted from primary to sensitivity per SAP §10 amendment 7 | Primary models omit eGFR |
| 2026-05-25 | TBD | Confounder 2 | BMI initially specified as time-varying; visit-level missingness reduced random-slope-estimable RIDs from 296 to 144; demoted to baseline-only per SAP §10 amendment 8 | BMI now baseline-only |
| 2026-05-25 | TBD | Inference rule | Pre-specified inference rule for Phase III: ρ "supported" if 95% CrI excludes 0 AND P(ρ has hypothesised sign) ≥ 0.975 | Applied in Tables 4 and Discussion §4.1 |
| 2026-05-25 | TBD | Hypothesis check | H1 (APOE → GFAP slope) and H2 (APOE not → p-tau217 slope) directionally falsified | Reported transparently in Abstract, Intro §1.6, Results §3.3, Discussion §4.1 (no retrofitting) |
| 2026-05-25 | TBD | Authorship | Single-author paper finalised after invitation/withdrawal cycle with Dr. Feng-Ping Chu; courtesy withdrawal sent and acknowledged | SAP §10 amendment 5 marked REVERTED; Dr. Chu thanked in Acknowledgements with his permission |
| 2026-05-25 | TBD | First public commit | Repository initialised; .gitignore excludes raw ADNI CSVs per DUA; brms `.rds` fit objects excluded (size); to be archived on Zenodo separately | v1.0.0 |
| 2026-09-11 | Claude (Anthropic), reviewed by author | R2 mediation bug | `phase3_me_correction_v1.R` reconstructed BLUP slopes as `random_deviation + fixed_main_effect_of_time`, omitting the `APOE4 x fixed_interaction_effect` term — silently stripped most of the APOE4 signal from the mediation outcome before `mediate()` ever saw it. Confirmed: pre-fix, APOE4 carriers'/non-carriers' mean reconstructed slopes were numerically identical (both 0.0578); post-fix, carriers rose to 0.0885 (+0.0307, exactly matching the primary model's time:APOE4 coefficient) | `phase3_me_correction_v2_fixed.R`; mediation re-run in `mediation_analysis_v2_fixed.R` |
| 2026-09-11 | Claude, reviewed by author | R2 mediation sensitivities | Collider check (dx3 excluded): `mediation_analysis_v3_no_dx3.R`. Continuous-amyloid mediator: `mediation_analysis_v4_continuous_amyloid.R` | PM 6.4% (primary) -> 7.7% (no dx3) -> 0.9%, p=0.088 (continuous); reported as exploratory given sensitivity to specification |
| 2026-09-11 | Claude, reviewed by author | R2 numerical inconsistencies | Abstract/Results/Fig.2 legend/Sec 3.3 reported the eGFR-restricted sensitivity's beta=0.0392/0.493 as if it were the primary model's (688/296) beta=0.0307/0.2877 | All instances corrected and cross-checked against Table 2 |
| 2026-09-11 | Claude, reviewed by author | R2 amyloid subgroup artifact | Naive A+/A- subgroup contrast (one CI excludes zero, one doesn't) is not itself evidence of effect modification. Formal 3-way interaction (time x APOE4 x A_pos, full cohort) found no significant modification (p=0.23) | `amyloid_APOE_time_interaction_v1.R`; Sec 3.5/3.5a/4.1/4.3 rewritten |
| 2026-09-11 | Claude, reviewed by author | R2 hippocampal subgroup artifact | Same artifact as above, applied to the p-tau217-hippocampal-ratio rho: by-group pooled brms model found A-/A+ slope-slope correlations statistically indistinguishable (P(diff>0)=0.73 vs pre-specified 0.975 threshold) | `hippocampal_amyloid_bygroup_joint_model_v1.R`, `hippocampal_amyloid_bygroup_extract_diff_v1.R`; Sec 3.5 rewritten |
| 2026-09-11 | Claude, reviewed by author | R2 diagnosis-stage sensitivity | time x diagnosis (additive) and time x APOE4 x diagnosis (3-way) tested; primary APOE4 x time estimate robust (beta=0.0272 vs 0.0307); no individual 3-way term significant | `diagnosis_stage_interaction_v1.R`; Sec 3.5c |
| 2026-09-11 | Claude, reviewed by author | R2 APOE allele dose | Real ADNI genotype (APGEN1/APGEN2) pulled; excluding e4 homozygotes did not attenuate the effect (beta=0.0332 vs 0.0307); continuous 0/1/2 dose model significant (p=0.011) | `apoe_allele_dose_sensitivity_v1.R`; Sec 3.3 |
| 2026-09-11 | Claude, reviewed by author | R2 trajectory linearity | Quadratic-time and 2-df spline sensitivity on primary p-tau217 model; no evidence of nonlinearity (LRT chi-sq(2)=2.62, p=0.27) | `trajectory_linearity_check_v1.R`; Sec 3.5c |
| 2026-09-11 | Claude, reviewed by author | R2 Bayesian reporting completeness | Full accounting of priors (all brms defaults) and per-model Rhat/ESS/divergences/PPC revealed the ms's "all Rhat<=1.05" claim was false: 3 GFAP-paired models did not fully converge (max Rhat 1.36, min bulk ESS 10) | `brms_rhat_ess_divergence_v1.R`, `brms_full_ess_diagnostics_v1.R`, `brms_posterior_predictive_check_v1.R`, `hippocampal_bygroup_full_diagnostics_v1.R`; Supplementary Tables S6-S7; Sec 3.4/4.6 |
| 2026-09-11 | TBD | R2 archive update | v1.1.0: all scripts above added; README/CITATION.cff/title updated to match the corrected (direct-pathway-dominant) narrative | v1.1.0 |

---

## Pre-commit verification (per submission)

Before every git tag (e.g., `v1.0-submission`), verify:

- [ ] All numbers in `Mauuscripts/Main_Text_v1.docx` Tables and inline text match the corresponding cell in `Analysis/outputs/*.csv`
- [ ] All figures in `figures/v4/` are regenerated from the latest pipeline (check file mtime)
- [ ] `sessionInfo()` log saved to `Analysis/logs/sessionInfo_<date>.txt`
- [ ] `.gitignore` is intact and `/ADNI_csv/` is not staged (`git ls-files | grep ADNI_csv/` should be empty)
- [ ] No `.rds` larger than 100 MB staged (`git ls-files | xargs -I {} du -h {} | sort -rh | head`)
- [ ] SAP §10 amendments dated correctly; any new amendments added since last submission are flagged in the cover letter

---

## Known limitations of this self-review approach

This log is the author's own audit trail, not external peer review. It does not substitute for journal peer review or replication by independent investigators. Readers are encouraged to:

- Re-run the full pipeline from the public scripts
- Inspect the SAP amendment history for any post-hoc modifications
- Compare numerical outputs in `Analysis/outputs/` against the manuscript tables

Issues, bugs, or replication failures should be reported via GitHub Issues at the repository URL.
