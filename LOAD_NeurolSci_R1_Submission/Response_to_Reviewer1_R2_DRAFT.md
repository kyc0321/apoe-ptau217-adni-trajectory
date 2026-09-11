# Response to Reviewer #1 — Neurological Sciences, decision 2026-09-10

Manuscript: "APOE ε4 Accelerates Plasma p-tau217 Trajectory: A Bayesian Joint
Longitudinal Analysis of ADNI" (revised title; see response to comment 1/4
below for why the title changed).

Status key: **DONE** = fixed/re-analyzed and reflected in the revised manuscript.
**IN PROGRESS** = not yet completed, planned before resubmission (deadline 28 Sep 2026).

---

## Comment 1 — Numerical inconsistencies (abstract/text vs. Table 1/Table 2)

**DONE.** We thank the reviewer for this careful check. Two distinct problems were found and fixed:

1. Table 1 vs. body text (age, sex, hypertension, diabetes, stroke, measurement counts): the body text has been corrected to match Table 1 exactly (mean age 73.9y, 43.7% women, 60.3% hypertension, 15.7% diabetes, 7.7% stroke, GFAP/p-tau217 measurement medians 2(2-2)/2(2-3)).
2. The more serious issue: the abstract and Results text reported β=0.0392 (interaction) and β=0.493 (APOE4 main effect) as the primary p-tau217 findings. On tracing these through the underlying results table, we confirmed these values belong to the eGFR-restricted **sensitivity** analysis (n=227/109), not the primary model (n=688/296, β=0.0307 interaction, β=0.2877 main effect). This has been corrected throughout the abstract, Results text (Section 3.3 also had the same mix-up, caught during the Section 7 dose-response work below), Figure 2 legend, and all other locations where these values were cited. We are grateful this was caught, as it affects the headline effect size.

## Comment 2 — 325 → ~296 sample-size accounting

**DONE.** The 325 participants are the full eligibility-screen cohort (age≥65, baseline diagnosis, APOE available, ≥2 GFAP and ≥2 p-tau217 measurements). The regression models (Table 2 onward) additionally require complete-case data on all Tier 1/Tier 2 covariates, which restricts N from 325 to 296. We traced this exactly: **all 29 excluded participants were missing exactly one variable, baseline BMI** — no other covariate (age, sex, education, diagnosis, APOE4, hypertension, diabetes) had any missingness in the full 325. We checked whether this BMI missingness was informative: it did not differ by APOE4 status (9.2% vs 8.6% in non-carriers/carriers, p=1.0, so no exposure-related selection bias), but was concentrated in the CN group (14.4% missing vs 2.8% MCI, 0% AD; p<0.001). Since diagnosis is a covariate retained in all models rather than the exposure of interest, we do not believe this introduces exposure-related bias, but we now state the pattern explicitly (Section 2.1) rather than leaving the 325→296 transition unexplained.

## Comment 3 — Mediation analysis / BLUP construction

**DONE — this was a real, confirmed bug, not only a conceptual concern.** The reviewer's intuition was exactly right. We traced the code and found that the BLUP slope used as the mediation outcome was reconstructed as `random_deviation + fixed_main_effect_of_time`, omitting the `APOE4 × fixed_interaction_effect` term that the underlying model (`time_yrs*APOE4 + ... + (1+time_yrs|RID)`) actually estimates. Because `ranef()` returns each participant's deviation net of their *full* fixed-effects prediction (which, for an APOE4 carrier, includes the interaction term), omitting that term when reconstructing the slope silently stripped out most of the APOE4-associated signal before the mediation model ever saw it. Confirmation: before the fix, APOE4 carriers' and non-carriers' mean reconstructed slopes were numerically identical (both 0.0578); after the fix, carriers' mean slope rose to 0.0885, a difference of +0.0307 — exactly matching the primary model's time×APOE4 coefficient.

We re-ran the mediation analysis with the corrected BLUP. Results changed substantially: ADE (direct effect) = 0.0288 (95% CI 0.0236-0.0343, p<0.001), ACME (indirect, via amyloid) = 0.00196 (95% CI 0.00019-0.00411, p=0.045), Total effect = 0.0307 (95% CI 0.0259-0.0352, p<0.001) — the total now exactly reproduces the primary model's coefficient, as it should for a correctly decomposed effect. Proportion mediated = 6.4% (95% CI 0.6-13.9%). We have revised the manuscript's title, abstract, Results, and Discussion accordingly: the effect is now reported as predominantly **direct**, with a small, exploratory (see comment 5) amyloid-mediated component, rather than "primarily amyloid-mediated" as originally claimed. The manuscript title was changed to remove this now-unsupported claim.

## Comment 4 — Causal language / baseline diagnosis as a possible collider

**DONE for the collider check; causal language softened throughout.** We re-ran the corrected mediation analysis excluding baseline diagnosis (dx3) from both the mediator and outcome models. The proportion-mediated result was robust to this exclusion (7.7%, 95% CI 1.9-16.4%, p=0.012, vs. 6.4% with dx3 included) — the small indirect effect does not appear to be an artifact of conditioning on a variable that may lie downstream of APOE4/amyloid. We have also revised "causal mediation" language throughout to state explicit assumptions and to avoid asserting that APOE4 "works by first making amyloid worse" given the small and only partially robust size of the mediated component (see comment 5).

## Comment 5 — Amyloid measurement detail and continuous modeling

**DONE — both parts.** We re-ran the mediation analysis modeling baseline Aβ42/Aβ40 continuously instead of dichotomizing at 0.085. The direct effect remained essentially unchanged (ADE=0.0302, p<0.001), but the indirect pathway attenuated to marginal significance (proportion mediated 0.9%, 95% CI 0.4-15.6%, p=0.088). We now report the indirect/mediated finding as exploratory and sensitive to how amyloid is operationalized, rather than a robust result, throughout the manuscript.

We also added the requested Methods detail (Section 2.1): plasma Aβ42/Aβ40 is measured on the same Fujirebio Lumipulse G1200 platform as p-tau217 (singlicate; documented QC-run precision 5.1-5.4% CV for Aβ42, 3.7-4.2% CV for Aβ40, per the primary ADNI UPENN methods documentation), and the amyloid and p-tau217-baseline measurements come from the same blood draw in 99.1% of participants (median gap 0 days; one outlier with a 2,427-day gap, excluded from this timing check but not from the main analysis). Amyloid status is therefore not a delayed or discordant-visit proxy in this cohort.

## Comment 6 — Linearity of the biomarker trajectory (quadratic-time / spline sensitivity)

**DONE.** We added a quadratic time term and, separately, a 2-df restricted cubic spline (both fixed-effect only; a random quadratic slope was not identifiable given a median of 2-3 measurements per participant) to the primary p-tau217 model. A likelihood-ratio test found no evidence that the nonlinear specification fits better than the linear one (quadratic term: χ²(2)=2.62, p=0.27; spline curvature term p=0.11, its APOE4 interaction p=0.83). The reported APOE4×time interaction is therefore not an artifact of forcing a linear form on the data. As the reviewer anticipated, this check has limited power given the sparse per-participant sampling (71% of participants have only 2 measurements), which we now note explicitly as a limitation on this specific check rather than treating the result as definitive proof of linearity.

## Comment 7 — APOE allele dose (ε2/ε3/ε4 genotype counts, homozygote sensitivity)

**DONE.** We obtained genotype (APGEN1/APGEN2, ADNI APOERES file) for all 296 analytic participants; the genotype-derived carrier status agreed with the previously used binary flag for every observation. Distribution: ε3/ε3 n=139 (reference), ε3/ε4 n=104, ε4/ε4 n=28, ε2/ε3 n=19, ε2/ε4 n=6. Three sensitivity models were fit: (1) excluding the 28 ε4 homozygotes, the interaction persisted at similar magnitude in the remaining heterozygote/non-carrier participants (β=0.0332, p=0.002, vs. β=0.0307 in the full primary model) — the association is not disproportionately driven by homozygotes; (2) a continuous 0/1/2-allele dose model was itself significant (β=0.0205 per allele, p=0.011), consistent with a general dose-related effect; (3) modeled as a 3-level factor, the heterozygote effect was clearly significant (β=0.0331, p=0.002) while the smaller homozygote group (n=28) did not reach significance alone (β=0.0181, p=0.37) — given this subgroup's small size we interpret this as underpowered rather than a genuine non-monotonic pattern, and note it explicitly as such in the manuscript rather than over-interpreting it. All of this is now reported in Results Section 3.3.

## Comment 8 — Disease-stage heterogeneity (time×diagnosis, APOE×time×diagnosis)

**DONE.** Motivated by the reviewer's general point about testing interactions formally in the pooled model rather than relying on underpowered subgroup splits, we applied this logic twice.

First, we re-examined our own amyloid-stratified subgroup analysis (Section 3.5), which originally reported near-complete attenuation of the APOE×time effect within the amyloid-positive subgroup. A formal three-way time×APOE4×amyloid-status interaction in the full cohort found **no significant evidence of effect modification by amyloid status** (β=-0.025, SE=0.021, p=0.23) — the apparent subgroup attenuation was an artifact of informally contrasting a significant and a non-significant subgroup estimate, not a real interaction. This is now corrected in the manuscript (Section 3.5) and reinforces the comment-3/5 finding that the direct pathway dominates.

Second, we applied the same approach to diagnosis, as directly requested. Adding time×diagnosis as an additive (stage-dependent-slope) term left the primary APOE4×time estimate essentially unchanged (β=0.0272, p=0.007, vs. β=0.0307 in the primary model) — including diagnosis as a stage-dependent slope rather than only a covariate does not materially alter the APOE4 estimate, addressing the reviewer's alternative suggestion directly. A full three-way time×APOE4×diagnosis model showed a marginally significant improvement in overall fit (χ²(4)=11.38, p=0.023), but no individual three-way term reached significance (MCI: p=0.95; AD: p=0.46, underpowered at n=80 AD observations); we report this conservatively as inconclusive rather than as evidence that the APOE-time slope differs by disease stage. Both analyses are now in Results Section 3.5c.

Third, we applied the same logic once more to a finding the reviewer did not directly raise but that shared the identical statistical structure: the manuscript's original claim that the p-tau217–hippocampal-ratio correlation "depends substantially on the A+/A− contrast" (Section 3.5), based on the same kind of naive significant-vs-non-significant subgroup comparison. A single pooled joint model with group-specific random-effects covariance (rather than two separately-fit subgroup models) found the A− and A+ slope-slope correlations were statistically indistinguishable (−0.44 [−0.83, 0.08] vs. −0.26 [−0.63, 0.13]; posterior P(difference>0)=0.73, well short of our pre-specified 0.975 threshold). We no longer interpret the hippocampal finding as amyloid-status-dependent, for the same reason the APOE×time×amyloid finding was revised above.

## Comment 9 — Bayesian model reporting (priors, ESS, divergences, posterior predictive checks; BLUP-normality discrepancy)

**DONE, and this uncovered a real, previously undisclosed convergence problem.** We added two new supplementary tables: **Table S6** reports the priors actually used (brms defaults throughout — flat priors on fixed effects; weakly-informative Student-t(3, ...) on intercepts; half-Student-t(3,0,...) on random-effect/residual SDs; LKJ(1) on the random-effects correlation matrix; no informative priors were set). **Table S7** reports, for each of the six primary bivariate models: max Rhat across *all* parameters (not just fixed effects), number of parameters with Rhat>1.05, divergent transitions, and posterior-predictive Bayesian p-values for the mean and SD of each submodel's response.

This full accounting revealed that the manuscript's claim of "all Rhat ≤ 1.05" (Results Section 3.4, as originally written) was **incorrect**: the three p-tau217-paired models (which drive the primary conclusions) did converge cleanly (max Rhat 1.02-1.05, 0 divergent transitions), but the three **GFAP-paired models did not fully converge** — the cross-response correlation parameter and a number of participant-level random slopes had Rhat as high as 1.36 (ADAS-Cog 13 x GFAP). Posterior predictive checks showed no mean/SD misfit in any model. We have corrected the "all Rhat ≤ 1.05" claim (Section 3.4), added the convergence caveat wherever the GFAP null findings are reported or interpreted (Section 3.4, Section 4.6 Limitations), and flagged that the GFAP-outcome null coupling should be read as provisional rather than definitively established, since a convergence problem cannot be ruled out as a partial explanation.

We also resolved the BLUP-normality discrepancy directly in the Supplementary Figure S1 legend: the main-text Shapiro-Wilk normality test concerns the **frequentist LME-derived BLUP point estimates** (used in the mediation analysis and the two-stage sensitivities S-G/S-H), while Figure S1's Q-Q plots show the **brms model's own posterior-mean random effects** — a different object entirely, extracted from a model that (for the p-tau217-paired submodels) samples the full joint posterior rather than relying on a BLUP point estimate. The two need not agree, and the legend now says so explicitly, along with the terminology fix (six separate bivariate models, not one multivariate five-response model, consistent with our response to comment 10) and the GFAP convergence caveat.

## Comment 10 — Bivariate vs. multivariate terminology; Figure/Methods model mislabeling

**DONE.** The reviewer is correct that the manuscript was internally inconsistent. Methods actually describes six *separate bivariate* Bayesian joint models (one per biomarker-outcome pair), each with its own 2×2 random-effects covariance matrix — not a single multivariate model coupling all responses. The "Phase III — Multivariate..." heading has been corrected to "Phase III — Bivariate Bayesian joint longitudinal models," Methods now explicitly states "six independent models in total, not a single multivariate model," and the Figure 3 legend (which incorrectly said the covariance structure "couples the five submodels") has been corrected to describe six separate 2×2 matrices. Separately, the Figure 2 legend attributed the APOE×time estimate to "the brms model," but this quantity is estimated in the primary LME (Phase II, lmerTest), not brms (Phase III); this has been corrected, along with the number itself, which had the same eGFR-sensitivity-vs-primary-model mix-up identified in comment 1.

---

## Summary of status

| # | Comment | Status |
|---|---|---|
| 1 | Numerical inconsistencies | **DONE** |
| 2 | 325→296 sample accounting | **DONE** |
| 3 | Mediation/BLUP construction | **DONE** (real bug found and fixed) |
| 4 | Collider concern / causal language | **DONE** |
| 5 | Amyloid measurement detail | **DONE** |
| 6 | Linearity sensitivity | **DONE** |
| 7 | APOE allele dose | **DONE** |
| 8 | Diagnosis-stage heterogeneity | **DONE** (both amyloid-status and diagnosis versions) |
| 9 | Bayesian reporting completeness | **DONE** (uncovered a real GFAP-model convergence issue, now disclosed) |
| 10 | Bivariate/multivariate terminology | **DONE** |
