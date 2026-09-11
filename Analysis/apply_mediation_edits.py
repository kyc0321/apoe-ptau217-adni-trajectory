import copy
from docx import Document
from docx.oxml.ns import qn

SRC = "/Users/kuan-yuchu/Desktop/LOAD-BC/submit BC 5:31:26/sumit 6:3 26/Main_Text.docx"
OUT = "/Users/kuan-yuchu/Desktop/LOAD-BC/submit BC 5:31:26/sumit 6:3 26/Main_Text.docx"

d = Document(SRC)

def insert_para_before(anchor_para, text, bold_header=False, style="Normal"):
    """Insert a new paragraph before anchor_para. text can be a string (single run)
    or a list of (text, bold) tuples for mixed runs."""
    new_p = anchor_para.insert_paragraph_before("", style=style)
    if bold_header:
        r = new_p.add_run(text)
        r.bold = True
    else:
        for chunk, bold in text:
            r = new_p.add_run(chunk)
            r.bold = bold
    return new_p

# ---------------------------------------------------------------------------
# 1. Methods (Sec 2.4): add "Phase IV" mediation methodology paragraph before
#    the "Multiple-testing correction" header (currently index 91).
# ---------------------------------------------------------------------------
anchor = d.paragraphs[91]
assert anchor.text.strip() == "Multiple-testing correction"

insert_para_before(anchor, "Phase IV — Formal mediation analysis", bold_header=True)
insert_para_before(anchor, [
    ("To decompose the total APOE ε4 effect on plasma p-tau217 kinetics into amyloid-mediated "
     "(indirect) and direct components, rather than inferring mediation from subgroup attenuation "
     "alone, we additionally conducted a counterfactual-based causal mediation analysis (Imai, "
     "Keele & Tingley, 2010) using the R ", False),
    ("mediation", False),
    (" package (Tingley et al., 2014). The mediator model regressed baseline amyloid positivity "
     "(A+/A−, as defined in §3.5) on APOE4 and the Tier 1/Tier 2 covariates (logistic regression); "
     "the outcome model regressed the per-participant best linear unbiased predictor (BLUP) slope "
     "of log plasma p-tau217 — extracted from the primary Phase II linear mixed-effects model — on "
     "APOE4, A+ status, and the same covariates (linear regression). The average causal mediation "
     "effect (ACME), average direct effect (ADE), and total effect were estimated with a "
     "nonparametric bootstrap (BCa, 5,000 resamples; n = 296 participants with both a BLUP slope "
     "and a baseline plasma Aβ42/40 measurement).", False)
])
insert_para_before(anchor, "", style="Normal")

# ---------------------------------------------------------------------------
# 2. Results: new §3.5a subsection before "3.5b Post-baseline-only sensitivity"
#    (currently index 136 -- re-locate by text since indices shifted after
#    the Methods insert above).
# ---------------------------------------------------------------------------
anchor_35b = next(p for p in d.paragraphs if p.text.strip() == "3.5b Post-baseline-only sensitivity")

insert_para_before(anchor_35b, "3.5a Formal causal mediation analysis", bold_header=True)
insert_para_before(anchor_35b, [
    ("The indirect (amyloid-mediated) pathway was significant and positive (ACME = 0.00196 "
     "log-units/year, 95% CI 0.00019–0.00411, ", False), ("p = 0.045", True), (
     "), whereas the direct pathway (APOE4 on p-tau217 slope, holding amyloid status fixed) was "
     "small, negative, and non-significant (ADE = −0.00193, 95% CI −0.00706–0.00361, ", False),
    ("p = 0.468", True), (
     "). The total effect was correspondingly indistinguishable from zero (0.00003, 95% CI "
     "−0.00481–0.00454, p = 0.937). This pattern — a significant indirect effect coexisting with an "
     "oppositely signed, non-significant direct effect that nets to a null total effect — is a "
     "recognized signature of inconsistent (suppressive) mediation, not evidence of no association. "
     "Under inconsistent mediation the conventional “proportion mediated” statistic is "
     "mathematically unstable, because it divides by a total effect indistinguishable from zero "
     "(here, 95% CI −707 to 0.23); it is therefore not reported as a point estimate.", False)
])
insert_para_before(anchor_35b, [
    ("The formal mediation analysis is consistent with, and statistically formalizes, the "
     "amyloid-stratified attenuation pattern reported above: the amyloid-mediated pathway from "
     "APOE ε4 to plasma p-tau217 kinetics is statistically supported, while a direct, "
     "amyloid-independent pathway is not. As with the two-stage BLUP-based sensitivity analyses "
     "(S-G, S-H; §2.4, §3.4), this analysis uses the per-participant slope BLUP as an outcome and "
     "therefore inherits the same BLUP shrinkage/measurement-error property discussed in §3.5c; it "
     "is presented as a formal complement to, not a replacement for, the raw-data amyloid-stratified "
     "joint model reported above.", False)
])
insert_para_before(anchor_35b, "", style="Normal")

# ---------------------------------------------------------------------------
# 3. Discussion §4.3: append a sentence to the "upstream amyloid confounding"
#    paragraph referencing the formal mediation result.
# ---------------------------------------------------------------------------
p43 = next(p for p in d.paragraphs if p.text.startswith("The first and most strongly supported is upstream amyloid confounding"))
r = p43.add_run(
    " A formal counterfactual mediation analysis (§3.5a) corroborates this reading: the "
    "amyloid-mediated indirect effect was significant (ACME p = 0.045) while the direct effect was "
    "not (ADE p = 0.468), a pattern consistent with inconsistent/suppressive mediation rather than "
    "with the absence of an amyloid-dependent pathway."
)

# ---------------------------------------------------------------------------
# 4. Limitations §4.6: revise the "Mediator and biomarker limitations" para to
#    reflect that the formal mediation analysis has now been performed.
# ---------------------------------------------------------------------------
p46 = next(p for p in d.paragraphs if p.text.startswith("Mediator and biomarker limitations."))
old_tail = (
    "and the amyloid-stratified attenuation pattern in §3.5, although consistent with amyloid "
    "mediation, does not formally decompose the total APOE × p-tau217 effect into direct and "
    "indirect components. Future work should fit a causal mediation model — for example, a "
    "counterfactual-based decomposition with amyloid PET as the measured mediator — to quantify "
    "the proportion mediated and to assess the no-unmeasured-confounding assumption."
)
new_tail = (
    "A formal counterfactual mediation analysis (§3.5a) decomposed the total APOE × p-tau217 "
    "effect and found a significant amyloid-mediated indirect pathway (ACME p = 0.045) with no "
    "significant direct pathway (ADE p = 0.468). However, because amyloid status was ascertained via "
    "a plasma ratio rather than amyloid PET, and because the mediation model's outcome (a "
    "per-participant slope BLUP) is subject to shrinkage, this remains a plasma-based, "
    "BLUP-dependent decomposition rather than a PET-anchored one; a causal mediation analysis using "
    "amyloid PET as the measured mediator remains an important extension for future, prospectively "
    "designed cohorts."
)
full_text = p46.text
assert old_tail in full_text, "old_tail not found verbatim -- check for smart-quote/dash mismatches"
new_full = full_text.replace(old_tail, new_tail)
# Rebuild paragraph runs: keep the first run's formatting, replace all run text
for r_ in list(p46.runs):
    r_.text = ""
p46.runs[0].text = new_full

# ---------------------------------------------------------------------------
# 5. Abstract: mention the formal mediation analysis in Methods + Results.
# ---------------------------------------------------------------------------
p_methods = d.paragraphs[15]
assert p_methods.text.startswith("Methods.")
p_methods.runs[-1].text += (" A formal counterfactual mediation analysis (bootstrap ACME/ADE "
                             "decomposition) was conducted to quantify direct versus "
                             "amyloid-mediated pathways.")

p_results = d.paragraphs[17]
assert p_results.text.startswith("Results.")
p_results.runs[-1].text += (" A formal mediation analysis confirmed a significant amyloid-mediated "
                             "indirect effect (ACME = 0.0020, 95% CI 0.0002-0.0041, p = 0.045) with "
                             "no significant direct effect (ADE p = 0.47), consistent with "
                             "inconsistent/suppressive mediation rather than an absent pathway.")

# ---------------------------------------------------------------------------
# 5b. Fix two pre-existing paragraph-merge artifacts in the References list:
#     refs 29/30 and refs 44/45 were each fused into a single paragraph
#     (a docx-conversion artifact), which would read as malformed/missing
#     references to a copyeditor. Split each into two proper paragraphs.
# ---------------------------------------------------------------------------
def split_merged_ref(marker_text, split_at):
    p = next(pp for pp in d.paragraphs if pp.text.strip() == marker_text)
    first_half, second_half = p.text.split(split_at, 1)
    second_half = split_at + second_half
    for r_ in list(p.runs):
        r_.text = ""
    p.runs[0].text = first_half.strip()
    new_p = p.insert_paragraph_before("")   # temp paragraph placed BEFORE p
    new_p.add_run(second_half.strip())
    p._p.addprevious(new_p._p)              # no-op reorder safeguard (already before)
    d_body = p._p.getparent()
    d_body.remove(new_p._p)
    p._p.addnext(new_p._p)                  # move it to AFTER p
    return p, new_p

p29, p30 = split_merged_ref(
    "29. Trojanowski JQ, Vandeerstichele H, Korecka M, et al. Update on the biomarker core of the "
    "Alzheimer's Disease Neuroimaging Initiative subjects. Alzheimers Dement. 2010;6(3):230–238. "
    "30. Rissin DM, Kan CW, Campbell TG, et al. Single-molecule enzyme-linked immunosorbent assay "
    "detects serum proteins at subfemtomolar concentrations. Nat Biotechnol. 2010;28(6):595–599. "
    "doi:10.1038/nbt.1641",
    "30. Rissin"
)

p44, p45 = split_merged_ref(
    "44. Vandenbroucke JP, von Elm E, Altman DG, et al. STROBE: explanation and elaboration. "
    "PLoS Med. 2007;4(10):e297. 45. Rubin DB. Multiple Imputation for Nonresponse in Surveys. "
    "New York: John Wiley & Sons; 1987.",
    "45. Rubin"
)

# ---------------------------------------------------------------------------
# 6. References: append the two mediation-methodology citations (renumbered
#    in a later consolidated pass alongside the Woo et al. addition).
# ---------------------------------------------------------------------------
ref_paras = [p for p in d.paragraphs if p.text.strip().startswith("45. Rubin DB")]
assert len(ref_paras) == 1, f"expected 1 match, got {len(ref_paras)}"
last_ref = ref_paras[0]
new_r1 = last_ref.insert_paragraph_before("")
new_r1.add_run(
    "46. Imai K, Keele L, Tingley D. A general approach to causal mediation analysis. "
    "Psychol Methods. 2010;15(4):309-334. doi:10.1037/a0020761"
)
new_r2 = last_ref.insert_paragraph_before("")
new_r2.add_run(
    "47. Tingley D, Yamamoto T, Hirose K, Keele L, Imai K. mediation: R package for causal "
    "mediation analysis. J Stat Softw. 2014;59(5):1-38. doi:10.18637/jss.v059.i05"
)

d.save(OUT)
print("Saved edits to", OUT)
