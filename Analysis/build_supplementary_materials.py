# -*- coding: utf-8 -*-
"""
Build a single, unified, correctly-numbered Supplementary Materials document,
replacing two inconsistent pre-existing files (Tables_final_v1.docx's S2/S3
labels don't match what the main text actually calls them; a much older
LOAD_Trajectories_Supplementary_v4.docx uses an unrelated S1-S11 scheme from
a stale draft with a different title/author-credential string and is not
blinded). Numbering here is anchored to what Main_Text_NeurolSci_R1.docx
ALREADY correctly references (Table S1=hippocampal slope, Table S2=VIF,
Figure S1=BLUP QQ) and extends forward for newly-demoted content.
"""
import json
from docx import Document
from docx.shared import Pt, Cm
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

FOLDER = "/Users/kuan-yuchu/Desktop/LOAD-BC-NS/submit NeurolSci"

with open('/tmp/supp_table_data.json') as f:
    T = json.load(f)

doc = Document()
style = doc.styles["Normal"]
style.font.name = "Times New Roman"
style.font.size = Pt(10)
sec = doc.sections[0]
sec.left_margin = sec.right_margin = Cm(2.54)
sec.top_margin = sec.bottom_margin = Cm(2.54)


def h(text, size=13, bold=True):
    p = doc.add_paragraph()
    r = p.add_run(text)
    r.bold = bold
    r.font.size = Pt(size)
    return p


def body(text):
    p = doc.add_paragraph()
    p.add_run(text)
    return p


def caption(text):
    p = doc.add_paragraph()
    r = p.add_run(text)
    r.bold = True
    return p


def table(rows, header_bold=True):
    ncol = len(rows[0])
    t = doc.add_table(rows=1, cols=ncol)
    t.style = "Table Grid"
    hdr = t.rows[0].cells
    for i, val in enumerate(rows[0]):
        hdr[i].text = ""
        p = hdr[i].paragraphs[0]
        r = p.add_run(val)
        r.bold = header_bold
        r.font.size = Pt(8)
    for row in rows[1:]:
        cells = t.add_row().cells
        for i, val in enumerate(row):
            cells[i].text = ""
            p = cells[i].paragraphs[0]
            r = p.add_run(val)
            r.font.size = Pt(8)
    doc.add_paragraph()
    return t


# --- Header (blinded: title + running title only, no author identity) ---
h("Supplementary Materials", 14)
h("APOE ε4 Acts on Plasma p-tau217 Trajectory Primarily via Amyloid Burden: "
  "A Bayesian Joint Longitudinal Analysis of ADNI", 12)
doc.add_paragraph()

h("Contents", 12)
for line in [
    "Supplementary Table S1. Per-diagnosis hippocampal volume slope",
    "Supplementary Table S2. Collinearity check (VIF) for Phase II fixed effects",
    "Supplementary Table S3. LCGMM model selection (exploratory diagnostic, full reporting)",
    "Supplementary Table S4. Phase III sensitivity analyses (four variants)",
    "Supplementary Table S5. Cross-assay concordance for plasma p-tau217 at baseline",
    "Supplementary Figure S1. BLUP normal Q-Q plots",
    "Supplementary Figure S2. Cross-assay consistency of the APOE × time estimate",
    "Supplementary Figure S3. LCGMM two-class trajectory solution",
]:
    p = doc.add_paragraph(style="List Bullet")
    p.add_run(line)
doc.add_page_break()

# --- Table S1: hippocampal volume slope ---
caption("Supplementary Table S1. Per-diagnosis hippocampal volume slope "
        "(biological plausibility check).")
table(T['hippo'])
body("Published expectations: CN −50 to −150; MCI −200 to −400; AD −300 to −600 mm³/year. "
     "Referenced in the main text, Results 3.1.")
doc.add_paragraph()

# --- Table S2: VIF ---
caption("Supplementary Table S2. Collinearity check (variance inflation factors) "
        "for Phase II fixed effects.")
table(T['vif'])
body("All VIFs < 5 (conventional threshold). Max VIF = 1.31. Referenced in the main "
     "text, Results 3.5c.")
doc.add_paragraph()

# --- Table S3: LCGMM model selection ---
caption("Supplementary Table S3. LCGMM model selection (exploratory diagnostic, "
        "full reporting).")
table(T['lcgmm'])
body("Bootstrap LRT (k=2 vs k=1, B=200): observed LRT = 32.64, p = 0 (200 valid "
     "replicates). For plasma p-tau217, all four pre-specified strong-support criteria "
     "are met for the k=2 solution. Referenced in the main text, Results 3.2 (this table "
     "was moved out of the main text to keep total main-text illustrations within the "
     "journal's 4-6 limit).")
doc.add_paragraph()

# --- Table S4a-d: Phase III sensitivity variants (was "Table 4 -- Sensitivity...") ---
caption("Supplementary Table S4. Phase III sensitivity analyses (four variants).")
body("Referenced in the main text, Results 3.4 (moved out of the main text along with "
     "Table S3 to keep total main-text illustrations within the journal's 4-6 limit).")
sub_defs = [
    ("S4a", "Amyloid-positive subset (plasma Aβ42/40 < 0.085; n_id = 143-149).",
     "Mediation check addressing upstream amyloid confounding. Cognitive ρ values "
     "attenuate but persist; hippocampal ρ attenuates to non-significance.", 'sens_a'),
    ("S4b", "Post-baseline-only outcome window (time_yrs ≥ 0; n_id = 167-269).",
     "Addresses time-frame asymmetry between biomarker and outcome trajectories. "
     "Cognitive ρ values essentially unchanged; hippocampal ρ attenuates to "
     "non-significance.", 'sens_b'),
    ("S4c", "Sensitivity S-G: Naive two-stage (BLUP-as-predictor, no measurement-error "
     "correction).",
     "This sensitivity inflates significance because it ignores measurement error in "
     "the BLUP slope. Reported only for transparency.", 'sens_c'),
    ("S4d", "Sensitivity S-H: MI measurement-error correction of two-stage "
     "(M=100, Rubin's pooling).",
     "Conservative comparison; treats BLUP slope as fixed-but-noisy predictor and "
     "cannot share information through random-effects covariance.", 'sens_d'),
]
for label, subcap, note, key in sub_defs:
    caption(f"Table {label}. {subcap}")
    table(T[key])
    body(note)
    doc.add_paragraph()

# --- Table S5: cross-assay concordance ---
caption("Supplementary Table S5. Cross-assay concordance for plasma p-tau217 at "
        "baseline (analytic cohort).")
table(T['crossassay'])
body("Only the C2N comparison has adequate baseline overlap (n=208) within the "
     "analytic cohort; Lilly (n=23) and Janssen (n=14) are reported but underpowered. "
     "Referenced in the main text, Results 3.6 (moved out of the main text to keep "
     "total main-text illustrations within the journal's 4-6 limit).")
doc.add_paragraph()
doc.add_page_break()

# --- Figure legends S1-S3 (images are separate uploaded files, per the journal's
#     accepted fallback for large figures; filenames noted below) ---
h("Supplementary Figure Legends", 12)

caption("Supplementary Figure S1. Best linear unbiased predictor (BLUP) normal Q-Q "
        "plots. [Image file: SuppFig_BLUP_QQ.tiff]")
body("Q-Q plots of the BLUP-extracted random intercepts and random slopes for each of "
     "the five longitudinal outcomes in the primary brms joint model (plasma GFAP, "
     "plasma p-tau217, ADAS-Cog 13, MMSE, hippocampal volume). Points follow the "
     "reference line closely across all submodels, supporting the multivariate "
     "normality assumption of the random-effects covariance and the validity of the "
     "multivariate Gaussian likelihood.")
doc.add_paragraph()

caption("Supplementary Figure S2. Cross-assay consistency of plasma p-tau217 "
        "trajectory APOE × time estimates. [Image file: Figure4_crossassay.tiff -- "
        "was main-text Figure 4 in the original submission]")
body("Forest plot of the APOE × time interaction coefficient on plasma p-tau217 from "
     "the primary brms joint model alongside two pre-specified sensitivity "
     "comparators: (i) a frequentist linear mixed-effects model and (ii) a "
     "multiple-imputation-on-BLUP analysis. All three approaches return directionally "
     "concordant estimates, with overlapping confidence/credible intervals supporting "
     "the robustness of the primary brms result. Moved out of the main text to keep "
     "total main-text illustrations within the journal's 4-6 limit.")
doc.add_paragraph()

caption("Supplementary Figure S3. Latent Class Growth Mixed Model (LCGMM) of plasma "
        "p-tau217 -- two-class solution. [Image file: Figure5_LCGMM_pT217_classes.tiff "
        "-- was main-text Figure 5 in the original submission]")
body("Group-mean trajectories for the two-class LCGMM solution selected by BIC and "
     "bootstrap likelihood-ratio test. The steep-slope class (n = 102) had a mean "
     "estimated slope of +0.072 log-units/year and contained 78.4% amyloid-positive "
     "participants at baseline; the gentle-slope class (n = 223) had a mean slope of "
     "+0.011 log-units/year and contained 29.4% amyloid-positive participants. The "
     "unsupervised latent-class structure thus recapitulates the amyloid-mediation "
     "pattern documented in main-text §3.5. Entropy = 0.82. Moved out of the main text "
     "to keep total main-text illustrations within the journal's 4-6 limit.")

out = f"{FOLDER}/Supplementary_Materials_R1.docx"
doc.save(out)
print("Saved:", out)

# verify blinding
full = "\n".join(p.text for p in doc.paragraphs)
for leak in ["Kuan-Yu Chu", "kyc0321", "MPH", "MHA"]:
    print(f"  {leak!r} present:", leak in full)
print("Tables:", len(doc.tables))
