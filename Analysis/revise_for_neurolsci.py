# -*- coding: utf-8 -*-
"""
Revise Main_Text_NeurolSci.docx + Title Page.docx for Neurological Sciences
(Springer) resubmission, per the Editor-in-Chief's "does not comply with
submission guidelines" decision. See REVISION_NOTES.md for the full audit.
"""
import json
import copy
from docx import Document
from docx.shared import Pt, Cm
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from docx.enum.text import WD_ALIGN_PARAGRAPH

FOLDER = "/Users/kuan-yuchu/Desktop/LOAD-BC-NS/submit NeurolSci"
SRC_MAIN = f"{FOLDER}/Main_Text_NeurolSci.docx"
TABLES_SRC = "/Users/kuan-yuchu/Desktop/LOAD-BC-NS/Mauuscripts/Tables_final_v1.docx"

# =============================================================================
# 1. Load original reference list programmatically (avoid hand-transcription
#    errors) -- keyed by ORIGINAL alphabetical position, 1-47.
# =============================================================================
with open('/tmp/orig_refs.json') as f:
    orig_ref_list = json.load(f)
assert len(orig_ref_list) == 47
OLD_FULL_TEXT = {i + 1: txt for i, txt in enumerate(orig_ref_list)}

# citation-key -> OLD reference number (1-47)
OLD_REF_NUM = {
    "GBD 2019 Dementia Forecasting Collaborators, 2022": 17, "Veitch et al., 2019": 44,
    "Jack et al., 2018": 21, "Liddelow et al., 2017": 27, "Habib et al., 2020": 18,
    "Jack et al., 2024": 22, "Cicognola et al., 2021": 11, "Pereira et al., 2021": 32,
    "Benedet et al., 2021": 6, "Chatterjee et al., 2021": 10, "Mielke et al., 2022": 29,
    "Bellaver et al., 2023": 5, "Janelidze et al., 2020": 23, "Palmqvist et al., 2020": 31,
    "Barthélemy et al., 2020": 3, "Ashton et al., 2021": 1, "Mattsson-Carlgren et al., 2020": 28,
    "Ashton et al., 2024": 2, "Therriault et al., 2023": 40, "Yamazaki et al., 2019": 47,
    "Serrano-Pozo et al., 2021": 39, "Janelidze et al., 2021": 24, "Janelidze et al., 2023b": 26,
    "Cullen et al., 2021": 12, "Brum et al., 2023": 8, "Petersen et al., 2010": 33,
    "Weiner et al., 2017": 45, "Janelidze et al., 2023a": 25, "Rissin et al., 2010": 35,
    "Mohs et al., 1997": 30, "Folstein et al., 1975": 15, "Fischl, 2012": 14,
    "Inker et al., 2021": 20, "Bates et al., 2015": 4, "Proust-Lima et al., 2017": 34,
    "Bürkner, 2017": 9, "Rizopoulos, 2012": 37, "Imai et al., 2010": 19,
    "Tingley et al., 2014": 41, "Benjamini and Hochberg, 1995": 7, "Woo et al., 2025": 46,
    "Fortea et al., 2024": 16, "Vandenbroucke et al., 2007": 43, "Rubin, 1987": 38,
}
REMOVED_OLD_NUMS = {13, 36, 42}  # Farrer 1997, Rizopoulos 2024/JMbayes2, Trojanowski 2010 -- never cited

CITATION_ORDER_FULL = [
    "GBD 2019 Dementia Forecasting Collaborators, 2022", "Veitch et al., 2019",
    "Jack et al., 2018", "Liddelow et al., 2017", "Habib et al., 2020",
    "Jack et al., 2024", "Cicognola et al., 2021", "Pereira et al., 2021",
    "Benedet et al., 2021", "Chatterjee et al., 2021", "Mielke et al., 2022",
    "Bellaver et al., 2023", "Janelidze et al., 2020", "Palmqvist et al., 2020",
    "Barthélemy et al., 2020", "Ashton et al., 2021", "Mattsson-Carlgren et al., 2020",
    "Ashton et al., 2024", "Therriault et al., 2023", "Yamazaki et al., 2019",
    "Serrano-Pozo et al., 2021", "Janelidze et al., 2021", "Janelidze et al., 2023b",
    "Cullen et al., 2021", "Brum et al., 2023", "Petersen et al., 2010",
    "Weiner et al., 2017", "Janelidze et al., 2023a", "Rissin et al., 2010",
    "Mohs et al., 1997", "Folstein et al., 1975", "Fischl, 2012",
    "Inker et al., 2021", "Proust-Lima et al., 2017", "Bürkner, 2017",
    "Rizopoulos, 2012", "Rubin, 1987", "Imai et al., 2010", "Tingley et al., 2014",
    "Benjamini and Hochberg, 1995", "Bates et al., 2015", "Vandenbroucke et al., 2007",
    "Woo et al., 2025", "Fortea et al., 2024",
]
assert len(CITATION_ORDER_FULL) == 44
assert set(CITATION_ORDER_FULL) == set(OLD_REF_NUM.keys())

# Trimmed to 40 (journal limit) per user's confirmed picks 2026-08-09: cut the
# 4 least load-bearing citations (each part of a multi-citation support group,
# vs. Cullen et al. 2021 which was in the smallest group and kept).
CUT_FOR_40_LIMIT = {"Janelidze et al., 2020", "Yamazaki et al., 2019",
                     "Rubin, 1987", "Tingley et al., 2014"}
CITATION_ORDER = [k for k in CITATION_ORDER_FULL if k not in CUT_FOR_40_LIMIT]
assert len(CITATION_ORDER) == 40
assert set(OLD_FULL_TEXT.keys()) - REMOVED_OLD_NUMS == set(OLD_REF_NUM.values())

NEW_NUM = {key: i + 1 for i, key in enumerate(CITATION_ORDER)}
# Final reference list content, keyed by NEW number 1-44, sourced verbatim from original doc
FULL_REF_TEXT = {NEW_NUM[key]: OLD_FULL_TEXT[OLD_REF_NUM[key]] for key in CITATION_ORDER}

# Sanity check: print a few to confirm correct alignment
for n in [1, 16, 33, 34, 38, 40]:
    print(n, ':', FULL_REF_TEXT[n][:70])

# =============================================================================
# 2. In-text citation replacements (exact substring -> bracket numbers),
#    hand-verified against NEW_NUM above and exact original text.
# =============================================================================
REPLACEMENTS = [
    ("(GBD 2019 Dementia Forecasting Collaborators, 2022)", "[1]"),
    ("(Veitch et al., 2019; Jack et al., 2018)", "[2,3]"),
    ("(Liddelow et al., 2017; Habib et al., 2020)", "[4,5]"),
    ("(Jack et al., 2024; Cicognola et al., 2021)", "[6,7]"),
    ("(Pereira et al., 2021; Benedet et al., 2021; Chatterjee et al., 2021; Mielke et al., 2022; Bellaver et al., 2023; Janelidze et al., 2020)", "[8-12]"),  # Janelidze 2020 cut (40-ref limit)
    ("(Palmqvist et al., 2020; Barthélemy et al., 2020; Ashton et al., 2021)", "[13-15]"),
    ("(Mattsson-Carlgren et al., 2020)", "[16]"),
    ("(Ashton et al., 2024; Therriault et al., 2023; Yamazaki et al., 2019; Serrano-Pozo et al., 2021)", "[17,18,19]"),  # Yamazaki 2019 cut
    ("(Janelidze et al., 2021, 2023b; Cullen et al., 2021)", "[20,21,22]"),
    ("(Brum et al., 2023)", "[23]"),
    ("(Petersen et al., 2010; Weiner et al., 2017; Janelidze et al., 2023a)", "[24-26]"),
    ("(Rissin et al., 2010)", "[27]"),
    ("(Mohs et al., 1997)", "[28]"),
    ("(Folstein et al., 1975)", "[29]"),
    ("(Fischl, 2012)", "[30]"),
    ("(Inker et al., 2021)", "[31]"),
    ("lcmm v1.9.2; Proust-Lima et al., 2017)", "lcmm v1.9.2 [32])"),
    ("(Bürkner, 2017; Rizopoulos, 2012)", "[33,34]"),
    ("(Imai et al., 2010; Tingley et al., 2014)", "[35]"),  # Tingley 2014 cut, Imai alone
    ("q < 0.05; Benjamini and Hochberg, 1995)", "q < 0.05) [36]"),
    ("(Janelidze et al., 2021; Brum et al., 2024)", "[20,23]"),  # typo fix: text said 2024, ref list has 2023
    ("(Woo et al., 2025)", "[39]"),
    ("(Mattsson-Carlgren et al., 2020; Ashton et al., 2024; Therriault et al., 2023; Fortea et al., 2024)", "[16-18,40]"),
    ("lme4 1.1-37, Bates et al., 2015; lmerTest 3.1-3; lcmm 1.9.2; brms 2.21.0, Bürkner 2017; nlme 3.1-164; JMbayes2 0.5-1 was attempted but abandoned for lacking a survival component",
     "lme4 1.1-37 [37]; lmerTest 3.1-3; lcmm 1.9.2; brms 2.21.0 [33]; nlme 3.1-164; JMbayes2 0.5-1 was attempted but abandoned for lacking a survival component"),
    # NOTE: "Rubin's rules" citation addition reverted -- Rubin 1987 was cut for
    # the 40-ref limit, so this mention is left as an uncited named-method
    # reference (common/acceptable practice), same as the original.
    ("(Vandenbroucke et al., 2007)", "[38]"),
    ("(Cicognola et al., 2021; Benedet et al., 2021; p = 0.96 here)", "[7,9] (p = 0.96 here)"),
]

def apply_replacements(text):
    for old, new in REPLACEMENTS:
        text = text.replace(old, new)
    return text

print("\n", len(REPLACEMENTS), "replacement rules ready.")

# =============================================================================
# 3. Table data (extracted from Mauuscripts/Tables_final_v1.docx, the file
#    that was NEVER included in the actual journal submission package --
#    this is the fix for the missing-tables defect).
# =============================================================================
with open('/tmp/table_data.json') as f:
    TABLE_DATA = json.load(f)

TABLE_CAPTIONS = {
    'table1': "Table 1 Baseline characteristics of the analytic cohort (N = 325)",
    'table3': "Table 2 APOE ε4 × time interaction on plasma biomarker trajectories (primary linear mixed-effects model, N≈296 RIDs)",
    'table4': "Table 3 Multivariate Bayesian joint longitudinal model: cross-response random-slope correlations (primary analysis)",
}
# NOTE: original manuscript's "Table 1/3/4" numbering is renumbered to Table
# 1/2/3 here, since Table 2 (LCGMM) and Table 5 (cross-assay) move to
# Supplementary per the user's confirmed plan -- main text now has exactly
# 3 tables, consecutively numbered with no gaps.

def add_heading(doc, text, level):
    p = doc.add_heading(level=level)
    p.add_run(text)
    return p

def add_body_para(doc, text):
    p = doc.add_paragraph()
    p.add_run(text)
    return p

def build_table(doc, rows, caption):
    cap = doc.add_paragraph()
    r = cap.add_run(caption)
    r.bold = True
    ncol = len(rows[0])
    t = doc.add_table(rows=1, cols=ncol)
    t.style = "Table Grid"
    hdr = t.rows[0].cells
    for i, val in enumerate(rows[0]):
        hdr[i].text = ""
        p = hdr[i].paragraphs[0]
        rn = p.add_run(val)
        rn.bold = True
        rn.font.size = Pt(8.5)
    for row in rows[1:]:
        cells = t.add_row().cells
        for i, val in enumerate(row):
            cells[i].text = ""
            p = cells[i].paragraphs[0]
            rn = p.add_run(val)
            rn.font.size = Pt(8.5)
    doc.add_paragraph()
    return t

print("Table/heading helpers ready.")

# =============================================================================
# 4. Assemble the revised, blinded main text
# =============================================================================
with open('/tmp/full_paras.json') as f:
    P = json.load(f)  # original paragraphs, indexable exactly as audited

doc = Document()
style = doc.styles["Normal"]
style.font.name = "Times New Roman"
style.font.size = Pt(10)
sec = doc.sections[0]
sec.left_margin = sec.right_margin = Cm(2.54)
sec.top_margin = sec.bottom_margin = Cm(2.54)

# --- Title + running title (no author identity: blinded per double-blind review) ---
tp = doc.add_paragraph()
tr = tp.add_run(P[0])
tr.bold = True
tr.font.size = Pt(13)
rt = doc.add_paragraph()
rt.add_run(P[1]).italic = True
doc.add_paragraph()  # spacer

# --- Abstract / Plain Language Summary / Abbreviations (unchanged, no citations here) ---
add_heading(doc, "Abstract", 1)
for i in [21, 22, 23, 24]:
    add_body_para(doc, P[i])
add_body_para(doc, P[25])  # Keywords

add_heading(doc, "Plain Language Summary", 1)
for i in [27, 28, 29]:
    add_body_para(doc, P[i])

add_heading(doc, "Abbreviations", 1)
add_body_para(doc, P[32])

# --- 1. Introduction ---
add_heading(doc, "1. Introduction", 1)
for i in [36, 37, 38]:
    add_body_para(doc, apply_replacements(P[i]))

# --- 2. Methods ---
add_heading(doc, "2. Methods", 1)
add_heading(doc, "2.1 Data source, population, and measurements", 2)
for i in [43, 44]:
    add_body_para(doc, apply_replacements(P[i]))

add_heading(doc, "2.2 Covariates", 2)
add_body_para(doc, apply_replacements(P[47]))

add_heading(doc, "2.3 Statistical analysis", 2)
add_heading(doc, P[50], 3)  # "Phase I -- LCGMM exploratory diagnostic"
add_body_para(doc, apply_replacements(P[51]))
add_heading(doc, P[52], 3)  # "Phase II -- Primary linear mixed-effects model"
for i in [53, 54, 56]:
    add_body_para(doc, apply_replacements(P[i]))
add_heading(doc, P[57], 3)  # "Phase III -- Multivariate Bayesian joint longitudinal model"
for i in [58, 59, 60, 61]:
    add_body_para(doc, apply_replacements(P[i]))
add_heading(doc, P[62], 3)  # "Phase IV -- Formal mediation analysis"
add_body_para(doc, apply_replacements(P[63]))
add_heading(doc, P[64], 3)  # "Multiple-testing correction"
add_body_para(doc, apply_replacements(P[65]))
add_heading(doc, P[66], 3)  # "Reproducibility"
add_body_para(doc, apply_replacements(P[67]))
# Relocated from the original "Declarations" section (STROBE + SAP statements
# reference Vandenbroucke et al. 2007, which needs to stay a genuine in-text
# citation -- these do not fit the title page's defined "Statements and
# Declarations" categories, so they stay in Methods, not on the Title Page).
add_body_para(doc, apply_replacements(P[135]))
add_body_para(doc, P[136])

# --- 3. Results ---
add_heading(doc, "3. Results", 1)
add_heading(doc, "3.1 Study population", 2)
add_body_para(doc, P[72])
build_table(doc, TABLE_DATA['table1'], TABLE_CAPTIONS['table1'])
for i in [73, 74]:
    add_body_para(doc, P[i])

add_heading(doc, "3.2 Phase I — LCGMM exploratory trajectory class diagnostic", 2)
add_body_para(doc, "LCGMM model selection details are reported in Supplementary Table S3 (moved out of the main text to keep total main-text illustrations within the journal's 4-6 limit).")
for i in [77]:
    add_body_para(doc, P[i])
add_body_para(doc, P[78].replace("(Figure 5)", "(Supplementary Figure S3)"))

add_heading(doc, "3.3 Phase II — APOE ε4 modulation of plasma biomarker slopes", 2)
add_body_para(doc, P[80].replace("Table 3", "Table 2"))
build_table(doc, TABLE_DATA['table3'], TABLE_CAPTIONS['table3'])

add_heading(doc, "3.4 Phase III — Joint individual-level association between biomarker and outcome trajectories", 2)
add_body_para(doc, P[82].replace("Table 4", "Table 3"))
build_table(doc, TABLE_DATA['table4'], TABLE_CAPTIONS['table4'])
add_body_para(doc, P[83].replace("Supplementary Tables S5–S6", "Supplementary Table S4 (parts c-d)"))

add_heading(doc, "3.5 Amyloid-stratified sensitivity analysis", 2)
for i in [85, 86]:
    add_body_para(doc, apply_replacements(P[i]))
add_body_para(doc, apply_replacements(P[87]).replace("Supplementary Table S3", "Supplementary Table S4a"))
add_heading(doc, P[88], 3)  # "3.5a Formal causal mediation analysis"
for i in [89, 90]:
    add_body_para(doc, P[i])
add_heading(doc, P[91], 3)  # "3.5b Post-baseline-only sensitivity"
add_body_para(doc, P[92].replace("Supplementary Table S4)", "Supplementary Table S4b)"))
add_heading(doc, P[93], 3)  # "3.5c Statistical diagnostics"
add_body_para(doc, P[94])

add_heading(doc, "3.6 Cross-assay plasma p-tau217 concordance", 2)
add_body_para(doc, P[96].replace("Table 5 and Figure 4", "Supplementary Table S5 and Supplementary Figure S2"))
add_body_para(doc, P[97])

# --- 4. Discussion ---
add_heading(doc, "4. Discussion", 1)
for hidx, bidxs in [
    (101, [102, 103, 104, 105, 106]),
    (107, [108, 109, 110]),
    (111, [112, 113, 114]),
    (115, [116]),
    (117, [118]),
    (119, [120, 121, 122]),
    (123, [124]),
    (125, [126]),
]:
    add_heading(doc, P[hidx], 2)
    for b in bidxs:
        add_body_para(doc, apply_replacements(P[b]))

# --- Figure Legends (main-text figures only: 1-3; Figs 4-5 move to Supplementary) ---
add_heading(doc, "Figure Legends", 1)
for i in [139, 140, 141]:
    txt = P[i]
    m = txt.split(". ", 1)
    fig_no = m[0].replace("Figure ", "")
    rest = m[1] if len(m) > 1 else ""
    # strip a trailing period per guideline ("no punctuation... at the end of the caption")
    if rest.endswith("."):
        rest = rest[:-1]
    p = doc.add_paragraph()
    r = p.add_run(f"Fig. {fig_no} ")
    r.bold = True
    p.add_run(rest)

# --- References (renumbered, consecutive) ---
add_heading(doc, "References", 1)
for n in range(1, 41):
    p = doc.add_paragraph()
    p.add_run(f"[{n}] {FULL_REF_TEXT[n]}")

out_path = f"{FOLDER}/Main_Text_NeurolSci_R1.docx"
doc.save(out_path)
print("\nSaved:", out_path)
