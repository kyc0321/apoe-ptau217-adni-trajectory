# GitHub + Zenodo Publication — Step-by-Step

This document walks through publishing this repository on GitHub and minting a citable Zenodo DOI. Estimated total time: **30–45 minutes**.

The repository has been **initialised locally** (`.git/` exists) but **nothing has been committed yet**. The first commit, push, and Zenodo linkage are user-controlled to give you final review.

---

## Step 1 — Pre-commit verification (5 min)

Before the first commit, verify the safety net:

```bash
cd /Users/kuan-yuchu/Desktop/LOAD

# 1.1 Confirm .gitignore is excluding raw ADNI data
git status --ignored | grep ADNI_csv
# Expected: "Ignored: ADNI_csv/" — must appear under Ignored, NOT Untracked

# 1.2 Confirm no .rds fit objects will be committed (they're large; Zenodo gets them)
git status | grep '\.rds'
# Expected: empty output

# 1.3 Preview what WILL be committed
git status
# Review the list; should contain README.md, LICENSE.md, CITATION.cff,
# Analysis/*.R, Analysis/*.md, Analysis/outputs/*.csv, Mauuscripts/*.md,
# Mauuscripts/*.docx (the manuscript drafts), figures/v4/*.png, *.pdf
```

If anything looks wrong, fix `.gitignore` before proceeding.

---

## Step 2 — First commit (2 min)

```bash
# 2.1 Configure your git identity (one-time, if not already set globally)
git config user.name  "Kuan-Yu Chu"
git config user.email "kyc0321@gmail.com"

# 2.2 Stage everything not ignored
git add .

# 2.3 First commit
git commit -m "Initial release: ADNI plasma GFAP/p-tau217 trajectory analysis

- N = 325 analytic cohort; ≥2/≥2 measurement threshold
- Phase II LME: APOE × time on plasma p-tau217 (p = 0.0023), null on GFAP
- Phase III brms multivariate joint LME: p-tau217 slope ρ = +0.67/-0.69/-0.39
  with ADAS-Cog 13, MMSE, hippocampal ratio (all CrI exclude zero)
- GFAP slope null across all outcomes
- SAP v1 locked 2026-02-12 with 10 documented amendments dated 2026-05-25
- Single-author submission

H1 and H2 directionally falsified by completed analysis; reported transparently."

# 2.4 (Optional) Tag the submission version
git tag -a v1.0-submission -m "First journal submission, 2026-05-25"
```

---

## Step 3 — Create the GitHub repository (5 min)

### Option A — Using GitHub CLI (`gh`), recommended

```bash
# 3A.1 Install gh if not present
brew install gh        # macOS

# 3A.2 Authenticate (opens browser)
gh auth login

# 3A.3 Create the repo and push (substitute desired name)
gh repo create LOAD-ADNI-pT217 \
  --public \
  --description "APOE ε4 and plasma p-tau217 trajectory in ADNI older adults — reproducibility repository" \
  --source=. \
  --push
```

### Option B — Using the GitHub web UI

1. Go to https://github.com/new
2. Repository name: `LOAD-ADNI-pT217` (or your choice; the name affects the citation URL)
3. Description: *APOE ε4 and plasma p-tau217 trajectory in ADNI older adults — reproducibility repository*
4. Visibility: **Public**
5. **Do NOT** initialise with README, .gitignore, or LICENSE (we already have them).
6. Click "Create repository". GitHub will display the push instructions.
7. From your terminal:
   ```bash
   cd /Users/kuan-yuchu/Desktop/LOAD
   git remote add origin https://github.com/[your-username]/LOAD-ADNI-pT217.git
   git branch -M main
   git push -u origin main
   git push --tags         # if you created v1.0-submission
   ```

After this step, the repository is live at `https://github.com/[your-username]/LOAD-ADNI-pT217`.

---

## Step 4 — Link Zenodo (10 min)

Zenodo gives every GitHub release a **citable DOI**. This is the value you put in the manuscript Data Availability statement.

1. Go to https://zenodo.org and log in with your GitHub account (top-right "Log in").
2. Click your username (top-right) → **GitHub**.
3. Find `LOAD-ADNI-pT217` in the list. Toggle the switch to **ON**.
4. Back on GitHub, go to the repository page → **Releases** → **Create a new release**:
   - Tag: `v1.0` (or `v1.0-submission`)
   - Title: `v1.0 — Initial journal submission`
   - Description: (paste the same text as your first-commit message, or summarise the key findings)
   - Click **Publish release**
5. Zenodo automatically receives a webhook from GitHub, archives the release, and mints a DOI. **Wait 1–2 minutes**, then refresh the Zenodo dashboard.
6. Click the new entry to see the DOI (looks like `10.5281/zenodo.XXXXXXX`).

---

## Step 5 — Upload the large `.rds` brms fit objects to Zenodo (5 min)

These are excluded from the GitHub repo (see `.gitignore`) because each is 20–50 MB. Archive them on the same Zenodo record:

1. On the Zenodo record page (created in Step 4), click **Edit**.
2. Drag and drop all files in `Analysis/outputs/*.rds` (9 files total: 6 from `phase3_brms_joint_v1.R`, 3 from the GFAP refit).
3. Also upload `Analysis/outputs/*.RData` if any.
4. Click **Publish** to update the record. The DOI is unchanged; this is the same versioned record.

---

## Step 6 — Update the manuscript with the DOI (5 min)

Once the DOI is minted, update three places:

1. **README.md** — replace the badge `[Zenodo DOI pending]` and the citation block:
   ```bash
   sed -i '' 's/DOI: \[pending\]/DOI: 10.5281\/zenodo.XXXXXXX/g' README.md
   sed -i '' 's/Zenodo DOI: \[pending\]/Zenodo DOI: 10.5281\/zenodo.XXXXXXX/g' README.md
   ```

2. **CITATION.cff** — add the `doi:` field at the top:
   ```yaml
   doi: 10.5281/zenodo.XXXXXXX
   ```

3. **Manuscript Data Availability statement** — in `Mauuscripts/Main_Text_v1.md` and `Main_Text_v1.docx`, replace `[Zenodo DOI TBD]` with the real DOI.

4. Commit the changes and push:
   ```bash
   git add README.md CITATION.cff Mauuscripts/Main_Text_v1.md Mauuscripts/Main_Text_v1.docx
   git commit -m "Add Zenodo DOI 10.5281/zenodo.XXXXXXX to citation metadata"
   git push
   ```

---

## Step 7 — Final pre-submission audit (5 min)

```bash
# 7.1 Repository looks healthy
git log --oneline -5
git tag

# 7.2 Public URL renders correctly (open in browser)
open "https://github.com/[your-username]/LOAD-ADNI-pT217"

# 7.3 Zenodo record renders correctly
open "https://doi.org/10.5281/zenodo.XXXXXXX"

# 7.4 Manuscript draft references the DOI in three places (README, CITATION.cff, Main_Text)
grep -r "zenodo" README.md CITATION.cff Mauuscripts/
```

After this, the public repository and DOI are ready to be cited in the manuscript's Cover Letter and Data Availability statement.

---

## Troubleshooting

- **"Refusing to push: ADNI_csv files detected"** — your `.gitignore` is not working. Run `git ls-files | grep ADNI_csv` to find which files are tracked, then `git rm --cached <file>` for each, then re-commit.
- **Zenodo did not receive the release** — confirm the GitHub-Zenodo toggle is ON for the repo (Step 4.3); confirm the release is *published* (not a draft).
- **`.rds` files too large for Zenodo individual file upload** — Zenodo's per-file limit is 50 GB; should be fine for single brms fits.
- **Want to release a v1.1 after revision** — create a new GitHub release with a new tag (e.g., `v1.1-revised`); Zenodo will mint a new DOI for that version and link it to the parent record.

---

## What this repository is, and what it is not

**It is:**
- A reproducible analysis pipeline that, given ADNI data access, regenerates every number in the manuscript
- A pre-specified Statistical Analysis Plan with a fully transparent amendment history
- An audit trail of analytical decisions

**It is not:**
- A data repository — raw ADNI data is governed by the ADNI DUA and is not redistributed
- A substitute for peer review — the self-review log documents the author's own audit trail, not external review
- A claim to clinical or regulatory application — the findings require external replication before clinical translation
