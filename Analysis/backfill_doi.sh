#!/usr/bin/env bash
# Backfill Zenodo DOI into README, CITATION.cff, Main_Text after Zenodo mints it.
# Usage:  bash backfill_doi.sh  10.5281/zenodo.XXXXXXX  https://github.com/kyc0321/LOAD-ADNI-pT217
# Run from project root: /Users/kuan-yuchu/Desktop/LOAD/

set -euo pipefail

DOI="${1:-}"
REPO_URL="${2:-}"
if [[ -z "$DOI" ]]; then
  echo "Usage: bash backfill_doi.sh <DOI> [<GitHub repo URL>]"
  echo "Example A (Zenodo only): bash backfill_doi.sh 10.5281/zenodo.10000000"
  echo "Example B (GitHub+Zenodo): bash backfill_doi.sh 10.5281/zenodo.10000000 https://github.com/user/repo"
  exit 1
fi
# If REPO_URL not provided, use Zenodo record URL as the fallback
if [[ -z "$REPO_URL" ]]; then
  REPO_URL="https://doi.org/${DOI}"
  echo "No GitHub URL provided — using Zenodo URL as repository reference."
fi

cd "$(dirname "$0")/.."

# 1. README.md
sed -i.bak \
  -e "s|DOI-pending|DOI-${DOI//\//%2F}-blue|g" \
  -e "s|Zenodo DOI: \[pending\]|Zenodo DOI: ${DOI}|g" \
  -e "s|\[pending\]|${DOI}|g" \
  -e "s|https://github.com/\[username\]/LOAD|${REPO_URL}|g" \
  README.md
rm README.md.bak

# 2. CITATION.cff — insert doi key at the top if not present
if ! grep -q "^doi:" CITATION.cff; then
  awk -v doi="$DOI" '
    NR==1 {print; print "doi: " doi; next}
    {print}
  ' CITATION.cff > CITATION.cff.tmp && mv CITATION.cff.tmp CITATION.cff
fi
sed -i.bak "s|https://github.com/\[username\]/LOAD|${REPO_URL}|g" CITATION.cff
rm CITATION.cff.bak

# 3. Main_Text_v1.md — Data Availability statement
sed -i.bak \
  -e "s|\[GitHub URL TBD; Zenodo DOI TBD pending submission\]|${REPO_URL} (Zenodo DOI: ${DOI})|g" \
  -e "s|\[Zenodo DOI TBD\]|${DOI}|g" \
  Mauuscripts/Main_Text_v1.md
rm Mauuscripts/Main_Text_v1.md.bak

# 4. SAP_v1.md — references to GitHub/Zenodo
sed -i.bak \
  -e "s|https://github.com/\[username\]/LOAD|${REPO_URL}|g" \
  Mauuscripts/SAP_v1.md 2>/dev/null || true
rm -f Mauuscripts/SAP_v1.md.bak

# 5. Regenerate Main_Text_v1.docx with new DOI
python3 /tmp/make_main_text_docx.py 2>&1 || echo "Skipped DOCX regen (run python3 /tmp/make_main_text_docx.py manually)"

echo "Backfilled DOI: ${DOI}"
echo "Backfilled repo URL: ${REPO_URL}"
echo "Files updated: README.md, CITATION.cff, Mauuscripts/Main_Text_v1.md (+.docx), Mauuscripts/SAP_v1.md"
echo ""
echo "Next: git add -A && git commit -m 'Add Zenodo DOI ${DOI}' && git push"
