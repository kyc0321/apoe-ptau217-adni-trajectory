# Inclusion-criterion sensitivity: N at >=2 / >=3 / >=4 measurements
# Author: Kuan-Yu Chu | 2026-05-25
# Purpose: inform whether to relax >=3 measurement threshold to improve N

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
})
read_csv  <- function(f, ...) utils::read.csv(f, stringsAsFactors = FALSE, na.strings = c("", "NA"))
write_csv <- function(x, f)   utils::write.csv(x, f, row.names = FALSE)

set.seed(20260525)

DATA  <- "/Users/kuan-yuchu/Desktop/LOAD/ADNI_csv"
OUT   <- "/Users/kuan-yuchu/Desktop/LOAD/Analysis/outputs"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- Load ADNIMERGE (for age, dx, education at baseline) ----
load("/tmp/adnimerge_unpack/ADNIMERGE/data/adnimerge.rdata")  # provides 'adnimerge'
load("/tmp/adnimerge_unpack/ADNIMERGE/data/ptdemog.rdata")    # provides 'ptdemog'

# ---- Load plasma biomarkers ----
upenn <- read_csv(file.path(DATA, "UPENN_PLASMA_FUJIREBIO_QUANTERIX_24May2026.csv"))

# Coerce numeric (LOQ sentinels like "-4.0" should be dropped)
upenn <- upenn %>%
  mutate(
    pT217_F = suppressWarnings(as.numeric(pT217_F)),
    GFAP_Q  = suppressWarnings(as.numeric(GFAP_Q)),
    RID     = as.integer(RID),
    EXAMDATE = as.Date(EXAMDATE)
  ) %>%
  # ADNI uses -4.0 / -1 as missing sentinels; drop non-positive
  mutate(
    pT217_F = ifelse(pT217_F <= 0, NA, pT217_F),
    GFAP_Q  = ifelse(GFAP_Q  <= 0, NA, GFAP_Q)
  )

cat("UPENN plasma rows (raw):", nrow(upenn), "\n")
cat("Unique RIDs in UPENN:", length(unique(upenn$RID)), "\n")
cat("Rows with non-missing pT217_F:", sum(!is.na(upenn$pT217_F)), "\n")
cat("Rows with non-missing GFAP_Q :", sum(!is.na(upenn$GFAP_Q )), "\n\n")

# ---- APOE genotype ----
apoe <- read_csv(file.path(DATA, "APOERES_24May2026.csv"), ) %>%
  mutate(RID = as.integer(RID)) %>%
  filter(!is.na(GENOTYPE), GENOTYPE != "") %>%
  distinct(RID, .keep_all = TRUE) %>%
  mutate(
    APOE4_n = stringr::str_count(GENOTYPE, "4"),
    APOE4   = as.integer(APOE4_n >= 1)
  ) %>%
  select(RID, GENOTYPE, APOE4)

cat("RIDs with APOE genotype:", nrow(apoe), "\n\n")

# ---- Baseline age and dx from ADNIMERGE ----
adni_bl <- adnimerge %>%
  filter(VISCODE == "bl") %>%
  transmute(
    RID      = as.integer(RID),
    age_bl   = as.numeric(AGE),
    sex      = PTGENDER,
    edu      = as.numeric(PTEDUCAT),
    dx_bl    = DX.bl
  ) %>%
  distinct(RID, .keep_all = TRUE)

cat("Baseline rows in ADNIMERGE:", nrow(adni_bl), "\n")
cat("dx_bl distribution:\n"); print(table(adni_bl$dx_bl, useNA = "ifany")); cat("\n")

# Map ADNIMERGE dx_bl to {CN, MCI, AD}
adni_bl <- adni_bl %>%
  mutate(
    dx3 = case_when(
      dx_bl %in% c("CN", "SMC")               ~ "CN",
      dx_bl %in% c("EMCI", "LMCI", "MCI")     ~ "MCI",
      dx_bl %in% c("AD", "Dementia")          ~ "AD",
      TRUE                                    ~ NA_character_
    )
  )

# ---- Apply inclusion gates sequentially ----
count_per_rid <- upenn %>%
  group_by(RID) %>%
  summarise(
    n_gfap  = sum(!is.na(GFAP_Q)),
    n_ptau  = sum(!is.na(pT217_F)),
    .groups = "drop"
  )

elig <- count_per_rid %>%
  inner_join(adni_bl,  by = "RID") %>%
  inner_join(apoe,     by = "RID") %>%
  filter(!is.na(age_bl), age_bl >= 65) %>%
  filter(!is.na(dx3))

cat("=== After age>=65 + dx in {CN,MCI,AD} + APOE available ===\n")
cat("Eligible RIDs with at least 1 biomarker measurement:", nrow(elig), "\n\n")

# ---- N at each inclusion threshold ----
thresh_tbl <- tibble(
  threshold = c(">=1", ">=2", ">=3", ">=4", ">=5"),
  k = 1:5
) %>%
  rowwise() %>%
  mutate(
    N_both    = sum(elig$n_gfap >= k & elig$n_ptau >= k),
    N_gfap    = sum(elig$n_gfap >= k),
    N_ptau    = sum(elig$n_ptau >= k)
  ) %>%
  ungroup()

cat("=== Inclusion-criterion sensitivity ===\n")
print(thresh_tbl)
cat("\n")

# Diagnosis breakdown at >=2, >=3, >=4
for (k in 2:4) {
  sub <- elig %>% filter(n_gfap >= k, n_ptau >= k)
  cat(sprintf("--- Threshold >=%d (N=%d) ---\n", k, nrow(sub)))
  print(table(sub$dx3, useNA = "ifany"))
  cat(sprintf("  APOE4+ : %d (%.1f%%)\n", sum(sub$APOE4 == 1), 100*mean(sub$APOE4 == 1)))
  cat(sprintf("  Female : %d (%.1f%%)\n", sum(sub$sex == "Female"), 100*mean(sub$sex == "Female")))
  cat(sprintf("  Age mean (SD): %.1f (%.1f)\n\n",
              mean(sub$age_bl, na.rm = TRUE), sd(sub$age_bl, na.rm = TRUE)))
}

# Save
write_csv(thresh_tbl, file.path(OUT, "inclusion_sensitivity.csv"))
cat("Saved:", file.path(OUT, "inclusion_sensitivity.csv"), "\n")
