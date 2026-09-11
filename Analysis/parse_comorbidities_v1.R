# Parse HTN / DM / Stroke from ADNI medical history
# Sources: RECMHIST (free-text MHDESC), INITHEALTH (IHDESC, IHSYMPTOM)
# Output: comorbidities.csv with RID, HTN_base, DM_base, stroke_ever

suppressPackageStartupMessages({ library(dplyr); library(stringr) })

DATA <- "/Users/kuan-yuchu/Desktop/LOAD/ADNI_csv"
OUT  <- "/Users/kuan-yuchu/Desktop/LOAD/Analysis/outputs"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

rec <- utils::read.csv(file.path(DATA, "RECMHIST_24May2026.csv"),
                       stringsAsFactors = FALSE, na.strings = c("", "NA"))
ini <- utils::read.csv(file.path(DATA, "INITHEALTH_24May2026.csv"),
                       stringsAsFactors = FALSE, na.strings = c("", "NA"))

cat("RECMHIST rows:", nrow(rec), " unique RIDs:", length(unique(rec$RID)), "\n")
cat("INITHEALTH rows:", nrow(ini), " unique RIDs:", length(unique(ini$RID)), "\n\n")

# Combine free-text fields
rec <- rec %>% transmute(RID = as.integer(RID), txt = toupper(MHDESC),
                         dt = as.Date(VISDATE))
ini <- ini %>% transmute(RID = as.integer(RID),
                         txt = toupper(paste(IHSYMPTOM, IHDESC)),
                         dt  = as.Date(VISDATE))
all_mh <- bind_rows(rec, ini) %>% filter(!is.na(RID), !is.na(txt))

# Regex patterns
re_htn    <- "HYPERTEN|\\bHTN\\b|HIGH BLOOD PRESS|ELEVATED BP"
re_dm     <- "DIABETES|\\bDM\\b|TYPE 2 DIABETE|TYPE II DIABETE|NIDDM|T2DM|IDDM|HYPERGLYCEM"
re_stroke <- "\\bSTROKE(S|D)?\\b|\\bCVA\\b|CEREBROVASCULAR ACC|CEREBRAL INFARC|CEREBELLAR INFARC|ISCHEMIC STROKE|HEMORRHAGIC STROKE|HEMORRH.+CEREB|INTRACRANIAL HEMORR|\\bTIA\\b|\\bTIAS\\b|\\bTIA'S\\b|TRANSIENT ISCHEMIC"

mh_flag <- all_mh %>%
  mutate(
    is_htn    = str_detect(txt, re_htn),
    is_dm     = str_detect(txt, re_dm),
    is_stroke = str_detect(txt, re_stroke)
  )

# Aggregate per RID
comorb <- mh_flag %>% group_by(RID) %>%
  summarise(
    HTN_base    = as.integer(any(is_htn, na.rm = TRUE)),
    DM_base     = as.integer(any(is_dm, na.rm = TRUE)),
    stroke_ever = as.integer(any(is_stroke, na.rm = TRUE)),
    .groups = "drop"
  )

cat("=== Comorbidity prevalence (all RIDs with any medical history) ===\n")
cat("HTN  :", sum(comorb$HTN_base),    sprintf("(%.1f%%)\n", 100*mean(comorb$HTN_base)))
cat("DM   :", sum(comorb$DM_base),     sprintf("(%.1f%%)\n", 100*mean(comorb$DM_base)))
cat("Stroke:", sum(comorb$stroke_ever),sprintf("(%.1f%%)\n", 100*mean(comorb$stroke_ever)))
cat("Total RIDs flagged at least one:", nrow(comorb), "\n\n")

# Spot-check examples
cat("=== Example free-text matches (top 10 each) ===\n")
cat("\n-- HTN --\n"); print(head(unique(mh_flag$txt[mh_flag$is_htn]), 10))
cat("\n-- DM --\n");  print(head(unique(mh_flag$txt[mh_flag$is_dm]), 10))
cat("\n-- Stroke --\n"); print(head(unique(mh_flag$txt[mh_flag$is_stroke]), 10))

utils::write.csv(comorb, file.path(OUT, "comorbidities.csv"), row.names = FALSE)
cat("\nSaved:", file.path(OUT, "comorbidities.csv"), "\n")
