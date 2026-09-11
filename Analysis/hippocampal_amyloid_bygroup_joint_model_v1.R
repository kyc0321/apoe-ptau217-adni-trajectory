suppressPackageStartupMessages({library(dplyr); library(brms)})
options(brms.backend = "rstan")
set.seed(20260911)

BASE <- "/Users/kuan-yuchu/Desktop/#2 LOAD-BC-NS"
DATA <- file.path(BASE, "ADNI_csv")
OUT  <- file.path(BASE, "Analysis/outputs")

biom_long <- read.csv(file.path(OUT, "analytic_biom_long.csv"))
biom_long$dx3 <- factor(biom_long$dx3, levels = c("CN","MCI","AD"))

upenn <- read.csv(file.path(DATA, "UPENN_PLASMA_FUJIREBIO_QUANTERIX_24May2026.csv")) %>%
  mutate(RID=as.integer(RID), EXAMDATE=as.Date(EXAMDATE),
         AB42_AB40_F=suppressWarnings(as.numeric(AB42_AB40_F)),
         AB42_AB40_F=ifelse(AB42_AB40_F<=0|AB42_AB40_F>1, NA, AB42_AB40_F))
ab_bl <- upenn %>% filter(!is.na(AB42_AB40_F)) %>% group_by(RID) %>%
  slice_min(EXAMDATE, with_ties=FALSE) %>% ungroup() %>%
  transmute(RID, A_pos = as.integer(AB42_AB40_F < 0.085))

bl_dates <- biom_long %>% group_by(RID) %>% summarise(bl_date = as.Date(min(EXAMDATE)), .groups="drop")

fs <- read.csv(file.path(DATA, "ADSP_PHC_T1_FS_24May2026.csv"), stringsAsFactors=FALSE, na.strings=c("","NA"))
fs_long <- fs %>% transmute(
  RID = as.integer(RID), EXAMDATE = as.Date(PHC_SCANDATE),
  hippo = Left.Hippocampus_combat + Right.Hippocampus_combat,
  ICV = EstimatedTotalIntraCranialVol_combat,
  hippo_ratio = hippo/ICV
) %>% filter(!is.na(hippo_ratio), RID %in% unique(biom_long$RID)) %>%
  left_join(bl_dates, by="RID") %>%
  mutate(time_yrs = as.numeric(EXAMDATE - bl_date)/365.25) %>%
  filter(!is.na(time_yrs)) %>% transmute(RID, time_yrs, outc = hippo_ratio*1000)

attr_z <- biom_long %>% distinct(RID, .keep_all=TRUE) %>%
  transmute(RID, age_z, sex=factor(sex), edu_z, dx3, APOE4) %>%
  left_join(ab_bl, by="RID") %>% filter(!is.na(A_pos))

biom_red <- biom_long %>% filter(!is.na(BMI_z)) %>% transmute(RID, time_yrs, bm=log_ptau) %>% filter(!is.na(bm))

rid_keep <- intersect(intersect(unique(biom_red$RID), unique(fs_long$RID)), unique(attr_z$RID))
cat("n_id pooled (has p-tau217 + hippo + A_pos) =", length(rid_keep), "\n")
bm_df <- biom_red[biom_red$RID %in% rid_keep, ]
oc_df <- fs_long[fs_long$RID %in% rid_keep, ]
bm_rows <- bm_df %>% mutate(outc=0, has_bm=TRUE, has_oc=FALSE)
oc_rows <- oc_df %>% mutate(bm=0, has_bm=FALSE, has_oc=TRUE)
d <- bind_rows(bm_rows, oc_rows) %>% left_join(attr_z, by="RID")
d$RID <- factor(d$RID)
d$A_pos_f <- factor(d$A_pos, levels=c(0,1), labels=c("Aneg","Apos"))
cat("A+ :", sum(attr_z$RID %in% rid_keep & attr_z$A_pos==1), " A- :", sum(attr_z$RID %in% rid_keep & attr_z$A_pos==0), "\n")

# Single POOLED model, random-effects covariance allowed to vary BY amyloid group
# (gr(RID, by=A_pos_f)) -- this gives a SEPARATE correlation estimate per group
# within one model fit on the full sample, rather than two underpowered split fits.
f_bm <- bf(bm   | subset(has_bm) ~ time_yrs*A_pos_f + age_z + sex + dx3 + (1 + time_yrs | p | gr(RID, by = A_pos_f)))
f_oc <- bf(outc | subset(has_oc) ~ time_yrs*A_pos_f + age_z + sex + edu_z + dx3 + (1 + time_yrs | p | gr(RID, by = A_pos_f)))

cat("\nFitting pooled by-group joint model (this may take several minutes)...\n")
fit <- brm(f_bm + f_oc + set_rescor(FALSE), data = d, family = gaussian(),
           chains = 4, iter = 4000, warmup = 1500, cores = 4, seed = 20260911, refresh = 0,
           control = list(adapt_delta = 0.97, max_treedepth = 13))

saveRDS(fit, file.path(OUT, "brmsfit_Hippo_pT217_bygroup_Apos.rds"))

ps <- as_draws_df(fit)
cor_cols <- grep("^cor_RID", names(ps), value=TRUE)
cat("\nCorrelation parameter names found:\n"); print(cor_cols)

rh <- brms::rhat(fit)
cat("\nmax Rhat (all params):", round(max(rh, na.rm=TRUE),4), "\n")
np <- nuts_params(fit)
cat("divergent transitions:", sum(np$Parameter=="divergent__" & np$Value==1), "/", sum(np$Parameter=="divergent__"), "\n")

for (cc in cor_cols) {
  draws <- ps[[cc]]
  cat(sprintf("%-55s mean=%.3f [%.3f, %.3f]\n", cc, mean(draws), quantile(draws,.025), quantile(draws,.975)))
}

if (length(cor_cols) == 2) {
  diff <- ps[[cor_cols[1]]] - ps[[cor_cols[2]]]
  cat(sprintf("\nDifference in correlation (group1 - group2): mean=%.3f [%.3f, %.3f], P(diff>0)=%.3f\n",
              mean(diff), quantile(diff,.025), quantile(diff,.975), mean(diff>0)))
}
