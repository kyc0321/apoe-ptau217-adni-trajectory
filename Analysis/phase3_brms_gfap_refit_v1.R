# Refit 3 GFAP-paired brms joint models with longer chains for convergence
# Author: Kuan-Yu Chu | 2026-05-25
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(brms)
})
options(brms.backend = "rstan")
set.seed(20260525)

DATA <- "/Users/kuan-yuchu/Desktop/LOAD/ADNI_csv"
OUT  <- "/Users/kuan-yuchu/Desktop/LOAD/Analysis/outputs"

biom_long <- utils::read.csv(file.path(OUT, "analytic_biom_long.csv"), stringsAsFactors = FALSE)
attr_tbl  <- utils::read.csv(file.path(OUT, "analytic_baseline_attr.csv"), stringsAsFactors = FALSE)
biom_long$dx3 <- factor(biom_long$dx3, levels = c("CN","MCI","AD"))
attr_tbl$dx3  <- factor(attr_tbl$dx3,  levels = c("CN","MCI","AD"))
zsc <- function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
attr_tbl <- attr_tbl %>% mutate(age_z = zsc(age_bl), edu_z = zsc(edu))

adas <- utils::read.csv(file.path(DATA, "ADAS_24May2026.csv"), stringsAsFactors = FALSE, na.strings = c("","NA"))
mmse <- utils::read.csv(file.path(DATA, "MMSE_24May2026.csv"), stringsAsFactors = FALSE, na.strings = c("","NA"))
fs   <- utils::read.csv(file.path(DATA, "ADSP_PHC_T1_FS_24May2026.csv"), stringsAsFactors = FALSE, na.strings = c("","NA"))

bl_dates <- biom_long %>% group_by(RID) %>%
  summarise(bl_date = as.Date(min(EXAMDATE)), .groups = "drop")
build_outc <- function(df, ycol, datecol) {
  df %>% mutate(RID = as.integer(RID), EXAMDATE = as.Date(.data[[datecol]]),
                y = suppressWarnings(as.numeric(.data[[ycol]]))) %>%
    filter(!is.na(EXAMDATE), !is.na(y), RID %in% attr_tbl$RID) %>%
    left_join(bl_dates, by = "RID") %>%
    mutate(time_yrs = as.numeric(EXAMDATE - bl_date)/365.25) %>%
    filter(!is.na(time_yrs))
}
adas_long <- build_outc(adas, "TOTAL13", "VISDATE")  %>% transmute(RID, time_yrs, outc = y)
mmse_long <- build_outc(mmse, "MMSCORE", "VISDATE")  %>% transmute(RID, time_yrs, outc = y)
fs_long   <- fs %>% mutate(hippo_ratio = (Left.Hippocampus_combat + Right.Hippocampus_combat) /
                             EstimatedTotalIntraCranialVol_combat) %>%
  build_outc("hippo_ratio", "PHC_SCANDATE") %>% transmute(RID, time_yrs, outc = y * 1000)

biom_red <- biom_long %>% filter(!is.na(BMI_z)) %>% transmute(RID, time_yrs, bm = log_gfap) %>% filter(!is.na(bm))

fit_one <- function(outc_df, lab) {
  cat("\n=== Refit:", lab, "===  ", as.character(Sys.time()), "\n")
  rid_keep <- intersect(unique(biom_red$RID), unique(outc_df$RID))
  bm_df <- biom_red[biom_red$RID %in% rid_keep, ]
  oc_df <- outc_df[outc_df$RID %in% rid_keep, ]
  bm_rows <- bm_df %>% mutate(outc = 0, has_bm = TRUE, has_oc = FALSE)
  oc_rows <- oc_df %>% mutate(bm = 0, has_bm = FALSE, has_oc = TRUE)
  d <- bind_rows(bm_rows, oc_rows) %>%
    left_join(attr_tbl %>% select(RID, age_z, sex, edu_z, dx3, APOE4), by = "RID")
  d$RID <- factor(d$RID)
  f_bm <- bf(bm   | subset(has_bm) ~ time_yrs + APOE4 + age_z + sex + dx3 + (1 + time_yrs | p | RID))
  f_oc <- bf(outc | subset(has_oc) ~ time_yrs + APOE4 + age_z + sex + edu_z + dx3 + (1 + time_yrs | p | RID))
  fit <- brm(f_bm + f_oc + set_rescor(FALSE), data = d, family = gaussian(),
             chains = 4, iter = 6000, warmup = 2000, cores = 4,
             seed = 20260525, refresh = 0,
             control = list(adapt_delta = 0.99, max_treedepth = 14))
  ps <- as_draws_df(fit)
  cor_nm <- grep("^cor_RID__bm_time_yrs__outc_time_yrs$", names(ps), value = TRUE)
  draws <- ps[[cor_nm]]
  summ <- data.frame(label = lab, parameter = cor_nm,
                     mean = mean(draws), median = median(draws),
                     lo = quantile(draws, .025), hi = quantile(draws, .975),
                     p_dir = mean(draws < 0),
                     rhat = max(brms::rhat(fit)[grep("cor_RID", names(brms::rhat(fit)))]),
                     n_id = length(rid_keep), n_bm = nrow(bm_df), n_oc = nrow(oc_df))
  cat("  rho =", round(summ$mean,3), "[", round(summ$lo,3), ",", round(summ$hi,3), "]  Rhat =", round(summ$rhat,3), "\n")
  saveRDS(fit, file.path(OUT, paste0("brmsfit_", lab, "_refit.rds")))
  summ
}

results <- bind_rows(
  fit_one(adas_long, "ADAS_GFAP"),
  fit_one(mmse_long, "MMSE_GFAP"),
  fit_one(fs_long,   "Hippo_GFAP")
)
utils::write.csv(results, file.path(OUT, "table4_brms_GFAP_refit.csv"), row.names = FALSE)
cat("\n=== GFAP refit results ===\n"); print(results)
