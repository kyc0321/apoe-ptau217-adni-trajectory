suppressPackageStartupMessages({library(brms)})
fit <- readRDS("Analysis/outputs/brmsfit_Hippo_pT217_bygroup_Apos.rds")
ps <- as_draws_df(fit)

slope_slope_aneg <- "cor_RID__bm_time_yrs:A_pos_fAneg__outc_time_yrs:A_pos_fAneg"
slope_slope_apos <- "cor_RID__bm_time_yrs:A_pos_fApos__outc_time_yrs:A_pos_fApos"

d_aneg <- ps[[slope_slope_aneg]]
d_apos <- ps[[slope_slope_apos]]

cat("A- (amyloid-negative) slope-slope rho: mean=", round(mean(d_aneg),3),
    " [", round(quantile(d_aneg,.025),3), ",", round(quantile(d_aneg,.975),3), "]\n", sep="")
cat("A+ (amyloid-positive) slope-slope rho: mean=", round(mean(d_apos),3),
    " [", round(quantile(d_apos,.025),3), ",", round(quantile(d_apos,.975),3), "]\n", sep="")

diff <- d_apos - d_aneg
cat("\nDifference (A+ minus A-): mean=", round(mean(diff),3),
    " [", round(quantile(diff,.025),3), ",", round(quantile(diff,.975),3), "]",
    " P(diff>0)=", round(mean(diff>0),3), "\n", sep="")

# Also report bulk/tail ESS specifically for these two parameters (console flagged low ESS somewhere)
s <- summary(fit)
cat("\nFull-model max Rhat:", round(max(brms::rhat(fit), na.rm=TRUE), 4), "\n")
