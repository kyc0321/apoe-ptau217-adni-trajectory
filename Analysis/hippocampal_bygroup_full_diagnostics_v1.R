suppressPackageStartupMessages({library(brms)})
fit <- readRDS("Analysis/outputs/brmsfit_Hippo_pT217_bygroup_Apos.rds")

rh <- brms::rhat(fit)
cat("max Rhat (all params):", round(max(rh, na.rm=TRUE), 4), "\n")
cat("n params Rhat>1.05:", sum(rh > 1.05, na.rm=TRUE), "\n")

np <- nuts_params(fit)
n_div <- sum(np$Parameter=="divergent__" & np$Value==1)
n_tot <- sum(np$Parameter=="divergent__")
cat("divergent transitions:", n_div, "/", n_tot, "\n")

smry <- as.data.frame(summary(fit)$fixed)
cat("Bulk_ESS range (fixed effects):", paste(range(smry$Bulk_ESS, na.rm=TRUE), collapse="-"), "\n")
cat("Tail_ESS range (fixed effects):", paste(range(smry$Tail_ESS, na.rm=TRUE), collapse="-"), "\n")

# full ESS across ALL parameters (not just fixed effects) -- this is what the low-ESS warning was about
library(posterior)
ps <- as_draws_df(fit)
ess_all <- summarise_draws(as_draws_array(fit), "ess_bulk", "ess_tail")
cat("\nOverall min ess_bulk (all params):", round(min(ess_all$ess_bulk, na.rm=TRUE),0), "\n")
cat("Overall min ess_tail (all params):", round(min(ess_all$ess_tail, na.rm=TRUE),0), "\n")
worst <- ess_all[order(ess_all$ess_bulk),][1:5,]
cat("\nWorst 5 parameters by bulk ESS:\n")
print(worst[,c("variable","ess_bulk","ess_tail")])

# PPC Bayes-p for mean/SD
yrep_bm <- posterior_predict(fit, resp="bm", ndraws=500)
y_bm <- get_y(fit, resp="bm")
yrep_oc <- posterior_predict(fit, resp="outc", ndraws=500)
y_oc <- get_y(fit, resp="outc")
cat("\nPPC bm: Bayes-p(mean)=", round(mean(apply(yrep_bm,1,mean)>=mean(y_bm,na.rm=TRUE)),2),
    " Bayes-p(sd)=", round(mean(apply(yrep_bm,1,sd)>=sd(y_bm,na.rm=TRUE)),2), "\n", sep="")
cat("PPC outc: Bayes-p(mean)=", round(mean(apply(yrep_oc,1,mean)>=mean(y_oc,na.rm=TRUE)),2),
    " Bayes-p(sd)=", round(mean(apply(yrep_oc,1,sd)>=sd(y_oc,na.rm=TRUE)),2), "\n", sep="")
