#!/usr/bin/env Rscript
# =============================================================================
# 03_DEP/a_non_imputed  --  PRIMARY DEP on the NON-imputed normalized matrix.
#
# limma handles per-protein NAs directly; no imputation before testing (the
# gold-standard choice -- imputing before a DE test inflates false positives).
# One fit: ~ 0 + group_time + (1 | Subject_ID), group_time {H_pre,CR_pre,CR_post};
# contrasts CRvH_Baseline / CR_Training / Resid. See R/dep_model.R.
# =============================================================================

suppressPackageStartupMessages({ library(here) })
source(here("R", "dep_model.R"))

CONTRASTS <- c("CRvH_Baseline = CR_pre - H_pre",
               "CR_Training = CR_post - CR_pre",
               "Resid = CR_post - H_pre")

dal <- readRDS(here("02_Normalization", "c_data", "03_DAList_normalized.rds"))
cat(sprintf("NON-imputed DEP: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))
out <- run_dep_model(dal, "~ 0 + group_time + (1 | Subject_ID)", CONTRASTS,
                     here("03_DEP", "a_non_imputed", "c_data"))
cat(sprintf("duplicateCorrelation rho = %.3f\n", out$rho))
cat("\nDA summary (non-imputed):\n"); print(as.data.frame(out$summary))
