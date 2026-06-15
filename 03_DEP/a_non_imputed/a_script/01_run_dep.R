#!/usr/bin/env Rscript
# =============================================================================
# 03_DEP/a_non_imputed  --  PRIMARY DEP on the NON-imputed normalized matrix.
# limma handles per-protein NAs (no impute-before-test). Full proteoDA reporting
# (tables -> c_data, plots -> b_reports). See R/dep_model.R.
#   ~ 0 + group_time + (1 | Subject_ID), group_time {H_pre, CR_pre, CR_post}
#   CRvH_Baseline / CR_Training / Resid
# =============================================================================

suppressPackageStartupMessages({ library(here) })
source(here("R", "dep_model.R"))

CONTRASTS <- c("CRvH_Baseline = CR_pre - H_pre",
               "CR_Training = CR_post - CR_pre",
               "Resid = CR_post - H_pre")

clear_dir <- function(d) { dir.create(d, recursive = TRUE, showWarnings = FALSE)
  unlink(setdiff(list.files(d, full.names = TRUE), file.path(d, ".gitkeep")), recursive = TRUE) }
clear_dir(here("03_DEP", "a_non_imputed", "c_data")); clear_dir(here("03_DEP", "a_non_imputed", "b_reports"))

dal <- readRDS(here("02_Normalization", "c_data", "DAList_normalized.rds"))
cat(sprintf("NON-imputed DEP: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

out <- run_dep_model(dal, "~ 0 + group_time + (1 | Subject_ID)", CONTRASTS,
                     out_dir = here("03_DEP", "a_non_imputed", "c_data"),
                     grouping_column = "group_time",
                     report_dir = here("03_DEP", "a_non_imputed", "b_reports"))
cat(sprintf("duplicateCorrelation rho = %.3f | proteoDA tables + plots written\n", out$rho))
