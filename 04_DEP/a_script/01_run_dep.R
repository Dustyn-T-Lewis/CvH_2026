#!/usr/bin/env Rscript
# =============================================================================
# 01_run_dep.R  --  CvH Stage 04: Differential abundance (single 3-group model)
#
# One limma fit on the NON-imputed normalized matrix:
#   ~ 0 + group_time + (1 | Subject_ID)   group_time in {H_pre, CR_pre, CR_post}
# proteoDA routes the random term through duplicateCorrelation (Smyth 2005);
# eBayes(robust=TRUE). Contrasts:
#   CRvH_Baseline = CR_pre  - H_pre    (cross-sectional cancer vs healthy)
#   CR_Training   = CR_post - CR_pre   (within-subject, paired recovery effect)
#   Resid         = CR_post - H_pre    (residual deviation after training; reversal)
# Significance: nominal p<0.10, BH-FDR, Pi-score<0.05 (Xiao 2014).
# =============================================================================

suppressPackageStartupMessages({
  library(proteoDA); library(here); library(readr); library(dplyr); library(tibble); library(purrr)
})
set.seed(42)
source(here("R", "pi_score.R"))

data_dir <- here("04_DEP", "c_data")
dir.create(file.path(data_dir, "per_contrast"), recursive = TRUE, showWarnings = FALSE)

dal <- readRDS(here("02_Normalization", "c_data", "03_DAList_normalized.rds"))
cat(sprintf("Loaded normalized DAList: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

dal <- add_design(dal, "~ 0 + group_time + (1 | Subject_ID)")
dal <- add_contrasts(dal, contrasts_vector = c(
  "CRvH_Baseline = CR_pre - H_pre",
  "CR_Training = CR_post - CR_pre",
  "Resid = CR_post - H_pre"
))
dal <- fit_limma_model(dal)
rho <- dal$eBayes_fit$correlation %||% NA_real_
cat(sprintf("Within-subject duplicateCorrelation rho = %.3f\n", rho))

dal <- extract_DA_results(dal, pval_thresh = 0.10, lfc_thresh = 0, adj_method = "BH")
res <- compute_pi_scores(dal$results, pi_thresh = 0.05)

ann <- as_tibble(dal$annotation) |> select(uniprot_id, gene, protein, description)
res <- map(res, ~ left_join(.x, ann, by = "uniprot_id"))

# per-contrast + combined long table
iwalk(res, ~ write_csv(.x, file.path(data_dir, "per_contrast", paste0(.y, ".csv"))))
write_csv(bind_rows(res), file.path(data_dir, "03_combined_results.csv"))

summary_tbl <- map_dfr(res, ~ tibble(
  contrast  = .x$contrast[1],
  n         = nrow(.x),
  sig_p10   = sum(.x$P.Value  < 0.10, na.rm = TRUE),
  sig_fdr10 = sum(.x$adj.P.Val < 0.10, na.rm = TRUE),
  sig_fdr05 = sum(.x$adj.P.Val < 0.05, na.rm = TRUE),
  up_pi     = sum(.x$sig_pi ==  1),
  down_pi   = sum(.x$sig_pi == -1)
))
write_csv(summary_tbl, file.path(data_dir, "02_DA_summary.csv"))

saveRDS(dal, file.path(data_dir, "01_DAList_dep.rds"))
writeLines(capture.output(sessionInfo()), file.path(data_dir, "sessionInfo.txt"))
cat("\nDA summary:\n"); print(as.data.frame(summary_tbl))
cat(sprintf("\nDone -> %s/\n", data_dir))
