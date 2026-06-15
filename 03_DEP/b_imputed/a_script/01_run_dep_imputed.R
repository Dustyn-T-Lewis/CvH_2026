#!/usr/bin/env Rscript
# =============================================================================
# 03_DEP/b_imputed  --  COMPARISON DEP on the imputed matrices.
#
# Runs the SAME model on each imputed DAList (MsCoreUtils hybrid + imp4p mixture),
# writing per-method results, and compares logFC to the non-imputed primary fit.
#
# CAVEAT: this arm is exploratory/sensitivity ONLY. Imputing before a DE test can
# inflate false positives; a_non_imputed is the primary analysis. This arm shows
# how/whether imputation moves the contrasts.
# =============================================================================

suppressPackageStartupMessages({ library(here); library(readr); library(dplyr); library(purrr); library(tibble) })
source(here("R", "dep_model.R"))

CONTRASTS <- c("CRvH_Baseline = CR_pre - H_pre",
               "CR_Training = CR_post - CR_pre",
               "Resid = CR_post - H_pre")
FORMULA <- "~ 0 + group_time + (1 | Subject_ID)"
nd <- here("02_Normalization", "c_data")
methods <- c(mscoreutils = "04_DAList_imputed_mscoreutils.rds", imp4p = "05_DAList_imputed_imp4p.rds")

runs <- imap(methods, function(rds, m) {
  dal <- readRDS(file.path(nd, rds))
  cat(sprintf("\n[%s] imputed DEP: %d x %d\n", m, nrow(dal$data), ncol(dal$data)))
  out <- run_dep_model(dal, FORMULA, CONTRASTS, here("03_DEP", "b_imputed", "c_data", m))
  cat(sprintf("[%s] rho = %.3f\n", m, out$rho)); print(as.data.frame(out$summary))
  out
})

# --- Compare logFC to the non-imputed primary fit ----------------------------
nonimp_file <- here("03_DEP", "a_non_imputed", "c_data", "03_combined_results.csv")
if (file.exists(nonimp_file)) {
  ni <- read_csv(nonimp_file, show_col_types = FALSE) |> select(uniprot_id, contrast, logFC_ni = logFC)
  cmp <- imap_dfr(runs, function(out, m) {
    bind_rows(out$results) |> select(uniprot_id, contrast, logFC) |>
      inner_join(ni, by = c("uniprot_id", "contrast")) |>
      group_by(contrast) |>
      summarise(method = m,
                spearman_vs_nonimputed = cor(logFC, logFC_ni, method = "spearman", use = "complete.obs"),
                .groups = "drop")
  })
  write_csv(cmp, here("03_DEP", "b_imputed", "c_data", "00_logFC_vs_nonimputed.csv"))
  cat("\nlogFC concordance vs non-imputed (Spearman):\n"); print(as.data.frame(cmp))
}
