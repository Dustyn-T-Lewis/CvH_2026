#!/usr/bin/env Rscript
# =============================================================================
# 03_DEP/b_imputed  --  COMPARISON DEP on the imputed matrices.
# imp4p is the CANONICAL imputation (full proteoDA tables + plots); mscoreutils
# is a lighter comparison (tables only). logFC is compared to the non-imputed
# primary fit.
# CAVEAT: exploratory/sensitivity only -- imputing before a DE test can inflate
# false positives; a_non_imputed is the primary analysis.
# =============================================================================

suppressPackageStartupMessages({ library(here); library(readr); library(dplyr); library(purrr); library(tibble) })
source(here("R", "dep_model.R"))

CONTRASTS <- c("CRvH_Baseline = CR_pre - H_pre",
               "CR_Training = CR_post - CR_pre",
               "Resid = CR_post - H_pre")
FORMULA   <- "~ 0 + group_time + (1 | Subject_ID)"
CANONICAL <- "imp4p"
clear_dir <- function(d) { dir.create(d, recursive = TRUE, showWarnings = FALSE)
  unlink(setdiff(list.files(d, full.names = TRUE), file.path(d, ".gitkeep")), recursive = TRUE) }
clear_dir(here("03_DEP", "b_imputed", "c_data")); clear_dir(here("03_DEP", "b_imputed", "b_reports"))
nd <- here("02_Normalization", "c_data")
methods <- c(imp4p = "DAList_imputed_imp4p.rds", mscoreutils = "DAList_imputed_mscoreutils.rds")

runs <- imap(methods, function(rds, m) {
  dal <- readRDS(file.path(nd, rds))
  cat(sprintf("\n[%s] imputed DEP: %d x %d%s\n", m, nrow(dal$data), ncol(dal$data),
              if (m == CANONICAL) " (canonical)" else ""))
  run_dep_model(dal, FORMULA, CONTRASTS,
                out_dir   = here("03_DEP", "b_imputed", "c_data", m),
                grouping_column = "group_time",
                report_dir = if (m == CANONICAL) here("03_DEP", "b_imputed", "b_reports", m) else NULL)
})

# --- logFC concordance vs the non-imputed primary ----------------------------
ni_file <- here("03_DEP", "a_non_imputed", "c_data", "combined_results_pi.csv")
if (file.exists(ni_file)) {
  ni <- read_csv(ni_file, show_col_types = FALSE) |> select(uniprot_id, contrast, logFC_ni = logFC)
  cmp <- imap_dfr(runs, function(out, m) {
    bind_rows(out$results) |> select(uniprot_id, contrast, logFC) |>
      inner_join(ni, by = c("uniprot_id", "contrast")) |> group_by(contrast) |>
      summarise(method = m,
                spearman_vs_nonimputed = cor(logFC, logFC_ni, method = "spearman", use = "complete.obs"),
                .groups = "drop")
  })
  write_csv(cmp, here("03_DEP", "b_imputed", "c_data", "logFC_vs_nonimputed.csv"))
  cat("\nlogFC concordance vs non-imputed (Spearman):\n"); print(as.data.frame(cmp))
}
