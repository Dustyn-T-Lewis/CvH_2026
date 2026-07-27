#!/usr/bin/env Rscript
# CvH imputed DEP, one arm per imputation method (a_non_imputed stays the reported analysis).
# missForest is canonical (feeds PCA and WGCNA); imp4p and mscoreutils are comparison arms.
# Each method writes its own b_reports/ + c_data/ subdir; logFC is checked against the
# non-imputed fit, since impute-before-test can inflate false positives.

pacman::p_load(proteoDA, here, readr, dplyr, tibble, purrr)

CANONICAL <- "missforest"
nd <- here("02_Normalization", "imputation", "c_data")
methods <- c(
  missforest = "DAList_imputed_missforest.rds",
  imp4p = "DAList_imputed_imp4p.rds",
  mscoreutils = "DAList_imputed_mscoreutils.rds"
)

contrasts <- c(
  CRvH_Baseline          = "0.5*CRE_pre + 0.5*PLA_pre - H_pre",
  CR_Training            = "0.5*CRE_post + 0.5*PLA_post - 0.5*CRE_pre - 0.5*PLA_pre",
  Resid                  = "0.5*CRE_post + 0.5*PLA_post - H_pre",
  Baseline_Supplement    = "CRE_pre - PLA_pre",
  Training_CRE           = "CRE_post - CRE_pre",
  Training_PLA           = "PLA_post - PLA_pre",
  Supplement_Interaction = "CRE_post - CRE_pre - PLA_post + PLA_pre"
)

clear_dir <- function(d) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  unlink(setdiff(list.files(d, full.names = TRUE), file.path(d, ".gitkeep")), recursive = TRUE)
}

runs <- imap(methods, function(rds, m) {
  out_c <- here("03_DEP", "b_imputed", m, "c_data")
  out_b <- here("03_DEP", "b_imputed", m, "b_reports")
  clear_dir(out_c)
  clear_dir(out_b)

  dal <- readRDS(file.path(nd, rds))
  cat(sprintf(
    "[%s] imputed DEP: %d x %d%s\n", m, nrow(dal$data), ncol(dal$data),
    if (m == CANONICAL) " (canonical)" else ""
  ))

  dal$metadata$model_cell <- dplyr::recode(dal$metadata$Group_Time,
    CRE_T1 = "CRE_pre", CRE_T2 = "CRE_post",
    PLA_T1 = "PLA_pre", PLA_T2 = "PLA_post", H_T1 = "H_pre"
  )
  dal <- add_design(dal, "~ 0 + model_cell + (1 | Subject_ID)")
  cm <- limma::makeContrasts(contrasts = contrasts, levels = dal$design$design_matrix)
  colnames(cm) <- names(contrasts)
  dal$design$contrast_matrix <- cm
  dal$design$contrast_vector <- unname(contrasts)
  dal <- fit_limma_model(dal)
  dal <- extract_DA_results(dal, pval_thresh = 0.10, lfc_thresh = 0, adj_method = "BH")

  write_limma_tables(dal,
    output_dir = out_c, overwrite = TRUE,
    annot_cols = c("uniprot_id", "gene", "protein", "description")
  )
  unlink(c(
    file.path(out_c, "combined_results.csv"), file.path(out_c, "DA_summary.csv"),
    file.path(out_c, "per_contrast_results")
  ), recursive = TRUE)
  write_limma_plots(dal,
    grouping_column = "group_time",
    table_columns = c("uniprot_id", "gene"), output_dir = out_b, overwrite = TRUE
  )

  ann <- as_tibble(dal$annotation) |> select(any_of(c("uniprot_id", "gene", "protein", "description")))
  res <- imap(dal$results, function(r, cname) {
    as_tibble(r, rownames = "uniprot_id") |>
      mutate(
        pi_score = P.Value^abs(logFC),
        sig_pi = case_when(
          is.na(pi_score) ~ NA_integer_,
          pi_score < 0.05 & logFC > 0 ~ 1L,
          pi_score < 0.05 & logFC < 0 ~ -1L, TRUE ~ 0L
        ),
        contrast = cname
      ) |>
      left_join(ann, by = "uniprot_id")
  })
  write_csv(bind_rows(res), file.path(out_c, "combined_results_pi.csv"))
  res
})

# Spearman logFC concordance of each arm against the non-imputed fit.
ni_file <- here("03_DEP", "a_non_imputed", "c_data", "combined_results_pi.csv")
if (file.exists(ni_file)) {
  ni <- read_csv(ni_file, show_col_types = FALSE) |> select(uniprot_id, contrast, logFC_ni = logFC)
  cmp <- imap_dfr(runs, function(res, m) {
    bind_rows(res) |>
      select(uniprot_id, contrast, logFC) |>
      inner_join(ni, by = c("uniprot_id", "contrast")) |>
      group_by(contrast) |>
      summarise(
        method = m,
        spearman_vs_nonimputed = cor(logFC, logFC_ni, method = "spearman", use = "complete.obs"),
        .groups = "drop"
      )
  })
  write_csv(cmp, here("03_DEP", "b_imputed", "logFC_vs_nonimputed.csv"))
  cat("\nlogFC concordance vs non-imputed (Spearman):\n")
  print(as.data.frame(cmp))
}
