#!/usr/bin/env Rscript
# CvH PRIMARY DEP on the NON-imputed normalized matrix (limma handles per-protein NAs).
#   ~ 0 + model_cell + (1 | Subject_ID) over 5 cells (H_pre + CRE/PLA x pre/post).
#   Supplement is kept in the model so pooled estimates are supplement-adjusted;
#   the reported contrasts average the CRE/PLA arms 50:50 and never test CRE vs PLA.
#   Contrasts: CRvH_Baseline / CR_Training / Resid.
#   Significance: Pi-score (Xiao 2014, Pi = P.Value^|logFC|) < 0.05 + BH-FDR.

pacman::p_load(proteoDA, here, readr, dplyr, tibble, purrr)

out_dir <- here("03_DEP", "a_non_imputed", "c_data")
report_dir <- here("03_DEP", "a_non_imputed", "b_reports")
clear_dir <- function(d) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  unlink(setdiff(list.files(d, full.names = TRUE), file.path(d, ".gitkeep")), recursive = TRUE)
}
clear_dir(out_dir)
clear_dir(report_dir)

dal <- readRDS(here("02_Normalization", "c_data", "DAList_normalized.rds"))
cat(sprintf("NON-imputed DEP: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

# Fit limma model
# (1 | Subject_ID) carries the within-subject pre->post pairing via duplicateCorrelation.

dal$metadata$model_cell <- dplyr::recode(dal$metadata$Group_Time,
  CRE_T1 = "CRE_pre", CRE_T2 = "CRE_post",
  PLA_T1 = "PLA_pre", PLA_T2 = "PLA_post", H_T1 = "H_pre"
)

dal <- add_design(dal, "~ 0 + model_cell + (1 | Subject_ID)")

# proteoDA::add_contrasts rejects weighted terms, so attach the contrast matrix
# with limma directly. Reported: pooled D / T / R = D + T. Supplement contrasts
# (CRE vs PLA) feed the supplementary figure only — never the main analysis.
contrasts <- c(
  CRvH_Baseline          = "0.5*CRE_pre + 0.5*PLA_pre - H_pre",
  CR_Training            = "0.5*CRE_post + 0.5*PLA_post - 0.5*CRE_pre - 0.5*PLA_pre",
  Resid                  = "0.5*CRE_post + 0.5*PLA_post - H_pre",
  Baseline_Supplement    = "CRE_pre - PLA_pre",
  Training_CRE           = "CRE_post - CRE_pre",
  Training_PLA           = "PLA_post - PLA_pre",
  Supplement_Interaction = "CRE_post - CRE_pre - PLA_post + PLA_pre"
)
cm <- limma::makeContrasts(contrasts = contrasts, levels = dal$design$design_matrix)
colnames(cm) <- names(contrasts)
dal$design$contrast_matrix <- cm
dal$design$contrast_vector <- unname(contrasts)

dal <- fit_limma_model(dal)
dal <- extract_DA_results(dal, pval_thresh = 0.10, lfc_thresh = 0, adj_method = "BH")

# Tables + plots
# Native proteoDA outputs; drop the loose CSVs and keep the consolidated workbook.

write_limma_tables(dal,
  output_dir = out_dir, overwrite = TRUE,
  annot_cols = c("uniprot_id", "gene", "protein", "description")
)
unlink(c(
  file.path(out_dir, "combined_results.csv"), file.path(out_dir, "DA_summary.csv"),
  file.path(out_dir, "per_contrast_results")
), recursive = TRUE)
write_limma_plots(dal,
  grouping_column = "group_time", table_columns = c("uniprot_id", "gene"),
  output_dir = report_dir, overwrite = TRUE
)

# Pi-score
# Xiao 2014: Pi = P.Value^|logFC|; sig_pi = +1 up / -1 down / 0 ns at Pi < 0.05.

ann <- as_tibble(dal$annotation) |> select(any_of(c("uniprot_id", "gene", "protein", "description")))
res <- imap(dal$results, function(r, cname) {
  as_tibble(r, rownames = "uniprot_id") |>
    mutate(
      pi_score = P.Value^abs(logFC),
      sig_pi = case_when(
        pi_score < 0.05 & logFC > 0 ~ 1L,
        pi_score < 0.05 & logFC < 0 ~ -1L, TRUE ~ 0L
      ),
      contrast = cname
    ) |>
    left_join(ann, by = "uniprot_id")
})
write_csv(bind_rows(res), file.path(out_dir, "combined_results_pi.csv"))

cat(sprintf(
  "duplicateCorrelation rho = %.3f | proteoDA tables + plots + pi written\n",
  dal$eBayes_fit$correlation %||% NA_real_
))
