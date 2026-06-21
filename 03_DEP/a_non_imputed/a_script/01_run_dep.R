#!/usr/bin/env Rscript
# CvH PRIMARY DEP on the NON-imputed normalized matrix (limma handles per-protein NAs).
#   ~ 0 + group_time + (1 | Subject_ID), group_time {H_pre, CR_pre, CR_post}.
#   Contrasts: CRvH_Baseline / CR_Training / Resid.
#   Significance: Pi-score (Xiao 2014, Pi = P.Value^|logFC|) < 0.05 + BH-FDR.

pacman::p_load(proteoDA, here, readr, dplyr, tibble, purrr)

out_dir    <- here("03_DEP", "a_non_imputed", "c_data")
report_dir <- here("03_DEP", "a_non_imputed", "b_reports")
clear_dir <- function(d) { dir.create(d, recursive = TRUE, showWarnings = FALSE)
  unlink(setdiff(list.files(d, full.names = TRUE), file.path(d, ".gitkeep")), recursive = TRUE) }
clear_dir(out_dir); clear_dir(report_dir)

dal <- readRDS(here("02_Normalization", "c_data", "DAList_normalized.rds"))
cat(sprintf("NON-imputed DEP: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

#### Fit limma model ####
# (1 | Subject_ID) carries the within-subject pre->post pairing via duplicateCorrelation.

dal <- add_design(dal, "~ 0 + group_time + (1 | Subject_ID)")
dal <- add_contrasts(dal, contrasts_vector = c(
  "CRvH_Baseline = CR_pre - H_pre",     # D: disease deviation from healthy
  "CR_Training   = CR_post - CR_pre",   # T: effect of training
  "Resid         = CR_post - H_pre"))   # R: what remains after training (R = D + T)
dal <- fit_limma_model(dal)
dal <- extract_DA_results(dal, pval_thresh = 0.10, lfc_thresh = 0, adj_method = "BH")

#### Tables + plots ####
# Native proteoDA outputs; drop the loose CSVs and keep the consolidated workbook.

write_limma_tables(dal, output_dir = out_dir, overwrite = TRUE,
                   annot_cols = c("uniprot_id", "gene", "protein", "description"))
unlink(c(file.path(out_dir, "combined_results.csv"), file.path(out_dir, "DA_summary.csv"),
         file.path(out_dir, "per_contrast_results")), recursive = TRUE)
write_limma_plots(dal, grouping_column = "group_time", table_columns = c("uniprot_id", "gene"),
                  output_dir = report_dir, overwrite = TRUE)

#### Pi-score ####
# Xiao 2014: Pi = P.Value^|logFC|; sig_pi = +1 up / -1 down / 0 ns at Pi < 0.05.

ann <- as_tibble(dal$annotation) |> select(any_of(c("uniprot_id", "gene", "protein", "description")))
res <- imap(dal$results, function(r, cname)
  as_tibble(r, rownames = "uniprot_id") |>
    mutate(pi_score = P.Value ^ abs(logFC),
           sig_pi = case_when(pi_score < 0.05 & logFC > 0 ~ 1L,
                              pi_score < 0.05 & logFC < 0 ~ -1L, TRUE ~ 0L),
           contrast = cname) |>
    left_join(ann, by = "uniprot_id"))
write_csv(bind_rows(res), file.path(out_dir, "combined_results_pi.csv"))

cat(sprintf("duplicateCorrelation rho = %.3f | proteoDA tables + plots + pi written\n",
            dal$eBayes_fit$correlation %||% NA_real_))
