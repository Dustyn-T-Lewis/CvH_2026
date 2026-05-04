# Run our two-model proteoDA/limma DEP on Sam's normalized DAList.
# Skips our normalization stage (Sam used cycloess, matches us).
# No imputation (matches Sam — limma handles NAs per-protein).
# Outputs go to 02-03_Sam's_Results/our_rerun/c_data/ to keep them
# isolated from our main pipeline outputs in 03_DEP/.

library(dplyr)
library(readr)
library(tibble)
library(proteoDA)
library(openxlsx)
library(stringr)

setwd(rprojroot::find_rstudio_root_file())
source("R/cvh_design.R")
set.seed(42)

cfg <- list(
  sam_rds     = "02-03_Sam's_Results/01_normalized_DAList_SURV_stringent_muscle.RDS",
  out_dir     = "02-03_Sam's_Results/our_rerun/c_data",
  per_dir     = "02-03_Sam's_Results/our_rerun/c_data/04_per_contrast_results",
  pval_thresh = 0.10,
  lfc_thresh  = 0,
  adj_method  = "BH",
  pi_thresh   = 0.05
)
dir.create(cfg$per_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load Sam's DAList and reconcile metadata to our schema -----------------

sam <- readRDS(cfg$sam_rds)
mat_all <- sam$data

ann_sam <- as.data.frame(sam$annotation)[, c("uniprot_id", "protein", "gene", "description")]

# Build CvH-schema metadata from Sam's columns and validate the
# repeated-measures structure before fitting. derive_cvh_analysis_metadata
# defaults to rejecting T2-only singletons; Sam's data has them (CR007_T2
# etc), so we replicate the derivation and pass the flag explicitly.
md <- as.data.frame(sam$metadata)
cvh_meta <- tibble::tibble(
  Col_ID     = normalize_cvh_col_id(md$sample_id),
  Subject_ID = normalize_cvh_subject_id(md$pid),
  Timepoint  = stringr::str_trim(md$timepoint),
  Supplement = dplyr::na_if(stringr::str_trim(md$supp), ""),
  Cancer     = stringr::str_trim(md$cancer),
  Group = dplyr::case_when(
    md$cancer == "SURV" & md$supp == "CRE" ~ "CR_CRE",
    md$cancer == "SURV" & md$supp == "PLA" ~ "CR_PLA",
    md$cancer == "CTL"                     ~ "PPS",
    TRUE ~ NA_character_
  ),
  Group_Time = dplyr::case_when(
    md$cancer == "SURV" & md$supp == "CRE" ~ paste0("CRE_", md$timepoint),
    md$cancer == "SURV" & md$supp == "PLA" ~ paste0("PLA_", md$timepoint),
    md$cancer == "CTL"  & md$timepoint == "T1" ~ "H_T1",
    TRUE ~ NA_character_
  )
)
stopifnot(all(!is.na(cvh_meta$Group)), all(!is.na(cvh_meta$Group_Time)))
assert_cvh_design_rules(
  cvh_meta[, c("Col_ID", "Subject_ID", "Group", "Timepoint", "Group_Time", "Supplement")],
  context = "Sam's reconciled metadata",
  allow_t2_only_singletons = TRUE
)

meta_all <- tibble::tibble(
  sample_id  = cvh_meta$Col_ID,
  subject    = cvh_meta$Subject_ID,
  group      = cvh_meta$Group_Time,
  timepoint  = cvh_meta$Timepoint,
  supplement = cvh_meta$Supplement
)

colnames(mat_all) <- normalize_cvh_col_id(colnames(mat_all))
stopifnot(setequal(meta_all$sample_id, colnames(mat_all)))
mat_all <- mat_all[, meta_all$sample_id]

cat(sprintf("Sam's data: %d proteins x %d samples\n", nrow(mat_all), ncol(mat_all)))
cat("Group counts:\n"); print(table(meta_all$group))

# --- Pi-score helper (idempotent) -------------------------------------------

compute_pi_scores <- function(results_list, pi_thresh = 0.05) {
  lapply(names(results_list), function(cname) {
    results_list[[cname]] |>
      tibble::rownames_to_column("uniprot_id") |>
      dplyr::mutate(
        pi_score = P.Value ^ abs(logFC),
        sig_pi = dplyr::case_when(
          pi_score < pi_thresh & logFC > 0 ~  1L,
          pi_score < pi_thresh & logFC < 0 ~ -1L,
          TRUE ~ 0L),
        contrast = cname)
  }) |> setNames(names(results_list))
}

# --- Two-model DEP (mirrors 03_DEP/a_script/01_run_dep.R) -------------------

models <- list(
  CRvH = list(
    model_tag = "CRvH",
    rds_name  = "01_limma_DAList_CRvH.rds",
    subset_fn = function(meta, mat) list(meta = meta, mat = mat),
    levels    = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1"),
    contrasts = c(
      "Cancer_vs_Healthy = (CRE_T1 + PLA_T1)/2 - H_T1",
      "Training_CR = (CRE_T2 + PLA_T2)/2 - (CRE_T1 + PLA_T1)/2")
  ),
  CR = list(
    model_tag = "CR",
    rds_name  = "01_limma_DAList_CR.rds",
    subset_fn = function(meta, mat) {
      idx <- which(meta$group != "H_T1")
      list(meta = meta[idx, ], mat = mat[, idx])
    },
    levels    = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2"),
    contrasts = c(
      "Baseline_Supplement = CRE_T1 - PLA_T1",
      "Training_CRE = CRE_T2 - CRE_T1",
      "Training_PLA = PLA_T2 - PLA_T1",
      "Supplement_Interaction = (CRE_T2 - CRE_T1) - (PLA_T2 - PLA_T1)")
  )
)

for (m_name in names(models)) {
  m <- models[[m_name]]
  sub <- m$subset_fn(meta_all, mat_all)
  meta <- sub$meta
  mat  <- sub$mat
  missing_levels <- setdiff(m$levels, unique(meta$group))
  if (length(missing_levels) > 0) {
    stop("Model ", m$model_tag, " missing required levels: ",
         paste(missing_levels, collapse = ", "))
  }
  meta$group <- factor(meta$group, levels = m$levels)

  cat(sprintf("\n=== Model: %s (%d proteins x %d samples) ===\n",
              m$model_tag, nrow(mat), ncol(mat)))
  print(table(meta$group))

  meta_df <- as.data.frame(meta)
  rownames(meta_df) <- meta$sample_id
  dal <- DAList(data = mat, annotation = ann_sam,
                metadata = meta_df, tags = list(norm_method = "cycloess"))

  dal <- add_design(dal, "~ 0 + group + (1 | subject)")
  colnames(dal$design$design_matrix) <- gsub("^group", "",
                                              colnames(dal$design$design_matrix))
  dal <- add_contrasts(dal, contrasts_vector = m$contrasts)
  dal <- fit_limma_model(dal)

  within_cor <- dal$eBayes_fit$correlation %||%
    dal$tags$duplicate_correlation %||% NA_real_
  if (!is.na(within_cor)) cat(sprintf("Within-subject correlation: %.3f\n", within_cor))

  dal <- extract_DA_results(dal, pval_thresh = cfg$pval_thresh,
                            lfc_thresh = cfg$lfc_thresh, adj_method = cfg$adj_method)
  saveRDS(dal, file.path(cfg$out_dir, m$rds_name))

  write_limma_tables(dal, output_dir = cfg$out_dir, overwrite = TRUE,
    contrasts_subdir = "04_per_contrast_results",
    summary_csv       = sprintf("02_DA_summary_%s.csv", m$model_tag),
    combined_file_csv = sprintf("03_combined_results_%s.csv", m$model_tag),
    spreadsheet_xlsx  = sprintf("05_results_%s.xlsx", m$model_tag))

  results_pi <- compute_pi_scores(dal$results, pi_thresh = cfg$pi_thresh)

  for (cname in names(results_pi)) {
    pc_path <- file.path(cfg$per_dir, paste0(cname, ".csv"))
    pc <- read_csv(pc_path, show_col_types = FALSE) |>
      dplyr::select(-dplyr::any_of(c("pi_score", "sig_pi")))
    pi_df <- results_pi[[cname]] |> dplyr::select(uniprot_id, pi_score, sig_pi)
    pc <- dplyr::left_join(pc, pi_df, by = "uniprot_id")
    write_csv(pc, pc_path)
  }

  cat(sprintf("Done: %s -- %d contrasts\n", m$model_tag, length(dal$results)))
}

writeLines(capture.output(sessionInfo()),
           file.path(cfg$out_dir, "sessionInfo.txt"))
cat("\nAll outputs in:", cfg$out_dir, "\n")
