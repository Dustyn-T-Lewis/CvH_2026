# CvH DEP — Both models: CRvH (all subjects) and CR (supplement effects)
# Refs: Ritchie 2015 (limma), Smyth 2005 (dupCor), Law 2020 (simple effects),
#       Xiao 2014 (Capital Pi: P^|logFC|, Bioinformatics 30:801)

library(dplyr)
library(readr)
library(tibble)
library(proteoDA)
library(openxlsx)

setwd(rprojroot::find_rstudio_root_file())
source("R/cvh_design.R")
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

build_da_summary <- function(results_list, cfg) {
  lapply(names(results_list), function(cname) {
    r <- results_list[[cname]]
    up <- r$logFC > 0; dn <- r$logFC < 0
    s <- function(x) sum(x, na.rm = TRUE)
    tibble::tibble(
      contrast   = cname,
      type       = c("up", "down", "nonsig"),
      sig.PVal   = c(s(r$P.Value < cfg$pval_thresh & up),
                     s(r$P.Value < cfg$pval_thresh & dn),
                     s(r$P.Value >= cfg$pval_thresh)),
      sig.FDR    = c(s(r$adj.P.Val < cfg$pval_thresh & up),
                     s(r$adj.P.Val < cfg$pval_thresh & dn),
                     s(r$adj.P.Val >= cfg$pval_thresh)),
      sig.Pi     = c(s(r$sig_pi == 1L), s(r$sig_pi == -1L), s(r$sig_pi == 0L)),
      sig.FDR.05 = c(s(r$adj.P.Val < 0.05 & up),
                     s(r$adj.P.Val < 0.05 & dn),
                     s(r$adj.P.Val >= 0.05)),
      sig.FDR.10 = c(s(r$adj.P.Val < 0.10 & up),
                     s(r$adj.P.Val < 0.10 & dn),
                     s(r$adj.P.Val >= 0.10)),
      pval_thresh = cfg$pval_thresh,
      lfc_thresh  = cfg$lfc_thresh,
      p_adj_method = cfg$adj_method)
  }) |> dplyr::bind_rows()
}
set.seed(42)

shared_cfg <- list(
  norm_csv    = "01_normalization/c_data/02_normalized.csv",
  norm_rds    = "01_normalization/c_data/03_DAList_normalized.rds",
  data_dir    = "03_DEP/c_data",
  per_dir     = "03_DEP/c_data/04_per_contrast_results",
  pval_thresh = 0.10,
  lfc_thresh  = 0,
  adj_method  = "BH",
  pi_thresh   = 0.05
)

models <- list(
  CRvH = list(
    model_tag    = "CRvH",
    proteoDA_dir = "03_DEP/b_reports/01_proteoDA_CRvH",
    rds_name     = "01_limma_DAList_CRvH.rds",
    subset_fn    = function(meta, mat) list(meta = meta, mat = mat),
    levels       = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1"),
    contrasts    = c(
      "Cancer_vs_Healthy = (CRE_T1 + PLA_T1)/2 - H_T1",
      "Training_CR = (CRE_T2 + PLA_T2)/2 - (CRE_T1 + PLA_T1)/2")
  ),
  CR = list(
    model_tag    = "CR",
    proteoDA_dir = "03_DEP/b_reports/02_proteoDA_CR",
    rds_name     = "01_limma_DAList_CR.rds",
    subset_fn    = function(meta, mat) {
      idx <- which(meta$group != "H_T1")
      list(meta = meta[idx, ], mat = mat[, idx])
    },
    levels       = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2"),
    # Law et al. 2020, F1000Res 9:1444
    contrasts    = c(
      "Baseline_Supplement = CRE_T1 - PLA_T1",
      "Training_CRE = CRE_T2 - CRE_T1",
      "Training_PLA = PLA_T2 - PLA_T1",
      "Supplement_Interaction = (CRE_T2 - CRE_T1) - (PLA_T2 - PLA_T1)")
  )
)

dir.create(shared_cfg$data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(shared_cfg$per_dir,  recursive = TRUE, showWarnings = FALSE)

# --- Load shared data once ---
df <- read_csv(shared_cfg$norm_csv, show_col_types = FALSE)
ann_cols   <- c("uniprot_id", "protein", "gene", "description")
ann        <- df[, ann_cols]
mat_all    <- as.matrix(df[, setdiff(names(df), ann_cols)])
rownames(mat_all) <- ann$uniprot_id

dal_norm <- readRDS(shared_cfg$norm_rds)
dal_meta <- as.data.frame(dal_norm$metadata)
assert_cvh_design_rules(
  dal_meta[, c("Col_ID", "Subject_ID", "Group", "Timepoint", "Group_Time", "Supplement")],
  context = "03_DEP normalized metadata",
  allow_t2_only_singletons = TRUE
)
meta_all <- tibble::tibble(
  sample_id  = dal_meta$Col_ID,
  group      = dal_meta$Group_Time,
  subject    = dal_meta$Subject_ID,
  timepoint  = dal_meta$Timepoint,
  supplement = dal_meta$Supplement
)
stopifnot(setequal(colnames(mat_all), meta_all$sample_id))

# --- Run each model ---
for (m_name in names(models)) {
  m   <- models[[m_name]]
  cfg <- c(shared_cfg, m[c("model_tag", "proteoDA_dir")])
  dir.create(cfg$proteoDA_dir, recursive = TRUE, showWarnings = FALSE)

  sub  <- m$subset_fn(meta_all, mat_all)
  meta <- sub$meta
  mat  <- sub$mat
  missing_levels <- setdiff(m$levels, unique(meta$group))
  if (length(missing_levels) > 0) {
    stop(
      "Model ", m$model_tag, " is missing required Group_Time levels: ",
      paste(missing_levels, collapse = ", ")
    )
  }
  meta$group <- factor(meta$group, levels = m$levels)

  cat(sprintf("\n=== Model: %s (%d proteins x %d samples) ===\n",
              m$model_tag, nrow(mat), ncol(mat)))
  print(table(meta$group))

  meta_df <- as.data.frame(meta)
  rownames(meta_df) <- meta$sample_id
  dal <- DAList(data = mat, annotation = as.data.frame(ann),
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
  saveRDS(dal, file.path(cfg$data_dir, m$rds_name))

  tryCatch(
    write_limma_plots(dal, grouping_column = "group", output_dir = cfg$proteoDA_dir,
                      table_columns = c("uniprot_id", "gene", "protein"),
                      title_column = "gene", overwrite = TRUE),
    error = function(e) cat(sprintf("write_limma_plots: %s\n", conditionMessage(e))))

  write_limma_tables(dal, output_dir = cfg$data_dir, overwrite = TRUE,
    contrasts_subdir = "04_per_contrast_results",
    summary_csv       = sprintf("02_DA_summary_%s.csv", cfg$model_tag),
    combined_file_csv = sprintf("03_combined_results_%s.csv", cfg$model_tag),
    spreadsheet_xlsx  = sprintf("05_results_%s.xlsx", cfg$model_tag))

  # Pi-score enrichment (idempotent: drops existing pi columns before join)
  results_pi <- compute_pi_scores(dal$results, pi_thresh = cfg$pi_thresh)

  combined <- read_csv(file.path(cfg$data_dir,
    sprintf("03_combined_results_%s.csv", cfg$model_tag)), show_col_types = FALSE) |>
    select(-matches("^pi_score_|^sig_pi_"))
  for (cname in names(results_pi)) {
    pi_df <- results_pi[[cname]] |> select(uniprot_id, pi_score, sig_pi)
    names(pi_df)[2:3] <- paste0(names(pi_df)[2:3], "_", cname)
    combined <- left_join(combined, pi_df, by = "uniprot_id")
  }
  write_csv(combined, file.path(cfg$data_dir,
    sprintf("03_combined_results_%s.csv", cfg$model_tag)))

  for (cname in names(results_pi)) {
    pc_path <- file.path(cfg$per_dir, paste0(cname, ".csv"))
    pc <- read_csv(pc_path, show_col_types = FALSE) |> select(-any_of(c("pi_score", "sig_pi")))
    pi_df <- results_pi[[cname]] |> select(uniprot_id, pi_score, sig_pi)
    pc <- left_join(pc, pi_df, by = "uniprot_id")
    write_csv(pc, pc_path)
  }

  da_summary <- build_da_summary(results_pi, cfg)
  write_csv(da_summary, file.path(cfg$data_dir,
    sprintf("02_DA_summary_%s.csv", cfg$model_tag)))

  print(dal$design$contrast_matrix)
  print(da_summary)
  cat(sprintf("Done: %s — %d contrasts\n", m$model_tag, length(dal$results)))
}

writeLines(capture.output(sessionInfo()), file.path(shared_cfg$data_dir, "sessionInfo.txt"))
