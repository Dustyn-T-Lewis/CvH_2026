# CvH DEP Overview — Multi-threshold significance table with/without outlier removal
# Compares FDR<0.05, FDR<0.10, Pi<0.05, P<0.05, P<0.01 across both models
# Outputs: 03_DEP/c_data/13_DEP_overview.csv, 14_DEP_overview.xlsx

library(proteoDA)
library(readxl)
library(readr)
library(dplyr)
library(MsCoreUtils)
library(openxlsx)

set.seed(42)
setwd(rprojroot::find_rstudio_root_file())
source("R/cvh_design.R")

`%||%` <- function(x, y) if (!is.null(x)) x else y

# BLOOD_CONTAMINANTS sourced from R/cvh_design.R (Geyer 2016 + HPA Ig)

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

add_sheet <- function(wb, name, df, title = NULL, notes = NULL) {
  openxlsx::addWorksheet(wb, name)
  start_row <- 1L
  if (!is.null(title)) {
    openxlsx::writeData(wb, name, title, startRow = 1)
    openxlsx::addStyle(wb, name,
      openxlsx::createStyle(textDecoration = "bold", fontSize = 12),
      rows = 1, cols = 1)
    start_row <- start_row + 1L
  }
  if (!is.null(notes)) {
    for (i in seq_along(notes)) {
      openxlsx::writeData(wb, name, notes[i], startRow = start_row)
      openxlsx::addStyle(wb, name,
        openxlsx::createStyle(fontSize = 10, fontColour = "#555555", wrapText = TRUE),
        rows = start_row, cols = 1)
      start_row <- start_row + 1L
    }
    start_row <- start_row + 1L
  }
  hs <- openxlsx::createStyle(textDecoration = "bold", border = "Bottom",
                               fgFill = "#DCE6F1")
  openxlsx::writeData(wb, name, df, startRow = start_row, headerStyle = hs)
  openxlsx::freezePane(wb, name, firstActiveRow = start_row + 1L,
                        firstActiveCol = 2)
  openxlsx::setColWidths(wb, name, cols = seq_len(ncol(df)), widths = "auto")
}

# =============================================================================
# SHARED CONFIG
# =============================================================================

cfg <- list(
  raw_file    = "00_input/CvH_raw.xlsx",
  meta_file   = "00_input/CvH_meta.csv",
  pheno_file  = "00_input/CRm_meta.csv",
  hpa_file    = "00_input/HPA_skeletal_muscle_annotations.tsv",
  norm_rds    = "01_normalization/c_data/03_DAList_normalized.rds",
  norm_int_rds = "01_normalization/c_data/00_report_intermediates.rds",
  min_reps    = 5L,
  min_groups  = 1L,
  outlier_k   = 3,
  mahal_p     = 0.01,
  norm_method = "cycloess",
  pi_thresh   = 0.05,
  pval_thresh = 0.10,
  lfc_thresh  = 0,
  adj_method  = "BH"
)

models <- list(
  CRvH = list(
    subset_fn = function(meta, mat) list(meta = meta, mat = mat),
    levels    = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1"),
    contrasts = c("Cancer_vs_Healthy = (CRE_T1 + PLA_T1)/2 - H_T1",
                  "Training_CR = (CRE_T2 + PLA_T2)/2 - (CRE_T1 + PLA_T1)/2")
  ),
  CR = list(
    subset_fn = function(meta, mat) {
      idx <- which(meta$group != "H_T1")
      list(meta = meta[idx, ], mat = mat[, idx])
    },
    levels    = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2"),
    contrasts = c("Baseline_Supplement = CRE_T1 - PLA_T1",
                  "Training_CRE = CRE_T2 - CRE_T1",
                  "Training_PLA = PLA_T2 - PLA_T1",
                  "Supplement_Interaction = (CRE_T2 - CRE_T1) - (PLA_T2 - PLA_T1)")
  )
)

# =============================================================================
# HELPER: Build full pipeline from raw data to DEP results
# =============================================================================

run_pca_simple <- function(mat, metadata) {
  if (anyNA(mat)) {
    stop("run_pca_simple expects a complete-case matrix")
  }
  mat <- log2(mat)
  pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  pc  <- as.data.frame(pca$x[, 1:3]) |> mutate(Col_ID = rownames(pca$x))
  list(pca = pca, scores = pc)
}

build_pipeline <- function(cfg, remove_outliers = TRUE) {
  if (remove_outliers && file.exists(cfg$norm_rds) && file.exists(cfg$norm_int_rds)) {
    dal <- readRDS(cfg$norm_rds)
    norm_int <- readRDS(cfg$norm_int_rds)
    outlier_ids <- norm_int$outlier_ids %||% character(0)
    assert_cvh_design_rules(
      dal$metadata[, c("Col_ID", "Subject_ID", "Group", "Timepoint", "Group_Time", "Supplement")],
      context = "04_dep_overview normalized branch",
      allow_t2_only_singletons = TRUE
    )
    dal$data <- impute_matrix(dal$data, method = "bpca")
    n_prot <- nrow(dal$data)
    n_samp <- ncol(dal$data)
    outlier_label <- if (length(outlier_ids) > 0) paste(outlier_ids, collapse = ", ") else "none"
    cat(sprintf("  Pipeline: %d proteins x %d samples (outliers removed: %s)\n",
                n_prot, n_samp, outlier_label))
    return(list(dal = dal, n_prot = n_prot, n_samp = n_samp, outlier_ids = outlier_ids))
  }

  # 1. Load raw
  raw <- read_excel(cfg$raw_file)
  annot_cols <- c("uniprot_id", "protein", "gene", "description", "n_seq")
  annotation <- raw[, annot_cols]
  intensity  <- raw[, setdiff(names(raw), annot_cols)]
  metadata   <- as.data.frame(load_cvh_analysis_metadata(
    meta_file = cfg$meta_file,
    pheno_file = cfg$pheno_file,
    raw_sample_ids = colnames(intensity)
  ))
  rownames(metadata) <- metadata$Col_ID
  intensity <- intensity[, metadata$Col_ID]

  # 2. HPA tissue filter
  hpa <- read_tsv(cfg$hpa_file, show_col_types = FALSE) |>
    select(Gene, Ensembl, Evidence,
           Protein_class    = `Protein class`,
           Subcellular_main = `Subcellular main location`,
           Interactions) |>
    distinct(Gene, .keep_all = TRUE)
  keep_hpa   <- annotation$gene %in% hpa$Gene
  intensity  <- intensity[keep_hpa, ]
  annotation <- annotation[keep_hpa, ] |> left_join(hpa, by = c("gene" = "Gene"))

  hpa_ig <- hpa$Gene[grepl("Immunoglobulin genes", hpa$Protein_class, fixed = TRUE)]
  blood_genes <- unique(c(BLOOD_CONTAMINANTS, hpa_ig))
  keep_blood <- !annotation$gene %in% blood_genes
  intensity  <- intensity[keep_blood, ]
  annotation <- annotation[keep_blood, ]

  # 3. Dedup
  if (any(duplicated(annotation$uniprot_id))) {
    annotation$row_mean <- rowMeans(data.matrix(intensity), na.rm = TRUE)
    keep_idx <- annotation |>
      mutate(row_idx = row_number()) |> group_by(uniprot_id) |>
      slice_max(row_mean, n = 1, with_ties = FALSE) |> pull(row_idx)
    annotation <- annotation[keep_idx, ]; intensity <- intensity[keep_idx, ]
    annotation$row_mean <- NULL
  }

  # 4. DAList
  int_mat <- as.data.frame(data.matrix(intensity))
  rownames(int_mat) <- annotation$uniprot_id
  annot_df <- as.data.frame(annotation); rownames(annot_df) <- annotation$uniprot_id
  meta_df  <- as.data.frame(metadata);   rownames(meta_df)  <- metadata$Col_ID
  dal <- DAList(data = int_mat, annotation = annot_df, metadata = meta_df)
  dal <- zero_to_missing(dal)

  # 5. Outlier removal (conditional)
  outlier_ids <- character(0)
  if (remove_outliers) {
    pct_missing <- colMeans(is.na(dal$data)) * 100
    miss_info <- dal$metadata |> select(Col_ID, Subject_ID, Group, Timepoint, Group_Time)
    paired_subjects <- miss_info |> filter(Group != "PPS") |>
      count(Subject_ID) |> filter(n == 2) |> pull(Subject_ID)
    miss_info$pct_missing <- pct_missing[miss_info$Col_ID]
    miss_info$delta_missing <- NA_real_
    for (subj in paired_subjects) {
      rows <- miss_info |> filter(Subject_ID == subj)
      t1 <- rows$Col_ID[rows$Timepoint == "T1"]
      t2 <- rows$Col_ID[rows$Timepoint == "T2"]
      if (length(t1) == 1 && length(t2) == 1) {
        d <- abs(pct_missing[t2] - pct_missing[t1])
        miss_info$delta_missing[miss_info$Col_ID == t1] <- d
        miss_info$delta_missing[miss_info$Col_ID == t2] <- d
      }
    }
    miss_thresh  <- quantile(pct_missing, 0.75) + 1.5 * IQR(pct_missing)
    delta_vals   <- miss_info$delta_missing[!is.na(miss_info$delta_missing)]
    delta_thresh <- if (length(delta_vals) > 2)
      quantile(delta_vals, 0.75) + 1.5 * IQR(delta_vals) else Inf
    miss_info$miss_flag <- miss_info$pct_missing > miss_thresh |
      (!is.na(miss_info$delta_missing) & miss_info$delta_missing > delta_thresh)

    samp_med   <- apply(log2(dal$data), 2, median, na.rm = TRUE)
    global_med <- median(samp_med); mad_val <- mad(samp_med)
    mad_flags  <- tibble(Col_ID = names(samp_med), sample_median = samp_med,
                         mad_flag = abs(samp_med - global_med) > cfg$outlier_k * mad_val)

    dal_filt <- filter_proteins_by_group(dal, min_reps = cfg$min_reps,
                                          min_groups = cfg$min_groups,
                                          grouping_column = "Group_Time")
    complete_mat <- dal$data[rowSums(is.na(dal$data)) == 0, , drop = FALSE]
    pca_obj   <- run_pca_simple(complete_mat, dal$metadata)
    pc3       <- pca_obj$pca$x[, 1:3]
    mahal     <- mahalanobis(pc3, colMeans(pc3), cov(pc3))
    pca_flags <- tibble(Col_ID = colnames(dal$data), mahal_dist = mahal,
                        pca_flag = mahal > qchisq(1 - cfg$mahal_p, df = 3))

    cor_mat     <- cor(log2(dal_filt$data), use = "pairwise.complete.obs")
    med_cor     <- apply(cor_mat, 2, function(x) median(x[x < 1], na.rm = TRUE))
    cor_flags   <- tibble(Col_ID = names(med_cor), median_cor = med_cor,
                          cor_flag = med_cor < median(med_cor) - cfg$outlier_k * mad(med_cor))

    outlier_diag <- miss_info |>
      left_join(mad_flags, by = "Col_ID") |>
      left_join(pca_flags |> select(Col_ID, pca_flag), by = "Col_ID") |>
      left_join(cor_flags |> select(Col_ID, cor_flag), by = "Col_ID") |>
      mutate(n_flags = as.integer(miss_flag) + as.integer(mad_flag) +
               as.integer(pca_flag) + as.integer(cor_flag),
             consensus_outlier = n_flags >= 3)
    outlier_ids <- outlier_diag |> filter(consensus_outlier) |> pull(Col_ID)
    if (length(outlier_ids) > 0) {
      keep_cols <- setdiff(colnames(dal$data), outlier_ids)
      dal$data <- dal$data[, keep_cols, drop = FALSE]
      dal$metadata <- dal$metadata[match(keep_cols, dal$metadata$Col_ID), , drop = FALSE]
    }
  }

  # 6. Protein filter
  dal <- filter_proteins_by_group(dal, min_reps = cfg$min_reps,
                                   min_groups = cfg$min_groups,
                                   grouping_column = "Group_Time")

  # 7. Normalize
  dal <- normalize_data(dal, norm_method = cfg$norm_method)

  # 8. Impute with bpca (faster than missForest; used for sensitivity analysis only)
  mat <- dal$data
  mat_imp <- impute_matrix(mat, method = "bpca")
  dal$data <- mat_imp

  n_prot <- nrow(dal$data)
  n_samp <- ncol(dal$data)
  outlier_label <- if (length(outlier_ids) > 0) paste(outlier_ids, collapse = ", ") else "none"
  cat(sprintf("  Pipeline: %d proteins x %d samples (outliers removed: %s)\n",
              n_prot, n_samp, outlier_label))

  list(dal = dal, n_prot = n_prot, n_samp = n_samp, outlier_ids = outlier_ids)
}

# =============================================================================
# HELPER: Fit models and extract multi-threshold counts
# =============================================================================

count_sig <- function(r, pi_r) {
  up <- r$logFC > 0; dn <- r$logFC < 0
  s <- function(x) sum(x, na.rm = TRUE)
  tibble(
    FDR_0.05_up   = s(r$adj.P.Val < 0.05 & up),
    FDR_0.05_down = s(r$adj.P.Val < 0.05 & dn),
    FDR_0.10_up   = s(r$adj.P.Val < 0.10 & up),
    FDR_0.10_down = s(r$adj.P.Val < 0.10 & dn),
    Pi_0.05_up    = s(pi_r$sig_pi == 1L),
    Pi_0.05_down  = s(pi_r$sig_pi == -1L),
    P_0.05_up     = s(r$P.Value < 0.05 & up),
    P_0.05_down   = s(r$P.Value < 0.05 & dn),
    P_0.01_up     = s(r$P.Value < 0.01 & up),
    P_0.01_down   = s(r$P.Value < 0.01 & dn)
  )
}

run_models <- function(pipeline_out, models, cfg, condition_label) {
  dal_full <- pipeline_out$dal
  assert_cvh_design_rules(
    dal_full$metadata[, c("Col_ID", "Subject_ID", "Group", "Timepoint", "Group_Time", "Supplement")],
    context = paste("04_dep_overview", condition_label),
    allow_t2_only_singletons = TRUE
  )
  meta_all <- tibble::tibble(
    sample_id  = dal_full$metadata$Col_ID,
    group      = dal_full$metadata$Group_Time,
    subject    = dal_full$metadata$Subject_ID,
    timepoint  = dal_full$metadata$Timepoint,
    supplement = dal_full$metadata$Supplement
  )
  mat_all <- dal_full$data
  ann     <- dal_full$annotation

  all_results <- list()

  for (m_name in names(models)) {
    m   <- models[[m_name]]
    sub <- m$subset_fn(meta_all, mat_all)
    meta <- sub$meta; mat <- sub$mat
    missing_levels <- setdiff(m$levels, unique(meta$group))
    if (length(missing_levels) > 0) {
      stop(
        "Overview model ", m_name, " is missing required Group_Time levels: ",
        paste(missing_levels, collapse = ", ")
      )
    }
    meta$group <- factor(meta$group, levels = m$levels)

    meta_df <- as.data.frame(meta); rownames(meta_df) <- meta$sample_id
    dal <- DAList(data = mat, annotation = as.data.frame(ann),
                  metadata = meta_df, tags = list(norm_method = "cycloess"))
    dal <- add_design(dal, "~ 0 + group + (1 | subject)")
    colnames(dal$design$design_matrix) <- gsub("^group", "",
                                                colnames(dal$design$design_matrix))
    dal <- add_contrasts(dal, contrasts_vector = m$contrasts)
    dal <- fit_limma_model(dal)
    dal <- extract_DA_results(dal, pval_thresh = cfg$pval_thresh,
                              lfc_thresh = cfg$lfc_thresh, adj_method = cfg$adj_method)

    results_pi <- compute_pi_scores(dal$results, pi_thresh = cfg$pi_thresh)

    for (cname in names(dal$results)) {
      row <- count_sig(dal$results[[cname]], results_pi[[cname]])
      row$model     <- m_name
      row$contrast  <- cname
      row$condition <- condition_label
      row$n_proteins <- nrow(mat)
      row$n_samples  <- ncol(mat)
      cor_val <- dal$eBayes_fit$correlation %||%
        dal$tags$duplicate_correlation %||% NA_real_
      row$within_cor <- round(cor_val, 4)
      all_results[[length(all_results) + 1]] <- row
    }
    cat(sprintf("  %s/%s: done (%d contrasts)\n", condition_label, m_name, length(dal$results)))
  }

  bind_rows(all_results)
}

# =============================================================================
# RUN BOTH CONDITIONS
# =============================================================================

cat("=== WITH outlier removal ===\n")
pipe_with <- build_pipeline(cfg, remove_outliers = TRUE)
res_with  <- run_models(pipe_with, models, cfg, "Outliers_removed")

cat("\n=== WITHOUT outlier removal ===\n")
pipe_without <- build_pipeline(cfg, remove_outliers = FALSE)
res_without  <- run_models(pipe_without, models, cfg, "All_samples")

# =============================================================================
# ASSEMBLE OVERVIEW TABLE
# =============================================================================

overview <- bind_rows(res_with, res_without) |>
  select(condition, model, contrast, n_proteins, n_samples, within_cor,
         FDR_0.05_up, FDR_0.05_down,
         FDR_0.10_up, FDR_0.10_down,
         Pi_0.05_up,  Pi_0.05_down,
         P_0.05_up,   P_0.05_down,
         P_0.01_up,   P_0.01_down) |>
  arrange(model, contrast, condition)

out_dir <- "03_DEP/c_data"
write_csv(overview, file.path(out_dir, "13_DEP_overview.csv"))

# Excel with formatting
wb <- createWorkbook()
add_sheet(wb, "Overview", overview,
  title = "DEP Significance Overview: With vs Without Outlier Removal",
  notes = c(
    "FDR = Benjamini-Hochberg adjusted p-value | Capital Pi = P.Value^|logFC| (Xiao 2014, Bioinformatics 30:801)",
    "P = unadjusted p-value | up/down = direction of logFC",
    sprintf("Outliers removed: %s",
            if (length(pipe_with$outlier_ids) > 0) paste(pipe_with$outlier_ids, collapse = ", ") else "none"),
    sprintf("Pipeline: HPA filter -> blood removal -> dedup -> protein filter (>=%d reps in >=%d group) -> cycloess -> bpca imputation -> limma+dupCor",
            cfg$min_reps, cfg$min_groups)))
saveWorkbook(wb, file.path(out_dir, "14_DEP_overview.xlsx"), overwrite = TRUE)

print(overview, n = Inf, width = Inf)
cat(sprintf("\nDone: %s/13_DEP_overview.csv + 14_DEP_overview.xlsx\n", out_dir))
