#!/usr/bin/env Rscript
# CvH_DEP.R — Consolidated DEP (auto-selects limma or DEqMS based on NAs)
# INPUT:  data_file, 00_input/CvH_meta.csv, 00_input/CvH_raw.xlsx (DEqMS only)
# OUTPUT: c_data/per_contrast_results/, combined_results.csv, DA_summary.csv
#         c_data/results.xlsx, results_list.rds, array_weights.csv
#         b_reports/ per-contrast volcanos + histograms

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(tidyverse, openxlsx, here)

base_dir   <- here::here()

# User config: set input path (non-imputed -> limma; imputed -> DEqMS)
data_file  <- file.path(base_dir, "01_normalization", "c_data", "01_normalized.csv")
# data_file <- file.path(base_dir, "02_Imputation", "c_data", "01_imputed.csv")

meta_file  <- file.path(base_dir, "00_input", "CvH_meta.csv")
raw_file   <- file.path(base_dir, "00_input", "CvH_raw.xlsx")
report_dir <- file.path(base_dir, "04_DEP", "b_reports")
data_dir   <- file.path(base_dir, "04_DEP", "c_data")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir,   recursive = TRUE, showWarnings = FALSE)

# --- 1. Load data & detect method ---
cat("Step 1: Load data\n")
stopifnot("Data file not found" = file.exists(data_file))
df <- read_csv(data_file, show_col_types = FALSE)

ann_cols   <- c("uniprot_id", "protein", "gene", "description", "n_seq")
ann        <- df[, intersect(ann_cols, names(df))]
samp_names <- setdiff(names(df), ann_cols)
mat        <- as.matrix(df[, samp_names])
rownames(mat) <- ann$uniprot_id

n_missing <- sum(is.na(mat))
HAS_NA    <- n_missing > 0
METHOD    <- if (HAS_NA) "limma" else "DEqMS"

cat(sprintf("  %d proteins x %d samples, %d NAs (%.1f%%)\n",
            nrow(mat), ncol(mat), n_missing, 100 * n_missing / length(mat)))
cat(sprintf("  Method: %s\n", METHOD))

meta <- tibble(sample_id = samp_names) |>
  left_join(read_csv(meta_file, show_col_types = FALSE), by = c("sample_id" = "Col_ID")) |>
  mutate(group = case_when(Group == "PPS" ~ "H_T1", TRUE ~ paste0(Supplement, "_", Timepoint)))
meta$group <- factor(meta$group, levels = c("CRE_T1","CRE_T2","PLA_T1","PLA_T2","H_T1"))

print(table(meta$group))
stopifnot(identical(colnames(mat), meta$sample_id))

# --- 2. Conditional library loading ---
if (METHOD == "limma") {
  pacman::p_load(proteoDA, limma)
} else {
  pacman::p_load(limma, DEqMS, readxl, ggrepel, patchwork)

  if ("n_seq" %in% names(ann)) {
    pep_count <- setNames(ann$n_seq, ann$uniprot_id)
  } else {
    ann <- left_join(ann,
      read_xlsx(raw_file) |> dplyr::select(uniprot_id, n_seq) |>
        distinct(uniprot_id, .keep_all = TRUE),
      by = "uniprot_id")
    pep_count <- setNames(ann$n_seq, ann$uniprot_id)
  }
  stopifnot(all(!is.na(pep_count)), all(pep_count >= 1))
  cat(sprintf("  Peptides: %d-%d (median %d)\n",
              min(pep_count), max(pep_count), median(pep_count)))
}

# --- 3. Contrast definitions ---
GROUP_LEVELS  <- c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1")
CONTRAST_SPECS <- c(
  "Training_CR = (CRE_T2 + PLA_T2)/2 - (CRE_T1 + PLA_T1)/2",
  "Cancer_vs_Healthy = (CRE_T1 + PLA_T1)/2 - H_T1",
  "Supplement_Interaction = (CRE_T2 - CRE_T1) - (PLA_T2 - PLA_T1)"
)

# --- 4. Model fitting ---
add_pi_flags <- function(df, pval_col = "P.Value", lfc_col = "logFC") {
  p <- df[[pval_col]]; lfc <- df[[lfc_col]]
  df$pi_score       <- -log10(p) * abs(lfc)
  df$pi_pval        <- rank(-df$pi_score) / nrow(df)
  df$sig_pi         <- ifelse(df$pi_pval < 0.05, sign(lfc), 0L)
  df$sig_nominal_01 <- ifelse(p < 0.01, sign(lfc), 0L)
  df
}

if (METHOD == "limma") {
  # --- 4A: limma via proteoDA ---
  cat("Step 2: DAList\n")
  meta_df <- as.data.frame(meta); rownames(meta_df) <- meta$sample_id
  dal <- DAList(data = mat, annotation = as.data.frame(ann), metadata = meta_df,
                tags = list(norm_method = "cycloess"))

  cat("Step 3: Design & contrasts\n")
  dal <- add_design(dal, "~ 0 + group + (1 | Subject_ID)")
  colnames(dal$design$design_matrix) <- gsub("^group", "", colnames(dal$design$design_matrix))
  dal <- add_contrasts(dal, contrasts_vector = CONTRAST_SPECS)

  design_mat   <- dal$design$design_matrix
  contrast_mat <- dal$design$contrast_matrix
  block_var    <- meta_df$Subject_ID

  cat("Step 4a: Array weights\n")
  aw <- arrayWeights(mat, design_mat)
  names(aw) <- colnames(mat)
  cat(sprintf("  Array weights: %.2f - %.2f (median %.2f)\n",
              min(aw), max(aw), median(aw)))
  write_csv(tibble(sample_id = names(aw), array_weight = round(aw, 4)),
            file.path(data_dir, "array_weights.csv"))

  cat("Step 4b: Fit model\n")
  dupcor <- duplicateCorrelation(mat, design_mat, block = block_var, weights = aw)
  cat(sprintf("  Within-subject correlation: %.3f\n", dupcor$consensus.correlation))

  fit <- lmFit(mat, design_mat, block = block_var,
               correlation = dupcor$consensus.correlation, weights = aw) |>
    contrasts.fit(contrast_mat) |>
    eBayes(robust = TRUE)

  dal$eBayes_fit <- fit
  dal$eBayes_fit$correlation <- dupcor$consensus.correlation

  cat("Step 4c: Extract results\n")
  dal <- extract_DA_results(dal, pval_thresh = 0.05, lfc_thresh = 0, adj_method = "BH")

  cat("Step 5: Tables\n")
  write_limma_tables(dal, output_dir = data_dir, overwrite = TRUE)

  contrast_names   <- names(dal$results)
  per_contrast_dir <- file.path(data_dir, "per_contrast_results")

  all_results <- list()
  for (cname in contrast_names) {
    fpath <- file.path(per_contrast_dir, paste0(cname, ".csv"))
    all_results[[cname]] <- read_csv(fpath, show_col_types = FALSE)
  }

  PVAL_COL     <- "P.Value"
  ADJ_PVAL_COL <- "adj.P.Val"

  cat("Step 5b: Plots\n")
  write_limma_plots(dal, grouping_column = "group", output_dir = report_dir,
                    table_columns = c("uniprot_id", "gene", "protein"),
                    title_column = "gene", overwrite = TRUE)

  static_dir <- file.path(report_dir, "static_plots")
  for (cname in contrast_names) {
    cdir <- file.path(report_dir, cname); dir.create(cdir, showWarnings = FALSE)
    html_f <- file.path(report_dir, paste0(cname, "_DA_report.html"))
    if (file.exists(html_f)) file.rename(html_f, file.path(cdir, basename(html_f)))
    pdfs <- list.files(static_dir, pattern = paste0("^", cname, "-"), full.names = TRUE)
    walk(pdfs, ~file.rename(.x, file.path(cdir, sub(paste0("^", cname, "-"), "", basename(.x)))))
  }
  if (length(list.files(static_dir)) == 0) unlink(static_dir, recursive = TRUE)

} else {
  # --- 4B: DEqMS ---
  cat("Step 2: Design & contrasts\n")
  design <- model.matrix(~ 0 + group, data = meta)
  colnames(design) <- gsub("^group", "", colnames(design))

  contrast_mat <- makeContrasts(
    Training_CR            = (CRE_T2 + PLA_T2)/2 - (CRE_T1 + PLA_T1)/2,
    Cancer_vs_Healthy      = (CRE_T1 + PLA_T1)/2 - H_T1,
    Supplement_Interaction = (CRE_T2 - CRE_T1) - (PLA_T2 - PLA_T1),
    levels = design
  )

  cat("Step 3a: Array weights\n")
  aw <- arrayWeights(mat, design)
  names(aw) <- colnames(mat)
  cat(sprintf("  Array weights: %.2f - %.2f (median %.2f)\n",
              min(aw), max(aw), median(aw)))
  write_csv(tibble(sample_id = names(aw), array_weight = round(aw, 4)),
            file.path(data_dir, "array_weights.csv"))

  cat("Step 3b: Fit model\n")
  dupcor <- duplicateCorrelation(mat, design, block = meta$Subject_ID, weights = aw)
  cat(sprintf("  Within-subject correlation: %.3f\n", dupcor$consensus.correlation))

  fit <- lmFit(mat, design, block = meta$Subject_ID,
               correlation = dupcor$consensus.correlation, weights = aw) |>
    contrasts.fit(contrast_mat) |> eBayes(robust = TRUE)
  fit$count <- pep_count[rownames(fit$coefficients)]
  fit <- spectraCounteBayes(fit)

  cat("Step 4: Results\n")
  contrast_names <- colnames(contrast_mat)
  per_contrast_dir <- file.path(data_dir, "per_contrast_results")
  dir.create(per_contrast_dir, showWarnings = FALSE, recursive = TRUE)

  all_results <- list()
  for (i in seq_along(contrast_names)) {
    cname <- contrast_names[i]
    res <- outputResult(fit, coef_col = i)
    res$gene <- NULL; res$uniprot_id <- rownames(res)

    res <- res |>
      left_join(as.data.frame(ann) |> dplyr::select(uniprot_id, gene, protein, description),
                by = "uniprot_id") |>
      arrange(sca.P.Value)

    all_results[[cname]] <- res
    cat(sprintf("  %s: %d DEqMS / %d limma sig (FDR<0.05)\n", cname,
                sum(res$sca.adj.pval < 0.05, na.rm = TRUE),
                sum(res$adj.P.Val < 0.05, na.rm = TRUE)))
  }

  PVAL_COL     <- "sca.P.Value"
  ADJ_PVAL_COL <- "sca.adj.pval"

  saveRDS(fit, file.path(data_dir, "DEqMS_fit.rds"))
}

# --- 5. Enrichment — pi-score flags, per-contrast CSVs, summary ---
cat("Step 6: Per-contrast output\n")

suffix <- if (METHOD == "DEqMS") "_DEqMS" else ""

for (cname in contrast_names) {
  all_results[[cname]] <- add_pi_flags(all_results[[cname]], pval_col = PVAL_COL)
  write_csv(all_results[[cname]],
            file.path(per_contrast_dir, paste0(cname, suffix, ".csv")))
}

combined <- map_dfr(names(all_results), ~mutate(all_results[[.x]], contrast = .x, .before = 1))
write_csv(combined, file.path(data_dir, paste0("combined_results", suffix, ".csv")))

if (METHOD == "limma") {
  da_path <- file.path(data_dir, "DA_summary.csv")
  da <- read_csv(da_path, show_col_types = FALSE)
  pi_rows <- map_dfr(contrast_names, function(cname) {
    res <- all_results[[cname]]
    tibble(contrast = cname, type = c("up", "down", "nonsig"),
           sig.Pi = c(sum(res$sig_pi == 1), sum(res$sig_pi == -1), sum(res$sig_pi == 0)),
           sig.Nominal.01 = c(sum(res$sig_nominal_01 == 1), sum(res$sig_nominal_01 == -1),
                              sum(res$sig_nominal_01 == 0)))
  })
  da <- left_join(da, pi_rows, by = c("contrast", "type"))
  write_csv(da, da_path)

} else {
  da_summary <- tibble(
    contrast = contrast_names,
    n_tested     = sapply(all_results, nrow),
    n_sig_deqms  = sapply(all_results, \(r) sum(r$sca.adj.pval < 0.05, na.rm = TRUE)),
    n_up_deqms   = sapply(all_results, \(r) sum(r$sca.adj.pval < 0.05 & r$logFC > 0, na.rm = TRUE)),
    n_down_deqms = sapply(all_results, \(r) sum(r$sca.adj.pval < 0.05 & r$logFC < 0, na.rm = TRUE)),
    n_sig_limma  = sapply(all_results, \(r) sum(r$adj.P.Val < 0.05, na.rm = TRUE)),
    n_sig_pi     = sapply(all_results, \(r) sum(r$pi_pval < 0.05, na.rm = TRUE)),
    n_sig_nom01  = sapply(all_results, \(r) sum(r[[PVAL_COL]] < 0.01, na.rm = TRUE))
  )
  write_csv(da_summary, file.path(data_dir, paste0("DA_summary", suffix, ".csv")))
}

saveRDS(all_results, file.path(data_dir, paste0("results_list", suffix, ".rds")))

# --- 6. Excel workbook ---
cat("Step 7: Excel workbook\n")
wb <- createWorkbook()
walk(contrast_names, ~{
  addWorksheet(wb, .x)
  writeData(wb, .x, all_results[[.x]])
})

da_sheet <- if (METHOD == "limma") {
  read_csv(file.path(data_dir, "DA_summary.csv"), show_col_types = FALSE)
} else {
  da_summary
}
addWorksheet(wb, "DA_Summary"); writeData(wb, "DA_Summary", da_sheet)
saveWorkbook(wb, file.path(data_dir, paste0("results", suffix, ".xlsx")), overwrite = TRUE)

# --- 7. DEqMS diagnostics (conditional) ---
if (METHOD == "DEqMS") {
  cat("Step 8: DEqMS diagnostics\n")

  pdf(file.path(report_dir, "07a_DEqMS_variance_diagnostic.pdf"), width = 10, height = 8)
  par(mfrow = c(1, 2))
  VarianceScatterplot(fit, xlab = "log2(Peptide count)")
  abline(h = log(fit$s2.prior), col = "green", lwd = 2)
  legend("topright", legend = c("DEqMS prior", "Limma prior"),
         col = c("red", "green"), lwd = 2, cex = 0.8)
  VarianceBoxplot(fit, n = 20, xlab = "Peptide count", main = "Variance by peptide count")
  dev.off()

  for (cname in contrast_names) {
    res <- all_results[[cname]] |>
      mutate(neg_log10_p = -log10(sca.P.Value),
             sig = case_when(sca.adj.pval < 0.05 & logFC > 0 ~ "Up",
                             sca.adj.pval < 0.05 & logFC < 0 ~ "Down", TRUE ~ "NS"),
             label = if_else(sca.adj.pval < 0.05, gene, NA_character_))

    cdir <- file.path(report_dir, cname); dir.create(cdir, showWarnings = FALSE)
    n_up <- sum(res$sig == "Up"); n_down <- sum(res$sig == "Down")

    p_vol <- ggplot(res, aes(logFC, neg_log10_p, color = sig)) +
      geom_point(alpha = 0.5, size = 1.5) +
      geom_text_repel(aes(label = label), size = 2.2, max.overlaps = 20,
                      show.legend = FALSE, na.rm = TRUE) +
      scale_color_manual(values = c(Up = "#B2182B", Down = "#2166AC", NS = "gray70"),
                         labels = c(Up = sprintf("Up (%d)", n_up),
                                    Down = sprintf("Down (%d)", n_down), NS = "NS")) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.4) +
      labs(x = "log2FC", y = "-log10(DEqMS p)",
           title = sprintf("%s (DEqMS)", cname), color = NULL) +
      theme_minimal(base_size = 12) + theme(legend.position = "bottom")
    ggsave(file.path(cdir, "volcano_DEqMS.pdf"), p_vol, width = 8, height = 6)

    p_hist <- ggplot(res, aes(sca.P.Value)) +
      geom_histogram(bins = 50, fill = "steelblue", color = "white", boundary = 0) +
      geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", alpha = 0.6) +
      labs(x = "DEqMS p-value", y = "Count",
           title = sprintf("%s — P-value Distribution", cname)) +
      theme_minimal(base_size = 12)
    ggsave(file.path(cdir, "pval_histogram_DEqMS.pdf"), p_hist, width = 7, height = 5)
  }

  p_compare <- wrap_plots(map(contrast_names, function(cname) {
    res <- all_results[[cname]]
    ggplot(res, aes(t, sca.t, color = log2(count + 1))) +
      geom_point(alpha = 0.5, size = 1) +
      geom_abline(slope = 1, linetype = "dashed", color = "red") +
      scale_color_viridis_c(name = "log2(n_seq)") +
      labs(x = "limma t", y = "DEqMS t (sca.t)", title = cname) +
      theme_minimal(base_size = 10)
  }), ncol = 2) +
    plot_annotation(title = "DEqMS vs limma t-statistics",
                    subtitle = "Colored by log2(peptide count); dashed = identity")
  ggsave(file.path(report_dir, "07c_DEqMS_vs_limma_tstat.pdf"), p_compare,
         width = 12, height = 10)
}

# --- 8. Summary ---
cat("Step 8: Summary\n")
cat(sprintf("  Method: %s | Data: %s\n", METHOD, basename(data_file)))

if (METHOD == "limma") {
  print(dal$design$contrast_matrix)
  cat("\n"); print(read_csv(file.path(data_dir, "DA_summary.csv"), show_col_types = FALSE))
} else {
  print(contrast_mat)
  cat("\n"); print(da_summary)
}

sessionInfo()
