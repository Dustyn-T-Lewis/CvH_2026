#!/usr/bin/env Rscript
# CvH Normalization Reports — diagnostic plots + supplementary workbook
#
# Reads: c_data/00_report_intermediates.rds (from 01_run_normalization.R)
#
# Outputs:
#   b_reports/04_diagnostics.pdf    — 4-page custom QC report
#   c_data/05_normalization_supp.xlsx — supplementary workbook (3 sheets)

library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(patchwork)
library(openxlsx)

setwd(rprojroot::find_rstudio_root_file())
PAL_GT <- c(CRE_T1 = "#2166AC", CRE_T2 = "#67A9CF",
            PLA_T1 = "#D6604D", PLA_T2 = "#F4A582",
            H_T1   = "#4DAF4A")

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

int <- readRDS("01_normalization/c_data/00_report_intermediates.rds")
list2env(int, envir = environment())

shape_tp <- c(T1 = 16, T2 = 17)

# --- PAGE 1: Filtering & missingness -----------------------------------------

p_filter <- ggplot(filter_bar_data, aes(step, n, fill = status)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = c(Retained = "#2166AC", Removed = "#B2182B")) +
  labs(x = NULL, y = "Proteins", fill = NULL, title = "Protein retention") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 25, hjust = 1))

p_miss <- ggplot(miss_bar_data, aes(reorder(Col_ID, -n * (status == "Detected")),
                                     n, fill = status)) +
  geom_col(aes(alpha = is_outlier), width = 0.8) +
  scale_fill_manual(values = c(Detected = "#2166AC", Missing = "#D6604D")) +
  scale_alpha_manual(values = c("FALSE" = 1, "TRUE" = 0.4), guide = "none") +
  facet_grid(~ Group_Time, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = "Proteins", fill = NULL,
       title = "Per-sample detection (all samples, outliers faded)") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
                   strip.text = element_text(face = "bold"))

# --- PAGE 2: Outlier diagnostics ---------------------------------------------
col_group <- c(CR_CRE = "#2166AC", CR_PLA = "#D6604D", PPS = "#4DAF4A")

p_out_miss <- ggplot(outlier_diag, aes(pct_missing, delta_missing,
                                        color = Group, shape = Timepoint)) +
  geom_point(size = 3) +
  geom_vline(xintercept = miss_thresh, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_hline(yintercept = delta_thresh, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_text_repel(data = \(d) filter(d, miss_flag),
                  aes(label = Subject_ID), size = 2.5, show.legend = FALSE) +
  scale_color_manual(values = col_group) +
  scale_shape_manual(values = shape_tp) +
  labs(x = "Sample missingness (%)", y = "Delta missingness (|T2 - T1|)",
       title = "A: Missingness",
       subtitle = sprintf("IQR threshold: %.1f%% miss, %.1f%% delta | %d flagged",
                           miss_thresh, delta_thresh, sum(outlier_diag$miss_flag))) +
  theme_minimal()

pca_outlier_df <- pca_pre$scores |>
  left_join(outlier_diag |> select(Col_ID, pca_flag), by = "Col_ID")

p_out_pca <- ggplot(pca_outlier_df, aes(PC1, PC2, color = Group, shape = Timepoint)) +
  geom_point(size = 3.5, alpha = 0.85) +
  geom_text_repel(data = \(d) filter(d, pca_flag),
                  aes(label = Subject_ID), size = 2.5, show.legend = FALSE) +
  scale_color_manual(values = col_group) +
  scale_shape_manual(values = shape_tp) +
  labs(x = sprintf("PC1 (%.1f%%)", pca_pre$var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", pca_pre$var_exp[2]),
       title = "B: PCA Mahalanobis",
       subtitle = sprintf("chi-sq cutoff p < %.2f | %d flagged",
                           cfg$mahal_p, sum(pca_outlier_df$pca_flag))) +
  coord_fixed() +
  theme_minimal()

p_out_mad <- ggplot(outlier_diag, aes(reorder(Subject_ID, sample_median),
                                       sample_median, color = Group, shape = Timepoint)) +
  geom_point(size = 2.5) +
  geom_text_repel(data = \(d) filter(d, mad_flag),
                  aes(label = Subject_ID), size = 2.5, show.legend = FALSE) +
  geom_hline(yintercept = global_med) +
  geom_hline(yintercept = global_med + c(-1, 1) * cfg$mad_k * mad_val,
             linetype = "dashed", color = "red", alpha = 0.5) +
  scale_color_manual(values = col_group) +
  scale_shape_manual(values = shape_tp) +
  labs(x = "Sample", y = "Median log2 intensity",
       title = "C: MAD median intensity",
       subtitle = sprintf("%dx MAD band | %d flagged",
                           cfg$mad_k, sum(outlier_diag$mad_flag))) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4))

p_out_cor <- ggplot(outlier_diag, aes(reorder(Subject_ID, median_cor),
                                       median_cor, color = Group, shape = Timepoint)) +
  geom_point(size = 2.5) +
  geom_text_repel(data = \(d) filter(d, cor_flag),
                  aes(label = Subject_ID), size = 2.5, show.legend = FALSE) +
  geom_hline(yintercept = median(outlier_diag$median_cor)) +
  geom_hline(yintercept = median(outlier_diag$median_cor) - cfg$mad_k * mad(outlier_diag$median_cor),
             linetype = "dashed", color = "red", alpha = 0.5) +
  scale_color_manual(values = col_group) +
  scale_shape_manual(values = shape_tp) +
  labs(x = "Sample", y = "Median pairwise correlation",
       title = "D: Inter-sample correlation",
       subtitle = sprintf("%dx MAD band | %d flagged",
                           cfg$mad_k, sum(outlier_diag$cor_flag))) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 4))

# --- PAGE 3: Post-normalization PCA ------------------------------------------

p_pca_post <- ggplot(pca_post$scores, aes(PC1, PC2,
                                            color = Group_Time, shape = Timepoint)) +
  geom_point(size = 3.5, alpha = 0.85) +
  stat_ellipse(aes(group = Group_Time), type = "norm", level = 0.68, linewidth = 0.7) +
  scale_color_manual(values = PAL_GT) + scale_shape_manual(values = shape_tp) +
  labs(x = sprintf("PC1 (%.1f%%)", pca_post$var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", pca_post$var_exp[2]),
       title = "Post-normalization PCA") +
  theme_minimal() + theme(legend.position = "bottom")

# --- PAGE 4: Variability -----------------------------------------------------

p_cv <- ggplot(subj_var, aes(reorder(Subject_ID, iqr), iqr)) +
  geom_line(aes(group = Subject_ID), color = "gray60", linewidth = 0.4) +
  geom_point(aes(color = Group_Time, shape = Timepoint), size = 2.5) +
  scale_color_manual(values = PAL_GT) +
  scale_shape_manual(values = shape_tp) +
  labs(x = "Subject", y = "IQR (log2 intensity)",
       title = "Per-subject variability",
       subtitle = sprintf("%d subjects | lines connect T1-T2 pairs",
                           length(unique(subj_var$Subject_ID)))) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
                   legend.position = "bottom")

p_eta2 <- ggplot(data.frame(eta2 = eta2_vals[!is.na(eta2_vals)]), aes(eta2)) +
  geom_histogram(bins = 50, fill = "#2166AC", color = "white", alpha = 0.8) +
  geom_vline(xintercept = median(eta2_vals, na.rm = TRUE),
             linetype = "dashed", color = "red") +
  annotate("text", x = median(eta2_vals, na.rm = TRUE) + 0.02, y = Inf,
           vjust = 2, size = 3.5, color = "red",
           label = sprintf("median = %.2f", median(eta2_vals, na.rm = TRUE))) +
  labs(x = expression(eta^2 ~ "(between-group / total)"),
       y = "Proteins", title = "Variance partition by group",
       subtitle = sprintf("eta-sq = SS_between / SS_total | median = %.2f | higher = more group-driven",
                           median(eta2_vals, na.rm = TRUE))) +
  theme_minimal()

# --- ASSEMBLE PDF (4 pages) --------------------------------------------------

pdf(file.path(cfg$report_dir, "04_diagnostics.pdf"), width = 20, height = 10)

print(
  p_filter / p_miss +
    plot_layout(heights = c(1, 1.2)) +
    plot_annotation(
      title = "Protein Filtering & Detection",
      subtitle = sprintf("%d raw -> %d retained | %d samples",
                         n_raw, dal_nrow, dal_ncol),
      theme = theme(plot.title = element_text(size = 18, face = "bold"),
                    plot.subtitle = element_text(size = 14)))
)

shared_legend <- get_legend(
  p_out_miss + theme(legend.position = "bottom",
                     legend.justification = "center")
)

print(
  ((p_out_miss + theme(legend.position = "none")) |
   (p_out_pca  + theme(legend.position = "none"))) /
  ((p_out_mad  + theme(legend.position = "none")) |
   (p_out_cor  + theme(legend.position = "none"))) /
  wrap_elements(shared_legend) +
    plot_layout(heights = c(1, 1, 0.08)) +
    plot_annotation(
      title = "Outlier Diagnostics (4-method consensus, >=3/4)",
      subtitle = sprintf("Consensus rule: sample removed if flagged by >=3 methods | %d removed",
                          n_outliers),
      theme = theme(plot.title = element_text(size = 18, face = "bold"),
                    plot.subtitle = element_text(size = 13)))
)

print(
  p_pca_post +
    plot_annotation(
      title = "Post-Normalization QC",
      theme = theme(plot.title = element_text(size = 18, face = "bold")))
)

print(
  (p_cv | p_eta2) +
    plot_annotation(
      title = "Variability Summary",
      theme = theme(plot.title = element_text(size = 18, face = "bold")))
)

dev.off()

# --- SUPPLEMENTARY WORKBOOK ---
wb <- createWorkbook()

add_sheet(wb, "Pipeline_Summary", filter_log,
  title = "Protein Filtering Pipeline",
  notes = "step: filter stage | n_before/n_after: counts | pct_of_raw: cumulative retention")

add_sheet(wb, "Outlier_Diagnostics", outlier_diag,
  title = "Per-Sample Outlier Diagnostics",
  notes = "4-method consensus: miss_flag/mad_flag/pca_flag/cor_flag | consensus_outlier: >= 3/4 flags")

add_sheet(wb, "Filtered_Proteins", filtered_proteins,
  title = "Proteins Removed by Filtering",
  notes = "removal_step: HPA or Missingness | identifiers: uniprot_id, gene, description")

saveWorkbook(wb, file.path(cfg$data_dir, "05_normalization_supp.xlsx"), overwrite = TRUE)

cat(sprintf("Done: 04_diagnostics.pdf + 05_normalization_supp.xlsx -> %s/, %s/\n",
            cfg$report_dir, cfg$data_dir))
