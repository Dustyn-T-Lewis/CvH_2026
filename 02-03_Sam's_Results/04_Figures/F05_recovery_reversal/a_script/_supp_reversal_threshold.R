#!/usr/bin/env Rscript
# Sourced by 02_supp_panels.R — expects style.R already loaded.
# F05 Supplementary Panel D: Reversal Classification Threshold Sensitivity
# Defends pattern classification — shows reversal % is stable across thresholds.
# Reads pre-computed threshold_sensitivity.csv if available, otherwise computes.
# Line plot: x = |logFC| threshold, y = percentage, color = category.

suppressPackageStartupMessages({
  library(tidyverse)
})

BASE    <- "02-03_Sam's_Results/04_Figures/F05_recovery_reversal"
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
DAT     <- file.path(BASE, "c_data", "panel_supp")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT,     recursive = TRUE, showWarnings = FALSE)
pdf_device <- grDevices::pdf  # force base pdf (cairo DLL fails)

# Data
pre_computed <- file.path(BASE, "c_data", "reversal_tests", "threshold_sensitivity.csv")
if (file.exists(pre_computed)) {
  thresh_df <- read_csv(pre_computed, show_col_types = FALSE)
  message("  Read pre-computed threshold sensitivity")
} else {
  dep <- read_csv("02-03_Sam's_Results/03_DEP/c_data/F05_combined_CvHvTCR.csv",
                  show_col_types = FALSE)
  fc_df <- dep |>
    transmute(gene,
              logFC_CvH = logFC_Cancer_vs_Healthy,
              logFC_TCR = logFC_Training_CR) |>
    filter(!is.na(logFC_CvH), !is.na(logFC_TCR))

  thresholds <- c(0.05, 0.1, 0.2, 0.3)
  res <- lapply(thresholds, function(thr) {
    # Threshold applied symmetrically to both axes: cancer effect must exceed
    # thr to be classifiable, and training response must also exceed thr.
    # This avoids over-counting reversal for near-zero cancer effects.
    classified <- fc_df |>
      mutate(category = case_when(
        abs(logFC_CvH) < thr               ~ "Negligible",
        (logFC_CvH > 0 & logFC_TCR < -thr) |
          (logFC_CvH < 0 & logFC_TCR > thr)  ~ "Reversed",
        (logFC_CvH > 0 & logFC_TCR > thr) |
          (logFC_CvH < 0 & logFC_TCR < -thr) ~ "Exacerbated",
        TRUE                                   ~ "Negligible"
      ))
    n_total <- nrow(classified)
    classified |>
      count(category) |>
      mutate(threshold = thr,
             pct       = 100 * n / n_total)
  })
  thresh_df <- bind_rows(res)
  message("  Computed threshold sensitivity from DEP results")
}

write_csv(thresh_df, file.path(DAT, "SUPP_reversal_threshold.csv"))

# Plot
cat_colors <- c("Reversed" = "#4393C3", "Exacerbated" = "#D6604D", "Negligible" = "grey60")
cat_levels <- c("Reversed", "Exacerbated", "Negligible")
thresh_df <- thresh_df |>
  mutate(category = factor(category, levels = cat_levels))

pS_thresh <- ggplot(thresh_df, aes(x = threshold, y = pct, color = category)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.5, shape = 16) +
  scale_color_manual(values = cat_colors, name = "Category") +
  scale_x_continuous(breaks = unique(thresh_df$threshold),
                     labels = sprintf("%.2f", unique(thresh_df$threshold))) +
  labs(title    = "Reversal Classification Threshold Sensitivity",
       subtitle = sprintf("logFC thresholds: %s | stable reversal %%",
                          paste(unique(thresh_df$threshold), collapse = ", ")),
       x = "|logFC| threshold",
       y = "Percentage of proteins") +
  FIG_THEME +
  theme(legend.position = "right")

PW <- 89; PH <- 70
ggsave(file.path(RPT_PNG, "SUPP_reversal_threshold.png"), pS_thresh,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_reversal_threshold.pdf"), pS_thresh,
       width = PW, height = PH, units = "mm", device = pdf_device)

message("F05 SUPP Panel D (reversal threshold sensitivity) saved")

# Expose for composite
pS_thresh_title    <- "Reversal Classification Threshold Sensitivity"
pS_thresh_subtitle <- sprintf("logFC thresholds: %s | stable reversal %%",
                               paste(unique(thresh_df$threshold), collapse = ", "))
pS_thresh          <- strip_for_composite(pS_thresh)
