#!/usr/bin/env Rscript
# 02_imputation_reports.R — Missingness and imputation quality reports (CvH)
#
# Reads:
#   c_data/00_report_intermediates.rds
#   c_data/benchmark/04_composite_ranking.csv  (preferred, composite ranking)
#   c_data/benchmark/03_benchmark_summary.csv  (fallback, old summary)
#
# Writes:
#   b_reports/01_missingness_report.pdf  (1 page, 16x10 in)
#   b_reports/02_imputation_report.pdf   (2 pages, 14x10 in)

library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(readr)
library(scales)
library(ggrepel)

setwd(rprojroot::find_rstudio_root_file())

# --- Load data ---------------------------------------------------------------
rpt <- readRDS("02_Imputation/c_data/00_report_intermediates.rds")

mat         <- rpt$mat
mat_imp     <- rpt$mat_imp
was_na      <- rpt$was_na
meta        <- rpt$meta
mc          <- rpt$miss_class
audit       <- rpt$mnar_audit
PAL_GT      <- rpt$PAL_GT
PAL_CLASS   <- rpt$PAL_CLASS
pct_miss    <- rpt$pct_miss
n_mar_prots <- rpt$n_mar_prots
n_mnar_prots <- rpt$n_mnar_prots
mar_miss_vals  <- rpt$mar_miss_vals
mnar_miss_vals <- rpt$mnar_miss_vals
total_miss_vals <- rpt$total_miss_vals
best        <- rpt$best
prot_pct    <- rpt$prot_pct
mnar_ids    <- rpt$mnar_ids
miss_by_group <- rpt$miss_by_group
classification_method <- rpt$classification_method

n_prot       <- nrow(mat)
n_samp       <- ncol(mat)
n_complete   <- sum(mc$classification == "Complete")
n_unreliable <- sum(!audit$imputation_reliable)

dir.create("02_Imputation/b_reports", showWarnings = FALSE, recursive = TRUE)

# --- Theme -------------------------------------------------------------------
THM <- theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        panel.grid.minor = element_blank())

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

# =============================================================================
# Report 1: Missingness (1 page)
# =============================================================================

# Panel A — Per-protein missingness histogram
pA1 <- mc |>
  filter(classification != "Complete") |>
  ggplot(aes(x = pct_miss, fill = classification)) +
  geom_histogram(binwidth = 5, boundary = 0, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = PAL_CLASS, name = NULL) +
  scale_x_continuous(labels = label_percent(scale = 1)) +
  labs(title = "A. Per-protein missingness", x = "% missing", y = "Count") +
  THM + theme(legend.position = "top")

# Panel B — Classification summary stacked bar
class_counts <- mc |>
  count(classification) |>
  mutate(classification = factor(classification, levels = c("Complete", "MAR", "MNAR")))

pB1 <- class_counts |>
  ggplot(aes(x = 1, y = n, fill = classification)) +
  geom_col(width = 0.6, color = "white") +
  geom_text(aes(label = paste0(classification, "\n", n)),
            position = position_stack(vjust = 0.5), size = 3.5, fontface = "bold") +
  scale_fill_manual(values = PAL_CLASS) +
  coord_flip() +
  labs(title = "B. Classification summary", x = NULL, y = "Proteins") +
  THM + theme(legend.position = "none", axis.text.y = element_blank())

# Panel C — Scatter: mean intensity vs % missing
pC1 <- mc |>
  filter(classification != "Complete") |>
  ggplot(aes(x = mean_intensity, y = pct_miss, color = classification)) +
  geom_point(alpha = 0.5, size = 1.2) +
  scale_color_manual(values = PAL_CLASS, name = NULL) +
  scale_y_continuous(labels = label_percent(scale = 1)) +
  labs(title = "C. Intensity vs missingness",
       x = "Mean log2 intensity", y = "% missing") +
  THM + theme(legend.position = "top")

# Panel D — Per-sample missingness bar
sample_miss <- data.frame(
  sample = colnames(mat),
  n_miss = colSums(is.na(mat))
) |>
  left_join(meta, by = c("sample" = "Col_ID")) |>
  mutate(Group_Time = factor(Group_Time, levels = names(PAL_GT)))

pD1 <- sample_miss |>
  ggplot(aes(x = reorder(sample, n_miss), y = n_miss, fill = Group_Time)) +
  geom_col() +
  scale_fill_manual(values = PAL_GT, name = NULL) +
  coord_flip() +
  labs(title = "D. Per-sample missingness", x = NULL, y = "Missing proteins") +
  THM + theme(legend.position = "top",
              axis.text.y = element_text(size = 6))

class_method <- classification_method %||% "k-means"
miss_sub <- sprintf("%s proteins \u00d7 %s samples | %.1f%% missing | %s classification",
                    comma(n_prot), n_samp, pct_miss, class_method)

fig1 <- (pA1 | pB1) / (pC1 | pD1) +
  plot_annotation(
    title = "Missingness Report",
    subtitle = miss_sub,
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

pdf("02_Imputation/b_reports/01_missingness_report.pdf", width = 16, height = 10)
print(fig1)
dev.off()
cat("Wrote: 02_Imputation/b_reports/01_missingness_report.pdf\n")

# =============================================================================
# Report 2: Imputation quality (2 pages)
# =============================================================================

# --- Page 1: Benchmark (if available) or OOB ---------------------------------

# Read OOB from summary file
oob_line <- grep("oob_error", readLines("02_Imputation/c_data/09_imputation_summary.txt"),
                 value = TRUE)
oob_val <- as.numeric(sub(".*=\\s*(?:c\\(NRMSE = )?([0-9.]+)\\)?", "\\1", oob_line, perl = TRUE))

# Try to load benchmark ranking — prefer new composite, fall back to old summary
composite_file <- "02_Imputation/c_data/benchmark/04_composite_ranking.csv"
bench_file <- "02_Imputation/c_data/benchmark/03_benchmark_summary.csv"
has_composite <- file.exists(composite_file)
has_bench <- has_composite || file.exists(bench_file)

if (has_composite) {
  comp_df <- read_csv(composite_file, show_col_types = FALSE)

  bench_classes <- sort(unique(comp_df$class))
  PAL_BENCH <- setNames(
    c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628",
      "#F781BF", "#999999")[seq_along(bench_classes)],
    bench_classes
  )

  top10 <- comp_df |> slice_min(rank, n = 10)

  pA2 <- top10 |>
    mutate(method = factor(method, levels = rev(method)),
           is_mf = method == "missForest") |>
    ggplot(aes(x = composite, y = method, fill = class, alpha = is_mf)) +
    geom_col() +
    geom_text(aes(label = sprintf("%.3f", composite)), hjust = -0.1, size = 3) +
    scale_fill_manual(values = PAL_BENCH, name = "Method class") +
    scale_alpha_manual(values = c(`TRUE` = 1.0, `FALSE` = 0.6), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(title = "A. Composite ranking (top 10)", x = "Composite score", y = NULL) +
    THM

  pB2 <- comp_df |>
    mutate(is_mf = method == "missForest") |>
    ggplot(aes(x = nrmse_mcar, y = fc_rho, color = class)) +
    geom_point(aes(size = is_mf), alpha = 0.7) +
    geom_text_repel(aes(label = method), size = 2.8, max.overlaps = 15) +
    scale_color_manual(values = PAL_BENCH, name = "Method class") +
    scale_size_manual(values = c(`TRUE` = 3.5, `FALSE` = 2), guide = "none") +
    labs(title = "B. Reconstruction vs fold-change fidelity",
         x = "NRMSE (MCAR)", y = expression("FC correlation ("*rho*")")) +
    THM

  bm_sub <- sprintf("#1 %s %.3f | #2 %s %.3f | OOB = %.3f",
                     comp_df$method[1], comp_df$composite[1],
                     comp_df$method[2], comp_df$composite[2], oob_val)

  page1 <- (pA2 | pB2) +
    plot_annotation(
      title = "Imputation Benchmark",
      subtitle = bm_sub,
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )

} else if (file.exists(bench_file)) {
  bench_sum <- read_csv(bench_file, show_col_types = FALSE)

  bench_plot_df <- bench_sum |>
    mutate(is_mf = method == "missForest" | method == "RF",
           method = factor(method, levels = rev(method)))

  pA2 <- bench_plot_df |>
    ggplot(aes(x = mean_nrmse, y = method, alpha = is_mf)) +
    geom_point(size = 3.5, color = "#377EB8") +
    geom_errorbarh(aes(xmin = mean_nrmse - sd_nrmse,
                       xmax = mean_nrmse + sd_nrmse),
                   height = 0.3, linewidth = 0.5, color = "#377EB8") +
    scale_alpha_manual(values = c(`TRUE` = 1.0, `FALSE` = 0.4), guide = "none") +
    labs(title = "A. Benchmark NRMSE ranking",
         subtitle = sprintf("Best: %s (%.4f) | missForest OOB: %.4f",
                            bench_sum$method[1], bench_sum$mean_nrmse[1], oob_val),
         x = "NRMSE (lower = better)", y = NULL) +
    THM

  page1 <- pA2 +
    plot_annotation(
      title = "Imputation Benchmark",
      subtitle = sprintf("12 methods | missForest | OOB = %.4f", oob_val),
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )

} else {
  # No benchmark file: show OOB error only
  oob_df <- data.frame(metric = "OOB NRMSE", value = oob_val)
  page1 <- ggplot(oob_df, aes(x = metric, y = value)) +
    geom_col(fill = "#377EB8", width = 0.4) +
    geom_text(aes(label = sprintf("%.4f", value)), vjust = -0.5, size = 5) +
    ylim(0, max(oob_val * 1.3, 0.5)) +
    labs(title = "missForest OOB Error",
         subtitle = "No benchmark data available",
         x = NULL, y = "NRMSE") +
    THM +
    plot_annotation(
      title = "Imputation Quality",
      subtitle = sprintf("missForest | OOB = %.4f", oob_val),
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )
}

# --- Page 2: Observed vs imputed + MNAR audit --------------------------------

# Panel A — Density: observed vs imputed
obs_vals <- as.numeric(mat[!was_na])
imp_vals <- as.numeric(mat_imp[was_na])

dens_df <- bind_rows(
  data.frame(value = obs_vals, type = "Observed"),
  data.frame(value = imp_vals, type = "Imputed")
)

pA3 <- dens_df |>
  ggplot(aes(x = value, fill = type, color = type)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = c(Observed = "#377EB8", Imputed = "#E41A1C"), name = NULL) +
  scale_color_manual(values = c(Observed = "#377EB8", Imputed = "#E41A1C"), name = NULL) +
  labs(title = "A. Observed vs imputed value distributions",
       x = "log2 intensity", y = "Density") +
  THM + theme(legend.position = "top")

# Panel B — MNAR audit scatter: pre vs post mean
pB3 <- audit |>
  ggplot(aes(x = pre_mean, y = post_mean, color = imputation_reliable)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c(`TRUE` = "#4DAF4A", `FALSE` = "#E41A1C"),
                     labels = c(`TRUE` = "Reliable", `FALSE` = "Unreliable"),
                     name = NULL) +
  labs(title = "B. MNAR imputation audit",
       x = "Pre-imputation mean", y = "Post-imputation mean") +
  THM + theme(legend.position = "top")

# Panel C — Cohen's d histogram for MNAR proteins
d_median <- median(audit$effect_d, na.rm = TRUE)
d_iqr    <- IQR(audit$effect_d, na.rm = TRUE)

pC3 <- audit |>
  ggplot(aes(x = effect_d)) +
  geom_histogram(binwidth = 0.1, fill = "#984EA3", color = "white", alpha = 0.7) +
  geom_vline(xintercept = d_median, linetype = "dashed", color = "black") +
  annotate("text", x = d_median, y = Inf, vjust = 2, hjust = -0.1, size = 3.2,
           label = sprintf("median = %.2f\nIQR = %.2f", d_median, d_iqr)) +
  labs(title = "C. Imputation shift (MNAR proteins)",
       x = "Cohen's d (pre vs post)", y = "Count") +
  THM

imp_sub <- sprintf("missForest | %s values imputed | %d unreliable (>50%% missing)",
                   comma(total_miss_vals), n_unreliable)

page2 <- pA3 / (pB3 | pC3) +
  plot_annotation(
    title = "Imputation Quality",
    subtitle = imp_sub,
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

pdf("02_Imputation/b_reports/02_imputation_report.pdf", width = 14, height = 10)
print(page1)
print(page2)
dev.off()
cat("Wrote: 02_Imputation/b_reports/02_imputation_report.pdf\n")
