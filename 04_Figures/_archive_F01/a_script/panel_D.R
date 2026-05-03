# Figure 1 — Panel D: logFC Density Histograms (5 CvH contrasts)
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(readr); library(ggplot2)
})

PD_W <- 140; PD_H <- 240
RPT_DIR <- "04_Figures/F01/b_reports"
DAT_DIR <- "04_Figures/F01/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# ── Read two DEP result files and merge logFC columns ──
dep_crvh <- read_csv("03_DEP/c_data/03_combined_results_CRvH.csv", show_col_types = FALSE) |>
  dplyr::select(uniprot_id, gene, starts_with("logFC_"))

dep_cr <- read_csv("03_DEP/c_data/03_combined_results_CR.csv", show_col_types = FALSE) |>
  dplyr::select(uniprot_id, starts_with("logFC_"))

dep_df <- inner_join(dep_crvh, dep_cr, by = "uniprot_id")

# ── Pivot to long format ──
contrast_order <- c("Cancer_vs_Healthy", "Training_CR",
                    "Baseline_Supplement", "Training_CRE",
                    "Training_PLA", "Supplement_Interaction")

lfc_long <- dep_df |>
  pivot_longer(starts_with("logFC_"), names_to = "contrast", values_to = "logFC") |>
  mutate(contrast = str_remove(contrast, "logFC_")) |>
  filter(!is.na(logFC), contrast %in% contrast_order) |>
  mutate(contrast = factor(contrast, levels = contrast_order))

# ── Bootstrap median |logFC| CI ──
set.seed(42)
boot_median_ci <- function(x, R = 2000, conf = 0.95) {
  meds <- replicate(R, median(sample(x, replace = TRUE)))
  qs   <- quantile(meds, c((1 - conf) / 2, (1 + conf) / 2))
  c(lower = unname(qs[1]), upper = unname(qs[2]))
}

lfc_stats <- lfc_long |>
  group_by(contrast) |>
  summarise(
    med_abs_lfc = median(abs(logFC)),
    ci_lo       = boot_median_ci(abs(logFC))[["lower"]],
    ci_hi       = boot_median_ci(abs(logFC))[["upper"]],
    n_above_05  = sum(abs(logFC) > 0.5),
    .groups     = "drop")

lfc_binwidth <- 4 / 50

lfc_stats$annotation <- sprintf(
  "Med.|logFC| = %.2f [%.2f, %.2f]\nn(>0.5) = %d",
  lfc_stats$med_abs_lfc, lfc_stats$ci_lo, lfc_stats$ci_hi, lfc_stats$n_above_05)

# ── Plot ──
pD <- ggplot(lfc_long, aes(x = logFC, fill = contrast)) +
  geom_histogram(bins = 50, color = "black", linewidth = 0.2, alpha = 0.85) +
  geom_density(aes(y = after_stat(count) * lfc_binwidth),
               alpha = 0.15, linewidth = 0.5, color = "grey20") +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.4) +
  geom_text(data = lfc_stats,
            aes(x = -0.90, y = Inf, label = annotation),
            inherit.aes = FALSE, hjust = -0.01, vjust = 1.15,
            size = scale_text(BASE_COUNT, PD_W),
            color = "grey20", fontface = "bold", lineheight = 0.9) +
  facet_wrap(~ contrast, ncol = 1, scales = "fixed",
             labeller = labeller(contrast = CTR_SHORT)) +
  coord_cartesian(xlim = c(-1, 1)) +
  scale_fill_manual(values = CONTRAST_COLORS) +
  labs(title = "Effect Size Distribution",
       x = expression(bold(log[2]~FC)), y = NULL, tag = "D") +
  FIG_THEME + theme(legend.position = "none",
                    strip.text = element_text(face = "bold", size = FIG_STRIP_SIZE,
                                              margin = margin(b = 2)),
                    panel.spacing.y = unit(4, "pt"),
                    axis.text.y = element_text(size = FIG_AXIS_TEXT - 1.5, color = "grey40"),
                    axis.ticks.y = element_blank())

# ── Save ──
ggsave(file.path(RPT_DIR, "panel_D_logfc_density.pdf"), pD,
       width = PD_W, height = PD_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_D_logfc_density.png"), pD,
       width = PD_W, height = PD_H, units = "mm", dpi = 300)

cat("Panel D saved:", PD_W, "x", PD_H, "mm\n")
cat("  Contrasts:", nlevels(lfc_long$contrast), "\n")
print(lfc_stats |> dplyr::select(contrast, med_abs_lfc, n_above_05))
