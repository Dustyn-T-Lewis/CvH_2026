# Figure 2 CR — Panel D: logFC Density Histograms (Effect Size Distributions)
# Contrasts: Baseline_Supplement, Training_CRE, Training_PLA, Supplement_Interaction
# Outputs: pD (ggplot object), panel_D_logfc_density_MAIN.pdf/.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F02/a_script/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
})

PD_W <- 120; PD_H <- 200

RPT_DIR <- "04_Figures/F02/CR/b_reports"
DAT_DIR <- "04_Figures/F02/CR/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

dep_df <- read_csv("03_DEP/c_data/03_combined_results_CR.csv",
                   show_col_types = FALSE)

pdf_device <- get_pdf_device()

cr_contrasts <- c("Baseline_Supplement", "Training_CRE",
                  "Training_PLA", "Supplement_Interaction")

lfc_long <- dep_df |>
  dplyr::select(gene, starts_with("logFC_")) |>
  pivot_longer(starts_with("logFC_"), names_to = "contrast", values_to = "logFC") |>
  mutate(contrast = str_remove(contrast, "logFC_")) |>
  filter(!is.na(logFC), contrast %in% cr_contrasts)

lfc_long$contrast <- factor(lfc_long$contrast, levels = cr_contrasts)

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
    .groups = "drop"
  )

lfc_binwidth <- 4 / 50

lfc_stats$annotation <- sprintf(
  "Med.|logFC| = %.2f [%.2f, %.2f]\nn(>0.5) = %d",
  lfc_stats$med_abs_lfc, lfc_stats$ci_lo, lfc_stats$ci_hi, lfc_stats$n_above_05
)

# ── KS test: Training_PLA vs Training_CRE (supplement effect on training) ──
lfc_cre <- lfc_long$logFC[lfc_long$contrast == "Training_CRE"]
lfc_pla <- lfc_long$logFC[lfc_long$contrast == "Training_PLA"]
ks_res  <- ks.test(abs(lfc_cre), abs(lfc_pla))

# Magnitude ratio: median |logFC| PLA / CRE with bootstrap CI
obs_ratio <- median(abs(lfc_pla)) / median(abs(lfc_cre))
set.seed(42)
boot_ratios <- replicate(2000, {
  median(sample(abs(lfc_pla), replace = TRUE)) /
    median(sample(abs(lfc_cre), replace = TRUE))
})
ratio_ci <- quantile(boot_ratios, c(0.025, 0.975))

n_prot <- n_distinct(lfc_long$gene)
dist_subtitle <- sprintf(
  "%s proteins (non-imputed) | KS D = %.2f, %s | Ratio PLA/CRE = %.2f [%.2f, %.2f]",
  format(n_prot, big.mark = ","),
  ks_res$statistic, fmt_p(ks_res$p.value),
  obs_ratio, ratio_ci[1], ratio_ci[2]
)

write.csv(as.data.frame(lfc_stats),
          file.path(DAT_DIR, "audit_panel_D_logfc_stats.csv"), row.names = FALSE)

pD <- ggplot(lfc_long, aes(x = logFC, fill = contrast)) +
  geom_histogram(bins = 50, color = "black", linewidth = 0.2, alpha = 0.85) +
  geom_density(aes(y = after_stat(count) * lfc_binwidth),
               alpha = 0.15, linewidth = 0.5, color = "grey20") +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey50", linewidth = 0.4) +
  geom_text(data = lfc_stats,
            aes(x = -1, y = Inf, label = annotation),
            inherit.aes = FALSE, hjust = -0.01, vjust = 1.15,
            size = scale_text(BASE_COUNT, PD_W),
            color = "grey20", fontface = "bold", lineheight = 0.9) +
  facet_wrap(~ contrast, ncol = 1, scales = "fixed",
             labeller = labeller(contrast = CTR_SHORT)) +
  coord_cartesian(xlim = c(-1, 1)) +
  scale_fill_manual(values = CONTRAST_COLORS[cr_contrasts]) +
  labs(title = "Effect Size Distribution (CR Model)",
       subtitle = dist_subtitle,
       x = expression(bold(log[2]~FC)), y = NULL,
       tag = "D") +
  FIG_THEME + theme(legend.position = "none",
                    plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE - 1.5,
                                                face = "bold.italic", color = "grey40"),
                    strip.text = element_text(face = "bold", size = FIG_STRIP_SIZE,
                                             margin = margin(b = 2)),
                    panel.spacing.y = unit(4, "pt"),
                    axis.text.y = element_text(size = FIG_AXIS_TEXT - 1.5,
                                               color = "grey40"),
                    axis.ticks.y = element_blank())

ggsave(file.path(RPT_DIR, "panel_D_logfc_density_MAIN.pdf"), pD,
       width = PD_W, height = PD_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_D_logfc_density_MAIN.png"), pD,
       width = PD_W, height = PD_H, units = "mm", dpi = 300)
