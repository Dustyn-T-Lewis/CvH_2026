# Figure 2 (CRvH) — Panel D: logFC Density Histograms (Effect Size Distributions)
# Contrasts: Cancer_vs_Healthy, Training_CR from CRvH model.
# Bootstrap median |logFC| CIs.
# Outputs: pD (ggplot object), panel_D_logfc_density_MAIN.pdf/.png

setwd(here::here())
source("04_Figures/F02/a_script/style.R")

pacman::p_load(dplyr, tidyr, stringr, readr, ggplot2)

PD_W <- 120; PD_H <- 120

RPT_DIR <- "04_Figures/F02/CRvH/b_reports"
DAT_DIR <- "04_Figures/F02/CRvH/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

# DEP results: new proteoDA long output -> wide per-contrast columns.
dep_df <- read_csv("03_DEP/a_non_imputed/c_data/combined_results_pi.csv",
                   show_col_types = FALSE) |>
  mutate(contrast = recode(contrast,
                           CRvH_Baseline = "Cancer_vs_Healthy",
                           CR_Training   = "Training_CR")) |>
  pivot_wider(id_cols = c(uniprot_id, gene, protein, description),
              names_from = contrast,
              values_from = c(logFC, t, P.Value, adj.P.Val, pi_score, sig_pi),
              names_glue = "{.value}_{contrast}")

pdf_device <- get_pdf_device()

# ── Reshape to long for CRvH contrasts ──
lfc_long <- dep_df |>
  dplyr::select(gene, starts_with("logFC_")) |>
  pivot_longer(starts_with("logFC_"), names_to = "contrast", values_to = "logFC") |>
  mutate(contrast = str_remove(contrast, "logFC_")) |>
  filter(!is.na(logFC), contrast %in% c("Cancer_vs_Healthy", "Training_CR"))

lfc_long$contrast <- factor(lfc_long$contrast,
                            levels = c("Cancer_vs_Healthy", "Training_CR"))

# ── Bootstrap 95% CI on median |logFC| per contrast ──
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

# ── KS test between contrasts ──
lfc_cvh <- lfc_long$logFC[lfc_long$contrast == "Cancer_vs_Healthy"]
lfc_tr  <- lfc_long$logFC[lfc_long$contrast == "Training_CR"]
ks_res  <- ks.test(abs(lfc_cvh), abs(lfc_tr))

# Magnitude ratio: median |logFC| Training_CR / Cancer_vs_Healthy with bootstrap CI
obs_ratio <- median(abs(lfc_tr)) / median(abs(lfc_cvh))
set.seed(42)
boot_ratios <- replicate(2000, {
  median(sample(abs(lfc_tr), replace = TRUE)) / median(sample(abs(lfc_cvh), replace = TRUE))
})
ratio_ci <- quantile(boot_ratios, c(0.025, 0.975))

dist_subtitle <- sprintf(
  "KS D = %.2f, %s | Ratio Tr/CvH = %.2f [%.2f, %.2f]",
  ks_res$statistic, fmt_p(ks_res$p.value),
  obs_ratio, ratio_ci[1], ratio_ci[2]
)

# ── Plot ──
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
  scale_fill_manual(values = CONTRAST_COLORS[c("Cancer_vs_Healthy", "Training_CR")]) +
  labs(title = "Effect Size Distribution",
       subtitle = dist_subtitle,
       x = expression(bold(log[2]~FC)), y = NULL,
       tag = "D") +
  FIG_THEME + theme(legend.position = "none",
                    strip.text = element_text(face = "bold", size = FIG_STRIP_SIZE,
                                              margin = margin(b = 2)),
                    panel.spacing.y = unit(4, "pt"),
                    plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE - 1.5,
                                                face = "bold.italic", color = "grey40"),
                    axis.text.y = element_text(size = FIG_AXIS_TEXT - 1.5, color = "grey40"),
                    axis.ticks.y = element_blank())

# ── Save ──
write.csv(as.data.frame(lfc_stats),
          file.path(DAT_DIR, "audit_panel_D_logfc_stats.csv"), row.names = FALSE)

ggsave(file.path(RPT_DIR, "panel_D_logfc_density_MAIN.pdf"), pD,
       width = PD_W, height = PD_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_D_logfc_density_MAIN.png"), pD,
       width = PD_W, height = PD_H, units = "mm", dpi = 300)
