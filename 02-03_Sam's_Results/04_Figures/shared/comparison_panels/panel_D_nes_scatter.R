# Shared NES Scatter Panel — Panel D in F04 (concordance) and F05 (reversal)

source("02-03_Sam's_Results/04_Figures/shared/style.R")
source("02-03_Sam's_Results/04_Figures/shared/print_scale_apply_380mm.R")

library(tidyverse)
library(ggrepel)

PG_W <- cfg$panel_w %||% 146
RPT_PNG <- cfg$rpt_png
RPT_PDF <- cfg$rpt_pdf
DAT     <- cfg$dat
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_D"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- grDevices::pdf  # force base pdf (cairo DLL fails)

fgsea_cache <- "02-03_Sam's_Results/04_Figures/shared/fgsea_tstat_all_sam.csv"
stopifnot("fGSEA cache missing" = file.exists(fgsea_cache))
fgsea_all <- read_csv(fgsea_cache, show_col_types = FALSE)
# leadingEdge stored as semicolon-delimited string — parse to list if needed
if (is.character(fgsea_all$leadingEdge)) {
  fgsea_all$leadingEdge <- strsplit(fgsea_all$leadingEdge, ";")
}

cx <- cfg$contrast_x
cy <- cfg$contrast_y
nes_x <- paste0("NES_", cx)
nes_y <- paste0("NES_", cy)
padj_x <- paste0("padj_", cx)
padj_y <- paste0("padj_", cy)
size_x <- paste0("size_", cx)
size_y <- paste0("size_", cy)

fgsea_hg <- fgsea_all |>
  filter(database %in% c("Hallmark", "GO Slim"),
         contrast %in% c(cx, cy))

fgsea_wide <- fgsea_hg |>
  dplyr::select(pathway, contrast, NES, padj, size, database) |>
  pivot_wider(id_cols = c(pathway, database), names_from = contrast,
              values_from = c(NES, padj, size)) |>
  filter(!is.na(.data[[nes_x]]), !is.na(.data[[nes_y]])) |>
  mutate(set_size = coalesce(.data[[size_x]], .data[[size_y]]))

fgsea_wide <- fgsea_wide |>
  mutate(
    sig_1 = !is.na(.data[[padj_x]]) & .data[[padj_x]] < 0.05,
    sig_2 = !is.na(.data[[padj_y]]) & .data[[padj_y]] < 0.05,
    significance = case_when(
      sig_1 & sig_2 ~ cfg$quadrant_defs$sig_both_label,
      sig_1         ~ cfg$quadrant_defs$sig_x_label,
      sig_2         ~ cfg$quadrant_defs$sig_y_label,
      TRUE          ~ "NS"
    ) |> factor(levels = names(cfg$sig_colors)),
    pathway_label = clean_pathway_name(pathway),
    db_shape = ifelse(database == "Hallmark", 24, 21)
  )

fgsea_sig <- fgsea_wide |> filter(significance != "NS")

message(sprintf("  %d total pathways (Hallmark: %d, GO Slim: %d) | %d significant",
                nrow(fgsea_wide),
                sum(fgsea_wide$database == "Hallmark"),
                sum(fgsea_wide$database == "GO Slim"),
                nrow(fgsea_sig)))

nes_cor_all <- cor.test(fgsea_wide[[nes_x]], fgsea_wide[[nes_y]], method = "spearman")
nes_ci_all  <- fisher_z_ci(nes_cor_all$estimate, nrow(fgsea_wide))
nes_cor_sig <- if (nrow(fgsea_sig) >= 3) {
  cor.test(fgsea_sig[[nes_x]], fgsea_sig[[nes_y]], method = "spearman")
} else NULL

nes_lim <- max(abs(c(fgsea_wide[[nes_x]], fgsea_wide[[nes_y]]))) * 1.15

qd <- cfg$quadrant_defs
n_q1 <- sum(fgsea_sig[[nes_x]] > 0 & fgsea_sig[[nes_y]] > 0)
n_q2 <- sum(fgsea_sig[[nes_x]] < 0 & fgsea_sig[[nes_y]] > 0)
n_q3 <- sum(fgsea_sig[[nes_x]] < 0 & fgsea_sig[[nes_y]] < 0)
n_q4 <- sum(fgsea_sig[[nes_x]] > 0 & fgsea_sig[[nes_y]] < 0)
n_metric <- qd$metric_count_fn(n_q1, n_q2, n_q3, n_q4)
n_total_sig <- nrow(fgsea_sig)
metric_frac <- if (n_total_sig > 0) n_metric / n_total_sig else 0

message(sprintf("  NES Spearman (all): rho = %.3f [%.3f, %.3f]",
                nes_cor_all$estimate, nes_ci_all[1], nes_ci_all[2]))

txt_pw   <- scale_text(BASE_PATHWAY, PG_W)          # no multiplier — prevents label overlap in composite
txt_quad <- scale_text(BASE_QUADRANT, PG_W) * 1.15

label_pw <- fgsea_sig |>
  mutate(nes_mag = abs(.data[[nes_x]]) + abs(.data[[nes_y]])) |>
  arrange(desc(nes_mag)) |>
  slice_head(n = 12) |>
  select(-nes_mag) |>
  mutate(
    label_fill     = cfg$sig_label_fill[as.character(significance)],
    label_text_col = cfg$sig_label_text[as.character(significance)],
    pathway_label  = pathway_label |>
      str_replace("Amino Acid Metabolic.*", "Amino Acid Metabolism") |>
      str_replace("Muscle System.*", "Muscle System") |>
      str_replace("Ketone Metabolic.*", "Ketone Metabolism") |>
      str_replace("^Trna Metabolic.*", "tRNA Metabolism") |>
      str_replace("^Establishment Or Maintenance Of Cell Polarity$", "Cell Polarity") |>
      str_replace("^Generation Of Precursor Metabolites And Energy$",
                  "Precursor Metabolites & Energy") |>
      str_replace("^Protein Localization To Plasma Membrane$",
                  "Plasma Membrane Protein Loc.") |>
      str_replace("^Extracellular Matrix Organization$", "ECM Organization") |>
      str_replace("^Epithelial Mesenchymal Transition$", "EMT") |>
      str_replace("^Microtubule-Based Movement$", "Microtubule Movement") |>
      str_replace("^Mitochondrion Organization$", "Mitochondrial Organization") |>
      str_replace("^UV Response Dn$", "UV Response (Down)") |>
      dplyr::recode(!!!cfg$display_overrides)
  )

ns_df  <- fgsea_wide |> filter(significance == "NS")
sig_df <- fgsea_wide |> filter(significance != "NS") |>
  mutate(draw_order = factor(significance, levels = cfg$sig_draw_order)) |>
  arrange(draw_order)

rho_sig_str <- if (!is.null(nes_cor_sig)) sprintf(", \u03c1(sig) = %.2f", nes_cor_sig$estimate) else ""
subtitle_str <- sprintf(
  "GO Slim + Hallmark | %d pathways (%d sig.) | \u03c1 = %.2f [%.2f, %.2f], %s%s\n%.0f%% %s | %s",
  nrow(fgsea_wide), n_total_sig,
  nes_cor_all$estimate, nes_ci_all[1], nes_ci_all[2],
  ifelse(nes_cor_all$p.value < 0.001, "p < 0.001", sprintf("p = %.3f", nes_cor_all$p.value)),
  rho_sig_str, metric_frac * 100,
  cfg$subtitle_metric, cfg$subtitle_interpretation
)

pD <- ggplot(mapping = aes(x = .data[[nes_x]], y = .data[[nes_y]])) +
  annotate("rect", xmin = qd$bg_red_1[1], xmax = qd$bg_red_1[2],
           ymin = qd$bg_red_1[3], ymax = qd$bg_red_1[4],
           fill = AGE_COLORS["Old"], alpha = 0.20, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = qd$bg_red_2[1], xmax = qd$bg_red_2[2],
           ymin = qd$bg_red_2[3], ymax = qd$bg_red_2[4],
           fill = AGE_COLORS["Old"], alpha = 0.20, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = qd$bg_blue_1[1], xmax = qd$bg_blue_1[2],
           ymin = qd$bg_blue_1[3], ymax = qd$bg_blue_1[4],
           fill = AGE_COLORS["Young"], alpha = 0.20, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = qd$bg_blue_2[1], xmax = qd$bg_blue_2[2],
           ymin = qd$bg_blue_2[3], ymax = qd$bg_blue_2[4],
           fill = AGE_COLORS["Young"], alpha = 0.20, color = "grey70", linewidth = 0.2) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = cfg$ref_slope, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  geom_point(data = ns_df, aes(shape = database),
             size = 1.0, fill = "grey70", color = "grey55", alpha = 0.40, stroke = 0.2) +
  geom_point(data = sig_df, aes(fill = significance, size = set_size, shape = database),
             color = ifelse(sig_df$database == "Hallmark", "black", "grey65"),
             alpha = 0.80, stroke = 0.4) +
  scale_fill_manual(values = cfg$sig_colors, name = "Significance") +
  scale_shape_manual(values = c("Hallmark" = 24, "GO Slim" = 21), name = "Database") +
  scale_size_continuous(range = c(1.5, 5), name = "Set size",
                        breaks = c(20, 50, 100, 200)) +
  geom_label_repel(data = label_pw, aes(label = pathway_label),
                   fill = label_pw$label_fill, color = label_pw$label_text_col,
                   size = txt_pw, fontface = "bold",
                   max.overlaps = 50,
                   segment.size = 0.5, segment.color = "grey20",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 1.0, point.padding = 0.45,
                   force = 35, force_pull = 0.15,
                   label.padding = unit(1, "pt"),
                   label.r = unit(0.5, "pt"),
                   label.size = cfg$label_border_size %||% 0.10, seed = 42) +
  annotate("label", x = nes_lim, y = nes_lim,
           label = sprintf("%s  n = %d", qd$label_tr, n_q1),
           hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
           color = qd$color_tr, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -nes_lim, y = nes_lim,
           label = sprintf("%s  n = %d", qd$label_tl, n_q2),
           hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
           color = qd$color_tl, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -nes_lim, y = -nes_lim,
           label = sprintf("%s  n = %d", qd$label_bl, n_q3),
           hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
           color = qd$color_bl, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = nes_lim, y = -nes_lim,
           label = sprintf("%s  n = %d", qd$label_br, n_q4),
           hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
           color = qd$color_br, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt")) +
  scale_x_continuous(expand = expansion(0, 0)) +
  scale_y_continuous(expand = expansion(0, 0)) +
  coord_fixed(ratio = 1, xlim = c(-nes_lim, nes_lim), ylim = c(-nes_lim, nes_lim)) +
  labs(title = cfg$title,
       subtitle = subtitle_str,
       x = sprintf("NES (%s)", cfg$axis_x_label),
       y = sprintf("NES (%s)", cfg$axis_y_label)) +
  FIG_THEME +
  theme(
    axis.text         = element_text(size = FIG_AXIS_TEXT, face = "bold", color = "grey30"),
    axis.title        = element_text(size = FIG_AXIS_TEXT, face = "bold"),
    legend.position   = "bottom",
    legend.title      = element_text(size = FIG_LEGEND_TITLE, face = "bold", color = "grey25"),
    legend.text       = element_text(size = FIG_LEGEND_TEXT, color = "grey20"),
    legend.key.size   = unit(2 * PRINT_SCALE, "mm"),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box        = "horizontal",
    legend.box.just   = "center",
    legend.spacing.x  = unit(3 * PRINT_SCALE, "mm"),
    legend.box.margin = margin(t = -2),
    plot.margin       = margin(0, 0, 0, 0)
  ) +
  guides(fill  = "none",
         shape = guide_legend(nrow = 1, order = 1,
                               keyheight = unit(4 * PRINT_SCALE, "mm"),
                               keywidth  = unit(4 * PRINT_SCALE, "mm"),
                               override.aes = list(size = 3 * PRINT_SCALE, fill = "grey50")),
         size  = guide_legend(nrow = 1, order = 2,
                               keyheight = unit(4 * PRINT_SCALE, "mm"),
                               keywidth  = unit(4 * PRINT_SCALE, "mm")))

ggsave(file.path(RPT_PNG, "MAIN_panel_D_nes_scatter.png"), pD,
       width = PG_W, height = PG_W, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "MAIN_panel_D_nes_scatter.pdf"), pD,
       width = PG_W, height = PG_W, units = "mm", device = pdf_device)

export_df <- fgsea_wide |>
  transmute(
    pathway, pathway_label, database,
    !!paste0("NES_", cx) := round(.data[[nes_x]], 3),
    !!paste0("NES_", cy) := round(.data[[nes_y]], 3),
    !!paste0("padj_", cx) := signif(.data[[padj_x]], 4),
    !!paste0("padj_", cy) := signif(.data[[padj_y]], 4),
    significance = as.character(significance),
    set_size
  ) |>
  arrange(significance, desc(abs(.data[[nes_x]]) + abs(.data[[nes_y]])))
write_csv(export_df, file.path(DAT, "panel_D", "nes_scatter.csv"))

pD_legend_grob <- cowplot::get_plot_component(pD, "guide-box-bottom", return_all = FALSE)
if (!is.null(pD_legend_grob)) {
  pD_legend_plot <- cowplot::ggdraw(pD_legend_grob)
  ggsave(file.path(RPT_PNG, "MAIN_panel_D_legend.png"), pD_legend_plot,
         width = 120, height = 14, units = "mm", dpi = 300)
}

pD_title    <- cfg$title
pD_subtitle <- subtitle_str
pD_legend   <- NULL
# Strip titles but KEEP legend (legend provides shape/size key for composite)
pD          <- pD + labs(title = NULL, subtitle = NULL, tag = NULL)

pw_conc_frac <- metric_frac  # F04 stitcher
pw_rev_frac  <- metric_frac  # F05 stitcher

cat(sprintf("%s Panel D done\n", cfg$fig_id))
