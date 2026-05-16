# F06_WGCNA — Main Panels (YvO-aligned structure)
# Panel A: Module-trait LMM heatmap (r_equiv, BH stars, CRvH + CR model sections)
# Panel B: Module-level NES scatters (fGSEA on module-member t-stat ranks per contrast)
#
# Run from A_CvH_2026/ root after 00_run_wgcna.R.

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(cowplot)
  library(ggrepel)
  library(fgsea)
})

source("02-03_Sam's_Results/04_Figures/shared/style.R")
source("02-03_Sam's_Results/04_Figures/shared/pathway_utils.R")

pdf_device <- grDevices::pdf  # force base pdf (cairo DLL fails on this system)

BASE      <- "02-03_Sam's_Results/04_Figures/F06_WGCNA"
PANEL_DIR <- file.path(BASE, "c_data")
WGCNA_DIR <- file.path(BASE, "c_data", "wgcna")
RPT_PDF   <- file.path(BASE, "b_reports", "main", "pdf", "panels")
RPT_PNG   <- file.path(BASE, "b_reports", "main", "png", "panels")
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)

# ── Load artifacts ─────────────────────────────────────────────────────────────

lmm_df         <- read_csv(file.path(WGCNA_DIR, "wgcna_lmm_contrasts.csv"),
                            show_col_types = FALSE)
mod_bio_labels  <- read_csv(file.path(PANEL_DIR, "mod_bio_labels.csv"),
                            show_col_types = FALSE)
module_df       <- read_csv(file.path(WGCNA_DIR, "wgcna_module_assignments.csv"),
                            show_col_types = FALSE)

# ── Panel A: Module-trait LMM heatmap ─────────────────────────────────────────

# Contrast ordering mirrors main CvH pipeline
crvh_contrasts <- c("Cancer_vs_Healthy", "Training_CR")
cr_contrasts   <- c("Baseline_Supplement", "Training_CRE",
                    "Training_PLA", "Supplement_Interaction")
all_contrasts  <- c(crvh_contrasts, cr_contrasts)

contrast_labels <- c(
  Cancer_vs_Healthy      = "CR vs H",
  Training_CR            = "Tr.(CR)",
  Baseline_Supplement    = "BL(CRE-PLA)",
  Training_CRE           = "Tr.(CRE)",
  Training_PLA           = "Tr.(PLA)",
  Supplement_Interaction = "CRExPLA"
)

# Module ordering: largest first
mod_order_colors <- mod_bio_labels |>
  arrange(desc(n_proteins)) |>
  pull(module_color)
mod_order <- paste0("ME", mod_order_colors)
mod_order <- intersect(mod_order, unique(lmm_df$module))

# Display labels and counts
mod_display_vec <- setNames(mod_bio_labels$display_label,
                             paste0("ME", mod_bio_labels$module_color))
gene_counts <- setNames(mod_bio_labels$n_proteins,
                         paste0("ME", mod_bio_labels$module_color))

# Build tile data frame
tile_df <- expand.grid(
  module   = mod_order,
  contrast = all_contrasts,
  stringsAsFactors = FALSE
) |> as_tibble() |>
  left_join(
    lmm_df |> select(module, contrast, r_equiv, p_raw, p_bh),
    by = c("module", "contrast")
  ) |>
  mutate(
    model    = if_else(contrast %in% crvh_contrasts, "CRvH", "CR"),
    sig_tier = case_when(
      p_bh  < 0.05 ~ "FDR",
      p_raw < 0.05 ~ "nominal",
      TRUE         ~ "ns"
    ),
    stars = case_when(
      p_bh < 0.001 ~ "***",
      p_bh < 0.01  ~ "**",
      p_bh < 0.05  ~ "*",
      TRUE         ~ ""
    ),
    label     = sprintf("%.2f", coalesce(r_equiv, NA_real_)),
    mod_label = mod_display_vec[module],
    ctr_label = contrast_labels[contrast]
  ) |>
  mutate(
    module   = factor(module, levels = rev(mod_order)),
    contrast = factor(contrast, levels = all_contrasts)
  )

# Protein count bar (module colours matching WGCNA assignment)
mod_color_raw <- gsub("^ME", "", as.character(rev(mod_order)))
count_df <- tibble(
  module   = factor(rev(mod_order), levels = rev(mod_order)),
  n        = gene_counts[rev(mod_order)],
  mod_col  = mod_color_raw
)
p_counts <- ggplot(count_df, aes(x = n, y = module)) +
  geom_col(fill = count_df$mod_col, color = "black",
           linewidth = 0.25, width = 0.65) +
  geom_text(aes(label = n), hjust = -0.1, size = 2.4, fontface = "bold") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.30))) +
  labs(x = "Proteins", y = NULL) +
  FIG_THEME +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid   = element_blank(),
    panel.border = element_blank(),
    axis.line.x  = element_line(color = "black", linewidth = 0.3),
    plot.margin  = margin(2, 4, 2, 0)
  )

# Main heatmap
n_mod <- length(mod_order)
p_heat <- ggplot(tile_df, aes(x = contrast, y = module)) +
  geom_tile(aes(fill = r_equiv), color = "white", linewidth = 0.5) +
  geom_tile(data = tile_df |> filter(sig_tier == "FDR"),
            aes(fill = r_equiv), color = "black", linewidth = 0.9) +
  geom_tile(data = tile_df |> filter(sig_tier == "nominal"),
            aes(fill = r_equiv), color = "grey30", linewidth = 0.5,
            linetype = "dashed") +
  geom_text(aes(label = label), size = 2.3, color = "grey20") +
  geom_text(data = tile_df |> filter(stars != ""),
            aes(label = stars), size = 2.3, vjust = -0.6, color = "grey10") +
  scale_fill_gradient2(
    low = "#4393C3", mid = "white", high = "#D6604D",
    midpoint = 0, limits = c(-1, 1), na.value = "grey92",
    name = "r-equiv",
    guide = guide_colorbar(barwidth = unit(3, "mm"), barheight = unit(35, "mm"))
  ) +
  scale_x_discrete(labels = contrast_labels, position = "top") +
  scale_y_discrete(labels = function(x) {
    lbl <- mod_display_vec[x]
    ifelse(is.na(lbl), x, lbl)
  }) +
  geom_vline(xintercept = length(crvh_contrasts) + 0.5,
             linewidth = 1.2, color = "black") +
  annotate("text", x = mean(seq_along(crvh_contrasts)), y = n_mod + 0.75,
           label = "CRvH Model", fontface = "bold", size = 3.5, hjust = 0.5) +
  annotate("text",
           x = length(crvh_contrasts) + mean(seq_along(cr_contrasts)),
           y = n_mod + 0.75,
           label = "CR Model", fontface = "bold", size = 3.5, hjust = 0.5) +
  labs(y = NULL, x = NULL,
       title = "WGCNA Module-Trait LMM Associations (Sam CvH, N=35)",
       subtitle = "Solid border = FDR < 0.05 (BH); dashed = nominal p < 0.05  |  * p<0.05  ** p<0.01  *** p<0.001") +
  coord_cartesian(clip = "off") +
  FIG_THEME +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 0, size = 8),
    axis.text.y  = element_text(size = 7.5),
    plot.margin  = margin(18, 5, 5, 5),
    legend.position = "right"
  )

panel_A <- p_heat + p_counts + plot_layout(widths = c(6, 1))

W_A <- 240; H_A <- max(120, n_mod * 14 + 25)
ggsave(file.path(RPT_PDF, "MAIN_panel_A_heatmap.pdf"), panel_A,
       width = W_A, height = H_A, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT_PNG, "MAIN_panel_A_heatmap.png"), panel_A,
       width = W_A, height = H_A, units = "mm", dpi = 300, limitsize = FALSE)

write_csv(
  tile_df |> select(module, contrast, model, r_equiv, p_raw, p_bh, sig_tier),
  file.path(PANEL_DIR, "01_panel_A_heatmap_data.csv"))
message("Panel A saved")

# ── Panel B: Module-level NES scatters ────────────────────────────────────────
# Mirrors YvO _panel_B_nes_scatters.R:
#   fGSEA run on per-module t-stat ranks (module gene sets as pathways),
#   scatter 1: NES(Cancer_vs_Healthy) vs NES(Training_CR) — concordance
#   scatter 2: NES(Cancer_vs_Healthy) vs NES(Training_CRE) — reversal

stopifnot(
  "CRvH DEP results missing" =
    file.exists("02-03_Sam's_Results/03_DEP/c_data/03_combined_results_CRvH.csv"),
  "CR DEP results missing" =
    file.exists("02-03_Sam's_Results/03_DEP/c_data/03_combined_results_CR.csv"),
  "WGCNA module assignments missing" =
    file.exists(file.path(WGCNA_DIR, "wgcna_module_assignments.csv"))
)

combined_crvh <- read_csv("02-03_Sam's_Results/03_DEP/c_data/03_combined_results_CRvH.csv",
                           show_col_types = FALSE)
combined_cr   <- read_csv("02-03_Sam's_Results/03_DEP/c_data/03_combined_results_CR.csv",
                           show_col_types = FALSE)

# Join all contrast t-stats on gene
combined <- combined_crvh |>
  select(gene, t_Cancer_vs_Healthy, t_Training_CR) |>
  left_join(
    combined_cr |> select(gene, t_Training_CRE, t_Training_PLA, t_Supplement_Interaction),
    by = "gene"
  ) |>
  filter(!is.na(gene), gene != "")

module_df_filt <- module_df |>
  filter(module_color != "grey", !is.na(gene), gene != "")
module_sets <- split(module_df_filt, module_df_filt$module_color) |>
  lapply(function(x) x$gene)
mod_sizes_n <- sapply(module_sets, length)

# Build t-stat rank vectors per contrast
build_ranks <- function(df, col) {
  vals <- df[[col]]; names(vals) <- df$gene
  vals <- vals[!is.na(vals)]; sort(vals, decreasing = TRUE)
}
ranks_CvH  <- build_ranks(combined, "t_Cancer_vs_Healthy")
ranks_TR   <- build_ranks(combined, "t_Training_CR")
ranks_CRE  <- build_ranks(combined, "t_Training_CRE")
ranks_PLA  <- build_ranks(combined, "t_Training_PLA")
ranks_INT  <- build_ranks(combined, "t_Supplement_Interaction")

run_module_fgsea <- function(ranks, module_sets) {
  res <- fgsea::fgseaMultilevel(pathways = module_sets, stats = ranks,
                                 minSize = 5, maxSize = 2000,
                                 nPermSimple = 10000, eps = 0)
  as.data.frame(res)
}

message("Panel B: running fGSEA on module gene sets for each contrast...")
fgsea_CvH <- run_module_fgsea(ranks_CvH,  module_sets)
fgsea_TR  <- run_module_fgsea(ranks_TR,   module_sets)
fgsea_CRE <- run_module_fgsea(ranks_CRE,  module_sets)
fgsea_PLA <- run_module_fgsea(ranks_PLA,  module_sets)
fgsea_INT <- run_module_fgsea(ranks_INT,  module_sets)

merge_fgsea <- function(res, suffix) {
  as_tibble(res) |>
    select(pathway, NES, padj, size) |>
    rename_with(~ paste0(., "_", suffix), c(NES, padj, size))
}

fgsea_wide <- merge_fgsea(fgsea_CvH, "CvH") |>
  left_join(merge_fgsea(fgsea_TR,  "TR"),  by = "pathway") |>
  left_join(merge_fgsea(fgsea_CRE, "CRE"), by = "pathway") |>
  left_join(merge_fgsea(fgsea_PLA, "PLA"), by = "pathway") |>
  left_join(merge_fgsea(fgsea_INT, "INT"), by = "pathway") |>
  mutate(
    module_color = pathway,
    n_proteins   = mod_sizes_n[pathway]
  ) |>
  left_join(mod_bio_labels |> select(module_color, display_label), by = "module_color") |>
  mutate(
    bio_label = ifelse(is.na(display_label), str_to_title(module_color), display_label),
    bio_label = str_wrap(bio_label, width = 14),
    # Concordance: both Cancer and Training same direction
    sig_conc = case_when(
      !is.na(padj_CvH) & padj_CvH < 0.05 & !is.na(padj_TR) & padj_TR < 0.05 ~ "Both sig",
      !is.na(padj_CvH) & padj_CvH < 0.05 ~ "Cancer only",
      !is.na(padj_TR)  & padj_TR  < 0.05 ~ "Training only",
      TRUE ~ "NS"),
    # Reversal: does Training_CRE oppose Cancer direction?
    sig_rev = case_when(
      !is.na(padj_CvH) & padj_CvH < 0.05 & !is.na(padj_CRE) & padj_CRE < 0.05 ~ "Both sig",
      !is.na(padj_CvH) & padj_CvH < 0.05 ~ "Cancer only",
      !is.na(padj_CRE) & padj_CRE < 0.05 ~ "Training only",
      TRUE ~ "NS")
  )

# Helper: is a color perceptually light?
is_light_hex <- function(hex) {
  rgb <- col2rgb(hex) / 255
  luminance <- 0.299 * rgb[1] + 0.587 * rgb[2] + 0.114 * rgb[3]
  luminance > 0.55
}

build_scatter <- function(df, x_col, y_col, x_lab, y_lab, quad_labels) {
  x_vals <- df[[x_col]]; y_vals <- df[[y_col]]
  nes_lim <- max(abs(c(x_vals, y_vals)), na.rm = TRUE) * 1.35

  sp <- cor.test(x_vals, y_vals, method = "spearman")
  rho_lbl <- sprintf("n = %d modules | rho = %.2f%s",
    nrow(df), sp$estimate,
    ifelse(sp$p.value < 0.001, ", p < 0.001",
           sprintf(", p = %.3f", sp$p.value)))

  q_tr <- sum(x_vals > 0 & y_vals > 0, na.rm = TRUE)
  q_bl <- sum(x_vals < 0 & y_vals < 0, na.rm = TRUE)
  q_tl <- sum(x_vals < 0 & y_vals > 0, na.rm = TRUE)
  q_br <- sum(x_vals > 0 & y_vals < 0, na.rm = TRUE)

  df$label_col <- sapply(df$module_color, function(mc) {
    if (is_light_hex(mc)) "black" else "white"
  })

  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
             fill = quad_labels$fill[1], alpha = 0.15) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
             fill = quad_labels$fill[2], alpha = 0.15) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
             fill = quad_labels$fill[3], alpha = 0.15) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
             fill = quad_labels$fill[4], alpha = 0.15) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                color = "black", linewidth = 0.3) +
    geom_point(aes(size = n_proteins), fill = df$module_color,
               color = "black", shape = 21, alpha = 0.85, stroke = 0.5) +
    geom_label_repel(aes(label = bio_label),
      fill = scales::alpha(df$module_color, 0.85), color = df$label_col,
      size = 4.0, fontface = "bold", lineheight = 0.75, max.overlaps = Inf,
      segment.size = 0.3, segment.color = "grey40", min.segment.length = 0,
      box.padding = 1.2, point.padding = 0.9, force = 80, force_pull = 0.05,
      label.padding = unit(1.5, "pt"), label.r = unit(1.5, "pt"),
      label.size = 0, seed = 7, show.legend = FALSE) +
    annotate("label", x = Inf,  y = Inf,  hjust = 1, vjust = 1,
             label = sprintf("%s  n=%d", quad_labels$label[1], q_tr),
             size = 3.5, fontface = "bold", color = quad_labels$text_col[1],
             fill = scales::alpha("white", 0.90), label.padding = unit(2, "pt")) +
    annotate("label", x = -Inf, y = -Inf, hjust = 0, vjust = 0,
             label = sprintf("%s  n=%d", quad_labels$label[2], q_bl),
             size = 3.5, fontface = "bold", color = quad_labels$text_col[2],
             fill = scales::alpha("white", 0.90), label.padding = unit(2, "pt")) +
    annotate("label", x = Inf,  y = -Inf, hjust = 1, vjust = 0,
             label = sprintf("%s  n=%d", quad_labels$label[3], q_br),
             size = 3.5, fontface = "bold", color = quad_labels$text_col[3],
             fill = scales::alpha("white", 0.90), label.padding = unit(2, "pt")) +
    annotate("label", x = -Inf, y = Inf,  hjust = 0, vjust = 1,
             label = sprintf("%s  n=%d", quad_labels$label[4], q_tl),
             size = 3.5, fontface = "bold", color = quad_labels$text_col[4],
             fill = scales::alpha("white", 0.90), label.padding = unit(2, "pt")) +
    scale_size_continuous(range = c(3, 9), name = "Proteins",
                          breaks = c(50, 150, 300, 500)) +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    scale_y_continuous(expand = expansion(mult = 0.02)) +
    coord_fixed(ratio = 1, xlim = c(-nes_lim, nes_lim),
                ylim = c(-nes_lim, nes_lim), clip = "off") +
    labs(subtitle = rho_lbl, x = x_lab, y = y_lab) +
    FIG_THEME +
    theme(
      axis.text    = element_text(size = 10, face = "bold", color = "grey30"),
      axis.title.x = element_text(size = 11, face = "bold"),
      axis.title.y = element_text(size = 11, face = "bold"),
      legend.position = "none",
      plot.subtitle   = element_text(size = 8),
      plot.margin  = margin(2, 4, 2, -4)
    )
}

quad_conc <- list(
  label    = c("Concordant Up", "Concordant Down", "Discordant", "Discordant"),
  fill     = c("#D6604D", "#D6604D", "#4393C3", "#4393C3"),
  text_col = c("#D6604D", "#D6604D", "#4393C3", "#4393C3"))

quad_rev <- list(
  label    = c("Exacerbated", "Exacerbated", "Reversed", "Reversed"),
  fill     = c("#D6604D", "#D6604D", "#2E7D32", "#2E7D32"),
  text_col = c("#D6604D", "#D6604D", "#2166AC", "#2166AC"))

p_top <- build_scatter(fgsea_wide,
  "NES_CvH", "NES_TR",
  "NES (Cancer vs Healthy)", "NES (Training CR)",
  quad_conc)

p_bottom <- build_scatter(fgsea_wide,
  "NES_CvH", "NES_CRE",
  "NES (Cancer vs Healthy)", "NES (Training CRE)",
  quad_rev)

scatters_panel <- (p_top / p_bottom) +
  plot_layout(heights = c(1, 1))

PB_W <- 220; PB_H <- 270

ggsave(file.path(RPT_PNG, "MAIN_panel_B_scatters.png"), scatters_panel,
       width = PB_W, height = PB_H, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "MAIN_panel_B_scatters.pdf"), scatters_panel,
       width = PB_W, height = PB_H, units = "mm", device = pdf_device)

# Separate legend PNG (for stitcher)
p_legend_src <- p_top +
  scale_size_continuous(range = c(3, 9), name = "Proteins",
    breaks = c(50, 150, 300, 500),
    guide = guide_legend(nrow = 1,
      override.aes = list(alpha = 0.7, fill = "grey60", stroke = 0))) +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text  = element_text(size = 10),
        legend.key   = element_rect(fill = NA, color = NA),
        legend.key.size = unit(5, "mm"),
        legend.background = element_rect(fill = NA, color = NA))
legend_grob <- cowplot::get_plot_component(p_legend_src,
                "guide-box-bottom", return_all = FALSE)
p_legend <- cowplot::ggdraw(legend_grob)
ggsave(file.path(RPT_PNG, "MAIN_panel_B_legend.png"), p_legend,
       width = 90, height = 16, units = "mm", dpi = 300)

# Save panel B data
fgsea_save <- fgsea_wide |>
  select(pathway, module_color, n_proteins, bio_label,
         NES_CvH, padj_CvH, NES_TR, padj_TR,
         NES_CRE, padj_CRE, NES_PLA, padj_PLA,
         NES_INT, padj_INT, sig_conc, sig_rev)
write_csv(fgsea_save, file.path(PANEL_DIR, "panel_B_module_fgsea.csv"))

message("Panel B NES scatters saved")
message("01_main_panels.R complete")
