# F06_WGCNA — Main Panels
# Panel A: Module-trait LMM heatmap (r_equiv, BH stars, model sections)
# Panel B: Top-module eigengene boxplots (pre/post by Group_Time, paired lines)
#
# Run from A_CvH_2026/ root after 00_run_wgcna.R.

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

source("04_Figures/shared/style.R")
# Override pdf device — cairo_pdf unavailable on this machine
pdf_device <- grDevices::pdf

BASE      <- "02-03_Sam's_Results/04_Figures/F06_WGCNA"
PANEL_DIR <- file.path(BASE, "c_data")
WGCNA_DIR <- file.path(BASE, "c_data", "wgcna")
RPT_PDF   <- file.path(BASE, "b_reports", "main", "pdf", "panels")
RPT_PNG   <- file.path(BASE, "b_reports", "main", "png", "panels")
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)

# pdf_device already set above

# ── Load artifacts ────────────────────────────────────────────────────────────

lmm_df         <- read_csv(file.path(WGCNA_DIR, "wgcna_lmm_contrasts.csv"),
                            show_col_types = FALSE)
mod_bio_labels  <- read_csv(file.path(PANEL_DIR, "mod_bio_labels.csv"),
                            show_col_types = FALSE)
MEs             <- readRDS(file.path(PANEL_DIR, "MEs.rds"))
meta            <- read_csv(file.path(PANEL_DIR, "meta.csv"), show_col_types = FALSE)

# ── Panel A: Module-trait LMM heatmap ────────────────────────────────────────

# Contrast ordering mirrors main CvH pipeline (CRvH first, CR second)
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

# Module ordering: largest first (matches heatmap row order top-to-bottom)
mod_order_colors <- mod_bio_labels |>
  arrange(desc(n_proteins)) |>
  dplyr::pull(module_color)
mod_order <- paste0("ME", mod_order_colors)
mod_order <- intersect(mod_order, unique(lmm_df$module))

# Display labels (M1: top_pathway)
mod_display_vec <- setNames(mod_bio_labels$display_label,
                             paste0("ME", mod_bio_labels$module_color))

# Gene count bar data
gene_counts <- setNames(mod_bio_labels$n_proteins,
                         paste0("ME", mod_bio_labels$module_color))

# Build tile data frame
tile_df <- expand.grid(
  module   = mod_order,
  contrast = all_contrasts,
  stringsAsFactors = FALSE
) |> as_tibble() |>
  left_join(
    lmm_df |> dplyr::select(module, contrast, r_equiv, p_raw, p_bh),
    by = c("module", "contrast")
  ) |>
  mutate(
    model    = if_else(contrast %in% crvh_contrasts, "CRvH", "CR"),
    sig_tier = case_when(
      p_bh  < 0.05 ~ "FDR",
      p_raw < 0.05 ~ "nominal",
      TRUE         ~ "ns"
    ),
    label  = paste0(
      sprintf("%.2f", coalesce(r_equiv, NA_real_)),
      ifelse(!is.na(p_bh) & p_bh < 0.001, "\n***",
      ifelse(!is.na(p_bh) & p_bh < 0.01,  "\n**",
      ifelse(!is.na(p_bh) & p_bh < 0.05,  "\n*", "")))),
    mod_label = mod_display_vec[module],
    ctr_label = contrast_labels[contrast]
  ) |>
  mutate(
    module   = factor(module, levels = rev(mod_order)),
    contrast = factor(contrast, levels = all_contrasts)
  )

# Gene count bar
count_df <- tibble(
  module = factor(mod_order, levels = rev(mod_order)),
  n = gene_counts[mod_order]
)
p_counts <- ggplot(count_df, aes(x = n, y = module)) +
  geom_col(fill = "grey70", width = 0.7) +
  geom_text(aes(label = n), hjust = -0.1, size = 2.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.35))) +
  labs(x = "Proteins", y = NULL) +
  FIG_THEME +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid   = element_blank(),
    panel.border = element_blank(),
    axis.line.x  = element_line(color = "black", linewidth = 0.3)
  )

# Main heatmap
n_mod <- length(mod_order)
p_heat <- ggplot(tile_df, aes(x = contrast, y = module)) +
  geom_tile(aes(fill = r_equiv), color = "white", linewidth = 0.5) +
  geom_tile(data = tile_df |> dplyr::filter(sig_tier == "FDR"),
            aes(fill = r_equiv), color = "black", linewidth = 0.8) +
  geom_tile(data = tile_df |> dplyr::filter(sig_tier == "nominal"),
            aes(fill = r_equiv), color = "grey30", linewidth = 0.5,
            linetype = "dashed") +
  geom_text(aes(label = label), size = 2.3, lineheight = 0.85) +
  scale_fill_gradient2(
    low = "#4393C3", mid = "white", high = "#D6604D",
    midpoint = 0, limits = c(-1, 1), na.value = "grey92",
    name = expression(r[equiv]),
    guide = guide_colorbar(barwidth = unit(3, "mm"), barheight = unit(40, "mm"))
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
       subtitle = "Solid border = FDR < 0.05 (BH); dashed = nominal p < 0.05") +
  coord_cartesian(clip = "off") +
  FIG_THEME +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 0, size = 8),
    axis.text.y  = element_text(size = 7.5),
    plot.margin  = margin(18, 5, 5, 5),
    legend.position = "right"
  )

panel_A <- p_heat + p_counts + plot_layout(widths = c(6, 1))

W_A <- 240; H_A <- max(120, n_mod * 12 + 20)
ggsave(file.path(RPT_PDF, "MAIN_panel_A_heatmap.pdf"), panel_A,
       width = W_A, height = H_A, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT_PNG, "MAIN_panel_A_heatmap.png"), panel_A,
       width = W_A, height = H_A, units = "mm", dpi = 300, limitsize = FALSE)

# Save panel data
write_csv(
  tile_df |> dplyr::select(module, contrast, model, r_equiv, p_raw, p_bh, sig_tier),
  file.path(PANEL_DIR, "01_panel_A_heatmap_data.csv"))

message("Panel A saved")

# ── Panel B: Top-module eigengene boxplots ────────────────────────────────────

# Select top 2 modules: strongest |r_equiv| for Cancer_vs_Healthy contrast
top2_mods <- lmm_df |>
  dplyr::filter(contrast == "Cancer_vs_Healthy", !is.na(r_equiv)) |>
  arrange(desc(abs(r_equiv))) |>
  head(2) |>
  dplyr::pull(module)

if (length(top2_mods) == 0) {
  top2_mods <- colnames(MEs)[1:min(2, ncol(MEs))]
}

me_long <- as.data.frame(MEs) |>
  rownames_to_column("sample_id") |>
  pivot_longer(-sample_id, names_to = "module", values_to = "eigengene") |>
  dplyr::filter(module %in% top2_mods) |>
  left_join(meta |> dplyr::select(sample_id, pid, group, cancer, timepoint),
            by = "sample_id") |>
  dplyr::filter(!is.na(group))

# Paired lines only for SURV (T1-T2 pairs)
surv_paired <- me_long |>
  dplyr::filter(cancer == "SURV") |>
  dplyr::group_by(module, pid) |>
  dplyr::filter(n() == 2) |>
  dplyr::ungroup()

# Module display labels for facet
mod_label_map <- setNames(
  mod_bio_labels$display_label,
  paste0("ME", mod_bio_labels$module_color)
)
me_long <- me_long |>
  mutate(mod_label = coalesce(mod_label_map[module], module))
surv_paired <- surv_paired |>
  mutate(mod_label = coalesce(mod_label_map[module], module))

group_levels <- c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1")
me_long    <- me_long    |> mutate(group = factor(group, levels = group_levels))
surv_paired <- surv_paired |> mutate(group = factor(group, levels = group_levels))

# p-value annotations from LMM
lmm_annot <- lmm_df |>
  dplyr::filter(module %in% top2_mods, contrast == "Cancer_vs_Healthy") |>
  mutate(
    mod_label  = coalesce(mod_label_map[module], module),
    star_label = paste0("r=", sprintf("%.2f", r_equiv),
                        ifelse(p_bh < 0.05, paste0(" (", sig_stars(p_bh), ")"), " (ns)"))
  )

panel_B <- ggplot(me_long, aes(x = group, y = eigengene)) +
  geom_boxplot(aes(fill = group), alpha = 0.7, width = 0.5,
               outlier.shape = NA, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1, alpha = 0.5, color = "grey30") +
  geom_line(data = surv_paired,
            aes(group = pid), color = "grey50", alpha = 0.4,
            linewidth = 0.35) +
  scale_fill_manual(values = GROUP_COLORS, guide = "none") +
  geom_text(data = lmm_annot, aes(label = star_label),
            x = 3, y = Inf, vjust = 1.5, hjust = 0.5,
            size = 2.8, fontface = "italic", color = "grey20") +
  facet_wrap(~ mod_label, scales = "free_y", ncol = 2) +
  labs(x = "Group_Time", y = "Module eigengene",
       title = "Top-Module Eigengenes by Group (Cancer_vs_Healthy)",
       subtitle = "Boxes: median+IQR; points: samples; lines: paired SURV subjects") +
  FIG_THEME +
  theme(
    axis.text.x  = element_text(angle = 35, hjust = 1, size = 7),
    strip.text   = element_text(size = 7)
  )

W_B <- 180; H_B <- 100
ggsave(file.path(RPT_PDF, "MAIN_panel_B_eigengene.pdf"), panel_B,
       width = W_B, height = H_B, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT_PNG, "MAIN_panel_B_eigengene.png"), panel_B,
       width = W_B, height = H_B, units = "mm", dpi = 300, limitsize = FALSE)

write_csv(me_long, file.path(PANEL_DIR, "02_panel_B_eigengene_long.csv"))

message("Panel B saved")
message("01_main_panels.R complete")
