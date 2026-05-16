#!/usr/bin/env Rscript
# F04 Supplementary: Comprehensive Enrichment — Training Blunting
# ComplexHeatmap showing pathway-level blunting split by response pattern
# Method: fGSEA on GO:BP + Reactome + Hallmark + KEGG_REF, reduced via
# collapsePathways() (Jaccard dedup disabled — collapsePathways sufficient).
#
# Sourced by 02_supp_panels.R — expects style.R already loaded.

source("02-03_Sam's_Results/04_Figures/shared/pathway_utils.R")

library(ComplexHeatmap)
library(circlize)
library(gridExtra)

BASE    <- "02-03_Sam's_Results/04_Figures/F04_training_concordance"
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
DAT     <- file.path(BASE, "c_data")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_supp"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- grDevices::pdf  # force base pdf (cairo DLL fails)

# Load data
dep <- read_csv("02-03_Sam's_Results/03_DEP/c_data/F04_combined_CREvPLA.csv",
                show_col_types = FALSE)

contrasts <- c("Training_CRE", "Training_PLA", "Supplement_Interaction")
stats_list <- setNames(lapply(contrasts, function(ctr) {
  col <- paste0("t_", ctr)
  if (!col %in% names(dep)) {
    message(sprintf("  Note: column %s not found, skipping contrast", col))
    return(NULL)
  }
  s <- setNames(dep[[col]], dep$gene)
  s[!is.na(s)]
}), contrasts)
# Remove NULLs
stats_list <- Filter(Negate(is.null), stats_list)

# Build pathway collection (no VARIANT, no WP)
pw_list <- build_pathway_collection(
  min_size = 15, max_size = 500,
  include_goslim = FALSE,
  exclude_variants = TRUE
)

# Run enrichment pipeline
set.seed(42)
ep <- run_enrichment_pipeline(
  stats_list     = stats_list,
  pw_list        = pw_list,
  jaccard_cutoff = 1,
  nperm          = 10000,
  min_size       = 15,
  max_size       = 500
)

long_df   <- ep$long_df
sig_union <- ep$sig_union

# Prepare wide data for ComplexHeatmap
pw_meta <- long_df |>
  group_by(pathway) |>
  summarise(size = max(size, na.rm = TRUE),
            database = first(database), .groups = "drop")

wide_df <- long_df |>
  dplyr::select(pathway, contrast, NES, padj) |>
  mutate(
    sig = !is.na(padj) & padj < 0.05,
    sig_label = case_when(
      !is.na(padj) & padj < 0.001 ~ "***",
      !is.na(padj) & padj < 0.01  ~ "**",
      !is.na(padj) & padj < 0.05  ~ "*",
      TRUE ~ NA_character_
    )
  ) |>
  pivot_wider(
    id_cols = pathway,
    names_from = contrast,
    values_from = c(NES, padj, sig, sig_label)
  ) |>
  left_join(pw_meta, by = "pathway") |>
  mutate(
    sig_Training_CRE = coalesce(sig_Training_CRE, FALSE),
    sig_Training_PLA = coalesce(sig_Training_PLA, FALSE)
  ) |>
  (\(df) {
    # Supplement_Interaction may be absent if that contrast had no t_ column
    if ("sig_Supplement_Interaction" %in% names(df)) {
      mutate(df, sig_Supplement_Interaction = coalesce(sig_Supplement_Interaction, FALSE))
    } else {
      mutate(df,
             sig_Supplement_Interaction       = FALSE,
             NES_Supplement_Interaction       = NA_real_,
             sig_label_Supplement_Interaction = NA_character_)
    }
  })()

# Pattern classification
PATTERN_COLS <- c(
  "Shared"           = "#457B9D",
  "Blunted"          = "#E63946",
  "Age-resistant"    = "#2A9D8F",
  "Interaction only" = "#9B5DE5"
)
pattern_levels <- c("Shared", "Blunted", "Age-resistant", "Interaction only")

pattern_df <- wide_df |>
  mutate(
    pattern = case_when(
      sig_Training_CRE & sig_Training_PLA          ~ "Shared",
      sig_Training_CRE & !sig_Training_PLA          ~ "Blunted",
      !sig_Training_CRE & sig_Training_PLA          ~ "Age-resistant",
      sig_Supplement_Interaction                    ~ "Interaction only",
      TRUE                                          ~ "Other"
    ),
    sort_key = case_when(
      pattern == "Age-resistant"    ~ NES_Training_PLA,
      pattern == "Interaction only" ~ NES_Supplement_Interaction,
      TRUE                         ~ NES_Training_CRE
    )
  ) |>
  filter(pattern %in% pattern_levels)

row_ord <- order(match(pattern_df$pattern, pattern_levels), -pattern_df$sort_key)
pattern_df <- pattern_df[row_ord, ]

for (p in pattern_levels) {
  message(sprintf("  %s: %d", p, sum(pattern_df$pattern == p)))
}

# Two-panel layout: Blunted (left) | Shared + Age-resistant + Interaction (right)
n_pw <- nrow(pattern_df)
left_df  <- pattern_df |> filter(pattern == "Blunted")
right_df <- pattern_df |> filter(pattern %in% c("Shared", "Age-resistant", "Interaction only"))
right_levels <- c("Shared", "Age-resistant", "Interaction only")

col_fun <- colorRamp2(
  c(-3, -1.5, 0, 1.5, 3),
  c("#08306B", "#4393C3", "white", "#D6604D", "#67000D")
)

nes_cols <- c("NES_Training_CRE", "NES_Training_PLA", "NES_Supplement_Interaction")
sig_cols <- c("sig_label_Training_CRE", "sig_label_Training_PLA", "sig_label_Supplement_Interaction")
col_labs <- c("Training (CRE)", "Training (PLA)", "Supplement Interaction")

make_layer_fun <- function(sig_mat, val_mat) {
  force(sig_mat); force(val_mat)
  function(j, i, x, y, w, h, fill) {
    lab <- pindex(sig_mat, i, j)
    sel <- !is.na(lab) & lab != ""
    if (any(sel)) {
      grid.text(lab[sel], x[sel], y[sel],
                gp = gpar(col = "white", fontsize = 9, fontface = "bold"))
    }
    v <- pindex(val_mat, i, j)
    na_sel <- is.na(v)
    if (any(na_sel)) {
      grid.text("N/A", x[na_sel], y[na_sel],
                gp = gpar(col = "grey50", fontsize = 6, fontface = "italic"))
    }
  }
}

# Left panel: Blunted
nes_L <- as.matrix(left_df[, nes_cols]); colnames(nes_L) <- col_labs
sig_L <- as.matrix(left_df[, sig_cols])

ht_L <- Heatmap(
  nes_L, name = "NES", col = col_fun, na_col = "grey95",
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_labels = clean_pathway_name(left_df$pathway),
  row_names_gp = gpar(fontsize = 6.5),
  column_names_gp = gpar(fontsize = 9, fontface = "bold"),
  column_names_rot = 30, column_names_side = "top",
  layer_fun = make_layer_fun(sig_L, nes_L),
  column_title = sprintf("Blunted (%d)", nrow(left_df)),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  heatmap_legend_param = list(direction = "horizontal"),
  width = unit(55, "mm")
)

# Right panel: Shared + Age-resistant + Interaction only
nes_R <- as.matrix(right_df[, nes_cols]); colnames(nes_R) <- col_labs
sig_R <- as.matrix(right_df[, sig_cols])
pat_R <- factor(right_df$pattern, levels = right_levels)

ha_R <- rowAnnotation(
  Pattern = pat_R,
  col = list(Pattern = PATTERN_COLS[right_levels]),
  annotation_name_gp = gpar(fontsize = 8),
  annotation_legend_param = list(
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7)
  )
)

ht_R <- Heatmap(
  nes_R, name = "NES_R", show_heatmap_legend = FALSE,
  col = col_fun, na_col = "grey95",
  row_split = pat_R, row_gap = unit(2, "mm"),
  cluster_rows = FALSE, cluster_columns = FALSE,
  row_labels = clean_pathway_name(right_df$pathway),
  row_names_gp = gpar(fontsize = 6.5),
  column_names_gp = gpar(fontsize = 9, fontface = "bold"),
  column_names_rot = 30, column_names_side = "top",
  layer_fun = make_layer_fun(sig_R, nes_R),
  right_annotation = ha_R,
  row_title_gp = gpar(fontsize = 8, fontface = "bold"),
  width = unit(55, "mm")
)

# Capture and arrange side by side
g_L <- grid.grabExpr(draw(ht_L, heatmap_legend_side = "bottom",
                          merge_legend = TRUE, padding = unit(c(2, 5, 2, 2), "mm")))
g_R <- grid.grabExpr(draw(ht_R, heatmap_legend_side = "bottom",
                          merge_legend = TRUE, padding = unit(c(2, 5, 2, 2), "mm")))

max_rows <- max(nrow(left_df), nrow(right_df))
fig_h_mm <- max(180, 4.5 * max_rows + 60)
fig_w_mm <- 360
title_grob <- textGrob(
  sprintf("Pathway-Level Training Response  |  %d pathways", n_pw),
  gp = gpar(fontsize = 12, fontface = "bold")
)

# Export figure
png(file.path(RPT_PNG, "SUPP_enrichment_heatmap.png"),
    width = fig_w_mm, height = fig_h_mm, units = "mm", res = 300)
grid.arrange(g_L, g_R, ncol = 2, widths = c(1, 1.2), top = title_grob)
dev.off()

pdf(file.path(RPT_PDF, "SUPP_enrichment_heatmap.pdf"),
    width = fig_w_mm / 25.4, height = fig_h_mm / 25.4)
grid.arrange(g_L, g_R, ncol = 2, widths = c(1, 1.2), top = title_grob)
dev.off()

# Export data
le_df <- long_df |>
  filter(pathway %in% sig_union) |>
  mutate(leading_edge = vapply(leadingEdge, paste, character(1), collapse = ";")) |>
  dplyr::select(pathway, contrast, leading_edge)

export_df <- long_df |>
  filter(pathway %in% pattern_df$pathway) |>
  dplyr::select(pathway, contrast, NES, padj, size, database) |>
  mutate(
    sig = !is.na(padj) & padj < 0.05,
    pathway_label = clean_pathway_name(pathway)
  ) |>
  left_join(
    pattern_df |> dplyr::select(pathway, pattern, sort_key),
    by = "pathway"
  ) |>
  left_join(le_df, by = c("pathway", "contrast")) |>
  mutate(NES = round(NES, 3), padj = signif(padj, 4)) |>
  arrange(match(pattern, pattern_levels), -sort_key, contrast)

write_csv(export_df, file.path(DAT, "panel_supp", "enrichment_blunting.csv"))

message(sprintf("F04 supplementary enrichment done: %d pathways, saved to %s",
                n_pw, RPT_PNG))
