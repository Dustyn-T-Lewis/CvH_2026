# F03/CRvH Supplementary: Cross-Contrast fGSEA Concordance Heatmap
# NES direction x significance across Cancer_vs_Healthy, Training_CR
# for pathways significant in at least one contrast.
# Reads cached fGSEA results from panel_C (no new analysis).

setwd(here::here())
source("04_Figures/shared/style.R")

pacman::p_load(tidyverse, ComplexHeatmap, circlize)

RPT <- "04_Figures/F03/CRvH/b_reports/supp"

RPT_PDF       <- file.path(RPT, "main", "pdf")

RPT_PNG       <- file.path(RPT, "main", "png")

RPT_SUPP_PDF  <- file.path(RPT, "supp", "pdf")

RPT_SUPP_PNG  <- file.path(RPT, "supp", "png")
DAT <- "04_Figures/F03/CRvH/c_data/supp"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# --- Load existing fGSEA results (no new analysis) ---
fgsea_df <- read_csv("04_Figures/F03/CRvH/c_data/01_panel_C_fgsea_results.csv",
                     show_col_types = FALSE)

# Filter to curated databases (exclude full GO:BP -- too many terms)
keep_db <- c("Hallmark", "KEGG", "Reactome", "GO Slim")
sig_df <- fgsea_df %>%
  filter(database %in% keep_db) %>%
  group_by(pathway) %>%
  filter(any(padj < 0.05)) %>%
  ungroup()

message(sprintf("  %d pathway-contrast pairs from %d unique pathways",
                nrow(sig_df), n_distinct(sig_df$pathway)))

# Contrast ordering
ctr_order <- c("Cancer_vs_Healthy", "Training_CR")
sig_df <- sig_df %>%
  filter(contrast %in% ctr_order) %>%
  mutate(contrast = factor(contrast, levels = ctr_order))

# --- Build NES matrix ---
nes_wide <- sig_df %>%
  select(pathway, contrast, NES) %>%
  pivot_wider(names_from = contrast, values_from = NES) %>%
  column_to_rownames("pathway")

padj_wide <- sig_df %>%
  select(pathway, contrast, padj) %>%
  pivot_wider(names_from = contrast, values_from = padj) %>%
  column_to_rownames("pathway")

# Ensure same column order
nes_mat  <- as.matrix(nes_wide[, ctr_order])
padj_mat <- as.matrix(padj_wide[, ctr_order])

# Replace NA NES with 0 (pathway not tested in that contrast)
nes_mat[is.na(nes_mat)]   <- 0
padj_mat[is.na(padj_mat)] <- 1

# Clean pathway labels
row_labels <- clean_pathway_name(rownames(nes_mat))

# Significance stars for cell annotation
sig_stars_mat <- matrix("", nrow = nrow(padj_mat), ncol = ncol(padj_mat))
sig_stars_mat[padj_mat < 0.05]  <- "*"
sig_stars_mat[padj_mat < 0.01]  <- "**"
sig_stars_mat[padj_mat < 0.001] <- "***"

# Column labels
col_labels <- c("CR vs H", "Tr.(CR)")

# --- Draw heatmap ---
max_abs <- max(abs(nes_mat), na.rm = TRUE)
col_fun <- colorRamp2(c(-max_abs, 0, max_abs), c("#2166AC", "white", "#B2182B"))

ht <- Heatmap(
  nes_mat,
  name = "NES",
  col  = col_fun,
  row_labels        = row_labels,
  column_labels     = col_labels,
  column_names_rot  = 0,
  cluster_columns   = FALSE,
  clustering_method_rows = "ward.D2",
  row_names_gp      = gpar(fontsize = 5),
  column_names_gp   = gpar(fontsize = 8, fontface = "bold"),
  row_names_max_width = unit(70, "mm"),
  layer_fun = function(j, i, x, y, w, h, fill) {
    stars <- pindex(sig_stars_mat, i, j)
    show <- stars != ""
    if (any(show)) {
      grid.text(stars[show], x[show], y[show],
                gp = gpar(fontsize = 5, col = "black"))
    }
  },
  column_title = "Cross-Contrast fGSEA Concordance (CRvH)",
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  heatmap_legend_param = list(
    title = "NES",
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7),
    legend_width = unit(30, "mm"),
    direction = "horizontal"
  ),
  show_row_dend = TRUE,
  row_dend_width = unit(10, "mm")
)

# Size based on number of pathways
n_pw <- nrow(nes_mat)
ht_h <- max(80, n_pw * 2.5 + 30)
ht_w <- 160

pdf(file.path(RPT_SUPP_PDF, "gsea_concordance_heatmap_SUPP.pdf"),
    width = ht_w / 25.4, height = ht_h / 25.4)
draw(ht, heatmap_legend_side = "bottom", padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()

png(file.path(RPT_SUPP_PNG, "gsea_concordance_heatmap_SUPP.png"),
    width = ht_w, height = ht_h, units = "mm", res = 300)
draw(ht, heatmap_legend_side = "bottom", padding = unit(c(2, 2, 2, 2), "mm"))
dev.off()

# --- Export NES matrix for audit ---
nes_export <- as_tibble(nes_mat, rownames = "pathway") %>%
  mutate(pathway_label = clean_pathway_name(pathway))
write_csv(nes_export, file.path(DAT, "gsea_concordance_matrix.csv"))

message(sprintf("  Heatmap: %d pathways x %d contrasts, height = %.0f mm",
                n_pw, ncol(nes_mat), ht_h))
cat("F03/CRvH Supplementary: GSEA concordance heatmap done.\n")
