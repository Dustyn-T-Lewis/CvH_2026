# Figure 5 — Panel A: Multi-Contrast NES Heatmap (Hallmark + GO Slim)
# Rows = pathways (clustered), Columns = 6 contrasts, Cells = NES
# Asterisks for padj < 0.05. Database color strip on left.
# Outputs: panel_A_nes_heatmap.pdf/png, nes_heatmap_data.csv

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

RPT <- "04_Figures/F05/b_reports"
DAT <- "04_Figures/F05/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# --- Build Hallmark + GO Slim pathway collection ---
pw_collection <- build_hallmark_goslim_collection(min_size = 10, max_size = 500)

# --- Load DEP results (t-statistics for ranking) ---
dep_crvh <- read_csv("03_DEP/c_data/03_combined_results_CRvH.csv",
                      show_col_types = FALSE)
dep_cr   <- read_csv("03_DEP/c_data/03_combined_results_CR.csv",
                      show_col_types = FALSE)

# --- Run fGSEA per contrast on Hallmark + GO Slim ---
set.seed(42)
fgsea_all_list <- list()

for (ctr in ALL_CONTRASTS) {
  tcol <- paste0("t_", ctr)
  dep_src <- if (ctr %in% ALL_CONTRASTS_CRVH) dep_crvh else dep_cr
  stats <- setNames(dep_src[[tcol]], dep_src$gene)
  stats <- stats[!is.na(stats) & is.finite(stats)]

  if (anyDuplicated(names(stats))) {
    dup_df <- tibble(gene = names(stats), t = unname(stats)) |>
      group_by(gene) |>
      slice_max(abs(t), n = 1, with_ties = FALSE) |>
      ungroup()
    stats <- setNames(dup_df$t, dup_df$gene)
  }
  stats <- sort(stats, decreasing = TRUE)

  res <- run_fgsea_perdb(
    ranks    = stats,
    pathways = pw_collection,
    nperm    = 10000,
    min_size = 10,
    max_size = 500
  )
  res$contrast <- ctr
  fgsea_all_list[[ctr]] <- res
}

fgsea_all <- bind_rows(fgsea_all_list) |>
  filter(contrast %in% ALL_CONTRASTS) |>
  select(pathway, database, contrast, NES, padj, size)

cat(sprintf("Hallmark + GO Slim fGSEA: %d rows, %d unique pathways, %d contrasts\n",
            nrow(fgsea_all), n_distinct(fgsea_all$pathway),
            n_distinct(fgsea_all$contrast)))

# --- Pivot to wide: rows = pathways, columns = contrasts ---
nes_wide <- fgsea_all |>
  pivot_wider(id_cols = c(pathway, database),
              names_from = contrast,
              values_from = c(NES, padj))

# Filter to pathways significant in at least 1 contrast (padj < 0.05)
padj_cols <- paste0("padj_", ALL_CONTRASTS)
padj_mat  <- as.matrix(nes_wide[, padj_cols])
any_sig   <- apply(padj_mat, 1, function(x) any(x < 0.05, na.rm = TRUE))
nes_sig   <- nes_wide[any_sig, ]

cat(sprintf("Pathways significant in >= 1 contrast: %d\n", nrow(nes_sig)))

# --- Build NES matrix for clustering ---
nes_cols <- paste0("NES_", ALL_CONTRASTS)
nes_mat  <- as.matrix(nes_sig[, nes_cols])
rownames(nes_mat) <- nes_sig$pathway

# Replace NA with 0 for clustering
nes_mat_clust <- nes_mat
nes_mat_clust[is.na(nes_mat_clust)] <- 0

# Hierarchical clustering on rows (pathways)
if (nrow(nes_mat_clust) > 2) {
  hc <- hclust(dist(nes_mat_clust), method = "ward.D2")
  row_order <- hc$order
} else {
  row_order <- seq_len(nrow(nes_mat_clust))
}

pathway_order <- nes_sig$pathway[row_order]

# --- Build long-form data for ggplot ---
plot_long <- nes_sig |>
  pivot_longer(cols = all_of(nes_cols),
               names_to = "contrast_nes", values_to = "NES_val") |>
  mutate(contrast = str_remove(contrast_nes, "^NES_"))

padj_long <- nes_sig |>
  select(pathway, all_of(padj_cols)) |>
  pivot_longer(cols = all_of(padj_cols),
               names_to = "contrast_padj", values_to = "padj_val") |>
  mutate(contrast = str_remove(contrast_padj, "^padj_"))

plot_df <- plot_long |>
  select(pathway, database, contrast, NES_val) |>
  left_join(padj_long |> select(pathway, contrast, padj_val),
            by = c("pathway", "contrast")) |>
  mutate(
    sig_label = case_when(
      is.na(padj_val) ~ "",
      padj_val < 0.001 ~ "***",
      padj_val < 0.01  ~ "**",
      padj_val < 0.05  ~ "*",
      TRUE ~ ""
    ),
    pathway = factor(pathway, levels = pathway_order),
    contrast = factor(contrast, levels = ALL_CONTRASTS),
    pathway_label = clean_pathway_name(as.character(pathway)),
    db = classify_database(as.character(pathway))
  )

# --- NES color limits (symmetric) ---
nes_lim <- max(abs(plot_df$NES_val), na.rm = TRUE)
nes_lim <- ceiling(nes_lim * 10) / 10

# --- Database annotation strip data ---
db_df <- plot_df |>
  distinct(pathway, db)

# --- Dynamic height ---
n_pathways <- length(pathway_order)
panel_h <- max(120, n_pathways * 3.5 + 40)
panel_w <- 280

cat(sprintf("Heatmap: %d pathways x %d contrasts, %d x %d mm\n",
            n_pathways, length(ALL_CONTRASTS), panel_w, round(panel_h)))

# --- Plot ---
pA <- ggplot(plot_df, aes(x = contrast, y = pathway, fill = NES_val)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sig_label), size = 3, color = "black", vjust = 0.75) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D6604D",
                       midpoint = 0, limits = c(-nes_lim, nes_lim),
                       name = "NES", na.value = "grey90") +
  scale_x_discrete(labels = CTR_SHORT[ALL_CONTRASTS], position = "top") +
  scale_y_discrete(labels = function(x) clean_pathway_name(x)) +
  labs(title = "Pathway Enrichment Across All Contrasts",
       subtitle = "Hallmark + GO Slim | fGSEA (per-database BH); * padj < 0.05",
       x = NULL, y = NULL, tag = "A") +
  FIG_THEME +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 0, size = 9, face = "bold"),
    axis.text.y = element_text(size = 7.5),
    legend.position = "right",
    panel.grid = element_blank(),
    panel.border = element_blank()
  )

# --- Database color strip (left annotation) ---
pA_db <- ggplot(db_df, aes(x = 1, y = pathway, fill = db)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = DB_COLORS, name = "Database") +
  scale_y_discrete(labels = NULL) +
  labs(x = NULL, y = NULL) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 9),
    legend.text  = element_text(size = 8),
    legend.key.size = unit(3, "mm"),
    plot.margin = margin(t = 5, r = 0, b = 5, l = 2)
  )

# Combine with patchwork
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  pA_combined <- pA_db + pA +
    plot_layout(widths = c(1, 30), guides = "collect") &
    theme(legend.position = "right")
} else {
  pA_combined <- pA
}

# --- Save ---
ggsave(file.path(RPT, "panel_A_nes_heatmap.pdf"), pA_combined,
       width = panel_w, height = panel_h, units = "mm", device = pdf_device,
       limitsize = FALSE)
ggsave(file.path(RPT, "panel_A_nes_heatmap.png"), pA_combined,
       width = panel_w, height = panel_h, units = "mm", dpi = 300,
       limitsize = FALSE)

# --- Export ---
export_df <- plot_df |>
  select(pathway, database = db, contrast, NES = NES_val,
         padj = padj_val, sig_label) |>
  arrange(pathway, contrast)

write.csv(as.data.frame(export_df),
          file.path(DAT, "nes_heatmap_data.csv"), row.names = FALSE)

cat(sprintf("Panel A (NES heatmap) done: %d pathways.\n", n_pathways))
