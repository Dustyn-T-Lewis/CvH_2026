# Figure 2 — Panel D (CRvH): fGSEA Grouped Bar Chart (Pathway Enrichment)
# Full collection: Hallmark + C2:CP (all subcollections) + GO:BP
# Per-database BH correction + Jaccard dedup.
# Reads CRvH result file for 2 contrasts.
# Outputs: pD_crvh (ggplot object)

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
})

# --- Config ---
DEP_FILE <- "03_DEP/c_data/03_combined_results_CRvH.csv"
RPT      <- "04_Figures/F02/b_reports"
DAT      <- "04_Figures/F02/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

CONTRASTS <- c("Cancer_vs_Healthy", "Training_CR")

pdf_device <- get_pdf_device()
PH_W <- 130

# --- Load DEP results ---
dep_df <- read_csv(DEP_FILE, show_col_types = FALSE) |>
  select(uniprot_id, gene, t_Cancer_vs_Healthy, t_Training_CR)

cat(sprintf("DEP (CRvH): %d proteins, %d contrasts\n", nrow(dep_df), length(CONTRASTS)))

# --- Build pathway collection ---
pw_collection <- build_pathway_collection(min_size = 10, max_size = 500)

# --- Run fGSEA per contrast ---
set.seed(42)
fgsea_perdb_all <- list()

for (ctr in CONTRASTS) {
  cat(sprintf("\n--- fGSEA: %s ---\n", ctr))
  tcol <- paste0("t_", ctr)
  stats <- setNames(dep_df[[tcol]], dep_df$gene)
  stats <- stats[!is.na(stats) & is.finite(stats)]

  # Deduplicate gene names: keep entry with highest |t|
  if (anyDuplicated(names(stats))) {
    dup_df <- tibble(gene = names(stats), t = unname(stats)) |>
      group_by(gene) |>
      slice_max(abs(t), n = 1, with_ties = FALSE) |>
      ungroup()
    stats <- setNames(dup_df$t, dup_df$gene)
    cat(sprintf("  Deduped genes: %d unique\n", length(stats)))
  }

  stats <- sort(stats, decreasing = TRUE)

  res <- run_fgsea_perdb(
    ranks          = stats,
    pathways       = pw_collection,
    jaccard_cutoff = 0.5,
    nperm          = 10000,
    min_size       = 10,
    max_size       = 500
  )
  res$contrast <- ctr
  fgsea_perdb_all[[ctr]] <- res
}

fgsea_combined <- bind_rows(fgsea_perdb_all)

# --- Export full results ---
fgsea_export <- fgsea_combined |>
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";")) |>
  arrange(database, contrast, padj)
write_csv(fgsea_export, file.path(DAT, "06_panel_H_CRvH_fgsea_results.csv"))
cat(sprintf("\nExported fGSEA results: %d rows\n", nrow(fgsea_export)))

# --- Build grouped bar chart ---
DISPLAY_DBS <- c("Hallmark", "KEGG", "Reactome", "GO:BP")

db_totals <- fgsea_combined |>
  filter(database %in% DISPLAY_DBS) |>
  distinct(pathway, database) |>
  count(database, name = "n_total")

count_df <- fgsea_combined |>
  filter(!is.na(padj), padj < 0.05, database %in% DISPLAY_DBS) |>
  group_by(contrast, database) |>
  summarise(
    Up   = sum(NES > 0),
    Down = sum(NES < 0),
    .groups = "drop"
  ) |>
  pivot_longer(cols = c(Up, Down), names_to = "direction",
               values_to = "count") |>
  left_join(db_totals, by = "database") |>
  mutate(fraction = count / n_total)

nonempty_dbs <- count_df |>
  group_by(database) |>
  filter(sum(count) > 0) |>
  pull(database) |>
  unique()
count_df <- count_df |> filter(database %in% nonempty_dbs)

db_labels <- setNames(
  sprintf("%s\n(n=%d)", db_totals$database, db_totals$n_total),
  db_totals$database
)[nonempty_dbs]

count_df$contrast  <- factor(count_df$contrast, levels = CONTRASTS)
count_df$database  <- factor(count_df$database, levels = intersect(DISPLAY_DBS, nonempty_dbs))
count_df$direction <- factor(count_df$direction, levels = c("Up", "Down"))

n_facets <- length(levels(count_df$database))
PH_H <- max(80, n_facets * 55)

lbl_sz <- scale_text(BASE_COUNT, PH_W)

# Background rectangles for 2 contrasts
bg_rects <- lapply(seq_along(CONTRASTS), function(i) {
  annotate("rect",
           xmin = i - 0.5, xmax = i + 0.5,
           ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS[CONTRASTS[i]],
           alpha = 0.20, color = "grey70", linewidth = 0.2)
})

pD_crvh <- ggplot(count_df, aes(x = contrast, y = fraction * 100, fill = direction)) +
  bg_rects +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.3) +
  geom_text(aes(y = fraction * 100 / 2,
                label = ifelse(count > 0, count, "")),
            position = position_dodge(width = 0.7),
            vjust = 0.5, hjust = 0.5, size = lbl_sz,
            color = "white", fontface = "bold", show.legend = FALSE) +
  facet_grid(database ~ ., scales = "free_y",
             labeller = as_labeller(db_labels)) +
  scale_x_discrete(labels = CTR_SHORT[CONTRASTS]) +
  scale_fill_manual(values = DIR_COLORS) +
  labs(title = "Pathway Enrichment (CRvH)",
       subtitle = "fGSEA (padj < 0.05)\nper-database BH correction",
       x = NULL, y = "% of database significant",
       tag = "D1") +
  FIG_THEME +
  theme(axis.text.x    = element_text(angle = 35, hjust = 1,
                                      size = FIG_AXIS_TEXT - 0.5),
        legend.position = "none",
        strip.text.y   = element_text(size = FIG_STRIP_SIZE, face = "bold", angle = 0))

# --- Audit summary ---
nes_summary <- fgsea_combined |>
  filter(!is.na(padj), padj < 0.05) |>
  group_by(contrast, database) |>
  summarise(
    n_sig = n(), n_up = sum(NES > 0), n_down = sum(NES < 0),
    median_NES = median(NES), mean_NES = mean(NES), sd_NES = sd(NES),
    min_padj = min(padj), median_padj = median(padj),
    .groups = "drop"
  )
write.csv(nes_summary, file.path(DAT, "panel_D_CRvH_nes_summary.csv"),
          row.names = FALSE)

# --- Save ---
ggsave(file.path(RPT, "panel_D_CRvH_fgsea.pdf"), pD_crvh,
       width = PH_W, height = PH_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_D_CRvH_fgsea.png"), pD_crvh,
       width = PH_W, height = PH_H, units = "mm", dpi = 300)

cat(sprintf("\nPanel D (CRvH) done: %d x %d mm, %d facets\n", PH_W, PH_H, n_facets))
