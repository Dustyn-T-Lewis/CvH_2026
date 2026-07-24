#!/usr/bin/env Rscript
# F03_Enrich_Volcanoes driver: per-contrast volcano-in-ring panels plus the
# three-ring composite, drawn by enrichVolcano from the reported DE (pi-score
# volcano) and the fgsea cache (NES arcs, deduplicated within database at cache
# build). The trio is the cancer / training / residual axis. Rings save to
# b_reports/main/*/panels; the composite to b_reports/main.

setwd(here::here())
source("04_Figures/shared/style.R")
source("04_Figures/shared/figure_supplement_helpers.R")
pacman::p_load(dplyr, readr, enrichVolcano, ggplot2, patchwork, ggtext)

BASE <- "04_Figures/F03_Enrich_Volcanoes"
RPT_PNG <- file.path(BASE, "b_reports/main/png")
RPT_PDF <- file.path(BASE, "b_reports/main/pdf")
PAN_PNG <- file.path(RPT_PNG, "panels")
PAN_PDF <- file.path(RPT_PDF, "panels")
DAT <- file.path(BASE, "c_data")
for (d in c(PAN_PNG, PAN_PDF, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

PI_THRESH <- 0.05
RING_N <- 12

contrasts <- tibble::tribble(
  ~ctr, ~role, ~tag, ~math,
  "CRvH_Baseline", "Cancer vs Healthy", "cancer", "Cancer − Healthy (baseline)",
  "CR_Training", "Training (CR)", "training", "CR post − CR pre",
  "Resid", "Residual", "residual", "CR post − Healthy (residual)"
)

dep <- read_csv(
  "03_DEP/a_non_imputed/c_data/combined_results_pi.csv",
  show_col_types = FALSE
)
fgsea_all <- read_csv("04_Figures/shared/fgsea_CRvH.csv", show_col_types = FALSE) |>
  filter(contrast %in% contrasts$ctr)

# Per contrast: a gene-level volcano (pi-score = significance, one point per
# gene) and the significant fgsea pathways (already deduplicated within database
# in the cache), capped at the top RING_N by FDR for a ring of even quarters.
prep <- lapply(seq_len(nrow(contrasts)), function(i) {
  ctr <- contrasts$ctr[i]
  volc <- dep |>
    filter(contrast == ctr, !is.na(logFC), !is.na(P.Value), !is.na(gene), gene != "") |>
    group_by(gene) |>
    slice_min(pi_score, n = 1, with_ties = FALSE) |>
    ungroup() |>
    transmute(gene, logFC, P.Value, padj = pi_score) |>
    as.data.frame()
  sig <- fgsea_all |>
    filter(contrast == ctr, !is.na(padj), padj < 0.05, size >= 15, size <= 500) |>
    arrange(padj)
  list(
    volc = volc,
    enrich = as.data.frame(slice_head(sig, n = RING_N)),
    n_dep = sum(volc$padj < PI_THRESH, na.rm = TRUE),
    n_sig = nrow(sig)
  )
})
names(prep) <- contrasts$role

# Volcano significance is the pi-score alone (it already embeds |logFC|); no
# extra fold gate, so the coloured points match the subtitle DEP count.
ring_args <- list(
  gene_col = "gene", logfc_col = "logFC", pval_col = "P.Value", padj_col = "padj",
  term_col = "pathway", nes_col = "NES", size_col = "size",
  genes_col = "leadingEdge", genes_sep = ";",
  p_threshold = PI_THRESH, logfc_threshold = 0,
  ring_radius = 5.5, volcano_radius = 5.2, arc_height_range = c(0.1, 3.2),
  label_size = 2.7, point_size = 2.0, point_alpha = 0.7,
  count_x_mult = 0.55, count_y_mult = 0.55
)

for (i in seq_len(nrow(contrasts))) {
  sub <- sprintf(
    "%s | %d DEPs, %d pathways",
    contrasts$math[i], prep[[i]]$n_dep, nrow(prep[[i]]$enrich)
  )
  p <- do.call(volcano_ring, c(
    list(prep[[i]]$volc, prep[[i]]$enrich,
      title = contrasts$role[i], subtitle = sub
    ),
    ring_args
  ))
  ggsave(file.path(PAN_PNG, sprintf("MAIN_F03_%s_ring.png", contrasts$tag[i])),
    p,
    width = 150, height = 150, units = "mm", dpi = 300
  )
  ggsave(file.path(PAN_PDF, sprintf("MAIN_F03_%s_ring.pdf", contrasts$tag[i])),
    p,
    width = 150, height = 150, units = "mm", device = pdf_device
  )
}

# Three-ring composite. The "shown / significant" ratio in each subtitle is the
# RING_N cap biting against the post-dedup total, so a reader can tell when a
# ring truncates vs. shows everything. Dedup rule + FDR ranges live in the caption.
composite_subtitle <- vapply(seq_len(nrow(contrasts)), function(i) {
  e <- prep[[i]]$enrich
  sprintf(
    "**%s**<br>%d DEPs · %d/%d pathways shown (%d↑ %d↓)",
    contrasts$math[i], prep[[i]]$n_dep, nrow(e), prep[[i]]$n_sig,
    sum(e$NES > 0), sum(e$NES < 0)
  )
}, character(1))

grid <- volcano_ring_grid(
  volc_dfs = lapply(prep, `[[`, "volc"),
  enrich_dfs = lapply(prep, `[[`, "enrich"),
  contrasts = contrasts$role,
  subtitles = composite_subtitle,
  gene_col = "gene", logfc_col = "logFC", pval_col = "P.Value", padj_col = "padj",
  term_col = "pathway", nes_col = "NES", size_col = "size",
  genes_col = "leadingEdge", genes_sep = ";",
  p_threshold = PI_THRESH, logfc_threshold = 0,
  x_scale = 0.9, y_scale = 0.93,
  label_size = 2.9, count_size = 2.8, ncol = 3,
  theme = volcano_ring_theme(base_size = 13)
)
fig <- grid$plot &
  theme(
    plot.subtitle = ggtext::element_markdown(
      hjust = 0.5, halign = 0.5, colour = "grey30",
      size = rel(0.7), lineheight = 1.15
    ),
    legend.position = "bottom",
    legend.key.width = unit(13, "mm"),
    legend.key.height = unit(2.5, "mm"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )
fig <- fig & guides(
  fill = guide_colorbar(direction = "horizontal", title.position = "top")
)
fig <- fig + plot_annotation(
  caption = paste0(
    "Volcano points: π-score (Xiao 2014) < 0.05 DEPs.  ",
    "Ring arcs: fgsea BH-FDR < 0.05, Jaccard-deduplicated within database; ",
    "arc height = −log10 FDR, colour = NES."
  ),
  theme = theme(
    plot.caption = element_text(size = 8, hjust = 0.5, colour = "grey35")
  )
)
ggsave(file.path(RPT_PNG, "MAIN_F03_enrich_volcanoes.png"), fig,
  width = 330, height = 175, units = "mm", dpi = 300, bg = "white", limitsize = FALSE
)
ggsave(file.path(RPT_PDF, "MAIN_F03_enrich_volcanoes.pdf"), fig,
  width = 330, height = 175, units = "mm", device = pdf_device,
  bg = "white", limitsize = FALSE
)

# Supplementary workbook: contrast key + every tested pathway per contrast.
ctr_key <- contrasts |>
  transmute(contrast = ctr, role, file_tag = tag, definition = math) |>
  as.data.frame()
pathway_sheets <- lapply(seq_len(nrow(contrasts)), function(i) {
  df <- fgsea_all |>
    filter(contrast == contrasts$ctr[i]) |>
    transmute(pathway, database, padj, pval, NES, size,
      shown = pathway %in% prep[[i]]$enrich$pathway
    ) |>
    arrange(padj) |>
    as.data.frame()
  list(name = contrasts$role[i], df = df)
})
overview_df <- data.frame(
  Sheet = c("contrast_map", contrasts$role),
  Description = c(
    "Contrast key: name, role, file tag, algebra",
    sprintf(
      "All fgsea pathways tested for %s (shown = drawn on the ring); sorted by BH-FDR",
      contrasts$role
    )
  )
)
build_workbook(
  file.path(DAT, "F03_supplementary.xlsx"),
  title = "F03: Enrichment volcano-in-ring panels",
  description = paste(
    "Per-contrast volcano-in-ring panels for the cancer / training / residual axis.",
    "Volcano significance is the Xiao pi-score (< 0.05); ring arcs are fgsea",
    "BH-FDR < 0.05 pathways after within-database Jaccard deduplication, top 12 by",
    "FDR, arc height = -log10 FDR, colour = NES."
  ),
  overview_df = overview_df,
  sheet_specs = c(list(list(name = "contrast_map", df = ctr_key)), pathway_sheets)
)

message(sprintf("F03: %d ring panels + composite saved", nrow(contrasts)))
