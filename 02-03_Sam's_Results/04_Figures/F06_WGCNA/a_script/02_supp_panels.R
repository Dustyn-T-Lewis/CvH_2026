# F06_WGCNA — Supplementary Panels (YvO-aligned, inline, no per-panel helpers)
#
# Supp A: Soft-power threshold (ggplot2, scale-independence + mean-connectivity)
# Supp B: Dendrogram + module color bar (WGCNA base-R via png, wrapped in ggplot)
# Supp C: Compartment GO enrichment per module (Fisher exact, BH)
# Supp D: Bicor vs Pearson sensitivity (Jaccard overlap heatmap)
# Supp E: Per-module triptych heatmaps (z-score | eigengene | top-5 ORA)
# Supp F: Per-module hub networks (igraph + TOM, x7 modules)
# Supp G: Inter-module hub chord diagram (circlize)
# Supp H: Module preservation stats (permutation-based quality, n=100 perms)
#
# Dropped vs old Sam 02_supp_panels.R:
#   - Module size bar chart (not in YvO)
#   - Module-trait simple Pearson heatmap (not in YvO)
#
# Run from A_CvH_2026/ root after 00_run_wgcna.R.

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(cowplot)
  library(png)
  library(grid)
  library(ggrepel)
  library(WGCNA)
  library(igraph)
  library(ggraph)
  library(tidygraph)
  library(ggforce)
  library(concaveman)
  library(graphlayouts)
  library(ggnewscale)
  library(circlize)
  library(fgsea)
})

source("02-03_Sam's_Results/04_Figures/shared/style.R")
source("02-03_Sam's_Results/04_Figures/shared/pathway_utils.R")

pdf_device <- get_pdf_device()

allowWGCNAThreads()
set.seed(42)

BASE      <- "02-03_Sam's_Results/04_Figures/F06_WGCNA"
PANEL_DIR <- file.path(BASE, "c_data")
WGCNA_DIR <- file.path(BASE, "c_data", "wgcna")
SUPP_PNG  <- file.path(BASE, "b_reports", "supp", "png", "panels")
SUPP_PDF  <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
MOD_PNG   <- file.path(BASE, "b_reports", "supp", "png", "modules")
MOD_PDF   <- file.path(BASE, "b_reports", "supp", "pdf", "modules")
for (d in c(SUPP_PNG, SUPP_PDF, MOD_PNG, MOD_PDF,
            file.path(PANEL_DIR, "supp"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ── Load shared objects ───────────────────────────────────────────────────────

mod_bio_labels  <- read_csv(file.path(PANEL_DIR, "mod_bio_labels.csv"),
                            show_col_types = FALSE)
module_df       <- read_csv(file.path(WGCNA_DIR, "wgcna_module_assignments.csv"),
                            show_col_types = FALSE)
hub_df_raw      <- read_csv(file.path(WGCNA_DIR, "wgcna_hub_proteins.csv"),
                            show_col_types = FALSE)
net             <- readRDS(file.path(WGCNA_DIR, "wgcna_network.rds"))
module_colors   <- readRDS(file.path(PANEL_DIR, "module_colors.rds"))
MEs             <- readRDS(file.path(PANEL_DIR, "MEs.rds"))
kME_all         <- readRDS(file.path(PANEL_DIR, "kME_all.rds"))
datExpr         <- readRDS(file.path(PANEL_DIR, "datExpr.rds"))
group_z         <- readRDS(file.path(PANEL_DIR, "group_z.rds"))
meta            <- read_csv(file.path(PANEL_DIR, "meta.csv"), show_col_types = FALSE)
ann             <- read_csv(file.path(PANEL_DIR, "imp_annotations.csv"), show_col_types = FALSE)
sft_csv         <- read_csv(file.path(WGCNA_DIR, "wgcna_sft_summary.csv"),
                            show_col_types = FALSE)
sft_fi          <- readRDS(file.path(WGCNA_DIR, "sft_fitIndices.rds"))
lmm_df          <- read_csv(file.path(WGCNA_DIR, "wgcna_lmm_contrasts.csv"),
                            show_col_types = FALSE)
enrich_df       <- read_csv(file.path(WGCNA_DIR, "wgcna_module_enrichment.csv"),
                            show_col_types = FALSE)
key_mods        <- readLines(file.path(WGCNA_DIR, "key_modules.txt")) |>
  trimws() |> (\(x) x[nzchar(x)])()

ALL_MODULES <- module_df |>
  filter(module_color != "grey") |>
  count(module_color, sort = TRUE) |>
  pull(module_color)

mod_display_vec <- setNames(mod_bio_labels$display_label,
                             mod_bio_labels$module_color)

soft_power <- sft_csv$selected_power[1]
n_proteins_total <- sft_csv$n_proteins[1]

uid2gene <- setNames(module_df$gene, module_df$uniprot_id)

# ── Supp A: Soft-power threshold (ggplot2) ────────────────────────────────────
message("Supp A: soft-power threshold...")

fit_df <- tibble(
  power    = sft_fi$Power,
  r2       = -sign(sft_fi$slope) * sft_fi$SFT.R.sq,
  mean_k   = sft_fi$mean.k.,
  slope    = sft_fi$slope
) |>
  mutate(selected = power == soft_power)

PA_W <- 240; PA_H <- 110
txt_label <- scale_text(BASE_GENE, PA_W) * 0.9

p_sft1 <- ggplot(fit_df, aes(power, r2)) +
  geom_hline(yintercept = 0.85, linetype = "dashed", color = "grey40",
             linewidth = 0.4) +
  geom_point(aes(fill = selected), shape = 21, size = 2.5, stroke = 0.4,
             color = "black") +
  geom_text_repel(aes(label = power), size = txt_label,
                  color = "grey25", fontface = "bold",
                  max.overlaps = 20, seed = 42,
                  min.segment.length = 0.3, segment.size = 0.25,
                  segment.color = "grey60") +
  scale_fill_manual(values = c("TRUE" = "#D6604D", "FALSE" = "grey70"),
                    guide = "none") +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  annotate("text", x = 19, y = 0.88,
           label = expression(R^2 == 0.85~threshold),
           size = txt_label * 0.85, color = "grey40", hjust = 1) +
  labs(x = "Soft Threshold (power)",
       y = expression(Scale~Free~Topology~R^2),
       title = "Scale Independence") +
  FIG_THEME +
  theme(plot.title = element_text(size = 10, face = "bold"))

p_sft2 <- ggplot(fit_df, aes(power, mean_k)) +
  geom_point(aes(fill = selected), shape = 21, size = 2.5, stroke = 0.4,
             color = "black") +
  geom_text_repel(aes(label = power), size = txt_label,
                  color = "grey25", fontface = "bold",
                  max.overlaps = 20, seed = 42,
                  min.segment.length = 0.3, segment.size = 0.25,
                  segment.color = "grey60") +
  scale_fill_manual(values = c("TRUE" = "#D6604D", "FALSE" = "grey70"),
                    guide = "none") +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  labs(x = "Soft Threshold (power)",
       y = "Mean Connectivity",
       title = "Mean Connectivity") +
  FIG_THEME +
  theme(plot.title = element_text(size = 10, face = "bold"))

p_sft <- (p_sft1 | p_sft2) +
  plot_annotation(
    title = "Scale-Free Topology Fit",
    subtitle = sprintf("Signed network | selected power = %d (R² = %.3f) | %s proteins",
                       soft_power, fit_df$r2[fit_df$power == soft_power],
                       format(n_proteins_total, big.mark = ",")),
    theme = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(face = "bold.italic", size = 10,
                                    color = "grey30"),
      plot.margin   = margin(4, 4, 4, 4)
    )
  )

ggsave(file.path(SUPP_PNG, "SUPP_soft_threshold.png"), p_sft,
       width = PA_W, height = PA_H, units = "mm", dpi = 300)
ggsave(file.path(SUPP_PDF, "SUPP_soft_threshold.pdf"), p_sft,
       width = PA_W, height = PA_H, units = "mm", device = pdf_device)

write_csv(fit_df, file.path(PANEL_DIR, "supp", "a01_sft_fit_indices.csv"))
message("  Supp A saved")

# ── Supp B: Dendrogram + module color bar ─────────────────────────────────────
message("Supp B: dendrogram...")

PB_W <- 240; PB_H <- 120

block_genes   <- net$blockGenes[[1]]
merged_cols   <- module_colors[block_genes]
unmerged_cols <- net$unmergedColors[block_genes]
color_matrix  <- cbind(unmerged_cols, merged_cols)
color_labels  <- c("Dynamic Tree Cut", "Merged Modules")

n_mods  <- length(unique(merged_cols[merged_cols != "grey"]))
n_genes <- length(merged_cols)
n_grey  <- sum(merged_cols == "grey")

dendro_tmp <- tempfile(fileext = ".png")
tryCatch({
  png(dendro_tmp, width = 3200, height = 1600, res = 300)
  par(mar = c(1, 4, 1, 0.5))
  plotDendroAndColors(net$dendrograms[[1]],
                      color_matrix, color_labels,
                      main = "",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      cex.colorLabels = 0.7, cex.axis = 0.8)
  dev.off()
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  message("Dendrogram render failed: ", e$message)
})

dendro_img <- if (file.exists(dendro_tmp) && file.size(dendro_tmp) > 0)
  readPNG(dendro_tmp) else NULL

sub_txt <- sprintf(
  "Signed network | power = %d | %d modules | %s proteins (%d unassigned)",
  soft_power, n_mods, format(n_genes, big.mark = ","), n_grey
)

p_dendro <- ggplot() +
  { if (!is.null(dendro_img))
      annotation_raster(dendro_img, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
    else
      annotate("text", x = 0.5, y = 0.5, label = "Dendrogram unavailable",
               size = 5, color = "grey50")
  } +
  labs(title = "Protein Dendrogram & Module Colors", subtitle = sub_txt) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void() +
  theme(plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 10, color = "grey30", face = "italic"),
        plot.margin   = margin(2, 2, 2, 2))

write_csv(
  tibble(uniprot_id = names(module_colors)[block_genes],
         unmerged_color = unmerged_cols, merged_color = merged_cols),
  file.path(PANEL_DIR, "supp", "a02_dendrogram_data.csv"))

ggsave(file.path(SUPP_PNG, "SUPP_dendrogram.png"), p_dendro,
       width = PB_W, height = PB_H, units = "mm", dpi = 300)
ggsave(file.path(SUPP_PDF, "SUPP_dendrogram.pdf"), p_dendro,
       width = PB_W, height = PB_H, units = "mm", device = pdf_device)
message("  Supp B saved")

# ── Supp C: Compartment GO enrichment per module ──────────────────────────────
message("Supp C: compartment enrichment...")

hpa_file <- "00_input/HPA_skeletal_muscle_annotations.tsv"
stopifnot("HPA file missing" = file.exists(hpa_file))
hpa <- read.delim(hpa_file, stringsAsFactors = FALSE)

hpa_sub <- hpa |>
  filter(Subcellular.main.location != "", Gene %in% module_df$gene) |>
  select(gene = Gene, location = Subcellular.main.location)

hpa_long <- hpa_sub |>
  separate_rows(location, sep = ",\\s*") |>
  filter(location != "")

COMPARTMENTS <- c("Mitochondria", "Cytosol", "Nucleoplasm",
                  "Plasma membrane", "Endoplasmic reticulum",
                  "Golgi apparatus", "Vesicles", "Cytoskeleton")

hpa_long <- hpa_long |>
  mutate(compartment = case_when(
    location %in% c("Microtubules", "Actin filaments",
                     "Intermediate filaments", "Cytoskeleton") ~ "Cytoskeleton",
    location %in% COMPARTMENTS ~ location,
    TRUE ~ NA_character_
  )) |>
  filter(!is.na(compartment))

compartment_sets <- hpa_long |>
  distinct(gene, compartment) |>
  group_by(compartment) |>
  summarise(genes = list(gene), n = n(), .groups = "drop")

universe <- unique(hpa_long$gene)
n_universe <- length(universe)
message(sprintf("  Compartment: %d compartments, universe %d genes",
                nrow(compartment_sets), n_universe))

results_comp <- list()
for (mc in ALL_MODULES) {
  mod_genes <- module_df$gene[module_df$module_color == mc]
  mod_genes_in_univ <- intersect(mod_genes, universe)
  n_mod <- length(mod_genes_in_univ)
  for (i in seq_len(nrow(compartment_sets))) {
    comp      <- compartment_sets$compartment[i]
    comp_genes <- compartment_sets$genes[[i]]
    in_both   <- length(intersect(mod_genes_in_univ, comp_genes))
    in_mod_only <- n_mod - in_both
    in_comp_only <- length(comp_genes) - in_both
    in_neither   <- n_universe - in_both - in_mod_only - in_comp_only
    ct <- matrix(c(in_both, in_mod_only, in_comp_only, in_neither), nrow = 2)
    ft <- fisher.test(ct, alternative = "greater")
    results_comp[[length(results_comp) + 1]] <- tibble(
      module = mc, compartment = comp,
      odds_ratio = ft$estimate, p_raw = ft$p.value,
      n_overlap = in_both, n_module = n_mod,
      n_compartment = length(comp_genes), n_universe = n_universe
    )
  }
}

enrich_comp <- bind_rows(results_comp) |>
  mutate(p_bh = p.adjust(p_raw, method = "BH")) |>
  arrange(p_bh)

write_csv(enrich_comp, file.path(PANEL_DIR, "supp", "a03_compartment_enrichment.csv"))
message(sprintf("  Compartment significant (FDR<0.05): %d / %d",
                sum(enrich_comp$p_bh < 0.05), nrow(enrich_comp)))

label_map <- setNames(mod_bio_labels$display_label, mod_bio_labels$module_color)
enrich_comp <- enrich_comp |>
  mutate(
    module_label  = ifelse(module %in% names(label_map), label_map[module], module),
    neg_log10_p   = -log10(pmax(p_bh, 1e-10)),
    sig_label     = case_when(
      p_bh < 0.001 ~ "***",
      p_bh < 0.01  ~ "**",
      p_bh < 0.05  ~ "*",
      TRUE          ~ ""
    )
  )

mod_ord_comp <- enrich_comp |>
  group_by(module_label) |>
  summarise(total_sig = -sum(log10(pmax(p_bh, 1e-10))), .groups = "drop") |>
  arrange(desc(total_sig)) |>
  pull(module_label)

enrich_comp <- enrich_comp |>
  mutate(
    module_label = factor(module_label, levels = rev(mod_ord_comp)),
    compartment  = factor(compartment, levels = COMPARTMENTS)
  )

p_comp <- ggplot(enrich_comp, aes(compartment, module_label, fill = neg_log10_p)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sig_label), size = 3.5, vjust = 0.75) +
  scale_fill_gradient(low = "grey95", high = "#1B5E20",
                      limits = c(0, max(enrich_comp$neg_log10_p, na.rm = TRUE)),
                      name = expression(-log[10](p[BH]))) +
  labs(x = NULL, y = NULL,
       title = "Module Subcellular Compartment Enrichment",
       subtitle = sprintf("One-sided Fisher exact, BH-corrected (%d tests)  |  * p<0.05  ** p<0.01  *** p<0.001",
                          nrow(enrich_comp))) +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8),
        axis.text.y = element_text(size = 7.5),
        plot.title    = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 9))

ggsave(file.path(SUPP_PNG, "SUPP_compartment_enrichment.png"), p_comp,
       width = 180, height = 130, units = "mm", dpi = 300)
ggsave(file.path(SUPP_PDF, "SUPP_compartment_enrichment.pdf"), p_comp,
       width = 180, height = 130, units = "mm", device = pdf_device)
message("  Supp C saved")

# ── Supp D: Bicor vs Pearson sensitivity ──────────────────────────────────────
message("Supp D: bicor sensitivity (rerunning blockwiseModules with bicor)...")

BICOR_RDS <- file.path(PANEL_DIR, "wgcna_bicor_result.rds")
PEARSON_MODS <- file.path(WGCNA_DIR, "wgcna_module_assignments.csv")

# Remove any columns with NAs before bicor (bicor requires complete data)
datExpr_bicor <- datExpr[, colSums(is.na(datExpr)) == 0]
message(sprintf("  Bicor datExpr: %d samples x %d proteins (complete columns only)",
                nrow(datExpr_bicor), ncol(datExpr_bicor)))

if (file.exists(BICOR_RDS)) {
  message("  Loading cached bicor result...")
  net_bicor <- readRDS(BICOR_RDS)
} else {
  t_bicor_start <- proc.time()

  powers <- 1:20
  sft_bicor <- pickSoftThreshold(datExpr_bicor, powerVector = powers,
                                  networkType = "signed",
                                  corFnc = "bicor",
                                  corOptions = list(maxPOutliers = 0.1,
                                                    use = "pairwise.complete.obs"),
                                  verbose = 2)
  r2_bicor  <- -sign(sft_bicor$fitIndices$slope) * sft_bicor$fitIndices$SFT.R.sq
  power_idx_b <- which(r2_bicor > 0.85)[1]
  bicor_power <- if (!is.na(power_idx_b)) powers[power_idx_b] else 6L
  message(sprintf("  Bicor soft power: %d (R^2 = %.3f)", bicor_power, r2_bicor[bicor_power]))

  cor <- WGCNA::cor   # ensure WGCNA cor
  net_bicor <- blockwiseModules(
    datExpr_bicor,
    power             = bicor_power,
    networkType       = "signed",
    TOMType           = "signed",
    corType           = "bicor",
    maxPOutliers      = 0.1,
    minModuleSize     = 30,
    mergeCutHeight    = 0.25,
    numericLabels     = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs          = FALSE,
    verbose           = 3
  )
  cor <- stats::cor   # restore

  t_elapsed <- (proc.time() - t_bicor_start)[3]
  message(sprintf("  Bicor blockwiseModules done (%.1f min)", t_elapsed / 60))
  saveRDS(net_bicor, BICOR_RDS)
}

bicor_colors <- labels2colors(net_bicor$colors)
n_bicor <- length(unique(bicor_colors[bicor_colors != "grey"]))
message(sprintf("  Bicor modules: %d (+ grey)", n_bicor))

pearson_df <- read_csv(PEARSON_MODS, show_col_types = FALSE)
bicor_df   <- tibble(uniprot_id = colnames(datExpr_bicor), bicor_module = bicor_colors)

compare <- pearson_df |>
  select(uniprot_id, pearson_module = module_color) |>
  inner_join(bicor_df, by = "uniprot_id")

pearson_mods <- setdiff(unique(compare$pearson_module), "grey")
bicor_mods   <- setdiff(unique(compare$bicor_module),   "grey")

jaccard_mat <- matrix(0, nrow = length(pearson_mods), ncol = length(bicor_mods),
                      dimnames = list(pearson_mods, bicor_mods))
for (pm in pearson_mods) {
  p_set <- compare$uniprot_id[compare$pearson_module == pm]
  for (bm in bicor_mods) {
    b_set <- compare$uniprot_id[compare$bicor_module == bm]
    inter <- length(intersect(p_set, b_set))
    uni   <- length(union(p_set, b_set))
    jaccard_mat[pm, bm] <- if (uni > 0) inter / uni else 0
  }
}

# Greedy best-match pairing
matched <- tibble(pearson_module = character(), bicor_module = character(),
                  jaccard = numeric(), n_shared = integer(),
                  n_pearson = integer(), n_bicor = integer())
used_p <- character(); used_b <- character()
repeat {
  rem <- jaccard_mat[!rownames(jaccard_mat) %in% used_p,
                     !colnames(jaccard_mat) %in% used_b, drop = FALSE]
  if (nrow(rem) == 0 || ncol(rem) == 0) break
  best_idx <- which(rem == max(rem), arr.ind = TRUE)[1, ]
  best_pm  <- rownames(rem)[best_idx[1]]
  best_bm  <- colnames(rem)[best_idx[2]]
  p_set    <- compare$uniprot_id[compare$pearson_module == best_pm]
  b_set    <- compare$uniprot_id[compare$bicor_module   == best_bm]
  matched  <- bind_rows(matched, tibble(
    pearson_module = best_pm, bicor_module = best_bm,
    jaccard  = rem[best_idx[1], best_idx[2]],
    n_shared = length(intersect(p_set, b_set)),
    n_pearson = length(p_set), n_bicor = length(b_set)
  ))
  used_p <- c(used_p, best_pm); used_b <- c(used_b, best_bm)
}

write_csv(matched, file.path(PANEL_DIR, "supp", "a04_bicor_sensitivity.csv"))
message(sprintf("  Matched modules: %d | Mean Jaccard: %.3f | Jaccard > 0.5: %d/%d",
                nrow(matched), mean(matched$jaccard),
                sum(matched$jaccard > 0.5), nrow(matched)))

jac_long <- as_tibble(jaccard_mat, rownames = "pearson_module") |>
  pivot_longer(-pearson_module, names_to = "bicor_module", values_to = "jaccard")

pm_order <- c(matched$pearson_module, setdiff(pearson_mods, matched$pearson_module))
bm_order <- c(matched$bicor_module,   setdiff(bicor_mods,   matched$bicor_module))

jac_long <- jac_long |>
  mutate(pearson_module = factor(pearson_module, levels = rev(pm_order)),
         bicor_module   = factor(bicor_module,   levels = bm_order))

p_bicor <- ggplot(jac_long, aes(bicor_module, pearson_module, fill = jaccard)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.2f", jaccard)),
            size = 2.8, colour = ifelse(jac_long$jaccard > 0.5, "white", "grey30")) +
  scale_fill_gradient2(low = "white", mid = "#B2DFDB", high = "#00695C",
                       midpoint = 0.35, limits = c(0, 1),
                       name = "Jaccard\nIndex") +
  labs(x = "Bicor module", y = "Pearson module",
       title = "Pearson vs Bicor Module Overlap",
       subtitle = sprintf("Mean matched Jaccard = %.3f  |  %d/%d modules > 0.5",
                          mean(matched$jaccard),
                          sum(matched$jaccard > 0.5), nrow(matched))) +
  FIG_THEME +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8),
        plot.title    = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10))

ggsave(file.path(SUPP_PNG, "SUPP_bicor_sensitivity.png"), p_bicor,
       width = 180, height = 140, units = "mm", dpi = 300)
ggsave(file.path(SUPP_PDF, "SUPP_bicor_sensitivity.pdf"), p_bicor,
       width = 180, height = 140, units = "mm", device = pdf_device)
message("  Supp D saved")

# ── Supp E: Per-module triptych heatmaps ──────────────────────────────────────
message("Supp E: per-module triptychs...")

gene_map <- setNames(ann$uniprot_id, ann$gene)
z_long   <- as.data.frame(group_z) |>
  rownames_to_column("gene") |>
  pivot_longer(-gene, names_to = "group", values_to = "z") |>
  mutate(uniprot_id = gene_map[gene]) |>
  inner_join(module_df |> select(uniprot_id, module_color), by = "uniprot_id") |>
  filter(module_color %in% ALL_MODULES) |>
  select(uniprot_id, gene, group, z, module = module_color)

me_long_e <- MEs |>
  as.data.frame() |>
  rownames_to_column("sample_id") |>
  pivot_longer(starts_with("ME"), names_to = "me_col", values_to = "eigengene") |>
  mutate(module = gsub("^ME", "", me_col)) |>
  filter(module %in% ALL_MODULES) |>
  inner_join(meta |> select(sample_id, pid, group, cancer, timepoint),
             by = "sample_id") |>
  select(sample_id, pid, group, cancer, timepoint, eigengene, module)

lmm_stats_e <- lmm_df |>
  filter(contrast %in% c("Cancer_vs_Healthy", "Training_CR", "Training_CRE")) |>
  mutate(module = gsub("^ME", "", module)) |>
  select(module, contrast, p_bh, r_equiv)

group_order_e <- c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1")
group_labels_e <- c(CRE_T1 = "CRE T1", CRE_T2 = "CRE T2",
                    PLA_T1 = "PLA T1", PLA_T2 = "PLA T2", H_T1 = "H T1")

PE_W <- 280; PE_H <- 110
txt_heat  <- 2.5; txt_title <- 3.5; txt_sig <- 2.8; txt_bar <- 2.5

pathway_slug_e <- function(mc) {
  lbl <- mod_display_vec[[mc]]
  if (is.na(lbl)) return(mc)
  gsub("[/ ]+", "_", tolower(lbl))
}

fmt_sig_e <- function(p) {
  if (length(p) == 0 || is.na(p)) return("ns")
  if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}

for (mod in ALL_MODULES) {
  label    <- coalesce(mod_display_vec[[mod]], str_to_title(mod))
  n_mod    <- z_long |> filter(module == mod) |> distinct(gene) |> nrow()
  title_t  <- paste0(label, " (n=", n_mod, ")")

  z_mod <- z_long |>
    filter(module == mod) |>
    mutate(group = factor(group, levels = group_order_e))
  gene_ord <- z_mod |> filter(group == "CRE_T1") |> arrange(z) |> pull(gene)
  z_mod$gene <- factor(z_mod$gene, levels = gene_ord)

  p_heat_e <- ggplot(z_mod, aes(x = group, y = gene, fill = z)) +
    geom_tile() +
    scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                         midpoint = 0, limits = c(-2, 2), oob = scales::squish,
                         guide = "none") +
    scale_x_discrete(labels = group_labels_e, position = "bottom") +
    labs(title = title_t, y = NULL, x = NULL) +
    FIG_THEME +
    theme(plot.title  = element_text(size = txt_title, face = "bold"),
          axis.text.x = element_text(size = 7, angle = 45, hjust = 1),
          axis.text.y = element_blank(),
          axis.ticks  = element_blank(),
          panel.border = element_blank(),
          plot.margin  = margin(2, 1, 2, 2))

  me_mod <- me_long_e |>
    filter(module == mod) |>
    mutate(timepoint = factor(timepoint, levels = c("T1", "T2")))

  p_cvh <- lmm_stats_e |>
    filter(module == mod, contrast == "Cancer_vs_Healthy") |>
    pull(p_bh)
  p_tr  <- lmm_stats_e |>
    filter(module == mod, contrast == "Training_CR") |>
    pull(p_bh)

  me_means <- me_mod |>
    group_by(cancer, timepoint) |>
    summarise(mean_me = mean(eigengene, na.rm = TRUE), .groups = "drop")

  p_eigen_e <- ggplot(me_mod, aes(x = group, y = eigengene, fill = group)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.shape = NA, linewidth = 0.3) +
    geom_jitter(width = 0.12, size = 0.8, alpha = 0.5, color = "grey30") +
    scale_fill_manual(values = GROUP_COLORS, guide = "none") +
    annotate("text", x = 4.5, y = max(me_mod$eigengene) * 0.95,
             label = paste0("CvH: ", fmt_sig_e(p_cvh)),
             size = txt_sig, fontface = "bold", color = GROUP_COLORS["H_T1"]) +
    annotate("text", x = 4.5, y = max(me_mod$eigengene) * 0.80,
             label = paste0("Tr.: ", fmt_sig_e(p_tr)),
             size = txt_sig, fontface = "bold", color = GROUP_COLORS["CRE_T1"]) +
    labs(y = "Eigengene", x = NULL) +
    FIG_THEME +
    theme(axis.text.x = element_text(size = 7, angle = 35, hjust = 1),
          axis.text.y = element_text(size = 7),
          panel.border = element_blank(),
          panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
          plot.margin  = margin(2, 1, 2, 1))

  bar_data_e <- enrich_df |>
    filter(module == mod, padj < 0.05) |>
    arrange(padj) |>
    head(5) |>
    mutate(
      neg_log10_p = -log10(padj),
      clean_name  = str_trunc(clean_pathway_name(pathway), 38, ellipsis = "..."),
      db_fill     = DB_COLORS[database]
    ) |>
    mutate(clean_name = make.unique(clean_name, sep = " ")) |>
    mutate(clean_name = factor(clean_name, levels = rev(clean_name)))

  if (nrow(bar_data_e) == 0) {
    p_bars_e <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No sig.\nenrichment",
               size = txt_bar, color = "grey50") +
      theme_void() + theme(plot.margin = margin(2, 2, 2, 1))
  } else {
    p_bars_e <- ggplot(bar_data_e, aes(x = neg_log10_p, y = clean_name)) +
      geom_col(aes(fill = db_fill), color = "black", linewidth = 0.3, width = 0.7) +
      geom_text(aes(label = clean_name, x = 0.3), hjust = 0, size = txt_bar,
                fontface = "bold",
                color = ifelse(bar_data_e$db_fill %in% c("#E41A1C", "#4DAF4A", "#377EB8"),
                               "white", "grey20")) +
      scale_fill_identity() +
      scale_x_continuous(expand = expansion(mult = c(0, 0)),
                         name = expression(-log[10](p[adj]))) +
      scale_y_discrete(labels = NULL) +
      labs(y = NULL) +
      FIG_THEME +
      theme(axis.text.x = element_text(size = 7),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            panel.border = element_blank(),
            axis.line.x  = element_line(color = "black", linewidth = 0.3),
            panel.grid   = element_blank(),
            plot.margin  = margin(2, 6, 2, 1))
  }

  z_legend_e <- ggplot(data.frame(z = seq(-2, 2, length.out = 100)),
                        aes(x = z, y = 1, fill = z)) +
    geom_tile() +
    scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                         midpoint = 0, limits = c(-2, 2),
                         name = "Z-score", guide = guide_colorbar(
                           barwidth = unit(40, "mm"), barheight = unit(3, "mm"))) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 7))

  row_e <- p_heat_e + p_eigen_e + p_bars_e + plot_layout(widths = c(3, 2, 3))
  single <- row_e / wrap_elements(z_legend_e) +
    plot_layout(heights = c(1, 0.12)) +
    plot_annotation(
      title = label,
      caption = "* FDR < .05  ** FDR < .01  *** FDR < .001 (LMM, BH-corrected)",
      theme = theme(
        plot.title   = element_text(face = "bold", size = 11),
        plot.caption = element_text(size = 7, color = "grey40", hjust = 0)
      )
    )

  slug <- pathway_slug_e(mod)
  fname <- sprintf("SUPP_triptych_%s_%s", mod, slug)
  ggsave(file.path(MOD_PNG, paste0(fname, ".png")), single,
         width = PE_W, height = PE_H, units = "mm", dpi = 300)
  ggsave(file.path(MOD_PDF, paste0(fname, ".pdf")), single,
         width = PE_W, height = PE_H, units = "mm", device = pdf_device)
  message(sprintf("  Triptych saved: %s", mod))
}
message("  Supp E triptychs complete")

# ── Supp F: Per-module hub networks ───────────────────────────────────────────
message("Supp F: hub networks...")

bg_genes_hub <- unique(module_df$gene[!is.na(module_df$gene) & module_df$gene != ""])
pw_full_hub  <- build_pathway_collection(min_size = 15, max_size = 500,
                                         include_goslim = FALSE)

cor <- WGCNA::cor   # needed for adjacency() inside hub network loop

PF_W <- 170; PF_H <- 170
txt_gene_hub  <- scale_text(BASE_GENE, PF_W)
txt_title_hub <- scale_text(BASE_STAT, PF_W) * 1.8
HULL_PALETTE  <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                   "#66A61E", "#E6AB02", "#A6761D", "#666666")

cor <- WGCNA::cor   # ensure WGCNA cor inside adjacency()

hub_network_list <- list()

for (mod in ALL_MODULES) {
  message(sprintf("  Hub network: %s", mod))
  mod_prots_all <- module_df$uniprot_id[module_df$module_color == mod]
  kme_col <- paste0("kME", mod)
  matched_prots <- intersect(mod_prots_all, rownames(kME_all))
  mod_kme_vals  <- setNames(kME_all[matched_prots, kme_col], matched_prots)
  mod_kme_vals  <- mod_kme_vals[!is.na(mod_kme_vals)]
  q90           <- quantile(mod_kme_vals, 0.90)
  hub_ids       <- names(mod_kme_vals[mod_kme_vals >= q90])

  hub_genes_hub <- uid2gene[hub_ids]
  hub_genes_hub <- hub_genes_hub[!is.na(hub_genes_hub) & hub_genes_hub != ""]
  hub_ids       <- hub_ids[hub_ids %in% names(hub_genes_hub)]

  n_mod_all <- sum(module_df$module_color == mod)
  message(sprintf("    %d hubs (Q90 of %d)", length(hub_ids), n_mod_all))

  mod_prots_expr <- intersect(mod_prots_all, colnames(datExpr))
  if (length(mod_prots_expr) < 5 || length(hub_ids) < 3) {
    message(sprintf("    Skipping %s: insufficient proteins", mod))
    hub_network_list[[mod]] <- NULL
    next
  }

  adj_mod <- adjacency(datExpr[, mod_prots_expr], power = soft_power,
                       type = "signed hybrid")
  tom_mod <- TOMsimilarity(adj_mod, TOMType = "signed")
  colnames(tom_mod) <- rownames(tom_mod) <- mod_prots_expr

  hub_ids_expr <- intersect(hub_ids, mod_prots_expr)
  if (length(hub_ids_expr) < 3) {
    hub_network_list[[mod]] <- NULL; next
  }

  tom_sub  <- tom_mod[hub_ids_expr, hub_ids_expr]
  tom_q90  <- quantile(tom_sub[upper.tri(tom_sub)], 0.90)

  g_hub <- graph_from_adjacency_matrix(tom_sub, mode = "undirected",
                                       weighted = TRUE, diag = FALSE)
  E(g_hub)$weight_orig <- E(g_hub)$weight
  g_hub <- delete_edges(g_hub, which(E(g_hub)$weight < tom_q90))
  iso <- which(igraph::degree(g_hub) == 0)
  if (length(iso) > 0) g_hub <- delete_vertices(g_hub, iso)
  if (vcount(g_hub) < 3) { hub_network_list[[mod]] <- NULL; next }

  node_uids  <- V(g_hub)$name
  node_genes <- uid2gene[node_uids]
  node_kme   <- setNames(kME_all[node_uids, kme_col], node_uids)
  top_label  <- names(sort(node_kme, decreasing = TRUE))[1:min(6, length(node_kme))]

  # ORA-based functional groups
  node_genes_clean <- node_genes[!is.na(node_genes) & node_genes != ""]
  if (length(node_genes_clean) >= 5) {
    ora_hub <- tryCatch(
      run_ora_deduplicated(node_genes_clean, bg_genes_hub, pw_full_hub,
                           jaccard_cutoff = 0.5, min_size = 10, max_size = 500,
                           padj_cutoff = 0.05),
      error = function(e) NULL
    )
  } else {
    ora_hub <- NULL
  }

  clean_pw <- function(name) {
    str_trunc(str_to_title(gsub("_", " ",
      gsub("^HALLMARK_|^GOBP_|^REACTOME_|^KEGG_MEDICUS_|^GOSLIM_", "", name))),
      35)
  }

  func_grp <- setNames(rep("Other", length(node_genes)), names(node_genes))
  if (!is.null(ora_hub) && nrow(ora_hub) > 0) {
    ora_hub  <- ora_hub[order(ora_hub$padj), ]
    gene_grp <- data.frame(gene = character(), grp = character())
    for (i in seq_len(min(nrow(ora_hub), 4))) {
      hits <- intersect(ora_hub$overlapGenes[[i]], node_genes_clean)
      if (length(hits) >= 3)
        gene_grp <- rbind(gene_grp,
          data.frame(gene = hits, grp = clean_pw(ora_hub$pathway[i])))
    }
    gene_grp <- gene_grp[!duplicated(gene_grp$gene), ]
    for (uid in node_uids) {
      g_name <- node_genes[uid]
      if (!is.na(g_name) && g_name %in% gene_grp$gene)
        func_grp[uid] <- gene_grp$grp[gene_grp$gene == g_name][1]
    }
  }

  V(g_hub)$gene     <- node_genes
  V(g_hub)$kME      <- node_kme
  V(g_hub)$func_grp <- func_grp[V(g_hub)$name]

  set.seed(42)
  lay_hub <- layout_with_stress(g_hub)
  nd_hub  <- data.frame(x = lay_hub[, 1], y = lay_hub[, 2],
                        name = node_uids, gene = node_genes,
                        kME = node_kme, func_grp = func_grp[node_uids],
                        stringsAsFactors = FALSE)
  nd_hub$n_in_grp <- as.integer(table(nd_hub$func_grp)[nd_hub$func_grp])

  grp_names <- setdiff(unique(nd_hub$func_grp[nd_hub$func_grp != "Other" &
                                               nd_hub$n_in_grp >= 3]), NA)
  hull_cols <- setNames(HULL_PALETTE[seq_along(grp_names)], grp_names)
  hull_nd   <- nd_hub |> filter(func_grp != "Other", n_in_grp >= 3)

  tg_hub <- as_tbl_graph(g_hub)
  disp_lbl <- coalesce(mod_display_vec[[mod]], str_to_title(mod))

  p_hub <- ggraph(tg_hub, layout = "manual", x = lay_hub[, 1], y = lay_hub[, 2])

  if (nrow(hull_nd) > 0 && length(grp_names) > 0) {
    p_hub <- p_hub +
      geom_mark_hull(data = hull_nd, aes(x = x, y = y, group = func_grp, fill = func_grp),
                     concavity = 2, expand = unit(2, "mm"), radius = unit(2, "mm"),
                     alpha = 0.15, linewidth = 0.6, show.legend = TRUE,
                     inherit.aes = FALSE) +
      scale_fill_manual(values = hull_cols, name = "Pathway")
  }

  p_hub <- p_hub +
    geom_edge_link(aes(width = weight_orig), alpha = 0.5,
                   color = "grey30", show.legend = FALSE) +
    scale_edge_width_continuous(range = c(0.5, 2.0))

  p_hub <- p_hub +
    new_scale_fill() +
    geom_node_point(aes(size = kME, fill = kME), shape = 21,
                    color = "black", stroke = 0.5) +
    scale_size_continuous(range = c(2.0, 7.0), guide = "none") +
    scale_fill_gradient(low = "grey95", high = mod,
                        name = "kME",
                        guide = guide_colorbar(barwidth = unit(3, "mm"),
                                               barheight = unit(18, "mm")))

  label_ids_hub <- top_label[top_label %in% nd_hub$name]
  p_hub <- p_hub +
    geom_label_repel(data = nd_hub[nd_hub$name %in% label_ids_hub, ],
                     aes(x = x, y = y, label = gene),
                     size = txt_gene_hub, fontface = "bold.italic",
                     fill = scales::alpha("white", 0.88), color = "grey10",
                     linewidth = 0.15, label.padding = unit(1.0, "mm"),
                     segment.size = 0.25, segment.color = "grey40",
                     box.padding = 0.45, point.padding = 0.2,
                     max.overlaps = 20, seed = 42, inherit.aes = FALSE)

  p_hub <- p_hub +
    labs(title    = disp_lbl,
         subtitle = sprintf("%d hubs (Q90, n=%d) | color: kME", length(hub_ids_expr), n_mod_all)) +
    theme_void() +
    theme(plot.title      = element_text(face = "bold", size = txt_title_hub, hjust = 0.5),
          plot.subtitle   = element_text(size = 7, hjust = 0.5, color = "grey40"),
          plot.background = element_rect(fill = "white", color = NA),
          plot.margin     = margin(4, 4, 4, 4),
          legend.position = "right",
          legend.title    = element_text(face = "bold", size = 8),
          legend.text     = element_text(size = 7))

  hub_network_list[[mod]] <- list(plot = p_hub, node_data = nd_hub)

  slug_f <- gsub("[/ ]+", "_", tolower(coalesce(mod_display_vec[[mod]], mod)))
  fname_f <- sprintf("SUPP_hub_%s_%s", mod, slug_f)
  ggsave(file.path(MOD_PNG, paste0(fname_f, ".png")), p_hub,
         width = PF_W, height = PF_H, units = "mm", dpi = 300)
  ggsave(file.path(MOD_PDF, paste0(fname_f, ".pdf")), p_hub,
         width = PF_W, height = PF_H, units = "mm", device = pdf_device)
  message(sprintf("  Hub network saved: %s", mod))
}

cor <- stats::cor   # restore

# Hub network composite
plots_hub <- Filter(Negate(is.null), lapply(hub_network_list, function(x) {
  if (is.null(x)) return(NULL)
  x$plot + theme(legend.position = "none", plot.margin = margin(2, 2, 2, 2))
}))

if (length(plots_hub) >= 2) {
  n_c <- min(3, length(plots_hub))
  n_r <- ceiling(length(plots_hub) / n_c)
  p_net_comp <- wrap_plots(plots_hub, ncol = n_c)
  ggsave(file.path(MOD_PNG, "SUPP_networks_composite.png"), p_net_comp,
         width = n_c * 170, height = n_r * 180, units = "mm",
         dpi = 300, limitsize = FALSE)
  ggsave(file.path(MOD_PDF, "SUPP_networks_composite.pdf"), p_net_comp,
         width = n_c * 170, height = n_r * 180, units = "mm",
         device = pdf_device, limitsize = FALSE)
  message("  Hub networks composite saved")
}
message("  Supp F hub networks complete")

# ── Supp G: Inter-module hub chord diagram ────────────────────────────────────
message("Supp G: hub chord diagram...")

chord_hubs <- hub_df_raw |>
  filter(module %in% key_mods) |>
  group_by(module) |>
  slice_max(kME, n = 8, with_ties = FALSE) |>
  ungroup() |>
  mutate(gene_label = coalesce(gene, uniprot_id),
         gene_label = ifelse(is.na(gene_label) | gene_label == "", uniprot_id, gene_label))

if (nrow(chord_hubs) >= 6) {
  link_df_chord <- chord_hubs |>
    transmute(from = str_to_title(module), to = gene_label, value = kME)

  mod_cols_chord  <- setNames(unique(chord_hubs$module), unique(chord_hubs$module))
  grid_cols_chord <- c(
    setNames(unname(mod_cols_chord), str_to_title(names(mod_cols_chord))),
    setNames(rep("grey80", nrow(link_df_chord)), link_df_chord$to)
  )

  pdf_chord <- file.path(SUPP_PDF, "SUPP_hub_chord.pdf")
  png_chord <- file.path(SUPP_PNG, "SUPP_hub_chord.png")

  draw_chord_g <- function() {
    circos.clear()
    circos.par(start.degree = 90, gap.degree = 2,
               track.margin = c(0.005, 0.005))
    chordDiagram(
      link_df_chord,
      grid.col        = grid_cols_chord,
      transparency    = 0.3,
      directional     = 1,
      direction.type  = "diffHeight",
      diffHeight      = -0.04,
      annotationTrack = "grid",
      preAllocateTracks = list(track.height = 0.08)
    )
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
      sector <- get.cell.meta.data("sector.index")
      xlim   <- get.cell.meta.data("xlim")
      ylim   <- get.cell.meta.data("ylim")
      is_mod <- sector %in% str_to_title(names(mod_cols_chord))
      circos.text(mean(xlim), ylim[1] + (if (is_mod) 0.8 else 0.4),
                  sector, facing = "clockwise", niceFacing = TRUE,
                  adj = c(0, 0.5),
                  cex = if (is_mod) 0.9 else 0.55,
                  font = if (is_mod) 2 else 1, col = "grey15")
    }, bg.border = NA)
    title("Key-Module Hub Proteins (top 8 by kME per module)",
          cex.main = 0.95, line = -1)
  }

  pdf(pdf_chord, width = 7.5, height = 7.5)
  draw_chord_g(); dev.off()

  png(png_chord, width = 7.5, height = 7.5, units = "in", res = 300)
  draw_chord_g(); dev.off()

  message(sprintf("  Supp G chord saved (%d modules, %d links)",
                  length(unique(link_df_chord$from)), nrow(link_df_chord)))
} else {
  message("  Supp G: insufficient hub data for chord diagram (< 6 rows)")
}

# ── Supp H: Module preservation (permutation-based quality) ───────────────────
# Substitution for YvO _supp_preservation.R:
# YvO splits Young vs Old (genuine biological reference network).
# Sam has no parallel cohort split. We substitute with permutation-based
# connectivity quality (n=100 random permutations of module labels), comparing
# actual within-module mean adjacency vs permuted null.
# NOTE: Zsummary-style preservation would require a second reference network.
#       This permutation quality metric is analogous but within a single network.
message("Supp H: module preservation (permutation quality)...")

PRES_RDS <- file.path(PANEL_DIR, "wgcna_perm_quality.rds")

if (file.exists(PRES_RDS)) {
  message("  Loading cached permutation quality...")
  perm_quality <- readRDS(PRES_RDS)
} else {
  t_pres_start <- proc.time()
  n_perm <- 100
  message(sprintf("  Running %d permutations per module...", n_perm))

  # Pre-compute full adjacency matrix (signed, same power as main run)
  adj_full <- adjacency(datExpr, power = soft_power, type = "signed hybrid")

  perm_quality <- lapply(ALL_MODULES, function(mod) {
    mod_ids <- intersect(module_df$uniprot_id[module_df$module_color == mod],
                         colnames(datExpr))
    if (length(mod_ids) < 5) return(NULL)

    # Actual within-module mean adjacency
    adj_sub  <- adj_full[mod_ids, mod_ids]
    observed <- mean(adj_sub[upper.tri(adj_sub)])

    # Permuted null: sample same number from all proteins
    all_ids  <- colnames(datExpr)
    n_mod_h  <- length(mod_ids)
    null_vals <- numeric(n_perm)
    for (i in seq_len(n_perm)) {
      samp       <- sample(all_ids, n_mod_h, replace = FALSE)
      adj_samp   <- adj_full[samp, samp]
      null_vals[i] <- mean(adj_samp[upper.tri(adj_samp)])
    }

    z_score <- (observed - mean(null_vals)) / sd(null_vals)
    tibble(
      module       = mod,
      observed_adj = observed,
      null_mean    = mean(null_vals),
      null_sd      = sd(null_vals),
      Zsummary     = z_score,
      module_size  = n_mod_h
    )
  })
  perm_quality <- bind_rows(perm_quality)

  t_elapsed <- (proc.time() - t_pres_start)[3]
  message(sprintf("  Permutation quality done (%.1f min)", t_elapsed / 60))
  saveRDS(perm_quality, PRES_RDS)
}

write_csv(perm_quality, file.path(PANEL_DIR, "supp", "a08_perm_quality.csv"))

perm_quality <- perm_quality |>
  mutate(
    bio_label = coalesce(mod_display_vec[module], str_to_title(module)),
    quality   = case_when(
      Zsummary > 10 ~ "Strong",
      Zsummary > 2  ~ "Moderate",
      TRUE          ~ "Weak"
    )
  ) |>
  arrange(Zsummary) |>
  mutate(bio_label = factor(bio_label, levels = bio_label))

p_pres <- ggplot(perm_quality, aes(x = Zsummary, y = bio_label)) +
  geom_vline(xintercept = c(2, 10), linetype = "dashed",
             color = "grey50", linewidth = 0.35) +
  geom_col(fill = perm_quality$module, color = "black",
           linewidth = 0.3, width = 0.65) +
  geom_text(aes(label = sprintf("%.1f", Zsummary), x = pmax(Zsummary / 2, 0.2)),
            size = 2.8, fontface = "bold", color = "white") +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(0, max(perm_quality$Zsummary, na.rm = TRUE) * 1.2),
                     breaks = c(0, 2, 10, 20, 30)) +
  labs(x = "Z-score (permutation quality)", y = NULL,
       title = "Module Connectivity Quality (Permutation Test)",
       subtitle = sprintf("n=%d permutations | Z>2: moderate | Z>10: strong | One-network alternative to Zsummary",
                          100)) +
  FIG_THEME +
  theme(axis.text.y = element_text(size = 7.5),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.border       = element_blank(),
        axis.line.x        = element_line(color = "black", linewidth = 0.3),
        legend.position    = "none")

ggsave(file.path(SUPP_PNG, "SUPP_preservation.png"), p_pres,
       width = 200, height = max(100, length(ALL_MODULES) * 14 + 30),
       units = "mm", dpi = 300, limitsize = FALSE)
ggsave(file.path(SUPP_PDF, "SUPP_preservation.pdf"), p_pres,
       width = 200, height = max(100, length(ALL_MODULES) * 14 + 30),
       units = "mm", device = pdf_device, limitsize = FALSE)
message("  Supp H preservation saved")

# ── Supp QC composite (A + B + C + D — mirrors YvO SUPP_F06_composite) ───────
message("Building Supp QC composite (panels A-D)...")

supp_panel_files <- c(
  "SUPP_soft_threshold.png",
  "SUPP_dendrogram.png",
  "SUPP_compartment_enrichment.png",
  "SUPP_bicor_sensitivity.png"
)
supp_exists <- file.exists(file.path(SUPP_PNG, supp_panel_files))

if (all(supp_exists)) {
  read_grob_supp <- function(f) rasterGrob(readPNG(file.path(SUPP_PNG, f)),
                                            interpolate = TRUE)
  pA_supp <- read_grob_supp(supp_panel_files[1])
  pB_supp <- read_grob_supp(supp_panel_files[2])
  pC_supp <- read_grob_supp(supp_panel_files[3])
  pD_supp <- read_grob_supp(supp_panel_files[4])

  bottom_row_supp <- wrap_elements(full = pC_supp) + wrap_elements(full = pD_supp) +
    plot_layout(widths = c(1, 1))
  composite_supp <- wrap_elements(full = pA_supp) /
                    wrap_elements(full = pB_supp) /
                    bottom_row_supp +
    plot_layout(heights = c(0.30, 0.35, 0.35)) +
    plot_annotation(theme = theme(plot.margin = margin(4, 6, 4, 6)))

  TAG_SZ_S <- 16
  composite_supp <- ggdraw(composite_supp) +
    draw_label("A", x = 0.02, y = 0.960, size = TAG_SZ_S, fontface = "bold",
               hjust = 0, vjust = 1) +
    draw_label("B", x = 0.02, y = 0.660, size = TAG_SZ_S, fontface = "bold",
               hjust = 0, vjust = 1) +
    draw_label("C", x = 0.02, y = 0.320, size = TAG_SZ_S, fontface = "bold",
               hjust = 0, vjust = 1) +
    draw_label("D", x = 0.52, y = 0.320, size = TAG_SZ_S, fontface = "bold",
               hjust = 0, vjust = 1)

  supp_pdf_dir <- file.path(BASE, "b_reports", "supp", "pdf")
  supp_png_dir <- file.path(BASE, "b_reports", "supp", "png")
  dir.create(supp_pdf_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(supp_png_dir, recursive = TRUE, showWarnings = FALSE)

  COMP_W_S <- 250; COMP_H_S <- 330
  ggsave(file.path(supp_pdf_dir, "SUPP_F06_composite.pdf"), composite_supp,
         width = COMP_W_S, height = COMP_H_S, units = "mm",
         device = pdf_device, limitsize = FALSE)
  ggsave(file.path(supp_png_dir, "SUPP_F06_composite.png"), composite_supp,
         width = COMP_W_S, height = COMP_H_S, units = "mm",
         dpi = 300, limitsize = FALSE)
  message("SUPP_F06_composite (A-D) saved")
} else {
  message("Cannot build supp composite — missing: ",
          paste(supp_panel_files[!supp_exists], collapse = ", "))
}

message("02_supp_panels.R complete")
