################################################################################
#   Figure 6 — Panel C: Hub Protein Networks (2x2 grid, top 4 modules)
#   Hub selection: kME >= module Q90 (data-driven, no cap)
#   Pathway DB: run_ora_deduplicated() with full multi-DB collection
#   Layout: stress (graphlayouts, deterministic)
#   Generates: panel_C_hub_network_MAIN.pdf/.png, c_data/04_panel_C_*.csv
################################################################################

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F06/a_script/style.R")
source("04_Figures/shared/pathway_utils.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(ggrepel)
  library(WGCNA)
  library(igraph)
  library(ggraph)
  library(ggforce)
  library(concaveman)
  library(graphlayouts)
  library(tidygraph)
  library(ggnewscale)
  library(fgsea)
  library(colorspace)
})

allowWGCNAThreads()
set.seed(42)

RPT <- "04_Figures/F06/b_reports"
DAT <- "04_Figures/F06/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

message("Panel C: Hub protein networks...")

# --- Load data ---
meta <- read_csv(file.path(DAT, "meta.csv"), show_col_types = FALSE)
meta$group_time <- factor(meta$group_time,
                          levels = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1"))
MEs     <- readRDS(file.path(DAT, "MEs.rds"))
kME_all <- readRDS(file.path(DAT, "kME_all.rds"))
datExpr <- readRDS(file.path(DAT, "datExpr.rds"))
module_df <- read_csv(file.path(DAT, "wgcna_module_assignments.csv"),
                      show_col_types = FALSE)
sft_csv   <- read.csv(file.path(DAT, "../../../05_WGCNA/c_data/wgcna/wgcna_sft_summary.csv"))
NET_POWER <- sft_csv$selected_power[1]

mod_bio_labels_df  <- read_csv(file.path(DAT, "mod_bio_labels.csv"), show_col_types = FALSE)
if (!"display_label" %in% colnames(mod_bio_labels_df)) {
  mod_bio_labels_df <- mod_bio_labels_df %>%
    mutate(display_label = paste0(module_id, ": ", module_color, " (n=", n_proteins, ")"))
}
display_label_vec <- setNames(mod_bio_labels_df$display_label, mod_bio_labels_df$module_color)

bg_genes    <- unique(module_df$gene)
KEY_MODULES <- readLines(file.path(DAT, "key_modules.txt"))
KEY_MODULES <- KEY_MODULES[nzchar(trimws(KEY_MODULES))]
uid2gene    <- setNames(module_df$gene, module_df$uniprot_id)

# --- Contrasting hull palette (Dark2-derived) ---
HULL_PALETTE <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                  "#66A61E", "#E6AB02", "#A6761D", "#666666")

# --- Pathway collection (full multi-DB) ---
pw_full <- build_pathway_collection(min_size = 15, max_size = 500, include_goslim = FALSE)

# --- Hub selection: Q90 per module ---
select_hubs_q90 <- function(mod) {
  mod_prots <- module_df$uniprot_id[module_df$module_color == mod]
  kme_col   <- paste0("kME", mod)
  matched   <- intersect(mod_prots, rownames(kME_all))
  mod_kme   <- setNames(kME_all[matched, kme_col], matched)
  mod_kme   <- mod_kme[!is.na(mod_kme)]
  q90       <- quantile(mod_kme, 0.90)
  names(mod_kme[mod_kme >= q90])
}

# --- Functional group assignment via ORA ---
assign_groups_ora <- function(gene_names, max_groups = 4, min_group_n = 3) {
  clean_pw_name <- function(name) {
    name %>%
      gsub("^HALLMARK_|^GOSLIM_|^GOBP_|^REACTOME_|^KEGG_MEDICUS_|^KEGG_", "", .) %>%
      gsub("_", " ", .) %>% str_to_title() %>% str_trunc(35)
  }

  ora_res <- tryCatch(
    run_ora_deduplicated(
      genes = gene_names, universe = bg_genes,
      pathways = pw_full, jaccard_cutoff = 0.5,
      min_size = 5, max_size = 500, padj_cutoff = 1
    ),
    error = function(e) NULL
  )

  if (is.null(ora_res) || nrow(ora_res) == 0) {
    return(setNames(rep("Other", length(gene_names)), gene_names))
  }

  top_pw <- ora_res %>%
    filter(padj < 0.1) %>%
    arrange(padj) %>%
    head(max_groups)

  gene_group <- setNames(rep("Other", length(gene_names)), gene_names)

  for (i in seq_len(nrow(top_pw))) {
    pw_genes <- top_pw$overlapGenes[[i]]
    shared <- intersect(pw_genes, gene_names)
    unassigned <- shared[gene_group[shared] == "Other"]
    if (length(unassigned) >= min_group_n) {
      gene_group[unassigned] <- clean_pw_name(top_pw$pathway[i])
    }
  }

  gene_group
}

# --- Build network for one module ---
build_hub_network <- function(mod) {
  hub_ids <- select_hubs_q90(mod)
  if (length(hub_ids) < 5) {
    message(sprintf("  Skipping %s: only %d hubs", mod, length(hub_ids)))
    return(NULL)
  }

  hub_expr <- datExpr[, hub_ids, drop = FALSE]
  cor_mat <- WGCNA::cor(hub_expr, use = "pairwise.complete.obs")

  adj_mat <- abs(cor_mat)^NET_POWER
  diag(adj_mat) <- 0

  # Threshold: keep edges > median adjacency
  threshold <- median(adj_mat[adj_mat > 0])
  adj_mat[adj_mat < threshold] <- 0

  g <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted = TRUE)

  # Gene names
  V(g)$gene <- uid2gene[V(g)$name]
  V(g)$gene[is.na(V(g)$gene)] <- V(g)$name[is.na(V(g)$gene)]

  # kME for node sizing
  kme_col <- paste0("kME", mod)
  V(g)$kME <- kME_all[V(g)$name, kme_col]

  # Functional groups via ORA
  hub_genes <- V(g)$gene
  func_groups <- assign_groups_ora(hub_genes)
  V(g)$func_group <- func_groups[hub_genes]

  # Layout: stress
  set.seed(42)
  layout <- layout_with_stress(g)

  tg <- as_tbl_graph(g)

  # Assign hull colors
  groups_present <- unique(V(g)$func_group)
  groups_present <- c(setdiff(groups_present, "Other"), "Other")
  group_colors <- setNames(
    c(HULL_PALETTE[seq_along(setdiff(groups_present, "Other"))], "grey80"),
    groups_present
  )

  mod_label <- if (mod %in% names(display_label_vec)) display_label_vec[mod] else mod

  p <- ggraph(tg, layout = layout) +
    geom_edge_link(aes(alpha = weight), edge_colour = "grey70", show.legend = FALSE) +
    scale_edge_alpha(range = c(0.1, 0.6)) +
    geom_node_point(aes(size = kME, fill = func_group),
                    shape = 21, color = "black", stroke = 0.4) +
    geom_node_text(aes(label = gene), size = 2.2, repel = TRUE,
                   max.overlaps = 20, segment.size = 0.2) +
    scale_size_continuous(range = c(2, 7), guide = "none") +
    scale_fill_manual(values = group_colors, name = "Function") +
    labs(title = mod_label) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      legend.position = "bottom",
      legend.text = element_text(size = 6.5),
      legend.title = element_text(size = 7.5, face = "bold"),
      legend.key.size = unit(3, "mm")
    ) +
    guides(fill = guide_legend(ncol = 2, override.aes = list(size = 3)))

  p
}

# --- Build all 4 networks ---
cor <- WGCNA::cor
plots <- lapply(KEY_MODULES, build_hub_network)
cor <- stats::cor

# Remove NULLs
plots <- Filter(Negate(is.null), plots)

if (length(plots) == 0) {
  message("  No hub networks generated — skipping Panel C")
} else {
  # 2x2 grid
  n_plots <- min(length(plots), 4)
  composite <- wrap_plots(plots[1:n_plots], ncol = 2) +
    plot_annotation(
      title = "Hub Protein Networks (kME >= Q90)",
      theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5))
    )

  W <- 260; H <- 260

  ggsave(file.path(RPT, "panel_C_hub_network_MAIN.pdf"), composite,
         width = W, height = H, units = "mm",
         device = pdf_device, limitsize = FALSE)
  ggsave(file.path(RPT, "panel_C_hub_network_MAIN.png"), composite,
         width = W, height = H, units = "mm",
         dpi = 300, limitsize = FALSE)

  # Save hub list
  hub_list <- lapply(KEY_MODULES, function(mod) {
    ids <- select_hubs_q90(mod)
    tibble(module = mod, uniprot_id = ids, gene = uid2gene[ids])
  })
  hub_export <- bind_rows(hub_list)
  write_csv(hub_export, file.path(DAT, "04_panel_C_hub_proteins.csv"))

  message("  Panel C (hub networks) saved")
}
