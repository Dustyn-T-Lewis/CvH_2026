# Reversal Panel E: RRHO2 Heatmap
# Stratified rank-rank hypergeometric overlap (Cahill et al. 2018)
# Uses the RRHO2 R package (Plaisier et al. 2010, NAR)
setwd(here::here())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")
pacman::p_load(tidyverse, RRHO2, fgsea)

PE_W <- 110

RPT_PNG <- "04_Figures/F04_Reversal/b_reports/supp/png"
RPT_PDF <- "04_Figures/F04_Reversal/b_reports/supp/pdf"
DAT <- "04_Figures/F04_Reversal/c_data"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_E"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# --- Load DEP data & build rank lists ---
source("04_Figures/F04_Reversal/a_script/f04_data.R") # dep_df (old column names)

rr_df <- dep_df %>%
  transmute(gene, t_1 = t_CRvH_Baseline, t_2 = t_CR_Training) %>%
  filter(!is.na(t_1) & !is.na(t_2)) %>%
  distinct(gene, .keep_all = TRUE)

n_shared <- nrow(rr_df)

# --- RRHO2 computation (Cahill et al. 2018) ---
list1 <- data.frame(gene = rr_df$gene, score = rr_df$t_1, stringsAsFactors = FALSE)
list2 <- data.frame(gene = rr_df$gene, score = rr_df$t_2, stringsAsFactors = FALSE)

rrho_obj <- RRHO2_initialize(
  list1, list2,
  labels = c("CRvH_Baseline", "CR_Training"),
  log10.ind = TRUE,
  multipleTesting = "none",
  boundary = 0.02,
  method = "hyper",
  stepsize = 20
)

hmat <- rrho_obj$hypermat
nr <- nrow(hmat)
nc <- ncol(hmat)
message(sprintf("  RRHO2 matrix: %d x %d", nr, nc))

# Locate NA boundary strip (zero-crossing)
na_rows <- which(apply(hmat, 1, function(r) all(is.na(r))))
na_cols <- which(apply(hmat, 2, function(c) all(is.na(c))))

if (length(na_rows) && length(na_cols)) {
  row_before <- 1:(min(na_rows) - 1)
  row_after <- (max(na_rows) + 1):nr
  col_before <- 1:(min(na_cols) - 1)
  col_after <- (max(na_cols) + 1):nc
} else {
  mid <- floor(nr / 2)
  row_before <- 1:mid
  row_after <- (mid + 1):nr
  col_before <- 1:mid
  col_after <- (mid + 1):nc
}

# Both lists sorted descending (Up first):
# [row_before, col_before] = UU (both up) = Exacerbated Up
# [row_after,  col_after]  = DD (both down) = Exacerbated Down
# [row_before, col_after]  = UD (list1 up, list2 down) = Reversed (Cancer↑ Tr↓)
# [row_after,  col_before] = DU (list1 down, list2 up) = Reversed (Cancer↓ Tr↑)
max_UU <- max(hmat[row_before, col_before], na.rm = TRUE)
max_DD <- max(hmat[row_after, col_after], na.rm = TRUE)
max_UD <- max(hmat[row_before, col_after], na.rm = TRUE)
max_DU <- max(hmat[row_after, col_before], na.rm = TRUE)

message(sprintf(
  "  Max -log10(p): UU=%.1f, DD=%.1f, UD=%.1f, DU=%.1f",
  max_UU, max_DD, max_UD, max_DU
))

# Extract hotspot genes from RRHO2 built-in gene lists
hotspot_genes <- list(
  UU = rrho_obj$genelist_uu$gene_list_overlap_uu,
  DD = rrho_obj$genelist_dd$gene_list_overlap_dd,
  UD = rrho_obj$genelist_ud$gene_list_overlap_ud,
  DU = rrho_obj$genelist_du$gene_list_overlap_du
)

n_UU <- length(hotspot_genes$UU)
n_DD <- length(hotspot_genes$DD)
n_UD <- length(hotspot_genes$UD)
n_DU <- length(hotspot_genes$DU)
message(sprintf("  Hotspot genes: UU=%d, DD=%d, UD=%d, DU=%d", n_UU, n_DD, n_UD, n_DU))

# --- Jet colormap heatmap ---
JET_COLORS <- c(
  "#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
  "yellow", "#FF7F00", "red", "#7F0000"
)

hmat_df <- expand.grid(row = 1:nr, col = 1:nc) %>%
  mutate(neg_log10_p = as.vector(hmat))

max_val <- max(hmat_df$neg_log10_p, na.rm = TRUE)
txt_quad <- 2.5

LABEL_FILL <- scales::alpha("white", 0.85)
LABEL_PADDING <- unit(1.5, "mm")

# Corner anchors
ann_x_left <- min(row_before)
ann_x_right <- max(row_after)
ann_y_bot <- min(col_before)
ann_y_top <- max(col_after)

pE_heat <- ggplot(hmat_df, aes(x = row, y = col, fill = neg_log10_p)) +
  geom_raster() +
  scale_fill_gradientn(
    colors = JET_COLORS,
    limits = c(0, max_val),
    na.value = "white",
    name = expression(-log[10](P)),
    guide = guide_colorbar(
      barwidth = unit(18, "mm"), barheight = unit(2.5, "mm"),
      title.position = "left", title.vjust = 0.5,
      title.theme = element_text(
        size = FIG_LEGEND_TITLE, face = "bold",
        color = "grey15"
      )
    )
  ) +
  # Bottom-left: UU (Exacerbated Up)
  annotate("label",
    x = ann_x_left, y = ann_y_bot,
    label = sprintf("Exacerbated Up\n(max %.0f, n=%d)", max_UU, n_UU),
    color = "grey15", fill = LABEL_FILL, linewidth = 0,
    label.padding = LABEL_PADDING, fontface = "bold", size = txt_quad,
    hjust = 0, vjust = 0
  ) +
  # Top-right: DD (Exacerbated Down)
  annotate("label",
    x = ann_x_right, y = ann_y_top,
    label = sprintf("Exacerbated Dn\n(max %.0f, n=%d)", max_DD, n_DD),
    color = "grey15", fill = LABEL_FILL, linewidth = 0,
    label.padding = LABEL_PADDING, fontface = "bold", size = txt_quad,
    hjust = 1, vjust = 1
  ) +
  # Top-left: DU (Reversed: Cancer Down, Training Up)
  annotate("label",
    x = ann_x_left, y = ann_y_top,
    label = sprintf("Reversed (C\u2193 T\u2191)\n(max %.0f, n=%d)", max_DU, n_DU),
    color = "grey15", fill = LABEL_FILL, linewidth = 0,
    label.padding = LABEL_PADDING, fontface = "bold", size = txt_quad,
    hjust = 0, vjust = 1
  ) +
  # Bottom-right: UD (Reversed: Cancer Up, Training Down)
  annotate("label",
    x = ann_x_right, y = ann_y_bot,
    label = sprintf("Reversed (C\u2191 T\u2193)\n(max %.0f, n=%d)", max_UD, n_UD),
    color = "grey15", fill = LABEL_FILL, linewidth = 0,
    label.padding = LABEL_PADDING, fontface = "bold", size = txt_quad,
    hjust = 1, vjust = 0
  ) +
  scale_x_continuous(expand = expansion(mult = 0.015)) +
  scale_y_continuous(expand = expansion(mult = 0.015)) +
  labs(
    title = "RRHO2: Cancer Recovery Reversal",
    subtitle = sprintf(
      "Stratified hypergeometric | %d genes | peak = %.0f",
      n_shared, max_val
    ),
    x = sprintf("Rank: Cancer vs Healthy (Up %s Down)", "\u2192"),
    y = sprintf("Rank: Training CR (Up %s Down)", "\u2192")
  ) +
  FIG_THEME +
  theme(
    axis.text        = element_blank(),
    axis.title.x     = element_text(size = FIG_AXIS_TEXT, face = "bold", margin = margin(t = 2)),
    axis.title.y     = element_text(size = FIG_AXIS_TEXT, face = "bold", margin = margin(r = 2)),
    axis.ticks       = element_blank(),
    panel.border     = element_blank(),
    panel.grid.major = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = FIG_LEGEND_TEXT, face = "bold"),
    legend.margin    = margin(2, 24, 0, 0, "mm"),
    plot.margin      = margin(0, 0, 0, 0, "mm")
  ) +
  coord_fixed(ratio = 1, clip = "off")

ggsave(file.path(RPT_PNG, "SUPP_F04_rrho2.png"), pE_heat,
  width = PE_W, height = PE_W, units = "mm", dpi = 300
)
ggsave(file.path(RPT_PDF, "SUPP_F04_rrho2.pdf"), pE_heat,
  width = PE_W, height = PE_W, units = "mm", device = pdf_device
)

# --- Export hotspot genes ---
hotspot_export <- bind_rows(
  tibble(quadrant = "Exacerbated Up", gene = hotspot_genes$UU),
  tibble(quadrant = "Exacerbated Down", gene = hotspot_genes$DD),
  tibble(quadrant = "Reversed (Cancer Up)", gene = hotspot_genes$UD),
  tibble(quadrant = "Reversed (Cancer Down)", gene = hotspot_genes$DU)
)
write_csv(hotspot_export, file.path(DAT, "panel_E", "rrho2_hotspot_genes.csv"))

# --- Per-quadrant ORA ---
pw_collection_E <- build_pathway_collection(
  min_size = 15, max_size = 500,
  include_goslim = FALSE,
  exclude_variants = TRUE
)
all_genes_E <- unique(rr_df$gene)

run_quadrant_ora <- function(gene_set, quadrant_name) {
  if (length(gene_set) < 5) {
    return(tibble())
  }
  tryCatch(
    run_ora_deduplicated(
      genes = unique(gene_set), universe = all_genes_E,
      pathways = pw_collection_E, jaccard_cutoff = 0.5,
      min_size = 15, max_size = 500, padj_cutoff = 0.05
    ) %>%
      mutate(
        quadrant = quadrant_name,
        pathway_label = clean_pathway_name(pathway)
      ),
    error = function(e) {
      message("  ORA error: ", e$message)
      tibble()
    }
  )
}

ora_rev_up <- run_quadrant_ora(hotspot_genes$UD, "Reversed (Cancer Up)")
ora_rev_down <- run_quadrant_ora(hotspot_genes$DU, "Reversed (Cancer Down)")
ora_exac_up <- run_quadrant_ora(hotspot_genes$UU, "Exacerbated Up")
ora_exac_down <- run_quadrant_ora(hotspot_genes$DD, "Exacerbated Down")

ora_concordant <- bind_rows(ora_exac_up, ora_exac_down)
ora_discordant <- bind_rows(ora_rev_up, ora_rev_down)

write_csv(ora_concordant, file.path(DAT, "panel_E", "rrho2_ora_concordant.csv"))
write_csv(ora_discordant, file.path(DAT, "panel_E", "rrho2_ora_discordant.csv"))

# --- Summary CSV ---
rrho2_meta <- tibble(
  quadrant = c(
    "Exacerbated Up", "Exacerbated Down",
    "Reversed (Cancer Up)", "Reversed (Cancer Down)"
  ),
  max_neg_log10_pvalue = round(c(max_UU, max_DD, max_UD, max_DU), 2),
  n_hotspot_genes = c(n_UU, n_DD, n_UD, n_DU),
  n_ora_pathways = c(
    nrow(ora_exac_up), nrow(ora_exac_down),
    nrow(ora_rev_up), nrow(ora_rev_down)
  ),
  matrix_rows = nr, matrix_cols = nc, n_shared_genes = n_shared
)
write_csv(rrho2_meta, file.path(DAT, "panel_E", "rrho2_summary.csv"))

# --- Export for composite ---
pE_title <- "RRHO2"
pE_subtitle <- sprintf("%d genes | peak = %.0f", n_shared, max_val)
pE_legend <- NULL
pE_heat <- pE_heat +
  labs(title = NULL, subtitle = NULL, tag = NULL) +
  coord_fixed(ratio = 1, clip = "off")

message("Reversal Panel E (RRHO2) done")
