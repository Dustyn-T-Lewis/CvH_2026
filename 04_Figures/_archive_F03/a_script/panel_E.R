# F03 Panel E: RRHO2 Concordance Map + Per-Quadrant ORA (CRE vs PLA)
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(msigdbr)
  library(fgsea)
  library(RRHO2)
})

PE_W <- 260

RPT <- "04_Figures/F03/b_reports"
DAT <- "04_Figures/F03/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_E"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

dep_df <- read_csv("03_DEP/c_data/03_combined_results_CR.csv", show_col_types = FALSE)

rr_df <- dep_df %>%
  transmute(gene, t_cre = t_Training_CRE, t_pla = t_Training_PLA) %>%
  filter(!is.na(t_cre) & !is.na(t_pla)) %>%
  distinct(gene, .keep_all = TRUE)

n_shared <- nrow(rr_df)

# --- RRHO2 computation ---
list1 <- data.frame(gene = rr_df$gene, score = rr_df$t_cre, stringsAsFactors = FALSE)
list2 <- data.frame(gene = rr_df$gene, score = rr_df$t_pla, stringsAsFactors = FALSE)

rrho_obj <- RRHO2_initialize(
  list1, list2,
  labels          = c("Training (CRE)", "Training (PLA)"),
  log10.ind       = TRUE,
  multipleTesting = "none",
  boundary        = 0.02,
  method          = "hyper",
  stepsize        = 20
)

hmat <- rrho_obj$hypermat
nr <- nrow(hmat); nc <- ncol(hmat)
message(sprintf("  RRHO2 matrix: %d x %d", nr, nc))

na_rows <- which(apply(hmat, 1, function(r) all(is.na(r))))
na_cols <- which(apply(hmat, 2, function(c) all(is.na(c))))

if (length(na_rows) > 0 && length(na_cols) > 0) {
  row_before <- 1:(min(na_rows) - 1)
  row_after  <- (max(na_rows) + 1):nr
  col_before <- 1:(min(na_cols) - 1)
  col_after  <- (max(na_cols) + 1):nc
} else {
  mid <- floor(nr / 2)
  row_before <- 1:mid
  row_after  <- (mid + 1):nr
  col_before <- 1:mid
  col_after  <- (mid + 1):nc
}

max_UU <- max(hmat[row_before, col_before], na.rm = TRUE)
max_DD <- max(hmat[row_after,  col_after],  na.rm = TRUE)
max_UD <- max(hmat[row_before, col_after],  na.rm = TRUE)
max_DU <- max(hmat[row_after,  col_before], na.rm = TRUE)

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

txt_quad <- scale_text(BASE_QUADRANT, PE_W)

# --- Jet colormap heatmap ---
JET_COLORS <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                "yellow", "#FF7F00", "red", "#7F0000")

hmat_df <- expand.grid(row = 1:nr, col = 1:nc) %>%
  mutate(neg_log10_p = as.vector(hmat))

max_val <- max(hmat_df$neg_log10_p, na.rm = TRUE)

ann_x_left  <- mean(row_before)
ann_x_right <- mean(row_after)
ann_y_bot   <- min(col_before) + 0.12 * diff(range(col_before))
ann_y_top   <- max(col_after)  - 0.08 * diff(range(col_after))

pE_heat <- ggplot(hmat_df, aes(x = row, y = col, fill = neg_log10_p)) +
  geom_raster() +
  scale_fill_gradientn(
    colors   = JET_COLORS,
    limits   = c(0, max_val),
    na.value = "white",
    name     = expression(-log[10](P)),
    guide    = guide_colorbar(
      barwidth = unit(30, "mm"), barheight = unit(3, "mm"),
      title.position = "top", title.hjust = 0.5,
      title.theme = element_text(size = 6, face = "bold"))
  ) +
  annotate("text", x = ann_x_left, y = ann_y_bot,
           label = sprintf("Concordant  CRE inc / PLA inc\nP < 1e-%d  |  n = %d genes",
                           round(max_UU), n_UU),
           color = "white", fontface = "bold", size = txt_quad * 0.85,
           lineheight = 0.85) +
  annotate("text", x = ann_x_right, y = ann_y_top,
           label = sprintf("Concordant  CRE dec / PLA dec\nP < 1e-%d  |  n = %d genes",
                           round(max_DD), n_DD),
           color = "white", fontface = "bold", size = txt_quad * 0.85,
           lineheight = 0.85) +
  annotate("text", x = ann_x_left, y = ann_y_top,
           label = sprintf("Discordant  CRE inc / PLA dec\nP < 1e-%d  |  n = %d genes",
                           round(max_UD), n_UD),
           color = "white", fontface = "bold", size = txt_quad * 0.85,
           lineheight = 0.85) +
  annotate("text", x = ann_x_right, y = ann_y_bot,
           label = sprintf("Discordant  CRE dec / PLA inc\nP < 1e-%d  |  n = %d genes",
                           round(max_DU), n_DU),
           color = "white", fontface = "bold", size = txt_quad * 0.85,
           lineheight = 0.85) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    title    = "Threshold-Free Concordance (RRHO2)",
    subtitle = sprintf("Stratified hypergeometric | %d shared genes", n_shared),
    x = expression("Training (CRE) rank"~(Up %->% Down)),
    y = expression("Training (PLA) rank"~(Up %->% Down))
  ) +
  FIG_THEME +
  theme(
    axis.text        = element_blank(),
    axis.title.x     = element_text(margin = margin(t = 2)),
    axis.title.y     = element_text(margin = margin(r = 2)),
    axis.ticks       = element_blank(),
    panel.border     = element_blank(),
    panel.grid.major = element_blank(),
    legend.position  = "bottom",
    legend.margin    = margin(0, 0, 0, 0),
    plot.margin = margin(2, 2, 2, 2, "mm")
  ) +
  coord_fixed(ratio = 1)

hotspot_export <- bind_rows(
  tibble(quadrant = "UU", gene = hotspot_genes$UU),
  tibble(quadrant = "DD", gene = hotspot_genes$DD),
  tibble(quadrant = "UD", gene = hotspot_genes$UD),
  tibble(quadrant = "DU", gene = hotspot_genes$DU)
)
write_csv(hotspot_export, file.path(DAT, "panel_E", "rrho2_hotspot_genes.csv"))

pw_collection_E <- build_pathway_collection(min_size = 15, max_size = 500)
all_genes_E <- rr_df$gene

run_quadrant_ora <- function(gene_set, quadrant_name) {
  if (length(gene_set) < 5) return(tibble())
  res <- tryCatch(
    run_ora_deduplicated(
      genes          = gene_set,
      universe       = all_genes_E,
      pathways       = pw_collection_E,
      jaccard_cutoff = 0.5,
      min_size       = 15,
      max_size       = 500,
      padj_cutoff    = 0.05
    ),
    error = function(e) { message("  ORA error: ", e$message); tibble() }
  )
  if (nrow(res) > 0) {
    res %>%
      mutate(quadrant = quadrant_name,
             pathway_label = clean_pathway_name(pathway)) %>%
      arrange(padj, size)
  } else {
    tibble()
  }
}

ora_UU <- run_quadrant_ora(hotspot_genes$UU, "Concordant Up")
ora_DD <- run_quadrant_ora(hotspot_genes$DD, "Concordant Down")
ora_UD <- run_quadrant_ora(hotspot_genes$UD, "Discordant (CRE Up / PLA Down)")
ora_DU <- run_quadrant_ora(hotspot_genes$DU, "Discordant (CRE Down / PLA Up)")

ora_all <- bind_rows(ora_UU, ora_DD, ora_UD, ora_DU)
write_csv(ora_all, file.path(DAT, "panel_E", "rrho2_ora_all.csv"))

ggsave(file.path(RPT, "panel_E_RRHO.pdf"), pE_heat,
       width = PE_W, height = PE_W, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_E_RRHO.png"), pE_heat,
       width = PE_W, height = PE_W, units = "mm", dpi = 300)

# --- Per-quadrant ORA bar charts ---
ORA_W      <- 260
BAR_H_MM   <- 10
MARGIN_MM  <- 35
MAX_PER_QUAD <- 12
txt_ora <- scale_text(BASE_STAT, PE_W)

make_ora_panel <- function(df, quad_name, quad_color, n_hotspot, txt_size) {
  n_bars <- nrow(df)
  if (n_bars == 0) return(NULL)
  df <- df %>%
    mutate(neg_log10_padj = -log10(padj)) %>%
    arrange(neg_log10_padj) %>%
    mutate(pathway_label = fct_inorder(pathway_label))
  ggplot(df, aes(x = neg_log10_padj, y = pathway_label)) +
    geom_col(fill = quad_color, width = 0.75) +
    geom_text(aes(x = 0.05, label = pathway_label),
              hjust = 0, size = txt_size * 0.85,
              color = "white", fontface = "bold") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
    scale_y_discrete(labels = NULL) +
    labs(
      title    = quad_name,
      subtitle = sprintf("%d pathways  |  %d hotspot genes", n_bars, n_hotspot),
      x = expression(-log[10](p[adj])),
      y = NULL
    ) +
    FIG_THEME +
    theme(
      plot.title         = element_text(size = 11, face = "bold", hjust = 0),
      plot.subtitle      = element_text(size = 9, hjust = 0, color = "grey40"),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor   = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin  = margin(4, 6, 4, 4, "mm")
    )
}

quad_meta <- list(
  list(name = "Concordant Up",                 slug = "concordant_up",     data = ora_UU, n_hot = n_UU),
  list(name = "Concordant Down",               slug = "concordant_down",   data = ora_DD, n_hot = n_DD),
  list(name = "Discordant (CRE Up / PLA Down)", slug = "discordant_cre_up", data = ora_UD, n_hot = n_UD),
  list(name = "Discordant (CRE Down / PLA Up)", slug = "discordant_cre_down", data = ora_DU, n_hot = n_DU)
)

for (qm in quad_meta) {
  if (nrow(qm$data) == 0) next
  q_df <- qm$data %>%
    mutate(neg_log10_padj = -log10(padj),
           pathway_label  = clean_pathway_name(pathway)) %>%
    arrange(desc(neg_log10_padj)) %>%
    slice_head(n = MAX_PER_QUAD)
  quad_color <- ORA_QUAD_COLORS_F3[[qm$name]]
  if (is.null(quad_color)) quad_color <- "grey50"
  p <- make_ora_panel(q_df, qm$name, quad_color, qm$n_hot, txt_ora)
  h <- nrow(q_df) * BAR_H_MM + MARGIN_MM
  ggsave(file.path(RPT, sprintf("panel_E_ORA_%s.pdf", qm$slug)), p,
         width = ORA_W, height = h, units = "mm", device = pdf_device)
  ggsave(file.path(RPT, sprintf("panel_E_ORA_%s.png", qm$slug)), p,
         width = ORA_W, height = h, units = "mm", dpi = 300, bg = "transparent")
}

rrho2_meta <- tibble(
  quadrant = c("Concordant_Up", "Concordant_Down",
               "Discordant_CREUp_PLADown", "Discordant_CREDown_PLAUp"),
  max_neg_log10_pvalue = round(c(max_UU, max_DD, max_UD, max_DU), 2),
  n_hotspot_genes = c(n_UU, n_DD, n_UD, n_DU),
  n_ora_pathways = c(nrow(ora_UU), nrow(ora_DD), nrow(ora_UD), nrow(ora_DU)),
  matrix_rows = nr, matrix_cols = nc, n_shared_genes = n_shared
)
write_csv(rrho2_meta, file.path(DAT, "panel_E", "rrho2_summary.csv"))

message("F03 Panel E done")
