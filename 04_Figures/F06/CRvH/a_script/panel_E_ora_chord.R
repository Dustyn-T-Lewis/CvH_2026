# F06 CRvH Panel E: Unified ORA-DEP Chord Diagram -- CRvH Model
# Single chord connecting ORA-enriched pathways to CRvH DEPs
# (Cancer_vs_Healthy + Training_CR, Pi < 0.05 in either).
#
# Ring design (outside -> inside):
#   Track 1: contrast identity color (proteins) | pathway name (bold on colored bg)
#   Track 2: contrast strip (proteins) | empty (pathways)
#   Track 3: logFC bar (proteins) | solid unique color (pathways)
#   Ribbons: protein Track 3 -> pathway Track 3, colored by pathway (alpha=0.25)
#
# ORA: Jaccard 0.3 dedup, rrvgo on GO:BP. All surviving pathways shown.
# Direction assigned per pathway via sign(mean logFC of overlap genes).
# Legend column on right side of panel.
# ---------------------------------------------------------------------------
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F06/a_script/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
  library(circlize)
})

set.seed(42)

# -- Figure-specific constants ------------------------------------------------

FIGURE_LABEL       <- "CRvH"
RELEVANT_CONTRASTS <- c("Cancer_vs_Healthy", "Training_CR")

RPT <- "04_Figures/F06/CRvH/b_reports"
DAT <- "04_Figures/F06/CRvH/c_data"
dir.create(file.path(DAT, "panel_E"), recursive = TRUE, showWarnings = FALSE)
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

PANEL_W <- 210   # mm -- total width including legend
PANEL_H <- 180   # mm
DPI     <- 300

PW_COLORS_10 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E",
                   "#E6AB02", "#A6761D", "#666666", "#1F78B4", "#FB9A99")

# -- Step 1: Load DEP results & build pool -----------------------------------

dep_df <- read_csv("03_DEP/c_data/03_combined_results_CRvH.csv", show_col_types = FALSE)
all_genes <- unique(dep_df$gene)

pi_cols   <- paste0("pi_score_", RELEVANT_CONTRASTS)
logfc_cols <- paste0("logFC_", RELEVANT_CONTRASTS)

# DEP membership per relevant contrast
dep_pool <- dep_df %>%
  mutate(across(all_of(pi_cols), ~ !is.na(.) & . < 0.05,
                .names = "sig_{.col}")) %>%
  filter(if_any(starts_with("sig_pi_score_"), ~ .))

# Rename sig columns
sig_cols <- paste0("sig_pi_score_", RELEVANT_CONTRASTS)
names(dep_pool)[match(sig_cols, names(dep_pool))] <- paste0("sig_", RELEVANT_CONTRASTS)

message(sprintf("DEP pool (%s): %d proteins", FIGURE_LABEL, nrow(dep_pool)))
for (ctr in RELEVANT_CONTRASTS) {
  n <- sum(dep_pool[[paste0("sig_", ctr)]])
  message(sprintf("  %s: %d", ctr, n))
}

# Primary contrast = lowest pi_score among relevant contrasts
primary_df <- dep_pool %>%
  select(gene, all_of(pi_cols)) %>%
  pivot_longer(-gene, names_to = "pi_col", values_to = "pi_val") %>%
  mutate(contrast = gsub("pi_score_", "", pi_col)) %>%
  filter(!is.na(pi_val), contrast %in% RELEVANT_CONTRASTS) %>%
  group_by(gene) %>%
  slice_min(pi_val, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(gene, primary_contrast = contrast)

# Display logFC = logFC from primary contrast
display_lfc <- dep_pool %>%
  select(gene, all_of(logfc_cols)) %>%
  left_join(primary_df, by = "gene") %>%
  pivot_longer(all_of(logfc_cols), names_to = "lfc_col", values_to = "lfc_val") %>%
  mutate(lfc_contrast = gsub("logFC_", "", lfc_col)) %>%
  filter(lfc_contrast == primary_contrast) %>%
  select(gene, primary_contrast, display_logFC = lfc_val)

message(sprintf("Primary contrast distribution:"))
print(table(display_lfc$primary_contrast))

# -- Step 2: ORA with Jaccard 0.3 -------------------------------------------

pw_collection <- build_pathway_collection(min_size = 15, max_size = 500,
                                           include_goslim = FALSE,
                                           exclude_variants = TRUE)

# Split DEPs by direction for separate Up/Down ORA
up_genes   <- display_lfc %>% filter(display_logFC > 0) %>% pull(gene)
down_genes <- display_lfc %>% filter(display_logFC <= 0) %>% pull(gene)
message(sprintf("DEP direction split: %d Up, %d Down", length(up_genes), length(down_genes)))

# Helper: run ORA + rrvgo dedup on a gene set
run_ora_with_rrvgo <- function(genes, label) {
  if (length(genes) < 5) {
    message(sprintf("  %s: too few genes (%d) -- skipping", label, length(genes)))
    return(tibble())
  }
  ora <- run_ora_deduplicated(
    genes = genes, universe = all_genes,
    pathways = pw_collection, jaccard_cutoff = 0.3,
    min_size = 15, max_size = 500, padj_cutoff = 0.05
  )
  message(sprintf("  %s ORA: %d sig pathways after Jaccard dedup", label, nrow(ora)))
  if (nrow(ora) == 0) return(tibble())

  # rrvgo on GO:BP
  gobp <- ora %>% filter(database == "GO:BP")
  other <- ora %>% filter(database != "GO:BP")
  if (nrow(gobp) >= 2) {
    tryCatch({
      gobp_msig <- msigdbr::msigdbr(species = "Homo sapiens",
                                     collection = "C5", subcollection = "GO:BP")
      go_map <- gobp_msig %>%
        distinct(gs_name, gs_exact_source) %>%
        filter(gs_name %in% gobp$pathway)
      if (nrow(go_map) > 0) {
        gobp <- gobp %>% left_join(go_map, by = c("pathway" = "gs_name"))
        has_go <- gobp %>% filter(!is.na(gs_exact_source))
        if (nrow(has_go) >= 2) {
          sim_mat <- rrvgo::calculateSimMatrix(
            has_go$gs_exact_source, orgdb = "org.Hs.eg.db",
            ont = "BP", method = "Wang")
          reduced <- rrvgo::reduceSimMatrix(
            sim_mat, setNames(has_go$padj, has_go$gs_exact_source),
            threshold = 0.7, orgdb = "org.Hs.eg.db")
          parent_ids <- reduced %>% filter(go == parent) %>% pull(go) %>% unique()
          keep_pw <- has_go %>% filter(gs_exact_source %in% parent_ids) %>% pull(pathway)
          gobp <- gobp %>% filter(pathway %in% keep_pw) %>% select(-gs_exact_source)
          message(sprintf("    rrvgo: %d GO:BP -> %d parents", nrow(has_go), length(keep_pw)))
        } else {
          gobp <- gobp %>% select(-any_of("gs_exact_source"))
        }
      }
    }, error = function(e) {
      message("    rrvgo failed: ", e$message)
      gobp <<- gobp %>% select(-any_of("gs_exact_source"))
    })
  }
  bind_rows(other, gobp) %>% arrange(padj)
}

ora_up   <- run_ora_with_rrvgo(up_genes, "Up")
ora_down <- run_ora_with_rrvgo(down_genes, "Down")

# Combine: remove duplicates (keep lower padj version)
ora_final <- bind_rows(
  if (nrow(ora_up) > 0) ora_up %>% mutate(ora_direction = "Up") else tibble(),
  if (nrow(ora_down) > 0) ora_down %>% mutate(ora_direction = "Down") else tibble()
) %>%
  group_by(pathway) %>%
  slice_min(padj, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(padj)

message(sprintf("Final ORA pathways: %d (%d from Up, %d from Down)",
                nrow(ora_final),
                sum(ora_final$ora_direction == "Up", na.rm = TRUE),
                sum(ora_final$ora_direction == "Down", na.rm = TRUE)))

if (nrow(ora_final) == 0) {
  message("No significant ORA pathways -- stopping")
  quit(save = "no", status = 0)
}

# -- Step 3: Build chord data -----------------------------------------------

chord_links <- list()
for (i in seq_len(nrow(ora_final))) {
  genes <- ora_final$overlapGenes[[i]]
  if (length(genes) == 0) next
  chord_links[[i]] <- tibble(
    gene          = genes,
    pathway       = ora_final$pathway[i],
    pathway_padj  = ora_final$padj[i],
    pathway_db    = ora_final$database[i],
    pathway_size  = ora_final$size[i]
  )
}
link_df <- bind_rows(chord_links)

# Clean pathway labels
shorten_pw <- function(x, max_chars = 28) {
  abbrev <- c(
    "Degradation Of The Extracellular Matrix" = "ECM Degradation",
    "Extracellular Matrix Organization"       = "ECM Organization",
    "Collagen Biosynthesis And Modifying Enzymes" = "Collagen Biosynthesis",
    "Assembly Of Collagen Fibrils And Other Multimeric Structures" = "Collagen Assembly",
    "Integrin Cell Surface Interactions"       = "Integrin Interactions",
    "Binding And Uptake Of Ligands By Scavenger Receptors" = "Scavenger Receptor Uptake",
    "Epithelial Mesenchymal Transition"        = "EMT",
    "Cellular Response To Acid Chemical"       = "Acid Response",
    "Interleukin 4 And Interleukin 13 Signaling" = "IL-4/IL-13 Signaling",
    "Collagen Fibril Organization"             = "Collagen Fibrils",
    "Muscle Cell Differentiation"              = "Muscle Diff.",
    "Sarcomere Organization"                   = "Sarcomere Org.",
    "Striated Muscle Contraction"              = "Striated Muscle",
    "Actomyosin Structure Organization"        = "Actomyosin Org.",
    "Membraneless Organelle Assembly"          = "Membraneless Org.",
    "Muscle Contraction"                       = "Muscle Contraction",
    "Actin Filament Based Movement"            = "Actin Movement",
    "Cellular Component Assembly Involved In Morphogenesis" = "Morphogenesis Assembly",
    "Muscle Cell Development"                  = "Muscle Cell Dev."
  )
  x <- ifelse(x %in% names(abbrev), abbrev[x], x)
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 2), ".."), x)
}

link_df <- link_df %>%
  mutate(
    pathway_full  = clean_pathway_name(pathway),
    pathway_label = shorten_pw(clean_pathway_name(pathway))
  )

# Merge with protein metadata
link_df <- link_df %>%
  left_join(display_lfc, by = "gene")

# Proteins in the chord = union of overlap genes
proteins <- sort(unique(link_df$gene))
n_prot   <- length(proteins)
message(sprintf("Proteins in chord: %d (from %d DEPs)", n_prot, nrow(dep_pool)))

# Drop proteins with missing contrast assignment
n_missing <- sum(is.na(link_df$primary_contrast))
if (n_missing > 0) {
  warning(sprintf("%d proteins in ORA overlap have no primary contrast -- dropping", n_missing))
  link_df <- link_df %>% filter(!is.na(primary_contrast))
}

# -- Step 4: Direction assignment per pathway --------------------------------

pw_direction <- link_df %>%
  group_by(pathway, pathway_full, pathway_label, pathway_padj, pathway_db,
           pathway_size) %>%
  summarize(
    mean_logFC = mean(display_logFC, na.rm = TRUE),
    n_overlap  = n_distinct(gene),
    .groups    = "drop"
  ) %>%
  mutate(direction = ifelse(mean_logFC > 0, "Up", "Down")) %>%
  arrange(pathway_padj) %>%
  group_by(pathway_full) %>%
  slice_min(pathway_padj, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(pathway_label = make.unique(pathway_label, sep = " ")) %>%
  arrange(pathway_padj)

message("\nPathway direction assignments:")
for (i in seq_len(nrow(pw_direction))) {
  cat(sprintf("  %2d. [%s] %-35s padj=%.1e  mean_lfc=%+.2f  %s  n=%d\n",
    i, pw_direction$pathway_db[i],
    pw_direction$pathway_label[i],
    pw_direction$pathway_padj[i],
    pw_direction$mean_logFC[i],
    pw_direction$direction[i],
    pw_direction$n_overlap[i]))
}

# -- Step 5: Order proteins & pathways ---------------------------------------

pw_direction <- pw_direction %>% arrange(direction, pathway_padj)
pathways_ordered <- pw_direction$pathway_label

# Assign each protein to its primary pathway (lowest padj among pathways it belongs to)
gene_primary_pw <- link_df %>%
  group_by(gene) %>%
  slice_min(pathway_padj, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(gene, primary_pw = pathway_label)

lfc_lookup <- setNames(display_lfc$display_logFC, display_lfc$gene)
ctr_lookup <- setNames(display_lfc$primary_contrast, display_lfc$gene)

proteins_ordered <- gene_primary_pw %>%
  mutate(
    pw_rank  = match(primary_pw, pathways_ordered),
    contrast = ctr_lookup[gene],
    abs_lfc  = abs(lfc_lookup[gene])
  ) %>%
  arrange(pw_rank, contrast, desc(abs_lfc)) %>%
  pull(gene) %>%
  unique()

# -- Step 6: Save chord data CSVs -------------------------------------------

write_csv(
  link_df %>% select(gene, pathway_label, pathway_padj, pathway_db,
                      primary_contrast, display_logFC),
  file.path(DAT, "panel_E", "chord_CRvH_combined.csv")
)

write_csv(
  pw_direction %>% select(pathway_label, pathway_padj, pathway_db,
                            n_overlap, mean_logFC, direction),
  file.path(DAT, "panel_E", "chord_CRvH_pathways.csv")
)
message("\nChord data saved to c_data/panel_E/")

# -- Step 7: Render chord diagram --------------------------------------------

logfc_color <- function(lfc) {
  lfc_clamped <- pmin(pmax(lfc, -2), 2)
  colorRamp2(c(-2, 0, 2), c("#4393C3", "white", "#D6604D"))(lfc_clamped)
}

direction_color <- function(mean_lfc) {
  lfc_clamped <- pmin(pmax(mean_lfc, -2), 2)
  colorRamp2(c(-2, 0, 2), c("#0D47A1", "#F5F5F5", "#B71C1C"))(lfc_clamped)
}

# Pathway identity colors
n_pw <- nrow(pw_direction)
pw_col_map <- setNames(PW_COLORS_10[seq_len(n_pw)], pathways_ordered)

# Pathway direction color map
pw_dir_col <- setNames(
  vapply(pw_direction$mean_logFC, function(m) direction_color(m), character(1)),
  pathways_ordered
)

# Sector sizes
pw_sizes <- setNames(
  pmax(pw_direction$n_overlap, 5),
  pathways_ordered
)
all_sectors <- c(proteins_ordered, pathways_ordered)
sector_sizes <- setNames(
  c(rep(3, n_prot), pw_sizes[pathways_ordered]),
  all_sectors
)

# Dynamic gaps
prot_gap <- max(0.2, 15 / n_prot)
gap_vec <- c(
  rep(prot_gap, n_prot - 1), 12,
  rep(2, max(n_pw - 1, 0)), 12
)

# Dynamic font
gene_cex <- max(0.28, min(0.55, 14 / n_prot))
pw_cex   <- 0.50

# Contrast color lookup
ctr_col_lookup <- CONTRAST_COLORS[ctr_lookup[proteins_ordered]]

# --- Render function ---
render_chord <- function(filepath, ext) {
  if (ext == "pdf") {
    pdf(filepath, width = PANEL_W / 25.4, height = PANEL_H / 25.4)
  } else {
    png(filepath, width = PANEL_W, height = PANEL_H, units = "mm", res = DPI,
        bg = "white")
  }

  # Layout: chord 72%, legend 28%
  layout(matrix(c(1, 2), ncol = 2), widths = c(0.72, 0.28))

  # --- Panel 1: Chord diagram ---
  par(mar = c(0.5, 0.5, 2.5, 0.5))

  circos.clear()
  circos.par(
    gap.after      = gap_vec,
    start.degree   = 90,
    clock.wise     = TRUE,
    cell.padding   = c(0, 0, 0, 0),
    track.margin   = c(0.005, 0.005)
  )
  circos.initialize(
    factors = factor(all_sectors, levels = all_sectors),
    xlim    = cbind(rep(0, length(all_sectors)), sector_sizes)
  )

  # TRACK 1 (outermost): gene names | direction-colored arc + name + black border
  circos.track(ylim = c(0, 1), track.height = 0.14, bg.border = NA,
    panel.fun = function(x, y) {
      s <- CELL_META$sector.index
      if (s %in% proteins_ordered) {
        circos.points(CELL_META$xcenter, 0.12, pch = 15,
                      col = ctr_col_lookup[s], cex = gene_cex * 2.0)
        circos.text(CELL_META$xcenter, 0.55, s,
                    cex = gene_cex, facing = "clockwise", niceFacing = TRUE,
                    font = 3)
      } else if (s %in% pathways_ordered) {
        circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1,
                    col = pw_dir_col[s], border = "black", lwd = 0.5)
        arc_w <- CELL_META$cell.xlim[2] - CELL_META$cell.xlim[1]
        label_len <- nchar(s)
        dyn_cex <- min(pw_cex, max(0.18, pw_cex * arc_w / max(label_len * 0.7, 8)))
        pw_idx <- match(s, pathways_ordered)
        mean_lfc <- pw_direction$mean_logFC[pw_idx]
        txt_col <- if (abs(mean_lfc) > 0.3) "white" else "grey20"
        circos.text(CELL_META$xcenter, 0.5, s,
                    cex = dyn_cex, font = 2, col = txt_col,
                    facing = "bending.inside", niceFacing = TRUE)
      }
    })

  # TRACK 2: contrast strip (proteins) | empty (pathways)
  circos.track(ylim = c(0, 1), track.height = 0.040, bg.border = NA,
    panel.fun = function(x, y) {
      s <- CELL_META$sector.index
      if (s %in% proteins_ordered) {
        col <- ctr_col_lookup[s]
        circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1,
                    col = col, border = NA)
      }
    })

  # TRACK 3 (innermost): logFC bar (proteins) | solid unique color (pathways)
  circos.track(ylim = c(0, 1), track.height = 0.055, bg.border = NA,
    panel.fun = function(x, y) {
      s <- CELL_META$sector.index
      if (s %in% proteins_ordered) {
        col <- logfc_color(lfc_lookup[s])
        circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1,
                    col = col, border = NA)
      } else if (s %in% pathways_ordered) {
        circos.rect(CELL_META$cell.xlim[1], 0, CELL_META$cell.xlim[2], 1,
                    col = pw_col_map[s], border = NA)
      }
    })

  # RIBBONS: protein -> pathway, colored by PATHWAY unique color
  link_unique <- link_df %>%
    filter(gene %in% proteins_ordered, pathway_label %in% pathways_ordered) %>%
    distinct(gene, pathway_label)

  for (i in seq_len(nrow(link_unique))) {
    gene_i <- link_unique$gene[i]
    pw_i   <- link_unique$pathway_label[i]

    pw_genes <- link_unique %>%
      filter(pathway_label == pw_i) %>%
      pull(gene)
    pw_genes <- pw_genes[pw_genes %in% proteins_ordered]
    idx <- which(pw_genes == gene_i)
    if (length(idx) == 0) next

    pw_xlim <- get.cell.meta.data("cell.xlim", sector.index = pw_i,
                                   track.index = 3)
    slot_w <- (pw_xlim[2] - pw_xlim[1]) / length(pw_genes)

    ribbon_col <- adjustcolor(pw_col_map[pw_i], alpha.f = 0.30)

    tryCatch(
      circos.link(
        gene_i, c(0, 1),
        pw_i, c(pw_xlim[1] + (idx - 1) * slot_w,
                pw_xlim[1] + idx * slot_w),
        col    = ribbon_col,
        border = NA,
        rou1   = get.cell.meta.data("cell.bottom.radius",
                                     sector.index = gene_i, track.index = 3),
        rou2   = get.cell.meta.data("cell.bottom.radius",
                                     sector.index = pw_i, track.index = 3)
      ),
      error = function(e) NULL
    )
  }

  title(main = paste(FIGURE_LABEL, "DEP - Pathway Chord"),
        cex.main = 1.3, font.main = 2, line = 0.5)
  circos.clear()


  # --- Panel 2: Legends ---
  par(mar = c(2, 0.5, 2.5, 1))
  plot.new()

  y_cursor <- 0.95

  # Legend helper
  draw_legend <- function(y_start, title, labels, colors, cex_title = 0.75,
                           cex_labels = 0.65, box_h = 0.022, spacing = 0.028) {
    text(0.05, y_start, title, adj = c(0, 1), font = 2, cex = cex_title)
    y <- y_start - 0.03
    for (i in seq_along(labels)) {
      rect(0.05, y - box_h, 0.15, y, col = colors[i], border = "grey50",
           lwd = 0.3)
      text(0.18, y - box_h / 2, labels[i], adj = c(0, 0.5), cex = cex_labels)
      y <- y - spacing
    }
    return(y - 0.01)
  }

  # 1. Primary Contrast
  y_cursor <- draw_legend(
    y_cursor, "Primary Contrast",
    labels = c("Cancer vs Healthy", "Training CR"),
    colors = CONTRAST_COLORS[RELEVANT_CONTRASTS]
  )

  # 2. Protein logFC gradient bar
  n_steps <- 50
  x_left <- 0.05; x_right <- 0.55
  bar_h <- 0.018

  text(0.05, y_cursor, "Protein logFC", adj = c(0, 1), font = 2, cex = 0.75)
  y_cursor <- y_cursor - 0.03
  lfc_seq <- seq(-2, 2, length.out = n_steps)
  lfc_cols <- logfc_color(lfc_seq)
  for (j in seq_len(n_steps)) {
    x0 <- x_left + (j - 1) / n_steps * (x_right - x_left)
    x1 <- x_left + j / n_steps * (x_right - x_left)
    rect(x0, y_cursor - bar_h, x1, y_cursor, col = lfc_cols[j], border = NA)
  }
  rect(x_left, y_cursor - bar_h, x_right, y_cursor, col = NA, border = "grey50",
       lwd = 0.3)
  text(x_left, y_cursor - bar_h - 0.012, "-2", cex = 0.55, adj = c(0, 1))
  text((x_left + x_right) / 2, y_cursor - bar_h - 0.012, "0",
       cex = 0.55, adj = c(0.5, 1))
  text(x_right, y_cursor - bar_h - 0.012, "+2", cex = 0.55, adj = c(1, 1))
  y_cursor <- y_cursor - bar_h - 0.04

  # 3. Pathway Direction -- single blue-white-red gradient bar
  text(0.05, y_cursor, "Pathway Direction", adj = c(0, 1), font = 2, cex = 0.75)
  y_cursor <- y_cursor - 0.03
  dir_seq <- seq(-2, 2, length.out = n_steps)
  dir_cols <- direction_color(dir_seq)
  for (j in seq_len(n_steps)) {
    x0 <- x_left + (j - 1) / n_steps * (x_right - x_left)
    x1 <- x_left + j / n_steps * (x_right - x_left)
    rect(x0, y_cursor - bar_h, x1, y_cursor, col = dir_cols[j], border = NA)
  }
  rect(x_left, y_cursor - bar_h, x_right, y_cursor, col = NA,
       border = "grey50", lwd = 0.3)
  text(x_left, y_cursor - bar_h - 0.012, "-2 (Down)", cex = 0.50, adj = c(0, 1))
  text(x_right, y_cursor - bar_h - 0.012, "+2 (Up)", cex = 0.50, adj = c(1, 1))
  y_cursor <- y_cursor - bar_h - 0.04

  # 4. Pathway Identity -- small color swatches
  legend_pw_labels <- pw_direction$pathway_label
  y_cursor <- draw_legend(
    y_cursor, "Pathway Identity",
    labels = legend_pw_labels,
    colors = pw_col_map[pathways_ordered],
    cex_labels = 0.50,
    box_h = 0.018,
    spacing = 0.024
  )

  dev.off()
}

# --- Render both formats ---
for (ext in c("pdf", "png")) {
  outfile <- file.path(RPT, paste0("panel_E_ora_chord_MAIN.", ext))
  render_chord(outfile, ext)
  message(sprintf("Rendered: %s", outfile))
}

message(sprintf("\nDone: F06 CRvH Panel E -- %s ORA Chord (%d proteins, %d pathways)",
                FIGURE_LABEL, n_prot, n_pw))
