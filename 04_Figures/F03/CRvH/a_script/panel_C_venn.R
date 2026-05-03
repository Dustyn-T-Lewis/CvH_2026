# F03/CRvH — Panel C (SUPP): DEP Direction Venn Diagrams
# Pi-score significant DEPs split by Up/Down across 2 CRvH contrasts
# Two Venn diagrams side by side: Up-regulated | Down-regulated
# Outputs: panel_C_venn_SUPP.pdf/.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F03/a_script/style.R")

library(dplyr)
library(readr)
library(VennDiagram)
library(grid)
library(png)
library(patchwork)
library(ggplot2)

DEP_FILE <- "03_DEP/c_data/03_combined_results_CRvH.csv"
RPT      <- "04_Figures/F03/CRvH/b_reports"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

CONTRASTS <- c("Cancer_vs_Healthy", "Training_CR")
dep_df    <- read_csv(DEP_FILE, show_col_types = FALSE)
pdf_device <- get_pdf_device()

PV_W <- 220
PV_H <- 120

# --- Extract Pi-score significant DEP gene lists per contrast x direction ---
get_dep_genes <- function(ctr, direction) {
  pi_col  <- paste0("pi_score_", ctr)
  lfc_col <- paste0("logFC_", ctr)
  dep_df %>%
    filter(!is.na(.data[[pi_col]]), .data[[pi_col]] < 0.05,
           if (direction == "Up") .data[[lfc_col]] > 0 else .data[[lfc_col]] <= 0) %>%
    pull(gene)
}

DISPLAY_NAMES <- unname(CTR_SHORT[CONTRASTS])
VENN_COLORS   <- unname(CONTRAST_COLORS[CONTRASTS])

up_lists <- lapply(CONTRASTS, function(ctr) get_dep_genes(ctr, "Up"))
names(up_lists) <- DISPLAY_NAMES

down_lists <- lapply(CONTRASTS, function(ctr) get_dep_genes(ctr, "Down"))
names(down_lists) <- DISPLAY_NAMES

for (i in seq_along(DISPLAY_NAMES)) {
  message(sprintf("  %s: %d Up, %d Down",
                  DISPLAY_NAMES[i], length(up_lists[[i]]), length(down_lists[[i]])))
}

# --- Render Venn to temp PNG ---
render_venn <- function(gene_lists, title_text, fill_alpha, filepath) {
  non_empty <- sapply(gene_lists, length) > 0
  if (sum(non_empty) < 2) {
    message(sprintf("  %s: fewer than 2 non-empty sets, skipping", title_text))
    return(NULL)
  }

  png(filepath, width = PV_W / 2, height = PV_H, units = "mm", res = 300)
  venn.plot <- venn.diagram(
    x = gene_lists[non_empty],
    filename = NULL,
    fill = VENN_COLORS[non_empty],
    alpha = fill_alpha,
    cat.cex = 0.8,
    cex = 0.7,
    cat.fontface = "bold",
    main = title_text,
    main.cex = 0.9,
    main.fontface = "bold",
    main.col = ifelse(grepl("Up", title_text),
                       unname(DIR_COLORS["Up"]),
                       unname(DIR_COLORS["Down"])),
    lwd = 0.5,
    cat.dist = rep(0.05, sum(non_empty))
  )
  grid.draw(venn.plot)
  dev.off()
  filepath
}

tmp_up   <- tempfile(fileext = ".png")
tmp_down <- tempfile(fileext = ".png")

render_venn(up_lists, "Up-regulated DEPs (Pi-score)", 0.35, tmp_up)
render_venn(down_lists, "Down-regulated DEPs (Pi-score)", 0.35, tmp_down)

img_up   <- readPNG(tmp_up)
img_down <- readPNG(tmp_down)

p_up   <- wrap_elements(rasterGrob(img_up, interpolate = TRUE))
p_down <- wrap_elements(rasterGrob(img_down, interpolate = TRUE))

p_venn <- p_up | p_down

ggsave(file.path(RPT, "panel_C_venn_SUPP.pdf"), p_venn,
       width = PV_W, height = PV_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_C_venn_SUPP.png"), p_venn,
       width = PV_W, height = PV_H, units = "mm", dpi = 300)

unlink(c(tmp_up, tmp_down))
message("F03/CRvH Panel C (Venn SUPP) done")
