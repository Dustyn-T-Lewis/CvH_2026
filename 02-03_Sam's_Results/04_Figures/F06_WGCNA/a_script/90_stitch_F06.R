# F06_WGCNA — Composite stitcher.
# Sources 01_main_panels.R (if not already run), then assembles the two main
# panel PNGs into the MAIN_F06_composite.{pdf,png}.
#
# Run from A_CvH_2026/ root after 00_run_wgcna.R and 01_main_panels.R.

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(cowplot)
  library(png)
  library(grid)
})

source("04_Figures/shared/style.R")
# Override pdf device — cairo_pdf unavailable on this machine
pdf_device <- grDevices::pdf

BASE     <- "02-03_Sam's_Results/04_Figures/F06_WGCNA"
RPT_PDF  <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG  <- file.path(BASE, "b_reports", "main", "png")
PANEL_PN <- file.path(RPT_PNG, "panels")
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)

# pdf_device already set above

read_panel <- function(file, dir = PANEL_PN) {
  path <- file.path(dir, file)
  if (!file.exists(path)) stop("Missing panel PNG: ", path)
  rasterGrob(readPNG(path), interpolate = TRUE)
}

pA_grob <- read_panel("MAIN_panel_A_heatmap.png")
pB_grob <- read_panel("MAIN_panel_B_eigengene.png")

# Layout constants (mm canvas)
COMP_W <- 470; COMP_H <- 290
TAG_SZ      <- 14
TITLE_SZ    <- 12
SUBTITLE_SZ <- 9

composite_final <- ggdraw(xlim = c(0, 1), ylim = c(0, 1)) +
  theme(plot.background = element_rect(fill = "white", color = NA)) +
  draw_grob(pA_grob,
            x = 0.00, y = 0.02, width = 0.56, height = 0.92,
            hjust = 0, vjust = 0) +
  draw_grob(pB_grob,
            x = 0.57, y = 0.10, width = 0.42, height = 0.78,
            hjust = 0, vjust = 0) +
  draw_label("A",
             x = 0.01, y = 0.97,
             size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("WGCNA Module–Trait Associations",
             x = 0.04, y = 0.97,
             size = TITLE_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Sam CvH N=35 | LMM (BH) | Signed Pearson, power auto-selected",
             x = 0.04, y = 0.935,
             size = SUBTITLE_SZ, fontface = "bold.italic", colour = "grey40",
             hjust = 0, vjust = 1) +
  draw_label("B",
             x = 0.58, y = 0.97,
             size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Top-Module Eigengene Profiles",
             x = 0.61, y = 0.97,
             size = TITLE_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Cancer_vs_Healthy — pre/post by Group",
             x = 0.61, y = 0.935,
             size = SUBTITLE_SZ, fontface = "bold.italic", colour = "grey40",
             hjust = 0, vjust = 1)

ggsave(file.path(RPT_PDF, "MAIN_F06_composite.pdf"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT_PNG, "MAIN_F06_composite.png"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm",
       dpi = 300, limitsize = FALSE)

message("F06_WGCNA composite saved: MAIN_F06_composite.{pdf,png}")
