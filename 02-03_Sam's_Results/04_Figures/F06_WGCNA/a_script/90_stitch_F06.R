# F06_WGCNA — End-to-end stitcher (YvO-aligned structure).
# Sources 01_main_panels.R then 02_supp_panels.R, then assembles the MAIN
# composite (A | B) at YvO dimensions: 470 x 300 mm canvas.
#
# Run from A_CvH_2026/ root:
#   Rscript "02-03_Sam's_Results/04_Figures/F06_WGCNA/a_script/90_stitch_F06.R"

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(cowplot)
  library(png)
  library(grid)
})

source("02-03_Sam's_Results/04_Figures/shared/style.R")
pdf_device <- get_pdf_device()

# ── Step 1: Main panels ───────────────────────────────────────────────────────
message("=== F06 Step 1: main panels ===")
source("02-03_Sam's_Results/04_Figures/F06_WGCNA/a_script/01_main_panels.R")

# ── Step 2: Supplementary panels ─────────────────────────────────────────────
message("=== F06 Step 2: supp panels ===")
source("02-03_Sam's_Results/04_Figures/F06_WGCNA/a_script/02_supp_panels.R")

# ── Step 3: Compose MAIN composite (A heatmap | B NES scatters) ───────────────
message("=== F06 Step 3: MAIN composite ===")

BASE    <- "02-03_Sam's_Results/04_Figures/F06_WGCNA"
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
PANEL_N <- file.path(RPT_PNG, "panels")
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)

read_panel_png <- function(file, dir = PANEL_N) {
  path <- file.path(dir, file)
  if (!file.exists(path)) stop("Missing panel PNG: ", path)
  rasterGrob(readPNG(path), interpolate = TRUE)
}

pA_grob     <- read_panel_png("MAIN_panel_A_heatmap.png")
pB_grob     <- read_panel_png("MAIN_panel_B_scatters.png")
pB_leg_grob <- read_panel_png("MAIN_panel_B_legend.png")

# Layout constants — match YvO 01_main_panels.R exactly
COMP_W <- 470; COMP_H <- 300
TAG_SZ      <- 14
TITLE_SZ    <- 12
SUBTITLE_SZ <- 9

mm2x <- function(mm) mm / COMP_W
mm2y <- function(mm) mm / COMP_H

# Coordinates in mm (origin = bottom-left)
CROP_L <- 0.01; CROP_R <- 0.80
CROP_B <- 0.17; CROP_T <- 0.85
SAVE_W <- COMP_W * (CROP_R - CROP_L)
SAVE_H <- COMP_H * (CROP_T - CROP_B)

composite_final <- ggdraw(xlim = c(CROP_L, CROP_R), ylim = c(CROP_B, CROP_T)) +
  theme(plot.background = element_rect(fill = "white", color = NA)) +
  # Panel B behind (left white margin hidden by Panel A)
  draw_grob(pB_grob, x = 0.43, y = 0.22, width = 0.55, height = 0.58,
            hjust = 0, vjust = 0) +
  draw_grob(pB_leg_grob, x = 0.70, y = 0.195, width = 0.24, height = 0.04,
            hjust = 0.5, vjust = 0.5) +
  draw_grob(pA_grob, x = 0.01, y = 0, width = 0.60, height = 0.96,
            hjust = 0, vjust = 0) +
  draw_label("A", x = mm2x(15), y = mm2y(248),
             size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("WGCNA Module-Trait Associations",
             x = mm2x(58), y = mm2y(248),
             size = TITLE_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("7 modules | LMM (BH) | CRvH + CR models | r-equiv",
             x = mm2x(58), y = mm2y(243),
             size = SUBTITLE_SZ, fontface = "bold.italic", colour = "grey40",
             hjust = 0, vjust = 1) +
  draw_label("B", x = mm2x(285), y = mm2y(248),
             size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Module-Level NES Scatters",
             x = mm2x(294), y = mm2y(248),
             size = TITLE_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("fGSEA on module-member t-stat ranks",
             x = mm2x(294), y = mm2y(243),
             size = SUBTITLE_SZ, fontface = "bold.italic", colour = "grey40",
             hjust = 0, vjust = 1)

ggsave(file.path(RPT_PDF, "MAIN_F06_composite.pdf"), composite_final,
       width = SAVE_W, height = SAVE_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT_PNG, "MAIN_F06_composite.png"), composite_final,
       width = SAVE_W, height = SAVE_H, units = "mm",
       dpi = 300, limitsize = FALSE)

message("F06 composite saved: MAIN_F06_composite.{pdf,png}")
message("=== F06 stitcher complete ===")
