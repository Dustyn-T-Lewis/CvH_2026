# F03 CRvH — Composite: A|B top, C|D bottom (2x2)
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(patchwork)
  library(ggplot2)
  library(png)
  library(grid)
})

RPT <- "04_Figures/F03/CRvH/b_reports"

RPT_PDF       <- file.path(RPT, "main", "pdf")

RPT_PNG       <- file.path(RPT, "main", "png")

RPT_SUPP_PDF  <- file.path(RPT, "supp", "pdf")

RPT_SUPP_PNG  <- file.path(RPT, "supp", "png")
dir.create(RPT_PDF,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PNG, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

read_panel_png <- function(filename) {
  path <- file.path(RPT_PNG, filename)
  if (!file.exists(path)) stop("Missing panel: ", path)
  rasterGrob(readPNG(path), interpolate = TRUE)
}

pA <- read_panel_png("panel_A_dep_counts_MAIN.png")
pB <- read_panel_png("panel_B_barcode_MAIN.png")
pC <- read_panel_png("panel_C_fgsea_MAIN.png")
pD <- read_panel_png("panel_D_fgsea_MAIN.png")

wrap_panel <- function(grob) {
  ggplot() + annotation_custom(grob) +
    theme_void() + theme(plot.margin = margin(2, 2, 2, 2))
}

composite <- (wrap_panel(pA) | wrap_panel(pB)) /
             (wrap_panel(pC) | wrap_panel(pD)) +
  plot_layout(heights = c(0.40, 0.60))

COMP_W <- 290; COMP_H <- 220

ggsave(file.path(RPT_PDF, "F03_CRvH_dep_pathway_MAIN.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT_PNG, "F03_CRvH_dep_pathway_MAIN.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm",
       dpi = 300, limitsize = FALSE)

cat("F03 CRvH composite saved\n")
