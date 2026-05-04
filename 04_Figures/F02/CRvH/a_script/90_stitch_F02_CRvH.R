# F02 CRvH — Composite: PCA + logFC Density
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(patchwork)
  library(ggplot2)
  library(png)
  library(grid)
})

RPT <- "04_Figures/F02/CRvH/b_reports"
pdf_device <- get_pdf_device()

read_panel_png <- function(filename) {
  path <- file.path(RPT, filename)
  if (!file.exists(path)) stop("Missing panel: ", path)
  rasterGrob(readPNG(path), interpolate = TRUE)
}

pC <- read_panel_png("panel_C_pca_MAIN.png")
pD <- read_panel_png("panel_D_logfc_density_MAIN.png")

wrap_panel <- function(grob) {
  ggplot() + annotation_custom(grob) +
    theme_void() + theme(plot.margin = margin(2, 2, 2, 2))
}

composite <- (wrap_panel(pC) | wrap_panel(pD)) +
  plot_layout(widths = c(0.55, 0.45))

COMP_W <- 260; COMP_H <- 120

ggsave(file.path(RPT, "F02_CRvH_proteome_MAIN.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "F02_CRvH_proteome_MAIN.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm",
       dpi = 300, limitsize = FALSE)

cat("F02 CRvH composite saved\n")
