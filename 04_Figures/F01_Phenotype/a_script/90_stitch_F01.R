# F01 — Composite Phenotype Figure Assembly
# Panels: A (Age), B (ALM), C (STS Max Power)
# Layout: A left | B/C stacked right
setwd(here::here())
source("04_Figures/shared/style.R")

pacman::p_load(patchwork, ggplot2, png, grid)

RPT <- "04_Figures/F01_Phenotype/b_reports"
pdf_device <- get_pdf_device()

read_panel_png <- function(filename) {
  path <- file.path(RPT, filename)
  if (!file.exists(path)) stop("Missing panel: ", path)
  img <- readPNG(path)
  rasterGrob(img, interpolate = TRUE)
}

pA <- read_panel_png("panel_A_age_MAIN.png")
pB <- read_panel_png("panel_B_alm_MAIN.png")
pC <- read_panel_png("panel_C_sts_MAIN.png")

wrap_panel <- function(grob) {
  ggplot() +
    annotation_custom(grob) +
    theme_void() +
    theme(plot.margin = margin(2, 2, 2, 2))
}

composite <- (wrap_panel(pA) | (wrap_panel(pB) / wrap_panel(pC))) +
  plot_layout(widths = c(0.35, 0.65))

COMP_W <- 260
COMP_H <- 180

ggsave(file.path(RPT, "F01_phenotype_MAIN.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "F01_phenotype_MAIN.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm",
       dpi = 300, limitsize = FALSE)

cat("F01 composite saved:", file.path(RPT, "F01_phenotype_MAIN.*"), "\n")
