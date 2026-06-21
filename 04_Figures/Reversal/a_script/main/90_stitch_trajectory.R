# Reversal — Companion composite: trajectory clusters (F) + residual volcano (G)
# A second reversal figure that pairs the pattern-first view (F: what reverses)
# with the residual view (G: what training did NOT fix). The 5-panel composite
# (90_stitch_main.R: scatter, NES, fry, RRHO2) stays separate to keep both readable.
# Reads the panel PNGs (run panel_F + panel_G first). Output: PNG only (cairo-less env).
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
suppressPackageStartupMessages({ library(patchwork); library(ggplot2); library(png); library(grid) })

PANELS <- "04_Figures/Reversal/b_reports/main/png/panels"
RPT    <- "04_Figures/Reversal/b_reports/main/png"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

# rebuild the two panels so the composite is always current
source("04_Figures/Reversal/a_script/main/panels/panel_F_trajectory.R")
source("04_Figures/Reversal/a_script/main/panels/panel_G_resid_volcano.R")

read_panel <- function(f) {
  p <- file.path(PANELS, f)
  if (!file.exists(p)) stop("missing panel PNG: ", p)
  ggplot() + annotation_custom(rasterGrob(readPNG(p), interpolate = TRUE)) +
    theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
}

fig <- read_panel("MAIN_panel_F_trajectory_composite.png") /
       read_panel("MAIN_panel_G_resid_volcano_composite.png") +
  plot_layout(heights = c(150, 120)) +                       # native panel heights (mm)
  plot_annotation(tag_levels = list(c("A", "B"))) &
  theme(plot.tag = element_text(size = 16, face = "bold"))

ggsave(file.path(RPT, "MAIN_Reversal_trajectory_composite.png"), fig,
       width = 220, height = 280, units = "mm", dpi = 300)
message("Reversal trajectory+residual companion composite saved")
