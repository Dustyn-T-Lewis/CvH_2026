# Figure 6 — Composite WGCNA Figure (CvH)
# Reads _MAIN PNGs from b_reports/ and stitches into single composite
# Layout: A (top-left, module-trait heatmap)
#         B (top-right, per-module triptych)
#         C (bottom, hub protein networks)
# Generates: F06_wgcna_MAIN.pdf/.png

setwd(here::here())
source("04_Figures/F06/a_script/style.R")

pacman::p_load(patchwork, png, grid)

RPT <- "04_Figures/F06/b_reports"

RPT_PDF       <- file.path(RPT, "main", "pdf")

RPT_PNG       <- file.path(RPT, "main", "png")

RPT_SUPP_PDF  <- file.path(RPT, "supp", "pdf")

RPT_SUPP_PNG  <- file.path(RPT, "supp", "png")
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

dir.create(RPT_PDF,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PNG, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

message("Stitching F06 composite figure...")

# --- Read PNGs ---
read_panel <- function(filename) {
  path <- file.path(RPT_PNG, filename)
  if (!file.exists(path)) {
    warning(sprintf("Missing panel: %s", path))
    return(grid::textGrob(sprintf("MISSING:\n%s", filename),
                          gp = gpar(col = "red", fontsize = 12)))
  }
  img <- readPNG(path)
  rasterGrob(img, interpolate = TRUE)
}

panel_A <- read_panel("panel_A_module_trait_MAIN.png")
panel_B <- read_panel("panel_B_triptych_MAIN.png")
panel_C <- read_panel("panel_C_hub_network_MAIN.png")

# --- Composite layout ---
# Top row: A (heatmap) + B (triptych), equal width
# Bottom row: C (hub networks), full width
top_row <- (wrap_elements(panel_A) + labs(tag = "A")) +
           (wrap_elements(panel_B) + labs(tag = "B")) +
           plot_layout(widths = c(1, 1.2))

bottom_row <- wrap_elements(panel_C) + labs(tag = "C")

composite <- top_row / bottom_row +
  plot_layout(heights = c(1, 0.9)) +
  plot_annotation(
    title    = "Figure 6: WGCNA Co-expression Network Analysis",
    subtitle = "Cancer Recovery vs Healthy",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(face = "italic", size = 11, hjust = 0.5,
                                   color = "grey40"),
      plot.tag      = element_text(face = "bold", size = 18)
    )
  )

# --- Save ---
W <- 400; H <- 500

ggsave(file.path(RPT_PDF, "F06_wgcna_MAIN.pdf"), composite,
       width = W, height = H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT_PNG, "F06_wgcna_MAIN.png"), composite,
       width = W, height = H, units = "mm",
       dpi = 300, limitsize = FALSE)

message("  F06 composite figure saved")
message(sprintf("  Outputs:\n    %s\n    %s",
                file.path(RPT_PDF, "F06_wgcna_MAIN.pdf"),
                file.path(RPT_PNG, "F06_wgcna_MAIN.png")))
