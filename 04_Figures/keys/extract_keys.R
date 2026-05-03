# Legend Extraction Utilities
# Helpers for extracting legends from panel plots into per-figure key composites.
# Usage: source this file, then call extract_key() on any ggplot object.

library(ggplot2)
library(patchwork)
library(cowplot)
library(grid)

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

KEYS_RPT <- "04_Figures/keys/b_reports"
dir.create(KEYS_RPT, recursive = TRUE, showWarnings = FALSE)

# Restore legend stripped by legend.position="none", extract as standalone grob
extract_key <- function(plot_obj, position = "bottom", ...) {
  p_with_legend <- plot_obj +
    theme(legend.position = position,
          legend.text     = element_text(size = 8),
          legend.title    = element_text(size = 9, face = "bold"),
          legend.key.size = unit(4, "mm"),
          ...)

  tryCatch(
    cowplot::get_legend(p_with_legend),
    error = function(e) {
      message("  Warning: could not extract legend: ", conditionMessage(e))
      NULL
    }
  )
}

wrap_key <- function(legend_grob) {
  if (is.null(legend_grob)) return(plot_spacer())
  wrap_elements(legend_grob)
}

save_key <- function(plot_obj, filename, width = 280, height = 100) {
  pdf_device <- get_pdf_device()

  ggsave(file.path(KEYS_RPT, paste0(filename, ".pdf")), plot_obj,
         width = width, height = height, units = "mm",
         device = pdf_device, limitsize = FALSE)
  ggsave(file.path(KEYS_RPT, paste0(filename, ".png")), plot_obj,
         width = width, height = height, units = "mm",
         dpi = 300, limitsize = FALSE)

  message(sprintf("  Key saved: %s.pdf/png (%d x %d mm)", filename, width, height))
}
