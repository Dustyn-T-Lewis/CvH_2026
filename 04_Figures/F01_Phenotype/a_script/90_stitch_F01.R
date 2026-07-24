# F01 Phenotype composite: Age (A, left) over ALM (B) and STS max power (C, stacked
# right). Composed from the live panel objects so the figure stays vector-sharp and
# inherits the shared theme, rather than re-importing rendered PNGs.
setwd(here::here())
source("04_Figures/shared/style.R")
pacman::p_load(patchwork, ggplot2)

A <- "04_Figures/F01_Phenotype/a_script"
source(file.path(A, "panel_A.R"))
source(file.path(A, "panel_B.R"))
source(file.path(A, "panel_C.R"))

RPT <- "04_Figures/F01_Phenotype/b_reports"
pdf_device <- get_pdf_device()

composite <- (pA | (pB / pC)) +
  plot_layout(widths = c(0.35, 0.65))

COMP_W <- 260
COMP_H <- 180
ggsave(file.path(RPT, "F01_phenotype_MAIN.pdf"), composite,
  width = COMP_W, height = COMP_H, units = "mm",
  device = pdf_device, limitsize = FALSE
)
ggsave(file.path(RPT, "F01_phenotype_MAIN.png"), composite,
  width = COMP_W, height = COMP_H, units = "mm", dpi = 300, limitsize = FALSE
)
message("F01 composite saved (live-grob)")
