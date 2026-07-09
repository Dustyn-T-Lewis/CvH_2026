# Reversal — Supplementary Composite (9-panel, two pages)
# Sources all 9 supp panel scripts and stitches into two composites:
#   Page 1 (diagnostics): ORA dedup, bootstrap r, circularity, threshold,
#                          GO Slim, fry leading
#   Page 2 (methods):     Melov proportion, CMap connectivity,
#                          directional asymmetry
# Panel tags added by patchwork; titles/subtitles kept on individual panels
# Output: SUPP_Reversal_diagnostics.{pdf,png}
#         SUPP_Reversal_methods.{pdf,png}
setwd(here::here())

source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(patchwork)
  library(cowplot)
})

RPT_PDF <- "04_Figures/Reversal/b_reports/supp/pdf"
RPT_PNG <- "04_Figures/Reversal/b_reports/supp/png"
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# ═══════════════════════════════════════════════════════════════════════════════
# Source panels — methods panels first (no AnnotationDbi), goslim last
# ═══════════════════════════════════════════════════════════════════════════════
message("=== Reversal SUPP: sourcing panels ===")

# Page 2 panels (methods) — source before AnnotationDbi-dependent panels
source("04_Figures/Reversal/a_script/supp/panels/SUPP_melov_proportion.R")
source("04_Figures/Reversal/a_script/supp/panels/SUPP_cmap_connectivity.R")
source("04_Figures/Reversal/a_script/supp/panels/SUPP_directional_asymmetry.R")

# Page 1 panels (diagnostics)
source("04_Figures/Reversal/a_script/supp/panels/SUPP_ora_dedup.R")
source("04_Figures/Reversal/a_script/supp/panels/SUPP_r_bootstrap.R")
source("04_Figures/Reversal/a_script/supp/panels/SUPP_fry_circularity.R")
source("04_Figures/Reversal/a_script/supp/panels/SUPP_reversal_threshold.R")
source("04_Figures/Reversal/a_script/supp/panels/SUPP_fry_leading.R")
# goslim LAST — loads AnnotationDbi which masks dplyr::select
source("04_Figures/Reversal/a_script/supp/panels/SUPP_goslim_bars.R")

# Reset paths (panel scripts overwrote them)
RPT_PDF <- "04_Figures/Reversal/b_reports/supp/pdf"
RPT_PNG <- "04_Figures/Reversal/b_reports/supp/png"

# ═══════════════════════════════════════════════════════════════════════════════
# PAGE 1: Diagnostics (original 6-panel 3x2 grid)
# ═══════════════════════════════════════════════════════════════════════════════
# Row 1: ORA dedup + Bootstrap r
# Row 2: Circularity + Threshold sensitivity
# Row 3: GO Slim bars (wider) + fry leading
supp_layout <- "
AABBBB
AABBBB
CCDDDD
CCDDDD
EEEEFF
EEEEFF
EEEEFF
"

supp_composite <- pS_ora_dedup + pS_r_boot +
                  pS_circ + pS_threshold +
                  pS_goslim + pS_fry_lead +
  plot_layout(design = supp_layout) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 14),
      plot.margin = margin(4, 4, 4, 4, "mm")))

SUPP_W <- 360; SUPP_H <- 400
ggsave(file.path(RPT_PNG, "SUPP_Reversal_diagnostics.png"), supp_composite,
       width = SUPP_W, height = SUPP_H, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(RPT_PDF, "SUPP_Reversal_diagnostics.pdf"), supp_composite,
       width = SUPP_W, height = SUPP_H, units = "mm", device = pdf_device)

message("Page 1: SUPP_Reversal_diagnostics saved")

# ═══════════════════════════════════════════════════════════════════════════════
# PAGE 2: Reversal Methods (3 stacked panels)
# ═══════════════════════════════════════════════════════════════════════════════
# Row 1: Melov proportion (stacked bar + permutation null)
# Row 2: CMap connectivity (running ES + permutation null)
# Row 3: Directional asymmetry (bars + bootstrap + threshold sensitivity)
methods_layout <- "
AAAA
BBBB
CCCC
"

methods_composite <- pS_melov + pS_cmap + pS_asym +
  plot_layout(design = methods_layout) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 14),
      plot.margin = margin(4, 4, 4, 4, "mm")))

METH_W <- 360; METH_H <- 340
ggsave(file.path(RPT_PNG, "SUPP_Reversal_methods.png"), methods_composite,
       width = METH_W, height = METH_H, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(RPT_PDF, "SUPP_Reversal_methods.pdf"), methods_composite,
       width = METH_W, height = METH_H, units = "mm", device = pdf_device)

message("Page 2: SUPP_Reversal_methods saved")
message("Reversal SUPP composites complete")
