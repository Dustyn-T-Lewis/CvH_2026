# F08 WGCNA — Figure-specific style
# Wraps shared/style.R + adds WGCNA-specific palettes for CvH
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

# ── WGCNA-specific palettes ──
# Light modules (need dark text on bar fills)
LIGHT_MODULES <- c("yellow", "pink", "white", "lightyellow", "lightgreen")

# Hull palette (Dark2-derived, for pathway hulls on hub networks)
HULL_PALETTE <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                  "#66A61E", "#E6AB02", "#A6761D", "#666666")

# ── GS phenotype mapping ──
# Default: Cancer_vs_Healthy is the primary phenotype for GS in CvH
# Override per-module via gs_phenotype_choices.csv if it exists
DEFAULT_GS_PHENO <- "cancer_binary"  # 1 = Cancer (CR), 0 = Healthy
DEFAULT_GS_LABEL <- "Cancer Status"
