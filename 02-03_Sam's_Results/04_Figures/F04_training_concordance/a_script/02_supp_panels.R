#!/usr/bin/env Rscript
# F04 Supp — 5-panel diagnostic composite (3+2 grid)
# A: ORA dedup sensitivity  B: Spearman rho bootstrap  C: Threshold sensitivity
# D: GO Slim distribution   E: fry leading-edge proteins
#
# Also sources the enrichment heatmap (standalone, not in composite grid).

setwd(rprojroot::find_rstudio_root_file())

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(ggplot2)
library(patchwork)
library(cowplot)

source("02-03_Sam's_Results/04_Figures/shared/style.R")
source("02-03_Sam's_Results/04_Figures/shared/pathway_utils.R")

pdf_device <- grDevices::pdf  # force base pdf (cairo DLL fails)

BASE <- "02-03_Sam's_Results/04_Figures/F04_training_concordance"

# Source enrichment heatmap first (standalone ComplexHeatmap, not in grid)
message("=== F04 SUPP: enrichment heatmap ===")
source("02-03_Sam's_Results/04_Figures/F04_training_concordance/a_script/_supp_enrichment_heatmap.R")

# Source diagnostic panels
message("=== F04 SUPP Composite: sourcing panels ===")
source("02-03_Sam's_Results/04_Figures/F04_training_concordance/a_script/_supp_ora_dedup.R")       # -> pS_ora_dedup
source("02-03_Sam's_Results/04_Figures/F04_training_concordance/a_script/_supp_rho_bootstrap.R")   # -> pS_rho_boot
source("02-03_Sam's_Results/04_Figures/F04_training_concordance/a_script/_supp_threshold_sens.R")  # -> pS_thresh
source("02-03_Sam's_Results/04_Figures/F04_training_concordance/a_script/_supp_goslim_bars.R")     # -> pS_goslim
source("02-03_Sam's_Results/04_Figures/F04_training_concordance/a_script/_supp_fry_leading.R")     # -> pS_fry_lead

# Output directories (set after panels so their vars don't persist)
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png")
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)

# Composite layout
COMP_W <- 280
COMP_H <- 275

cts <- composite_text_sizes(COMP_H)
TAG_SZ <- cts$tag
TTL_SZ <- cts$title + 1
SUB_SZ <- cts$subtitle + 0.5

# Fix individual panels before combining
pS_thresh <- pS_thresh +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = 7))

axis_fix <- theme(legend.position = "none",
                  axis.title.y = element_text(margin = margin(0, 2, 0, 0)),
                  axis.title.x = element_text(margin = margin(2, 0, 0, 0)))
pS_ora_dedup <- pS_ora_dedup + axis_fix +
  theme(plot.margin = margin(18, 10, 2, 5))
pS_rho_boot  <- pS_rho_boot + axis_fix +
  theme(plot.margin = margin(18, 10, 2, 5))
pS_thresh    <- pS_thresh + axis_fix +
  theme(plot.margin = margin(18, 14, 2, 5))
pS_goslim    <- pS_goslim + axis_fix +
  theme(plot.margin = margin(18, 10, 2, 5))
pS_fry_lead  <- pS_fry_lead + axis_fix +
  theme(plot.margin = margin(18, 10, 2, 5))

top_row <- (pS_ora_dedup | pS_rho_boot | pS_thresh) +
  plot_layout(widths = c(1, 1, 1.15), guides = "keep")
bot_row <- (pS_goslim | pS_fry_lead) +
  plot_layout(guides = "keep")

fig <- (top_row / bot_row) +
  plot_layout(heights = c(1, 1.05), guides = "keep")

# Panel tags, titles, subtitles via cowplot
X_A <- 0.080;  X_B <- 0.405;  X_C <- 0.675
X_D <- 0.020;  X_E <- 0.520
X_TTL      <- 0.022
SUB_OFFSET <- 0.016
Y_TOP <- 0.993;  Y_BOT <- 0.488

composite_final <- ggdraw(fig) +
  draw_label("A", x = X_A,         y = Y_TOP,               size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("ORA Dedup Sensitivity", x = X_A + X_TTL, y = Y_TOP, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Enriched pathways by Jaccard cutoff", x = X_A + X_TTL, y = Y_TOP - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("B", x = X_B,         y = Y_TOP,               size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Spearman rho Bootstrap", x = X_B + X_TTL, y = Y_TOP, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("1000 replicates, 95% CI", x = X_B + X_TTL, y = Y_TOP - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("C", x = X_C,         y = Y_TOP,               size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Threshold Sensitivity", x = X_C + X_TTL, y = Y_TOP, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Quadrant counts by significance criterion", x = X_C + X_TTL, y = Y_TOP - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("D", x = X_D,         y = Y_BOT,               size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("GO Slim Distribution", x = X_D + X_TTL, y = Y_BOT, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Proteins per GO Slim category by quadrant", x = X_D + X_TTL, y = Y_BOT - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("E", x = X_E,         y = Y_BOT,               size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("fry Leading-Edge Proteins", x = X_E + X_TTL, y = Y_BOT, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Top 20 by |t-stat| in Training PLA", x = X_E + X_TTL, y = Y_BOT - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40")

# Consolidated legend bar at bottom
legend_df <- data.frame(
  x = 1:3,
  fill_lab = factor(c("Concordant Up", "Concordant Down", "Discordant"),
                    levels = c("Concordant Up", "Concordant Down", "Discordant")))
legend_plot <- ggplot(legend_df, aes(x = x, y = 1, fill = fill_lab)) +
  geom_col() +
  scale_fill_manual(values = c("Concordant Up" = "#E57373",
                               "Concordant Down" = "#64B5F6",
                               "Discordant" = "#FFB74D"),
                    name = "Quadrant") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 7),
        legend.key.size = unit(4, "mm"),
        legend.spacing.x = unit(3, "mm"))

legend_grob <- cowplot::get_legend(legend_plot)
composite_final <- composite_final +
  draw_plot(legend_grob, x = 0.08, y = -0.003, width = 0.45, height = 0.035)

# Save
ggsave(file.path(RPT_PDF, "SUPP_F04_diagnostics.pdf"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_PNG, "SUPP_F04_diagnostics.png"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message(sprintf("F04 SUPP composite (5-panel) saved: %s x %s mm", COMP_W, COMP_H))
