#!/usr/bin/env Rscript
# F05 SUPP -- Composite stitch: Diagnostic Panels
# Sources 6 supplementary panels and composes a 3-row x 2-col grid:
#   Row 1: A (ORA dedup sensitivity)      | B (Pearson r bootstrap)
#   Row 2: C (Circularity diagnostic)     | D (Reversal threshold)
#   Row 3: E (GO Slim bars)               | F (fry leading edge)
# Output: b_reports/supp/pdf/SUPP_F05_diagnostics.{pdf,png}
# J Physiol double-column (178mm), ~270mm height.

setwd(rprojroot::find_rstudio_root_file())

source("02-03_Sam's_Results/04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(patchwork)
  library(ggplot2)
  library(cowplot)
})

# Source panels
source("02-03_Sam's_Results/04_Figures/F05_recovery_reversal/a_script/_supp_ora_dedup.R")          # -> pS_ora_dedup
source("02-03_Sam's_Results/04_Figures/F05_recovery_reversal/a_script/_supp_r_bootstrap.R")        # -> pS_r_boot
source("02-03_Sam's_Results/04_Figures/F05_recovery_reversal/a_script/_supp_fry_circularity.R")    # -> pS_circ
source("02-03_Sam's_Results/04_Figures/F05_recovery_reversal/a_script/_supp_reversal_threshold.R") # -> pS_thresh
source("02-03_Sam's_Results/04_Figures/F05_recovery_reversal/a_script/_supp_goslim_bars.R")        # -> pS_goslim
source("02-03_Sam's_Results/04_Figures/F05_recovery_reversal/a_script/_supp_fry_leading.R")        # -> pS_fry_lead

BASE    <- "02-03_Sam's_Results/04_Figures/F05_recovery_reversal"
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png")
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
pdf_device <- grDevices::pdf  # force base pdf (cairo DLL fails)

COMP_W <- 260    # wider — long titles need room across two columns
COMP_H <- 310    # taller for breathing room between rows
CTS    <- composite_text_sizes(COMP_H)
CTS$title    <- CTS$title + 1      # +1pt for readability
CTS$subtitle <- CTS$subtitle + 0.5 # slightly larger subtitles

# Fix Panel C left clipping
pS_circ <- pS_circ + theme(plot.margin = margin(18, 8, 6, 10))

# Tighten axis title spacing on all panels
axis_fix <- theme(axis.title.y = element_text(margin = margin(0, 2, 0, 0)),
                  axis.title.x = element_text(margin = margin(2, 0, 0, 0)))

# 3-row x 2-col grid
grid <- (pS_ora_dedup | pS_r_boot) /
        (pS_circ      | pS_thresh) /
        (pS_goslim    | pS_fry_lead) &
  theme(plot.margin = margin(18, 8, 6, 5),
        axis.title = element_text(size = 7, face = "bold")) &
  axis_fix

# Tag + title + subtitle placement via cowplot
X_LEFT     <- 0.015
X_RIGHT    <- 0.525
X_TTL      <- 0.028
SUB_OFFSET <- 0.013   # subtitle gap below title

# Row Y positions — titles sit above each row's plot region
Y_R1 <- 0.993
Y_R2 <- 0.663
Y_R3 <- 0.333

composite <- ggdraw(grid) +
  # Panel A (top-left): ORA dedup sensitivity
  draw_label("A", x = X_LEFT, y = Y_R1,
             fontface = "bold", size = CTS$tag, hjust = 0, vjust = 1) +
  draw_label(pS_ora_title, x = X_LEFT + X_TTL, y = Y_R1,
             fontface = "bold", size = CTS$title, hjust = 0, vjust = 1) +
  draw_label(pS_ora_subtitle, x = X_LEFT + X_TTL, y = Y_R1 - SUB_OFFSET,
             fontface = "bold.italic", size = CTS$subtitle, colour = "grey30",
             hjust = 0, vjust = 1) +
  # Panel B (top-right): Pearson r bootstrap
  draw_label("B", x = X_RIGHT, y = Y_R1,
             fontface = "bold", size = CTS$tag, hjust = 0, vjust = 1) +
  draw_label(pS_rboot_title, x = X_RIGHT + X_TTL, y = Y_R1,
             fontface = "bold", size = CTS$title, hjust = 0, vjust = 1) +
  draw_label(pS_rboot_subtitle, x = X_RIGHT + X_TTL, y = Y_R1 - SUB_OFFSET,
             fontface = "bold.italic", size = CTS$subtitle, colour = "grey30",
             hjust = 0, vjust = 1) +
  # Panel C (mid-left): Circularity diagnostic
  draw_label("C", x = X_LEFT, y = Y_R2,
             fontface = "bold", size = CTS$tag, hjust = 0, vjust = 1) +
  draw_label(pS_circ_title, x = X_LEFT + X_TTL, y = Y_R2,
             fontface = "bold", size = CTS$title, hjust = 0, vjust = 1) +
  draw_label(pS_circ_subtitle, x = X_LEFT + X_TTL, y = Y_R2 - SUB_OFFSET,
             fontface = "bold.italic", size = CTS$subtitle, colour = "grey30",
             hjust = 0, vjust = 1) +
  # Panel D (mid-right): Reversal threshold
  draw_label("D", x = X_RIGHT, y = Y_R2,
             fontface = "bold", size = CTS$tag, hjust = 0, vjust = 1) +
  draw_label(pS_thresh_title, x = X_RIGHT + X_TTL, y = Y_R2,
             fontface = "bold", size = CTS$title, hjust = 0, vjust = 1) +
  draw_label(pS_thresh_subtitle, x = X_RIGHT + X_TTL, y = Y_R2 - SUB_OFFSET,
             fontface = "bold.italic", size = CTS$subtitle, colour = "grey30",
             hjust = 0, vjust = 1) +
  # Panel E (bottom-left): GO Slim bars
  draw_label("E", x = X_LEFT, y = Y_R3,
             fontface = "bold", size = CTS$tag, hjust = 0, vjust = 1) +
  draw_label(pS_goslim_title, x = X_LEFT + X_TTL, y = Y_R3,
             fontface = "bold", size = CTS$title, hjust = 0, vjust = 1) +
  draw_label(pS_goslim_subtitle, x = X_LEFT + X_TTL, y = Y_R3 - SUB_OFFSET,
             fontface = "bold.italic", size = CTS$subtitle, colour = "grey30",
             hjust = 0, vjust = 1) +
  # Panel F (bottom-right): fry leading edge
  draw_label("F", x = X_RIGHT, y = Y_R3,
             fontface = "bold", size = CTS$tag, hjust = 0, vjust = 1) +
  draw_label(pS_lead_title, x = X_RIGHT + X_TTL, y = Y_R3,
             fontface = "bold", size = CTS$title, hjust = 0, vjust = 1) +
  draw_label(pS_lead_subtitle, x = X_RIGHT + X_TTL, y = Y_R3 - SUB_OFFSET,
             fontface = "bold.italic", size = CTS$subtitle, colour = "grey30",
             hjust = 0, vjust = 1)

# Save
ggsave(file.path(RPT_PDF, "SUPP_F05_diagnostics.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_PNG, "SUPP_F05_diagnostics.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message(sprintf("F05 SUPP composite saved -> %s / %s", RPT_PDF, RPT_PNG))
