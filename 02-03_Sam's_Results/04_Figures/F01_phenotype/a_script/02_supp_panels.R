#!/usr/bin/env Rscript
# F01 Supp — Chest Press (A) + Grip Strength (B) + Sit-to-Stand Power (C)
# Sam-parallel phenotype figure; mirrors YvO F01 supp layout.

setwd(rprojroot::find_rstudio_root_file())

library(patchwork)
library(cowplot)

source("02-03_Sam's_Results/04_Figures/shared/style.R")

# cairo_pdf requires X11 (unavailable on this system); fall back to base pdf().
get_pdf_device <- function() grDevices::pdf

BASE    <- "02-03_Sam's_Results/04_Figures/F01_phenotype"
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf")
PNL_PNG <- file.path(RPT_PNG, "panels")
PNL_PDF <- file.path(RPT_PDF, "panels")
DAT     <- file.path(BASE, "c_data", "supp")
for (d in c(PNL_PNG, PNL_PDF, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

TMPL <- "02-03_Sam's_Results/04_Figures/F01_phenotype/a_script/_prepost_template.R"

# ── Panel A: Chest Press 1RM ──────────────────────────────────────────────────

cfg <- list(
  pre_col = "pre_chest_press_lbs", post_col = "post_chest_press_lbs",
  y_label = "Chest Press 1RM (lbs)",
  delta_label = expression(bold(Delta ~ "Chest Press (lbs)")),
  title = "Chest Press 1RM", tag = "a", output_prefix = "pSA",
  file_tag = "panel_A_chest_press", audit_file = "panel_A_chest_press.csv",
  file_prefix = "SUPP", rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT)
source(TMPL)

# ── Panel B: Grip Strength ────────────────────────────────────────────────────

cfg <- list(
  pre_col = "pre_grip_lbs", post_col = "post_grip_lbs",
  y_label = "Grip Strength (lbs)",
  delta_label = expression(bold(Delta ~ "Grip Strength (lbs)")),
  title = "Grip Strength", tag = "b", output_prefix = "pSB",
  file_tag = "panel_B_grip", audit_file = "panel_B_grip.csv",
  file_prefix = "SUPP", rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT)
source(TMPL)

# ── Panel C: Sit-to-Stand Max Power ──────────────────────────────────────────

cfg <- list(
  pre_col = "pre_sts_max_pwr", post_col = "post_sts_max_pwr",
  y_label = "STS Max Power (W)",
  delta_label = expression(bold(Delta ~ "STS Power (W)")),
  title = "Sit-to-Stand Max Power", tag = "c", output_prefix = "pSC",
  file_tag = "panel_C_sts", audit_file = "panel_C_sts.csv",
  file_prefix = "SUPP", rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT)
source(TMPL)

# ── supp composite ────────────────────────────────────────────────────────────
composite <- (pSA_left + pSA_right +
              pSB_left + pSB_right +
              pSC_left + pSC_right) +
  plot_layout(ncol = 2, widths = c(0.65, 0.35), heights = c(1, 1, 1)) &
  theme(plot.margin = margin(6, 2, 2, 2))

COMP_W <- 85; COMP_H <- 88

.comp_txt <- function(h) list(
  title    = pmax(6, pmin(8, round(5 + h / 80))),
  subtitle = pmax(4, pmin(6, round(3 + h / 100))),
  tag      = 8
)

txt <- .comp_txt(COMP_H)
TAG_SZ <- txt$tag - 3; TTL_SZ <- txt$title - 2; SUB_SZ <- txt$subtitle - 1
Y_A <- 0.983; Y_B <- 0.660; Y_C <- 0.336
X_TAG <- 0.056; X_TTL <- 0.116; SUB_OFF <- 0.016

composite <- ggdraw(composite) +
  draw_label("A", x = X_TAG, y = Y_A, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSA_title, x = X_TTL, y = Y_A, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSA_subtitle, x = X_TTL, y = Y_A - SUB_OFF, size = SUB_SZ,
             fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30") +
  draw_label("B", x = X_TAG, y = Y_B, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSB_title, x = X_TTL, y = Y_B, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSB_subtitle, x = X_TTL, y = Y_B - SUB_OFF, size = SUB_SZ,
             fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30") +
  draw_label("C", x = X_TAG, y = Y_C, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSC_title, x = X_TTL, y = Y_C, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSC_subtitle, x = X_TTL, y = Y_C - SUB_OFF, size = SUB_SZ,
             fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30")

ggsave(file.path(RPT_PDF, "SUPP_F01_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT_PNG, "SUPP_F01_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message("F01 supp composite done")
