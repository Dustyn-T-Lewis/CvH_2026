#!/usr/bin/env Rscript
# F01 Main — Leg Extension (A) + DXA LBM (B) + ALM (C)
# Sam-parallel phenotype figure; mirrors YvO F01 main layout.

setwd(rprojroot::find_rstudio_root_file())

library(dplyr)
library(ggplot2)
library(ggsignif)
library(patchwork)
library(cowplot)

source("02-03_Sam's_Results/04_Figures/shared/style.R")
source("04_Figures/shared/figure_supplement_helpers.R")

# cairo_pdf requires X11 (unavailable on this system); fall back to base pdf().
get_pdf_device <- function() grDevices::pdf

BASE    <- "02-03_Sam's_Results/04_Figures/F01_phenotype"
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
PNL_PNG <- file.path(RPT_PNG, "panels")
PNL_PDF <- file.path(RPT_PDF, "panels")
DAT     <- file.path(BASE, "c_data")
for (d in c(PNL_PNG, PNL_PDF, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

TMPL <- "02-03_Sam's_Results/04_Figures/F01_phenotype/a_script/_prepost_template.R"

# ── Panel A: Leg Extension 1RM ────────────────────────────────────────────────

cfg <- list(
  pre_col = "pre_leg_ext_lbs", post_col = "post_leg_ext_lbs",
  y_label = "Leg Extension 1RM (lbs)",
  delta_label = expression(bold(Delta ~ "Leg Ext. (lbs)")),
  title = "Leg Extension 1RM", tag = "a", output_prefix = "pA",
  file_tag = "panel_A_leg_ext", audit_file = "panel_A_leg_ext.csv",
  file_prefix = "MAIN", rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT)
source(TMPL)

# ── Panel B: DXA Lean Body Mass ───────────────────────────────────────────────

cfg <- list(
  pre_col = "pre_LBM_kg", post_col = "post_LBM_kg",
  y_label = "DXA LBM (kg)",
  delta_label = expression(bold(Delta ~ "DXA LBM (kg)")),
  title = "DXA Lean Body Mass", tag = "b", output_prefix = "pB",
  file_tag = "panel_B_dxa_lbm", audit_file = "panel_B_dxa_lbm.csv",
  file_prefix = "MAIN", rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT)
source(TMPL)

# ── Panel C: Appendicular Lean Mass ──────────────────────────────────────────

cfg <- list(
  pre_col = "pre_ALM_kg", post_col = "post_ALM_kg",
  y_label = "ALM (kg)",
  delta_label = expression(bold(Delta ~ "ALM (kg)")),
  title = "Appendicular Lean Mass", tag = "c", output_prefix = "pC",
  file_tag = "panel_C_alm", audit_file = "panel_C_alm.csv",
  file_prefix = "MAIN", rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT)
source(TMPL)

# ── composite helper (mirrors YvO) ───────────────────────────────────────────
.strip <- function(p) p + labs(title = NULL, subtitle = NULL, tag = NULL) +
  theme(legend.position = "none")

.comp_txt <- function(h) list(
  title    = pmax(6, pmin(8, round(5 + h / 80))),
  subtitle = pmax(4, pmin(6, round(3 + h / 100))),
  tag      = 8
)

# ── single-column layout (85 × 125 mm) ───────────────────────────────────────
sc_cfg <- list(
  w = 85, h = 125,
  tag_x   = 0.02,
  ttl_x   = 0.08,
  sub_off = 0.022,
  y = c(A = 0.979, B = 0.614, C = 0.314)
)

pB_comp <- (pB_left | pB_right) + plot_layout(widths = c(0.65, 0.35))
pC_comp <- (pC_left | pC_right) + plot_layout(widths = c(0.65, 0.35))
pA_comp <- (pA_left | pA_right) + plot_layout(widths = c(0.65, 0.35))

sc <- (pA_comp / pB_comp / pC_comp) +
  plot_layout(heights = c(1.0, 0.8, 0.8)) &
  theme(plot.margin = margin(10, 2, 2, 2))

txt <- .comp_txt(sc_cfg$h)

sc <- ggdraw(sc) +
  draw_label("A", x = sc_cfg$tag_x, y = sc_cfg$y["A"],
             size = txt$tag, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pA_title, x = sc_cfg$ttl_x, y = sc_cfg$y["A"],
             size = txt$title, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pA_subtitle, x = sc_cfg$ttl_x, y = sc_cfg$y["A"] - sc_cfg$sub_off,
             size = txt$subtitle, fontface = "bold.italic",
             hjust = 0, vjust = 1, colour = "grey30") +
  draw_label("B", x = sc_cfg$tag_x, y = sc_cfg$y["B"],
             size = txt$tag, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pB_title, x = sc_cfg$ttl_x, y = sc_cfg$y["B"],
             size = txt$title, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pB_subtitle, x = sc_cfg$ttl_x, y = sc_cfg$y["B"] - sc_cfg$sub_off,
             size = txt$subtitle, fontface = "bold.italic",
             hjust = 0, vjust = 1, colour = "grey30") +
  draw_label("C", x = sc_cfg$tag_x, y = sc_cfg$y["C"],
             size = txt$tag, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pC_title, x = sc_cfg$ttl_x, y = sc_cfg$y["C"],
             size = txt$title, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pC_subtitle, x = sc_cfg$ttl_x, y = sc_cfg$y["C"] - sc_cfg$sub_off,
             size = txt$subtitle, fontface = "bold.italic",
             hjust = 0, vjust = 1, colour = "grey30")

ggsave(file.path(RPT_PDF, "MAIN_F01_composite_single_col.pdf"), sc,
       width = sc_cfg$w, height = sc_cfg$h, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT_PNG, "MAIN_F01_composite_single_col.png"), sc,
       width = sc_cfg$w, height = sc_cfg$h, units = "mm", dpi = 300)

# ── double-column layout (178 × 75 mm) ───────────────────────────────────────
dc_cfg <- list(
  w = 178, h = 75,
  tag_sz  = 7,
  ttl_sz  = 6,
  sub_sz  = 4,
  x_a     = 0.010,
  x_bc    = 0.360,
  ttl_off = 0.030,
  sub_off = 0.028,
  y_top   = 0.984,
  y_mid   = 0.516
)

pA_comp2 <- (pA_left | pA_right) + plot_layout(widths = c(0.65, 0.35))
pB_comp2 <- (pB_left | pB_right) + plot_layout(widths = c(0.65, 0.35))
pC_comp2 <- (pC_left | pC_right) + plot_layout(widths = c(0.65, 0.35))

dc <- wrap_elements(full = pA_comp2 & theme(plot.margin = margin(8, 2, 10, 2))) +
  wrap_elements(full = pB_comp2 & theme(plot.margin = margin(4, 2, 2, 2))) +
  wrap_elements(full = pC_comp2 & theme(plot.margin = margin(4, 2, 4, 2))) +
  plot_layout(design = "AB\nAC", widths = c(0.35, 0.65), heights = c(1, 1))

dc <- ggdraw(dc) +
  draw_label("A", x = dc_cfg$x_a, y = dc_cfg$y_top,
             size = dc_cfg$tag_sz, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pA_title, x = dc_cfg$x_a + dc_cfg$ttl_off, y = dc_cfg$y_top,
             size = dc_cfg$ttl_sz, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pA_subtitle, x = dc_cfg$x_a + dc_cfg$ttl_off,
             y = dc_cfg$y_top - dc_cfg$sub_off,
             size = dc_cfg$sub_sz, fontface = "bold.italic",
             hjust = 0, vjust = 1, colour = "grey30") +
  draw_label("B", x = dc_cfg$x_bc, y = dc_cfg$y_top,
             size = dc_cfg$tag_sz, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pB_title, x = dc_cfg$x_bc + dc_cfg$ttl_off, y = dc_cfg$y_top,
             size = dc_cfg$ttl_sz, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pB_subtitle, x = dc_cfg$x_bc + dc_cfg$ttl_off,
             y = dc_cfg$y_top - dc_cfg$sub_off,
             size = dc_cfg$sub_sz, fontface = "bold.italic",
             hjust = 0, vjust = 1, colour = "grey30") +
  draw_label("C", x = dc_cfg$x_bc, y = dc_cfg$y_mid,
             size = dc_cfg$tag_sz, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pC_title, x = dc_cfg$x_bc + dc_cfg$ttl_off, y = dc_cfg$y_mid,
             size = dc_cfg$ttl_sz, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pC_subtitle, x = dc_cfg$x_bc + dc_cfg$ttl_off,
             y = dc_cfg$y_mid - dc_cfg$sub_off,
             size = dc_cfg$sub_sz, fontface = "bold.italic",
             hjust = 0, vjust = 1, colour = "grey30")

ggsave(file.path(RPT_PDF, "MAIN_F01_composite.pdf"), dc,
       width = dc_cfg$w, height = dc_cfg$h, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT_PNG, "MAIN_F01_composite.png"), dc,
       width = dc_cfg$w, height = dc_cfg$h, units = "mm", dpi = 300)

message("F01 main composites done")
