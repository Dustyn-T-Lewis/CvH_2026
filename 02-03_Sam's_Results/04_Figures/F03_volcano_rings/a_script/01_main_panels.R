#!/usr/bin/env Rscript
# F03 Main — 4 volcano rings (Cancer_vs_Healthy, Training_CR, Training_CRE, Training_PLA)
# 2×2 composite + NES gradient legend

setwd(rprojroot::find_rstudio_root_file())

library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)

# Sam's-Results uses the main CvH style.R (via wrapper)
source("02-03_Sam's_Results/04_Figures/shared/style.R")

# Helpers present in YvO style.R but absent from CvH style.R
if (!exists("strip_for_composite")) {
  strip_for_composite <- function(p) {
    p + labs(title = NULL, subtitle = NULL, tag = NULL) +
      theme(legend.position = "none")
  }
}
if (!exists("composite_text_sizes")) {
  composite_text_sizes <- function(comp_h_mm) {
    list(title    = pmax(6, pmin(8, round(5 + comp_h_mm / 80))),
         subtitle = pmax(4, pmin(6, round(3 + comp_h_mm / 100))),
         tag      = 8)
  }
}

# cairo_pdf requires X11 (unavailable on this system); fall back to base pdf()
get_pdf_device <- function() grDevices::pdf

# Shared volcano-ring utility lives in YvO shared/ — source via absolute path
source(file.path(dirname(rprojroot::find_rstudio_root_file()),
                 "A_YvO_2026", "04_Figures", "shared", "volcano_ring.R"))

# ── Paths ──────────────────────────────────────────────────────────────────
SAM_DEP_DIR  <- "02-03_Sam's_Results/03_DEP/c_data/04_per_contrast_results"
FGSEA_CACHE  <- "02-03_Sam's_Results/04_Figures/shared/fgsea_cache"
BASE         <- "02-03_Sam's_Results/04_Figures/F03_volcano_rings"
RPT_PNL_PNG  <- file.path(BASE, "b_reports", "main", "png", "panels")
RPT_PNL_PDF  <- file.path(BASE, "b_reports", "main", "pdf", "panels")
RPT_PDF      <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG      <- file.path(BASE, "b_reports", "main", "png")
DAT          <- file.path(BASE, "c_data")
for (d in c(RPT_PNL_PNG, RPT_PNL_PDF, RPT_PDF, RPT_PNG, DAT))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Panel canvas (89 mm each, 2×2 into 178×180 mm composite)
VW <- 89; VH <- 89

CTRS <- c("Cancer_vs_Healthy", "Training_CR", "Training_CRE", "Training_PLA")

# ── Build dep_df (wide format with contrast-suffixed columns) ───────────────
# volcano_ring.R expects: logFC_<ctr>, P.Value_<ctr>, pi_score_<ctr>
dep_list <- lapply(CTRS, function(ctr) {
  d <- read_csv(file.path(SAM_DEP_DIR, paste0(ctr, ".csv")),
                show_col_types = FALSE)
  d |>
    select(gene, average_intensity,
           !!paste0("logFC_",    ctr) := logFC,
           !!paste0("P.Value_",  ctr) := P.Value,
           !!paste0("pi_score_", ctr) := pi_score)
})

dep_df <- dep_list[[1]]
for (i in 2:length(dep_list))
  dep_df <- dep_df |> left_join(dep_list[[i]], by = c("gene", "average_intensity"))

# ── Build fgsea_all (long format combining all per-contrast RDS) ────────────
fgsea_all <- bind_rows(lapply(CTRS, function(ctr) {
  rds <- file.path(FGSEA_CACHE, paste0(ctr, "_fgsea.rds"))
  stopifnot(file.exists(rds))
  readRDS(rds)
}))

# ── Helper: DEP + pathway stats subtitle string ─────────────────────────────
dbs_used <- c("Hallmark", "GO Slim", "GO:BP", "KEGG", "Reactome")

contrast_stats <- function(ctr) {
  pi_col <- paste0("pi_score_", ctr)
  n_dep  <- if (pi_col %in% names(dep_df)) sum(dep_df[[pi_col]] < 0.05, na.rm = TRUE) else 0
  ctr_rows <- fgsea_all[fgsea_all$contrast == ctr & fgsea_all$database %in% dbs_used, ]
  sprintf("%d DEPs (Π < 0.05)  |  %d / %d pathways (FDR < 0.05)",
          n_dep, sum(ctr_rows$padj < 0.05, na.rm = TRUE), sum(!is.na(ctr_rows$padj)))
}

# ── Shared panel builder ─────────────────────────────────────────────────────
build_panel <- function(ctr, title, subtitle_base, tag) {

  pi_col   <- paste0("pi_score_", ctr)
  n_dep    <- if (pi_col %in% names(dep_df)) sum(dep_df[[pi_col]] < 0.05, na.rm = TRUE) else 0
  ctr_rows <- fgsea_all[fgsea_all$contrast == ctr & fgsea_all$database %in% dbs_used, ]
  n_path   <- sum(ctr_rows$padj < 0.05, na.rm = TRUE)
  total_p  <- sum(!is.na(ctr_rows$padj))
  enriched_sub <- sprintf("%s | %d DEPs, %d pathways", subtitle_base, n_dep, n_path)

  top_terms <- select_ring_terms(fgsea_all, ctr)
  ring_data <- build_ring_with_gaps(top_terms, ctr, fgsea_all)

  # Handle sparse contrasts (Training_PLA: 1 sig pathway)
  if (nrow(ring_data) == 0) {
    message(sprintf("Panel %s (%s): no ring terms — building bare volcano", tag, ctr))
    ring_data <- tibble::tibble()
  }

  max_arc      <- if (nrow(ring_data) > 0) max(ring_data$arc_r1_var, na.rm = TRUE) else 4.8
  adaptive_gap <- 0.7 + 0.3 * (max_arc - 4.8) / 1.6

  .p <- make_volcano_ring(
    de_df              = dep_df,
    go_df              = fgsea_all,
    contrast           = ctr,
    contrast_title     = title,
    contrast_subtitle  = enriched_sub,
    ring_data_override = ring_data,
    label_size         = 2.7,
    label_gap          = adaptive_gap,
    title_size         = 5,
    subtitle_size      = 3.5,
    point_size         = 0.5,
    point_alpha        = 0.55,
    count_label_size   = scale_text(BASE_COUNT, VW) + 0.4,
    count_y_mult       = 0.75,
    count_x_mult       = 0.85,
    bg_color           = unname(CONTRAST_COLORS[ctr]),
    bg_alpha           = 0.20,
    show_legend        = FALSE
  ) + labs(tag = tag)

  # Save standalone
  fname <- tolower(gsub("[^a-zA-Z0-9]+", "_", ctr))
  ggsave(file.path(RPT_PNL_PNG, sprintf("MAIN_panel_%s_%s.png", tag, fname)),
         .p, width = VW, height = VH, units = "mm", dpi = 300)
  ggsave(file.path(RPT_PNL_PDF, sprintf("MAIN_panel_%s_%s.pdf", tag, fname)),
         .p, width = VW, height = VH, units = "mm", device = get_pdf_device())

  # Write ring terms CSV
  ring_out <- attr(.p, "ring_data")
  if (!is.null(ring_out) && nrow(ring_out) > 0) {
    dir.create(file.path(DAT, paste0("panel_", tag)), showWarnings = FALSE)
    write_csv(ring_out |> dplyr::select(-gene_list,
                -any_of(c("term_idx", "start_deg", "end_deg", "mid_deg",
                          "start_rad", "end_rad", "mid_rad", "arc_r1_var"))),
              file.path(DAT, paste0("panel_", tag), "ring_terms.csv"))
  }

  message(sprintf("F03 panel %s (%s) done  [DEPs: %d | sig pathways: %d / %d]",
                  tag, ctr, n_dep, n_path, total_p))

  strip_for_composite(.p)
}

# ── Build panels ──────────────────────────────────────────────────────────────

pA <- build_panel("Cancer_vs_Healthy", "Cancer vs Healthy",
                  "CR_Post − Healthy", "A")

pB <- build_panel("Training_CR", "Training Response (CR)",
                  "CR_Post − CR_Pre", "B")

pC <- build_panel("Training_CRE", "Training Response (CRE)",
                  "CRE_Post − CRE_Pre", "C")

pD <- build_panel("Training_PLA", "Training Response (PLA)",
                  "PLA_Post − PLA_Pre", "D")

# ── Composite (2×2 + NES legend) ─────────────────────────────────────────────

nes_legend <- build_nes_legend_bar(text_size = 5, title_size = 5,
                                   bar_margin = margin(0, 0, 0, 0, "mm"))

pA <- pA + theme(plot.margin = margin(4, -9, 0, 9, "mm"))
pB <- pB + theme(plot.margin = margin(4,  0, 0, 0, "mm"))
pC <- pC + theme(plot.margin = margin(0, -9, 0, 9, "mm"))
pD <- pD + theme(plot.margin = margin(0,  0, 0, 0, "mm"))

composite <- ((pA | pB) / (pC | pD)) + plot_layout(heights = c(1, 1))

layout_cfg <- list(
  w = 178, h = 180,
  tag_bump  =  4,
  ttl_bump  =  2,
  sub_bump  =  2,
  x_l   = 0.070,
  x_r   = 0.510,
  x_ttl = 0.040,
  y_top   = 0.960,
  y_bot   = 0.505,
  sub_off = 0.022,
  nes_x = 0.35, nes_y = 0.025, nes_w = 0.30, nes_h = 0.028
)

COMP_W <- layout_cfg$w; COMP_H <- layout_cfg$h
txt <- composite_text_sizes(COMP_H)
TAG_SZ <- txt$tag + layout_cfg$tag_bump
TTL_SZ <- txt$title + layout_cfg$ttl_bump
SUB_SZ <- txt$subtitle + layout_cfg$sub_bump
X_L <- layout_cfg$x_l; X_R <- layout_cfg$x_r; X_TTL <- layout_cfg$x_ttl
Y_TOP <- layout_cfg$y_top; Y_BOT <- layout_cfg$y_bot; SUB_OFF <- layout_cfg$sub_off

titles <- c("Cancer vs Healthy", "Training Response (CR)",
            "Training Response (CRE)", "Training Response (PLA)")
subs   <- c(contrast_stats("Cancer_vs_Healthy"), contrast_stats("Training_CR"),
            contrast_stats("Training_CRE"), contrast_stats("Training_PLA"))
tags   <- LETTERS[1:4]
xs     <- c(X_L, X_R, X_L, X_R)
ys     <- c(Y_TOP, Y_TOP, Y_BOT, Y_BOT)

composite <- ggdraw(composite)
for (i in 1:4) {
  composite <- composite +
    draw_label(tags[i], x = xs[i], y = ys[i] + 0.002, size = TAG_SZ,
               fontface = "bold", hjust = 0, vjust = 1) +
    draw_label(titles[i], x = xs[i] + X_TTL, y = ys[i], size = TTL_SZ,
               fontface = "bold", hjust = 0, vjust = 1) +
    draw_label(subs[i], x = xs[i] + X_TTL, y = ys[i] - SUB_OFF, size = SUB_SZ,
               fontface = "bold.italic", colour = "grey40", hjust = 0, vjust = 1)
}
composite <- composite +
  draw_plot(nes_legend, x = layout_cfg$nes_x, y = layout_cfg$nes_y,
            width = layout_cfg$nes_w, height = layout_cfg$nes_h)

ggsave(file.path(RPT_PDF, "MAIN_F03_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = get_pdf_device(),
       limitsize = FALSE)
ggsave(file.path(RPT_PNG, "MAIN_F03_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300, limitsize = FALSE)

message("F03 main composite done")
