#!/usr/bin/env Rscript
# F05 Main — Cancer Recovery Reversal (5-panel composite)
# A: Quadrant ORA  B: Pattern heatmap  C: fry barcode
# D: NES scatter   E: RRHO2

setwd(rprojroot::find_rstudio_root_file())

library(dplyr); library(tidyr); library(tibble)
library(stringr); library(readr); library(ggplot2)
library(patchwork); library(cowplot)

source("02-03_Sam's_Results/04_Figures/shared/style.R")

BASE    <- "02-03_Sam's_Results/04_Figures/F05_recovery_reversal"
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
PNL_PNG <- file.path(RPT_PNG, "panels")
PNL_PDF <- file.path(RPT_PDF, "panels")
DAT     <- file.path(BASE, "c_data")
for (d in c(RPT_PDF, RPT_PNG, PNL_PNG, PNL_PDF, DAT))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

pdf_device <- grDevices::pdf  # force base pdf (cairo DLL fails)

message("=== F05 Composite: sourcing panels ===")

source("02-03_Sam's_Results/04_Figures/F05_recovery_reversal/a_script/_panel_A_ORA.R")
n_total_A  <- nrow(scatter_df)
n_sig_A    <- n_sig
n_enrich_A <- n_enrich
r_spear_A  <- r_spear

# Panel D
cfg <- list(
  fig_id      = "F05",
  contrast_x  = "Cancer_vs_Healthy",
  contrast_y  = "Training_CR",
  title       = "Pathway-Level Reversal (fGSEA)",
  axis_x_label = "Cancer vs Healthy",
  axis_y_label = "Training CR",
  subtitle_metric = "reversed",
  subtitle_interpretation = "negative rho = training CR opposes cancer-associated changes",
  ref_slope   = -1,
  panel_w     = 146,
  label_border_size = 0.25,
  rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  sig_colors     = SIG_COLORS_F3,
  sig_label_fill = SIG_LABEL_FILL_F3,
  sig_label_text = SIG_LABEL_TEXT_F3,
  sig_draw_order = c("Sig CR only", "Sig CvH only", "Sig Both"),
  quadrant_defs = list(
    sig_both_label = "Sig Both",
    sig_x_label    = "Sig CvH only",
    sig_y_label    = "Sig CR only",
    bg_blue_1 = c(0, Inf, -Inf, 0),
    bg_blue_2 = c(-Inf, 0, 0, Inf),
    bg_red_1  = c(0, Inf, 0, Inf),
    bg_red_2  = c(-Inf, 0, -Inf, 0),
    label_tr = "Exacerbated",  color_tr = "#D6604D",
    label_tl = "Reversed",     color_tl = "#4393C3",
    label_bl = "Exacerbated",  color_bl = "#D6604D",
    label_br = "Reversed",     color_br = "#4393C3",
    metric_count_fn = function(q1, q2, q3, q4) q2 + q4
  ),
  display_overrides = c(
    "Unfolded Protein Response"         = "UPR",
    "Ribosome Biogenesis"               = "Ribo Bio",
    "Amino Acid Metabolism"             = "AA Metabolism",
    "Fatty Acid Metabolism"             = "FA Metabolism",
    "Fatty Acid Beta Oxidation"         = "FA Beta-Oxidation",
    "Plasma Membrane Protein Loc."      = "PM Protein Loc.",
    "Cytoplasmic Translation"           = "Cytoplasmic Transl.",
    "Mitochondrial Organization"        = "Mito Org.",
    "Precursor Metabolites & Energy"    = "Precursor Metab. & Energy",
    "Mitochondrial Transport"           = "Mito Transport",
    "Mitochondrial Protein Import"      = "Mito Protein Import",
    "Mitochondrial Protein Degradation" = "Mito Protein Deg.",
    "Extracellular Matrix Organization" = "ECM Org.",
    "Heme Metabolism"                   = "Heme Metab.",
    "Ketone Metabolism"                 = "Ketone Metab.",
    "Peroxisome"                        = "Peroxisome",
    "Muscle System"                     = "Muscle System"
  )
)
source("02-03_Sam's_Results/04_Figures/shared/comparison_panels/panel_D_nes_scatter.R")
n_pw_D    <- nrow(fgsea_wide)
n_sig_pw_D <- n_total_sig
rho_D     <- as.numeric(nes_cor_all$estimate)
rho_lo_D  <- nes_ci_all[1]
rho_hi_D  <- nes_ci_all[2]
pw_rev_D  <- pw_rev_frac

# Panel C
ROW_H_fry <- 0.078
cfg <- list(
  fig_id           = "F05",
  contrast_source  = "Cancer_vs_Healthy",
  contrast_test    = "Training_CR",
  test_contrast_expr = "supp_CRE_T2 - supp_CRE_T1",
  x_axis_label     = "Tr. CR",
  set_prefix       = "cvh",
  expected_up      = "Down",
  expected_down    = "Up",
  driving_up_sign  = "neg",
  driving_dn_sign  = "pos",
  has_circularity  = TRUE,
  dep_combined_csv = "02-03_Sam's_Results/03_DEP/c_data/F05_combined_CvHvTCR.csv",
  up_title_fmt = "CvH-Up DEPs (Pi < 0.05, n = %d) → Tr.(CR) ranked t",
  dn_title_fmt = "CvH-Down DEPs (Pi < 0.05, n = %d) → Tr.(CR) ranked t%s",
  fig_color    = unname(CONTRAST_COLORS["Cancer_vs_Healthy"] %||% "#4CAF50"),
  stat_corner_up = "bottomleft",
  stat_corner_dn = "topright",
  ora_flank_up_label = "Reversed (Up→Down)",
  ora_flank_dn_label = "Reversed (Down→Up)",
  ora_supp_up_label  = "CvH-Up (Reversed)",
  ora_supp_dn_label  = "CvH-Down (Reversed)",
  ora_supp_title     = "Leading-Edge ORA: fry Driving Proteins (Reversal)",
  ora_supp_subtitle  = "Hypergeometric ORA on reversal-driving proteins | top 3 per set",
  label_map = c(
    "Amino Acid Catabolic Process"                           = "AA Catabol.",
    "Fatty Acid Catabolic Process"                           = "FA Catabol.",
    "Amino Acid Metabolic Process"                           = "AA Metabolic Process",
    "Establishment Or Maintenance Of Cell Polarity"          = "Cell Polarity",
    "Generation Of Precursor Metabolites And Energy"         = "Precursor Metab. & Energy",
    "Protein Localization To Plasma Membrane"                = "Plasma Membrane Protein Loc.",
    "Organic Acid Catabolic Process"                         = "Organic Acid Catabolism",
    "Membraneless Organelle Assembly"                        = "Membraneless Org. Assembly"
  ),
  force_inside_labels = c("AA Catabol.", "FA Catabol."),
  long_label_mode     = "truncate",
  title        = "fry Gene-Set Rotation Test: Cancer Recovery Reversal",
  subtitle_fmt = "Rotation-based set test (exact GSEA analogue) | Circularity r = %.3f | dupCor = %.3f | n = %d proteins",
  panel_w     = 178,
  rpt_png     = PNL_PNG, rpt_pdf = PNL_PDF,
  rpt_sup_png = file.path(BASE, "b_reports", "supp", "png", "panels"),
  rpt_sup_pdf = file.path(BASE, "b_reports", "supp", "pdf", "panels"),
  dat         = DAT
)
source("02-03_Sam's_Results/04_Figures/shared/comparison_panels/panel_C_fry.R")
cor_imp_C <- cor_imp
n_all_C   <- n_all
circ_r_C  <- circ_r

# Panel E
cfg <- list(
  fig_id     = "F05",
  t_col_1    = "t_Cancer_vs_Healthy",
  t_col_2    = "t_Training_CR",
  dep_combined_csv = "02-03_Sam's_Results/03_DEP/c_data/F05_combined_CvHvTCR.csv",
  rrho_labels = c("Cancer vs Healthy", "Training CR"),
  title       = "Threshold-Free Reversal (RRHO2)",
  subtitle_fmt = "Stratified hypergeometric | %d shared genes | warm off-diagonal = training reverses cancer changes | No MTC (Cahill et al. 2018)",
  axis_label_1 = expression("Cancer vs Healthy rank"~(Up %->% Down)),
  axis_label_2 = expression("Training CR rank"~(Up %->% Down)),
  quadrant_labels = list(
    UU = "Exacerbated Up",
    DD = "Exacerbated Down",
    UD = "Reversed (CvHUp CRDn)",
    DU = "Reversed (CvHDn CRUp)"
  ),
  hotspot_export_names = list(
    UU = "Exacerbated_Up",
    DD = "Exacerbated_Down",
    UD = "Reversed_CvH_Up_CR_Down",
    DU = "Reversed_CvH_Down_CR_Up"
  ),
  ora_min_size = 10,
  ora_quadrant_names = list(
    UU = "Exacerbated Up",
    DD = "Exacerbated Down",
    UD = "Reversed (CvH Up / CR Down)",
    DU = "Reversed (CvH Down / CR Up)"
  ),
  ora_grouped = list(
    file_1_quads     = c("ora_UD", "ora_DU"),
    file_2_quads     = c("ora_UU", "ora_DD"),
    note_if_empty_2  = "No pathways enriched in exacerbation quadrants (padj<0.05)"
  ),
  ora_colors = ORA_QUAD_COLORS_F3,
  summary_quadrant_names = list(
    UU = "Exacerbated_Up",                  UU_slug = "exacerbated_up",
    DD = "Exacerbated_Down",                DD_slug = "exacerbated_down",
    UD = "Reversed_CvH_Up_CR_Down",         UD_slug = "reversed_cvh_up",
    DU = "Reversed_CvH_Down_CR_Up",         DU_slug = "reversed_cvh_down"
  ),
  rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  supp = list(
    rpt_png        = file.path(BASE, "b_reports", "supp", "png", "panels"),
    rpt_pdf        = file.path(BASE, "b_reports", "supp", "pdf", "panels"),
    ora_bar_title  = "Enriched Pathways by Reversal Quadrant",
    ora_quad_order = c("Reversed (CvH Up / CR Down)",
                       "Reversed (CvH Down / CR Up)",
                       "Exacerbated Up", "Exacerbated Down"),
    ora_quad_short = c(
      "Reversed (CvH Up / CR Down)"   = "Reversed\n(Up → Down)",
      "Reversed (CvH Down / CR Up)"   = "Reversed\n(Down → Up)",
      "Exacerbated Up"                 = "Exacerbated Up",
      "Exacerbated Down"               = "Exacerbated Down"
    )
  )
)
source("02-03_Sam's_Results/04_Figures/shared/comparison_panels/panel_E_rrho2.R")
n_shared_E <- n_shared
max_rev_E  <- max(max_UD, max_DU)
n_rev_E    <- if (max_UD >= max_DU) n_UD else n_DU

# Panel B
ROW_H <- 0.078
cfg <- list(
  fig_id     = "F05",
  contrast_x = "Cancer_vs_Healthy",
  contrast_y = "Training_CR",
  dep_combined_csv = "02-03_Sam's_Results/03_DEP/c_data/F05_combined_CvHvTCR.csv",
  title      = "Cancer Recovery Reversal Patterns",
  col_headers = c("CvH", "Tr.(CR)"),
  sort_col   = "logFC_Cancer_vs_Healthy",
  rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  classify_fn = function(dep_df) {
    dep_df |>
      dplyr::filter(!is.na(logFC_Cancer_vs_Healthy), !is.na(logFC_Training_CR)) |>
      dplyr::filter(pi_score_Cancer_vs_Healthy < 0.05 | pi_score_Training_CR < 0.05) |>
      dplyr::mutate(
        quadrant = dplyr::case_when(
          logFC_Cancer_vs_Healthy > 0 & logFC_Training_CR < 0 ~ "Reversed Up",
          logFC_Cancer_vs_Healthy < 0 & logFC_Training_CR > 0 ~ "Reversed Down",
          TRUE ~ "Non-reversed"
        ),
        sig_cat = dplyr::case_when(
          pi_score_Cancer_vs_Healthy < 0.05 & pi_score_Training_CR < 0.05 ~ "Both",
          pi_score_Cancer_vs_Healthy < 0.05 ~ "CvH",
          pi_score_Training_CR < 0.05       ~ "Tr.(CR)",
          TRUE ~ "NS"
        )
      )
  },
  QUAD_ORDER      = c("Reversed Up", "Reversed Down", "Non-reversed"),
  QUAD_COLORS     = c("Reversed Up" = "#B2182B", "Reversed Down" = "#2166AC",
                      "Non-reversed" = "#1B7837"),
  QUAD_BG         = c("Reversed Up" = "#F4D9D2", "Reversed Down" = "#D5DEEF",
                      "Non-reversed" = "#C8E0CD", "Tied" = "#EEEEEE"),
  ENDPOINT_COLORS = c("Reversed Up" = "#67001F", "Reversed Down" = "#053061",
                      "Non-reversed" = "#00441B"),
  SIG_COLORS      = c("Both" = "#2E7D32", "CvH" = "#E05A4E",
                      "Tr.(CR)" = "#5DA5DA", "NS" = "grey70"),
  display_labels = c(
    "Carbohydrate & Energy Metabolism" = "Carb. & Energy Metab.",
    "Amino Acid & Cofactor Metabolism" = "AA & Cofactor\nMetab."
  ),
  col_header_colors = c(
    CONTRAST_COLORS["Cancer_vs_Healthy"] %||% "#E05A4E",
    CONTRAST_COLORS["Training_CR"] %||% "#5DA5DA"
  ),
  bg_extend_right      = 0.05,
  bar_scale            = 0.20,
  bar_ref_width        = 32,
  key_y_base           = ROW_H * 15.5,
  key_dy               = ROW_H * 3.8,
  key_x_sig            = NULL,
  protein_count_x_mult = 15,
  count_tick_y_label   = ROW_H * 2.6,
  count_tick_filter    = function(df) dplyr::filter(df, val != 15),
  sig_cats       = c("Tr.(CR)", "CvH", "Both"),
  sig_cat_labels = c("Sig CR", "Sig CvH", "Sig Both")
)
source("02-03_Sam's_Results/04_Figures/shared/comparison_panels/panel_B_pattern_heatmap.R")

RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")

nudge_idx3 <- 0.15
quad_legend <- ggplot(inset_quad_df) +
  geom_rect(aes(xmin = (as.integer(quadrant) - 1) * 3.5 +
                       (as.integer(quadrant) == 3) * nudge_idx3,
                xmax = (as.integer(quadrant) - 1) * 3.5 + 0.7 +
                       (as.integer(quadrant) == 3) * nudge_idx3,
                ymin = -0.35, ymax = 0.35),
            fill = inset_quad_df$bg_color, color = "black", linewidth = 0.5) +
  geom_rect(aes(xmin = (as.integer(quadrant) - 1) * 3.5 + 0.10 +
                       (as.integer(quadrant) == 3) * nudge_idx3,
                xmax = (as.integer(quadrant) - 1) * 3.5 + 0.60 +
                       (as.integer(quadrant) == 3) * nudge_idx3,
                ymin = -0.15, ymax = 0.15),
            fill = inset_quad_df$bar_color, color = "black", linewidth = 0.3) +
  geom_text(aes(x = (as.integer(quadrant) - 1) * 3.5 + 0.85 +
                    (as.integer(quadrant) == 3) * nudge_idx3,
                y = 0, label = as.character(quadrant)),
            hjust = 0, size = 3.5, fontface = "bold", color = "grey15") +
  coord_cartesian(xlim = c(0, 10.5), ylim = c(-0.7, 0.7), clip = "off") +
  theme_void() +
  theme(plot.background = element_blank(), panel.background = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "mm"))

COMP_W <- 420; COMP_H <- 310
PRINT_SCALE <- 380 / 178
TAG_SZ <- round(10 * PRINT_SCALE * 0.85)
TTL_SZ <- round(10 * PRINT_SCALE * 0.85)
SUB_SZ <- round(7 * PRINT_SCALE * 0.85)

ttl_A <- "Quadrant ORA (Reversal)"
sub_A <- sprintf("N = %d | %d DEPs (Pi) | %d enriched (FDR) | rho = %.2f",
                 n_total_A, n_sig_A, n_enrich_A, r_spear_A)
ttl_B <- "Protein-to-Pathway"
sub_B <- sprintf("%d proteins | %d pathways", n_total, n_pw)
ttl_C <- "fry: Reversal"
sub_C <- sprintf("n = %d | dupCor = %.3f", n_all_C, cor_imp_C)
ttl_D <- "Pathway Reversal"
sub_D <- sprintf("rho = %.2f | %.0f%% reversed", rho_D, pw_rev_D * 100)
ttl_E <- "RRHO2 Reversal"
sub_E <- sprintf("%d genes | max %d", n_shared_E, n_rev_E)

layout <- paste(
  "##############",
  "AAAAAAAABBBBBB",
  "AAAAAAAABBBBBB",
  "AAAAAAAABBBBBB",
  "AAAAAAAABBBBBB",
  "AAAAAAAABBBBBB",
  "AAAAAAAABBBBBB",
  "##############",
  "##############",
  "CCCCCCDDDDEEEE",
  "CCCCCCDDDDEEEE",
  "CCCCCCDDDDEEEE",
  "CCCCCCDDDDEEEE",
  "CCCCCCDDDDEEEE",
  "CCCCCCDDDDEEEE",
  sep = "\n"
)

composite <- composite + plot_annotation(theme = theme(plot.margin = margin(-2.5, -1, -2.5, -1, "mm")))
pC_fry  <- pC_fry + plot_annotation(theme = theme(plot.margin = margin(3, 5, 0, 0, "mm")))
pD      <- pD + theme(plot.margin = margin(-2.8, 5, 2.8, -5, "mm"))
pE_heat <- pE_heat + theme(plot.margin = margin(-2.1, -0.2, 3.4, -3.5, "mm"),
                           axis.title = element_text(face = "bold", size = 8))
pB <- pB + coord_cartesian(xlim = c(-0.25, X_BAR_MAX + 1.75),
                           ylim = c(BAR_YMAX + ROW_H * 6.5, -ROW_H * 0.05),
                           expand = FALSE) +
           theme(plot.margin = margin(1, -28, 8, -14, "mm"))

fig <- wrap_elements(full = composite) +
       wrap_elements(full = pB) +
       wrap_elements(full = pC_fry) +
       wrap_elements(full = pD) +
       wrap_elements(full = pE_heat) +
       plot_layout(design = layout,
                   widths  = rep(1, 14),
                   heights = c(6.5, rep(10, 6), 4, 4.5, rep(12, 6)))

X_A <- 0.005;  X_B <- 0.549;  X_C <- 0.012;  X_D <- 0.406;  X_E <- 0.693
X_TTL <- 0.030; TAG_DY <- -0.002; SUB_OFFSET <- 0.020
Y_A <- 0.984; Y_B <- 0.984; Y_C <- 0.512; Y_D <- 0.511; Y_E <- 0.511

composite_final <- ggdraw(fig) +
  draw_label("A",   x = X_A,         y = Y_A - TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_A, x = X_A + X_TTL, y = Y_A,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_A, x = X_A + X_TTL, y = Y_A - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("B",   x = X_B,         y = Y_B - TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_B, x = X_B + X_TTL, y = Y_B,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_B, x = X_B + X_TTL, y = Y_B - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("C",   x = X_C,         y = Y_C - TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_C, x = X_C + X_TTL, y = Y_C,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_C, x = X_C + X_TTL, y = Y_C - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("D",   x = X_D,         y = Y_D - TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_D, x = X_D + X_TTL, y = Y_D,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_D, x = X_D + X_TTL, y = Y_D - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("E",   x = X_E,         y = Y_E - TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_E, x = X_E + X_TTL, y = Y_E,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_E, x = X_E + X_TTL, y = Y_E - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_plot(quad_legend, x = 0.64, y = 0.524, width = 0.30, height = 0.045)

ggsave(file.path(RPT_PDF, "MAIN_F05_composite.pdf"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_PNG, "MAIN_F05_composite.png"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message("F05 composite (5-panel, 3-column layout) saved")
