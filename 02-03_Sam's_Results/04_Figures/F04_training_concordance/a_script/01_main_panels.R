#!/usr/bin/env Rscript
# F04 Main — Training Concordance (5-panel composite)
# A: Quadrant ORA  B: Pattern heatmap  C: fry barcode
# D: NES scatter   E: RRHO2
#
# 3-column geometry-aware layout (panels labeled in visual reading order):
#   Top row:    A (Quadrant ORA)  | B (Pattern heatmap, full-height)
#   Bottom row: C (fry barcode)   | D (NES scatter) | E (RRHO2)
#
# Panel B (pattern heatmap) loads AnnotationDbi (via go_slim_categories.R);
# the S4 select() masking is repaired inside go_slim_categories.R.

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

BASE    <- "02-03_Sam's_Results/04_Figures/F04_training_concordance"
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")
PNL_PNG <- file.path(RPT_PNG, "panels")
PNL_PDF <- file.path(RPT_PDF, "panels")
DAT     <- file.path(BASE, "c_data")
for (d in c(RPT_PDF, RPT_PNG, PNL_PNG, PNL_PDF, DAT))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

pdf_device <- grDevices::pdf  # force base pdf (cairo DLL fails)

# Source panels (A first, then D/C/E, then B last for stat-snapshot ordering)

message("=== F04 Composite: sourcing panels ===")

# Panel A — ORA scatter + flanking bars (complex, extracted helper)
source("02-03_Sam's_Results/04_Figures/F04_training_concordance/a_script/_panel_A_ORA.R")
n_total_A    <- nrow(scatter_df)
n_sig_A      <- n_sig
n_enrich_A   <- n_enrich
r_spear_A    <- r_spear

# Panel D — fGSEA NES scatter (config wrapper → shared engine)
cfg <- list(
  fig_id      = "F04",
  contrast_x  = "Training_CRE",
  contrast_y  = "Training_PLA",
  title       = "Pathway-Level Concordance (fGSEA)",
  axis_x_label = "Training CRE",
  axis_y_label = "Training PLA",
  subtitle_metric = "concordant",
  subtitle_interpretation = "positive rho = shared pathway regulation across supplement groups",
  ref_slope   = 1,
  panel_w     = 146,
  label_border_size = 0.25,
  rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  sig_colors     = SIG_COLORS_F2,
  sig_label_fill = SIG_LABEL_FILL_F2,
  sig_label_text = SIG_LABEL_TEXT_F2,
  sig_draw_order = c("Sig PLA only", "Sig CRE only", "Sig Both", "Interaction"),
  quadrant_defs = list(
    sig_both_label = "Sig Both",
    sig_x_label    = "Sig CRE only",
    sig_y_label    = "Sig PLA only",
    bg_red_1  = c(0, Inf, 0, Inf),
    bg_red_2  = c(-Inf, 0, -Inf, 0),
    bg_blue_1 = c(0, Inf, -Inf, 0),
    bg_blue_2 = c(-Inf, 0, 0, Inf),
    label_tr = "Concordant Up",   color_tr = "#D6604D",
    label_tl = "Discordant",      color_tl = "#4393C3",
    label_bl = "Concordant Down", color_bl = "#D6604D",
    label_br = "Discordant",      color_br = "#4393C3",
    metric_count_fn = function(q1, q2, q3, q4) q1 + q3
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
    "Protein Folding"                   = "Protein Folding",
    "Spermatogenesis"                   = "Spermatogenesis",
    "Adipogenesis"                      = "Adipogenesis",
    "Heme Metabolism"                   = "Heme Metab.",
    "Ketone Metabolism"                 = "Ketone Metab.",
    "Detoxification"                    = "Detoxification"
  )
)
source("02-03_Sam's_Results/04_Figures/shared/comparison_panels/panel_D_nes_scatter.R")
n_pw_D       <- nrow(fgsea_wide)
n_sig_pw_D   <- n_total_sig
rho_D        <- as.numeric(nes_cor_all$estimate)
rho_lo_D     <- nes_ci_all[1]
rho_hi_D     <- nes_ci_all[2]
pw_conc_D    <- pw_conc_frac

# Panel C — fry rotation test (config wrapper → shared engine)
ROW_H_fry <- 0.078
cfg <- list(
  fig_id           = "F04",
  contrast_source  = "Training_CRE",
  contrast_test    = "Training_PLA",
  test_contrast_expr = "supp_PLA_T2 - supp_PLA_T1",
  x_axis_label     = "Tr. PLA",
  set_prefix       = "tcre",
  expected_up      = "Up",
  expected_down    = "Down",
  driving_up_sign  = "pos",
  driving_dn_sign  = "neg",
  has_circularity  = FALSE,
  dep_combined_csv = "02-03_Sam's_Results/03_DEP/c_data/F04_combined_CREvPLA.csv",
  up_title_fmt = "Tr.(CRE)-Up DEPs (Pi < 0.05, n = %d) → Tr.(PLA) ranked t",
  dn_title_fmt = "Tr.(CRE)-Down DEPs (Pi < 0.05, n = %d) → Tr.(PLA) ranked t%s",
  fig_color    = unname(DIR_COLORS["Up"]),
  stat_corner_up = "topright",
  stat_corner_dn = "bottomleft",
  ora_flank_up_label = "Concordant (Up)",
  ora_flank_dn_label = "Concordant (Down)",
  ora_supp_up_label  = "Up DEPs",
  ora_supp_dn_label  = "Down DEPs",
  ora_supp_title     = "Leading-Edge ORA: fry Driving Proteins",
  ora_supp_subtitle  = "Hypergeometric ORA on concordant driving proteins | top 3 per set",
  label_map = c(
    "Dynein Recruitment To The Kinetochore"                  = "Dynein-Kinetochore Rec.",
    "Microtubule-Based Movement"                             = "MT-Based Movement",
    "Regulation Of Mitotic Cell Cycle Phase Transition"       = "Reg. Mitotic Cell Cycle Trans.",
    "Negative Regulation Of Nuclear Division"                 = "Neg. Reg. Nuclear Division",
    "Integrin Cell Surface Interactions"                      = "Integrin Cell Surface Int.",
    "Scavenging By Class A Receptors"                        = "Class A Rec. Scav.",
    "Degradation Of The Extracellular Matrix"                = "ECM Degradation",
    "Collagen Chain Trimerization"                           = "Collagen Trimerization",
    "Respiratory Electron Transport"                          = "Respiratory ETC",
    "ATP Synthesis Coupled Electron Transport"                = "ATP Synth. (ETC)",
    "Non Integrin Membrane Ecm Interactions"                  = "Non-integrin ECM",
    "Formation Of The Dystrophin Glycoprotein Complex Dgc"   = "Dystrophin Complex",
    "Cargo Recognition For Clathrin Mediated Endocytosis"    = "Clathrin Cargo Rec.",
    "Class A Receptor Scavenging"                            = "Class A Receptor Scav.",
    "Collagen Trimerization"                                  = "Collagen Trimerization",
    "Hcmv Early Events"                                       = "HCMV Early Events"
  ),
  force_inside_labels = NULL,
  long_label_mode     = "wrap",
  title        = "fry Gene-Set Rotation Test: Training Concordance",
  subtitle_fmt = "Rotation-based set test (exact GSEA analogue, accounts for inter-gene correlation) | dupCor = %.3f | n = %d proteins",
  panel_w     = 178,
  rpt_png     = PNL_PNG, rpt_pdf = PNL_PDF,
  rpt_sup_png = file.path(BASE, "b_reports", "supp", "png", "panels"),
  rpt_sup_pdf = file.path(BASE, "b_reports", "supp", "pdf", "panels"),
  dat         = DAT
)
source("02-03_Sam's_Results/04_Figures/shared/comparison_panels/panel_C_fry.R")
cor_imp_C    <- cor_imp
n_all_C      <- n_all

# Panel E — RRHO2 (config wrapper → shared engine)
cfg <- list(
  fig_id     = "F04",
  t_col_1    = "t_Training_CRE",
  t_col_2    = "t_Training_PLA",
  dep_combined_csv = "02-03_Sam's_Results/03_DEP/c_data/F04_combined_CREvPLA.csv",
  rrho_labels = c("Training (CRE)", "Training (PLA)"),
  title       = "Threshold-Free Concordance (RRHO2)",
  subtitle_fmt = "Stratified hypergeometric | %d shared genes | warm corners = concordant gene regulation | No MTC (Cahill et al. 2018)",
  axis_label_1 = expression("Training (CRE) rank"~(Up %->% Down)),
  axis_label_2 = expression("Training (PLA) rank"~(Up %->% Down)),
  quadrant_labels = list(
    UU = "Concordant Up", DD = "Concordant Down",
    UD = "Discordant (CREUp PLADn)", DU = "Discordant (CREDn PLAUp)"
  ),
  hotspot_export_names = list(UU = "UU", DD = "DD", UD = "UD", DU = "DU"),
  ora_min_size = 15,
  ora_quadrant_names = list(
    UU = "Concordant Up", DD = "Concordant Down",
    UD = "Discordant (CRE Up / PLA Down)", DU = "Discordant (CRE Down / PLA Up)"
  ),
  ora_grouped = list(
    file_1_quads = c("ora_UU", "ora_DD"),
    file_2_quads = c("ora_UD", "ora_DU")
  ),
  ora_colors = ORA_QUAD_COLORS_F2,
  summary_quadrant_names = list(
    UU = "Concordant_Up",       UU_slug = "concordant_up",
    DD = "Concordant_Down",     DD_slug = "concordant_down",
    UD = "Discordant_CRE_Up_PLA_Down",  UD_slug = "discordant_cre_up",
    DU = "Discordant_CRE_Down_PLA_Up",  DU_slug = "discordant_cre_down"
  ),
  rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  supp = NULL
)
source("02-03_Sam's_Results/04_Figures/shared/comparison_panels/panel_E_rrho2.R")
n_shared_E   <- n_shared
max_UU_E     <- max_UU
n_UU_E       <- n_UU

# Panel B — Pattern heatmap (last: loads AnnotationDbi, clobbers select())
ROW_H <- 0.078
cfg <- list(
  fig_id     = "F04",
  contrast_x = "Training_CRE",
  contrast_y = "Training_PLA",
  dep_combined_csv = "02-03_Sam's_Results/03_DEP/c_data/F04_combined_CREvPLA.csv",
  title      = "Training Response Patterns",
  col_headers = c("Tr.(CRE)", "Tr.(PLA)"),
  sort_col   = "logFC_Training_CRE",
  rpt_png = PNL_PNG, rpt_pdf = PNL_PDF, dat = DAT,
  classify_fn = function(dep_df) {
    dep_df |>
      dplyr::filter(!is.na(logFC_Training_CRE), !is.na(logFC_Training_PLA)) |>
      dplyr::filter(pi_score_Training_CRE < 0.05 | pi_score_Training_PLA < 0.05) |>
      dplyr::mutate(
        quadrant = dplyr::case_when(
          logFC_Training_CRE > 0 & logFC_Training_PLA > 0 ~ "Concordant Up",
          logFC_Training_CRE < 0 & logFC_Training_PLA < 0 ~ "Concordant Down",
          TRUE ~ "Discordant"
        ),
        sig_cat = dplyr::case_when(
          pi_score_Training_CRE < 0.05 & pi_score_Training_PLA < 0.05 ~ "Both",
          pi_score_Training_CRE < 0.05 ~ "Tr.(CRE)",
          pi_score_Training_PLA < 0.05 ~ "Tr.(PLA)",
          TRUE ~ "NS"
        )
      )
  },
  QUAD_ORDER      = c("Concordant Up", "Concordant Down", "Discordant"),
  QUAD_COLORS     = c("Concordant Up" = "#B2182B", "Concordant Down" = "#2166AC",
                      "Discordant" = "#1B7837"),
  QUAD_BG         = c("Concordant Up" = "#F4D9D2", "Concordant Down" = "#D5DEEF",
                      "Discordant" = "#C8E0CD", "Tied" = "#EEEEEE"),
  ENDPOINT_COLORS = c("Concordant Up" = "#67001F", "Concordant Down" = "#053061",
                      "Discordant" = "#00441B"),
  SIG_COLORS      = c("Both" = "#2E7D32", "Tr.(CRE)" = "#E05A4E",
                      "Tr.(PLA)" = "#5DA5DA", "NS" = "grey70"),
  display_labels = c(
    "Carbohydrate & Energy Metabolism" = "Carb. & Energy Metab.",
    "Amino Acid & Cofactor Metabolism" = "AA & Cofactor\nMetab."
  ),
  col_header_colors = c(
    CONTRAST_COLORS["Training_CRE"] %||% "#E05A4E",
    CONTRAST_COLORS["Training_PLA"] %||% "#5DA5DA"
  ),
  x_sig                = 0.27,
  tile_w               = 0.80,
  bar_scale            = 0.46,
  bar_ref_width        = 14,
  key_y_base           = ROW_H * 15.5,
  key_dy               = ROW_H * 3.8,
  key_x_sig            = NULL,
  protein_count_x_mult = 7.5,
  count_tick_y_label   = ROW_H * 2.6,
  count_ticks_max      = 16,
  count_tick_filter    = function(df) dplyr::filter(df, !val %in% c(15, 25)),
  sig_cats       = c("Tr.(CRE)", "Tr.(PLA)", "Both"),
  sig_cat_labels = c("Sig CRE", "Sig PLA", "Sig Both"),
  inset_legend   = FALSE   # rendered separately below panel B in composite (see quad_legend)
)
source("02-03_Sam's_Results/04_Figures/shared/comparison_panels/panel_B_pattern_heatmap.R")

# Restore RPT paths (shared engines clobber RPT_PDF/RPT_PNG to panels subdir)
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")

# Quadrant legend (horizontal) — overlaid below panel B in the composite via draw_plot.
nudge_idx3 <- 0.30   # shift entire 3rd item ("Discordant") right (~3.6mm physical)
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
  theme(plot.background = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "mm"))

# Composite layout (3-column, geometry-aware)
COMP_W     <- 420
COMP_H     <- 310
PRINT_SCALE <- 380 / 178
TAG_SZ     <- round(10 * PRINT_SCALE * 0.85)
TTL_SZ     <- round(10 * PRINT_SCALE * 0.85)
SUB_SZ     <- round(7 * PRINT_SCALE * 0.85)

# Stat-snapshot title strings
ttl_A <- "Quadrant ORA (Concordance)"
sub_A <- sprintf("N = %d | %d DEPs (Pi) | %d enriched (FDR) | rho = %.2f",
                 n_total_A, n_sig_A, n_enrich_A, r_spear_A)
ttl_B <- "Protein-to-Pathway"
sub_B <- sprintf("%d proteins | %d pathways", n_total, n_pw)
ttl_C <- "fry: Concordance"
sub_C <- sprintf("n = %d | dupCor = %.3f", n_all_C, cor_imp_C)
ttl_D <- "Pathway Concordance"
sub_D <- sprintf("rho = %.2f | %.0f%% concordant",
                 rho_D, pw_conc_D * 100)
ttl_E <- "RRHO2 Concordance"
sub_E <- sprintf("%d genes | max %d", n_shared_E, n_UU_E)

# Patchwork layout
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
           theme(plot.margin = margin(1, 0, 8, -1, "mm"))   # 181mm panel width

fig <- wrap_elements(full = composite) +
       wrap_elements(full = pB) +
       wrap_elements(full = pC_fry) +
       wrap_elements(full = pD) +
       wrap_elements(full = pE_heat) +
       plot_layout(design = layout,
                   widths  = rep(1, 14),
                   heights = c(6.5, rep(10, 6), 4, 4.5, rep(12, 6)))

# Panel tags, titles, subtitles via cowplot
X_A <- 0.005;  X_B <- 0.549;  X_C <- 0.012;  X_D <- 0.406;  X_E <- 0.693
X_TTL      <- 0.030
TAG_DY     <- -0.002
SUB_OFFSET <- 0.020

Y_A <- 0.984;  Y_B <- 0.984
Y_C <- 0.512;  Y_D <- 0.511;  Y_E <- 0.511

composite_final <- ggdraw(fig) +
  draw_label("A",    x = X_A,           y = Y_A - TAG_DY,       size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_A,  x = X_A + X_TTL,   y = Y_A,                size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_A,  x = X_A + X_TTL,   y = Y_A - SUB_OFFSET,   size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("B",    x = X_B,           y = Y_B - TAG_DY,       size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_B,  x = X_B + X_TTL,   y = Y_B,                size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_B,  x = X_B + X_TTL,   y = Y_B - SUB_OFFSET,   size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("C",    x = X_C,           y = Y_C - TAG_DY,       size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_C,  x = X_C + X_TTL,   y = Y_C,                size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_C,  x = X_C + X_TTL,   y = Y_C - SUB_OFFSET,   size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("D",    x = X_D,           y = Y_D - TAG_DY,       size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_D,  x = X_D + X_TTL,   y = Y_D,                size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_D,  x = X_D + X_TTL,   y = Y_D - SUB_OFFSET,   size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("E",    x = X_E,           y = Y_E - TAG_DY,       size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_E,  x = X_E + X_TTL,   y = Y_E,                size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_E,  x = X_E + X_TTL,   y = Y_E - SUB_OFFSET,   size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_plot(quad_legend, x = 0.64, y = 0.524, width = 0.30, height = 0.045)

ggsave(file.path(RPT_PDF, "MAIN_F04_composite.pdf"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_PNG, "MAIN_F04_composite.png"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message("F04 composite (5-panel, 3-column layout) saved")
