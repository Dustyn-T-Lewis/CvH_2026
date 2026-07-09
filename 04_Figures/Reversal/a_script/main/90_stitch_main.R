# Reversal — Composite: Cancer Recovery Reversal (Panels A–E)
# 3-column geometry-aware layout:
#   Top row:    A (Quadrant ORA, 8 cols) | B (Pattern heatmap, 6 cols)
#   Bottom row: C (fry barcode, 6 cols)  | D (NES scatter, 4 cols) | E (RRHO2, 4 cols)
# Source order: A, B, D, E, then C (C last — AnnotationDbi masking)
# Outputs: MAIN_Reversal_composite.{pdf,png} + Reversal_supplementary.xlsx

setwd(here::here())

pacman::p_load(patchwork, cowplot, ggplot2)

# -- Source panels (C last for AnnotationDbi masking) -------------------------
message("=== Reversal Composite: sourcing panels ===")

source("04_Figures/Reversal/a_script/main/panels/panel_A_ORA.R")
n_total_A  <- nrow(scatter_df)
n_sig_A    <- sum(scatter_df$is_sig)
n_enrich_A <- if (exists("all_quad_ora") && nrow(all_quad_ora) > 0) sum(all_quad_ora$significant) else 0L
r_spear_A  <- cor(scatter_df$logFC_CvH, scatter_df$logFC_TR,
                   use = "complete.obs", method = "spearman")
r_pear_A   <- cor(scatter_df$logFC_CvH[scatter_df$is_sig],
                   scatter_df$logFC_TR[scatter_df$is_sig],
                   use = "complete.obs")
# Quadrant counts
q_counts_A <- table(scatter_df$quadrant)

source("04_Figures/Reversal/a_script/main/panels/panel_B_nes_scatter.R")
n_pw_B     <- nrow(fgsea_hg)
n_sig_pw_B <- nrow(fgsea_sig)
rho_B      <- as.numeric(nes_cor_all$estimate)
rho_lo_B   <- nes_ci_all[1]
rho_hi_B   <- nes_ci_all[2]
pw_rev_B   <- rev_frac

source("04_Figures/Reversal/a_script/main/panels/panel_D_fry.R")
cor_imp_D  <- cor_imp
n_all_D    <- n_all
circ_r_D   <- circ_r
fry_up_D   <- fry_up
fry_dn_D   <- fry_dn

source("04_Figures/Reversal/a_script/main/panels/panel_E_rrho2.R")
n_shared_E <- n_shared
max_rev_E  <- max(max_UD, max_DU)
n_rev_E    <- if (max_UD >= max_DU) n_UD else n_DU
rrho_max   <- c(UU = max_UU, DD = max_DD, UD = max_UD, DU = max_DU)
rrho_n     <- c(UU = n_UU, DD = n_DD, UD = n_UD, DU = n_DU)

source("04_Figures/Reversal/a_script/main/panels/panel_C_pattern_heatmap.R")
n_total_C  <- n_total
n_pw_C     <- n_pw

RPT_PDF <- "04_Figures/Reversal/b_reports/main/pdf"
RPT_PNG <- "04_Figures/Reversal/b_reports/main/png"
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

# -- Layout constants ---------------------------------------------------------
COMP_W  <- 420
COMP_H  <- 300
TAG_SZ  <- 16
TTL_SZ  <- 13
SUB_SZ  <- 9

# -- Stat-snapshot title strings -----------------------------------------------
ttl_A <- "Quadrant ORA"
sub_A <- sprintf("N = %d | %d DEPs (\u03a0) | %d enriched | \u03c1 = %.2f",
                 n_total_A, n_sig_A, n_enrich_A, r_spear_A)
ttl_B <- "Protein-to-Pathway"
sub_B <- sprintf("%d proteins | %d pathways", n_total_C, n_pw_C)
ttl_C <- "fry Rotation Test"
sub_C <- sprintf("n = %d | dupCor = %.3f | circ r = %.3f", n_all_D, cor_imp_D, circ_r_D)
ttl_D <- "Pathway NES"
sub_D <- sprintf("\u03c1 = %.2f [%.2f, %.2f] | %.0f%% reversed",
                 rho_B, rho_lo_B, rho_hi_B, pw_rev_B * 100)
ttl_E <- "RRHO2"
sub_E <- sprintf("%d genes | peak = %.0f", n_shared_E, max_rev_E)

# -- Compose: 2-row layout with title spacers --------------------------------
layout <- paste(
  "##############",     # row 1: title spacer (A + B)
  "AAAAAAAABBBBBB",     # rows 2-7: A (8 cols) + B (6 cols)
  "AAAAAAAABBBBBB",
  "AAAAAAAABBBBBB",
  "AAAAAAAABBBBBB",
  "AAAAAAAABBBBBB",
  "AAAAAAAABBBBBB",
  "##############",     # row 8: spacer
  "##############",     # row 9: title spacer (C + D + E)
  "CCCCCCDDDDEEEE",     # rows 10-15: C=6, D=4, E=4
  "CCCCCCDDDDEEEE",
  "CCCCCCDDDDEEEE",
  "CCCCCCDDDDEEEE",
  "CCCCCCDDDDEEEE",
  "CCCCCCDDDDEEEE",
  sep = "\n"
)

# Strip all panel titles/margins for composite embedding
composite <- composite +
  plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0, "mm")))
pD_fry <- pD_fry +
  plot_annotation(theme = theme(plot.margin = margin(0, 2, 0, 0, "mm")))
pB <- pB + theme(plot.margin = margin(0, 2, 0, 0, "mm"))
pE_heat <- pE_heat +
  theme(plot.margin = margin(0, 0, 0, 0, "mm"),
        axis.title = element_text(face = "bold", size = 7))

fig <- wrap_elements(full = composite) +     # A
       wrap_elements(full = p) +              # B
       wrap_elements(full = pD_fry) +         # C
       wrap_elements(full = pB) +             # D
       wrap_elements(full = pE_heat) +        # E
       plot_layout(design = layout,
                   widths  = rep(1, 14),
                   heights = c(5, rep(10, 6), 3, 5, rep(10, 6)))

# -- Panel tags, titles, subtitles via cowplot --------------------------------
X_A <- 0.005;  X_B <- 0.555;  X_C <- 0.005;  X_D <- 0.425;  X_E <- 0.715
X_TTL      <- 0.028
TAG_DY     <- -0.002
SUB_OFFSET <- 0.018
Y_A <- 0.985;  Y_B <- 0.985
Y_C <- 0.515;  Y_D <- 0.515;  Y_E <- 0.515

composite_final <- ggdraw(fig) +
  draw_label("A",    x = X_A,         y = Y_A + TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_A,  x = X_A + X_TTL, y = Y_A,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_A,  x = X_A + X_TTL, y = Y_A - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("B",    x = X_B,         y = Y_B + TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_B,  x = X_B + X_TTL, y = Y_B,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_B,  x = X_B + X_TTL, y = Y_B - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("C",    x = X_C,         y = Y_C + TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_C,  x = X_C + X_TTL, y = Y_C,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_C,  x = X_C + X_TTL, y = Y_C - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("D",    x = X_D,         y = Y_D + TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_D,  x = X_D + X_TTL, y = Y_D,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_D,  x = X_D + X_TTL, y = Y_D - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("E",    x = X_E,         y = Y_E + TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_E,  x = X_E + X_TTL, y = Y_E,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_E,  x = X_E + X_TTL, y = Y_E - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40")

# -- Save composite ------------------------------------------------------------
ggsave(file.path(RPT_PDF, "MAIN_Reversal_composite.pdf"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_PNG, "MAIN_Reversal_composite.png"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)
message("Reversal composite (5-panel, 3-column layout) saved")

# --- Supplementary Excel workbook with summary stats -------------------------
source("04_Figures/shared/figure_supplement_helpers.R")

cat("=== Reversal supplementary workbook ===\n")

# Build summary stats sheet
summary_stats <- data.frame(
  Panel = c("A", "A", "A", "A", "A",
            "B", "B", "B",
            "C", "C", "C", "C", "C",
            "D", "D", "D", "D", "D",
            "E", "E", "E", "E"),
  Metric = c(
    "Total proteins", "DEPs (Pi < 0.05)", "Enriched pathways (FDR < 0.05)",
    "Spearman rho (all)", "Pearson r (sig)",
    "Proteins classified", "GO Slim categories", "Pathways (n >= 2)",
    "fry cancer_up direction", "fry cancer_up p-value",
    "fry cancer_down direction", "fry cancer_down p-value",
    "Circularity r(t_CvH, t_TR)",
    "Pathways (Hallmark + GO Slim)", "Significant pathways",
    "Spearman rho (all)", "95% CI", "Reversed fraction",
    "Shared genes", "Max -log10(p) reversed",
    "Hotspot genes (reversed)", "ORA pathways (reversed)"),
  Value = c(
    as.character(n_total_A), as.character(n_sig_A), as.character(n_enrich_A),
    sprintf("%.3f", r_spear_A), sprintf("%.3f", r_pear_A),
    as.character(n_total_C), as.character(n_pw_C), as.character(n_pw_C),
    as.character(fry_up_D$direction), sprintf("%.4f", fry_up_D$PValue),
    as.character(fry_dn_D$direction), sprintf("%.4f", fry_dn_D$PValue),
    sprintf("%.3f", circ_r_D),
    as.character(n_pw_B), as.character(n_sig_pw_B),
    sprintf("%.3f", rho_B), sprintf("[%.3f, %.3f]", rho_lo_B, rho_hi_B),
    sprintf("%.1f%%", pw_rev_B * 100),
    as.character(n_shared_E), sprintf("%.1f", max_rev_E),
    as.character(n_rev_E),
    as.character(nrow(ora_concordant) + nrow(ora_discordant))),
  stringsAsFactors = FALSE
)

rev_specs <- list(
  list(name = "summary_stats",               df = summary_stats),
  list(name = "panel_A_ora_quadrant",         path = "04_Figures/Reversal/c_data/panel_A/ora_quadrant.csv"),
  list(name = "panel_B_pattern_class",        path = "04_Figures/Reversal/c_data/panel_C_heatmap/pattern_classification.csv"),
  list(name = "panel_B_sankey",               path = "04_Figures/Reversal/c_data/panel_C_heatmap/sankey_links.csv"),
  list(name = "panel_B_bar",                  path = "04_Figures/Reversal/c_data/panel_C_heatmap/bar_data.csv"),
  list(name = "panel_C_fry_results",          path = "04_Figures/Reversal/c_data/panel_D_fry/fry_results_all.csv"),
  list(name = "panel_C_fry_driving",          path = "04_Figures/Reversal/c_data/panel_D_fry/driving_proteins.csv"),
  list(name = "panel_D_nes_scatter",          path = "04_Figures/Reversal/c_data/panel_B/nes_scatter.csv"),
  list(name = "panel_E_rrho2_summary",        path = "04_Figures/Reversal/c_data/panel_E/rrho2_summary.csv"),
  list(name = "panel_E_rrho2_hotspot",        path = "04_Figures/Reversal/c_data/panel_E/rrho2_hotspot_genes.csv"),
  list(name = "panel_E_rrho2_ora_concord",    path = "04_Figures/Reversal/c_data/panel_E/rrho2_ora_concordant.csv"),
  list(name = "panel_E_rrho2_ora_discord",    path = "04_Figures/Reversal/c_data/panel_E/rrho2_ora_discordant.csv")
)

build_workbook(
  "04_Figures/Reversal/c_data/Reversal_supplementary.xlsx",
  title = "Reversal Figure \u2014 Source Data",
  description = "Cancer recovery reversal diagnostics: quadrant ORA, pathway NES scatter, per-protein pattern classification, fry rotation test, RRHO2.",
  overview_df = data.frame(
    Sheet = vapply(rev_specs, `[[`, character(1), "name"),
    Description = c(
      "Summary statistics for all panels (key metrics)",
      "Panel A: ORA by reversal-quadrant scatter",
      "Panel B: per-protein reversal pattern classification",
      "Panel B: pathway-protein sankey links",
      "Panel B: per-pattern bar chart counts",
      "Panel C: fry rotation test for reversal",
      "Panel C: reversal driving proteins",
      "Panel D: NES scatter per pathway with Spearman + Fisher Z CI",
      "Panel E: RRHO2 quadrant summary (max -log10p per quadrant)",
      "Panel E: RRHO2 hotspot genes per quadrant",
      "Panel E: ORA on exacerbated quadrant genes",
      "Panel E: ORA on reversed quadrant genes"),
    stringsAsFactors = FALSE),
  sheet_specs = rev_specs
)
cleanup_after_workbook(rev_specs,
  extra_subdirs = c("04_Figures/Reversal/c_data/panel_A",
                     "04_Figures/Reversal/c_data/panel_B",
                     "04_Figures/Reversal/c_data/panel_C_heatmap",
                     "04_Figures/Reversal/c_data/panel_D_fry",
                     "04_Figures/Reversal/c_data/panel_E"))

message("Reversal main stitcher complete")
