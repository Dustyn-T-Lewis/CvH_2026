# F04 driver: cancer-recovery reversal figure. Builds the 5-panel main composite
# (A quadrant ORA, B protein-to-pathway, C fry rotation, D pathway NES, E RRHO2)
# plus its source-data workbook, the trajectory/residual companion (F, G), and
# the two-page supplementary composite. Panels read native pipeline outputs
# through f04_data.R.

setwd(here::here())
source("04_Figures/shared/style.R")
pacman::p_load(patchwork, cowplot, ggplot2, png, grid, dplyr)

# Several panels pull in AnnotationDbi (GO.db / RRHO2), whose S4 select()/filter()
# mask dplyr's. Bind the dplyr versions in the global env so bare calls in every
# sourced panel resolve to dplyr regardless of package load order.
select <- dplyr::select
filter <- dplyr::filter

A <- "04_Figures/F04_Reversal/a_script"

# ── Main composite (A–E). Panel C sourced last: AnnotationDbi masks select(). ──
source(file.path(A, "panel_A_ORA.R"))
n_total_A <- nrow(scatter_df)
n_sig_A <- sum(scatter_df$is_sig)
n_enrich_A <- if (exists("all_quad_ora") && nrow(all_quad_ora) > 0) sum(all_quad_ora$significant) else 0L
r_spear_A <- cor(scatter_df$logFC_CvH, scatter_df$logFC_TR, use = "complete.obs", method = "spearman")
r_pear_A <- cor(scatter_df$logFC_CvH[scatter_df$is_sig], scatter_df$logFC_TR[scatter_df$is_sig], use = "complete.obs")

source(file.path(A, "panel_B_nes_scatter.R"))
n_pw_B <- nrow(fgsea_hg)
n_sig_pw_B <- nrow(fgsea_sig)
rho_B <- as.numeric(nes_cor_all$estimate)
rho_lo_B <- nes_ci_all[1]
rho_hi_B <- nes_ci_all[2]
pw_rev_B <- rev_frac

source(file.path(A, "panel_D_fry.R"))
cor_imp_D <- cor_imp
n_all_D <- n_all
circ_r_D <- circ_r
fry_up_D <- fry_up
fry_dn_D <- fry_dn

source(file.path(A, "panel_E_rrho2.R"))
n_shared_E <- n_shared
max_rev_E <- max(max_UD, max_DU)
n_rev_E <- if (max_UD >= max_DU) n_UD else n_DU

source(file.path(A, "panel_C_pattern_heatmap.R"))
n_total_C <- n_total
n_pw_C <- n_pw

RPT_PDF <- "04_Figures/F04_Reversal/b_reports/main/pdf"
RPT_PNG <- "04_Figures/F04_Reversal/b_reports/main/png"
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

COMP_W <- 420
COMP_H <- 320
TAG_SZ <- 16
TTL_SZ <- 13
SUB_SZ <- 9

ttl_A <- "Quadrant ORA"
sub_A <- sprintf(
  "N = %d | %d DEPs (Π) | %d enriched | ρ = %.2f",
  n_total_A, n_sig_A, n_enrich_A, r_spear_A
)
ttl_B <- "Pathway NES"
sub_B <- sprintf(
  "ρ = %.2f [%.2f, %.2f] | %.0f%% reversed",
  rho_B, rho_lo_B, rho_hi_B, pw_rev_B * 100
)
ttl_C <- "fry Rotation Test"
sub_C <- sprintf("n = %d | dupCor = %.3f | circ r = %.3f", n_all_D, cor_imp_D, circ_r_D)

# ORA composite spans the top; NES and fry share the bottom row. Blank rows
# reserve space for the drawn tags/titles above each region.
layout <- paste(
  "############",
  "AAAAAAAAAAAA", "AAAAAAAAAAAA", "AAAAAAAAAAAA",
  "AAAAAAAAAAAA", "AAAAAAAAAAAA", "AAAAAAAAAAAA",
  "############",
  "BBBBBBCCCCCC", "BBBBBBCCCCCC", "BBBBBBCCCCCC",
  "BBBBBBCCCCCC", "BBBBBBCCCCCC", "BBBBBBCCCCCC",
  sep = "\n"
)
SPACER_TOP <- 6
SPACER_MID <- 10

composite <- composite +
  plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0, "mm")))
pB <- pB + theme(plot.margin = margin(0, 2, 0, 0, "mm"))
pD_fry <- pD_fry +
  plot_annotation(theme = theme(plot.margin = margin(0, 2, 0, 0, "mm")))

fig <- wrap_elements(full = composite) +
  wrap_elements(full = pB) +
  wrap_elements(full = pD_fry) +
  plot_layout(
    design = layout, widths = rep(1, 12),
    heights = c(SPACER_TOP, rep(10, 6), SPACER_MID, rep(10, 6))
  )

X_A <- 0.005
X_B <- 0.005
X_C <- 0.505
X_TTL <- 0.030
TAG_DY <- -0.002
SUB_OFFSET <- 0.017
Y_top <- 0.985
Y_bot <- 0.512

tag_block <- function(d, tag, ttl, sub, x, y) {
  d +
    draw_label(tag, x = x, y = y + TAG_DY, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
    draw_label(ttl, x = x + X_TTL, y = y, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
    draw_label(sub, x = x + X_TTL, y = y - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40")
}

composite_final <- ggdraw(fig) |>
  tag_block("A", ttl_A, sub_A, X_A, Y_top) |>
  tag_block("B", ttl_B, sub_B, X_B, Y_bot) |>
  tag_block("C", ttl_C, sub_C, X_C, Y_bot)

ggsave(file.path(RPT_PDF, "MAIN_F04_composite.pdf"), composite_final,
  width = COMP_W, height = COMP_H, units = "mm", device = pdf_device
)
ggsave(file.path(RPT_PNG, "MAIN_F04_composite.png"), composite_final,
  width = COMP_W, height = COMP_H, units = "mm", dpi = 300
)
message("F04 main composite saved")

# ── Source-data workbook ──
source("04_Figures/shared/figure_supplement_helpers.R")

summary_stats <- data.frame(
  Panel = c("A", "A", "A", "A", "B", "B", "C", "C", "C", "D", "D", "D", "E", "E"),
  Metric = c(
    "Total proteins", "DEPs (Pi < 0.05)", "Enriched pathways (FDR < 0.05)", "Spearman rho (all)",
    "Proteins classified", "GO Slim categories",
    "fry cancer_up p-value", "fry cancer_down p-value", "Circularity r(t_CvH, t_TR)",
    "Significant pathways", "Spearman rho (all)", "Reversed fraction",
    "Shared genes", "Max -log10(p) reversed"
  ),
  Value = c(
    as.character(n_total_A), as.character(n_sig_A), as.character(n_enrich_A), sprintf("%.3f", r_spear_A),
    as.character(n_total_C), as.character(n_pw_C),
    sprintf("%.4f", fry_up_D$PValue), sprintf("%.4f", fry_dn_D$PValue), sprintf("%.3f", circ_r_D),
    as.character(n_sig_pw_B), sprintf("%.3f", rho_B), sprintf("%.1f%%", pw_rev_B * 100),
    as.character(n_shared_E), sprintf("%.1f", max_rev_E)
  ),
  stringsAsFactors = FALSE
)

rev_specs <- list(
  list(name = "summary_stats", df = summary_stats),
  list(name = "panel_A_ora_quadrant", path = "04_Figures/F04_Reversal/c_data/panel_A/ora_quadrant.csv"),
  list(name = "panel_B_pattern_class", path = "04_Figures/F04_Reversal/c_data/panel_C_heatmap/pattern_classification.csv"),
  list(name = "panel_B_sankey", path = "04_Figures/F04_Reversal/c_data/panel_C_heatmap/sankey_links.csv"),
  list(name = "panel_B_bar", path = "04_Figures/F04_Reversal/c_data/panel_C_heatmap/bar_data.csv"),
  list(name = "panel_C_fry_results", path = "04_Figures/F04_Reversal/c_data/panel_D_fry/fry_results_all.csv"),
  list(name = "panel_C_fry_driving", path = "04_Figures/F04_Reversal/c_data/panel_D_fry/driving_proteins.csv"),
  list(name = "panel_D_nes_scatter", path = "04_Figures/F04_Reversal/c_data/panel_B/nes_scatter.csv"),
  list(name = "panel_E_rrho2_summary", path = "04_Figures/F04_Reversal/c_data/panel_E/rrho2_summary.csv"),
  list(name = "panel_E_rrho2_hotspot", path = "04_Figures/F04_Reversal/c_data/panel_E/rrho2_hotspot_genes.csv"),
  list(name = "panel_E_rrho2_ora_concord", path = "04_Figures/F04_Reversal/c_data/panel_E/rrho2_ora_concordant.csv"),
  list(name = "panel_E_rrho2_ora_discord", path = "04_Figures/F04_Reversal/c_data/panel_E/rrho2_ora_discordant.csv")
)

build_workbook(
  "04_Figures/F04_Reversal/c_data/F04_supplementary.xlsx",
  title = "F04 Reversal — Source Data",
  description = "Cancer recovery reversal diagnostics: quadrant ORA, pathway NES scatter, per-protein pattern classification, fry rotation test, RRHO2.",
  overview_df = data.frame(
    Sheet = vapply(rev_specs, `[[`, character(1), "name"),
    Description = c(
      "Summary statistics for all panels",
      "Panel A: ORA by reversal-quadrant scatter",
      "Panel B: per-protein reversal pattern classification",
      "Panel B: pathway-protein sankey links",
      "Panel B: per-pattern bar chart counts",
      "Panel C: fry rotation test for reversal",
      "Panel C: reversal driving proteins",
      "Panel D: NES scatter per pathway with Spearman + Fisher Z CI",
      "Panel E: RRHO2 quadrant summary",
      "Panel E: RRHO2 hotspot genes per quadrant",
      "Panel E: ORA on concordant quadrant genes",
      "Panel E: ORA on discordant quadrant genes"
    ),
    stringsAsFactors = FALSE
  ),
  sheet_specs = rev_specs
)
cleanup_after_workbook(rev_specs,
  extra_subdirs = c(
    "04_Figures/F04_Reversal/c_data/panel_A", "04_Figures/F04_Reversal/c_data/panel_B",
    "04_Figures/F04_Reversal/c_data/panel_C_heatmap", "04_Figures/F04_Reversal/c_data/panel_D_fry",
    "04_Figures/F04_Reversal/c_data/panel_E"
  )
)

# ── Trajectory / residual companion (F, G) ──
source(file.path(A, "panel_F_trajectory.R"))
source(file.path(A, "panel_G_resid_volcano.R"))

PANELS <- "04_Figures/F04_Reversal/b_reports/main/png/panels"
MAIN_PNG <- "04_Figures/F04_Reversal/b_reports/main/png" # panel scripts overwrite RPT_PNG
read_panel <- function(f) {
  path <- file.path(PANELS, f)
  if (!file.exists(path)) stop("missing panel PNG: ", path)
  ggplot() +
    annotation_custom(rasterGrob(readPNG(path), interpolate = TRUE)) +
    theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0))
}
companion <- read_panel("MAIN_panel_F_trajectory_composite.png") /
  read_panel("MAIN_panel_G_resid_volcano_composite.png") +
  plot_layout(heights = c(150, 120)) +
  plot_annotation(tag_levels = list(c("A", "B"))) &
  theme(plot.tag = element_text(size = 16, face = "bold"))
ggsave(file.path(MAIN_PNG, "MAIN_F04_trajectory_companion.png"), companion,
  width = 220, height = 280, units = "mm", dpi = 300
)
message("F04 trajectory companion saved")

# ── Supplementary composites (two pages) ──
S <- file.path(A, "supp")
source(file.path(S, "SUPP_melov_proportion.R"))
source(file.path(S, "SUPP_cmap_connectivity.R"))
source(file.path(S, "SUPP_directional_asymmetry.R"))
source(file.path(S, "SUPP_ora_dedup.R"))
source(file.path(S, "SUPP_r_bootstrap.R"))
source(file.path(S, "SUPP_fry_circularity.R"))
source(file.path(S, "SUPP_reversal_threshold.R"))
source(file.path(S, "SUPP_fry_leading.R"))
source(file.path(S, "SUPP_goslim_bars.R")) # last: loads AnnotationDbi

RPT_SUPP_PDF <- "04_Figures/F04_Reversal/b_reports/supp/pdf"
RPT_SUPP_PNG <- "04_Figures/F04_Reversal/b_reports/supp/png"
dir.create(RPT_SUPP_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PNG, recursive = TRUE, showWarnings = FALSE)

diag_layout <- "AABBBB\nAABBBB\nCCDDDD\nCCDDDD\nEEEEFF\nEEEEFF\nEEEEFF"
diagnostics <- pS_ora_dedup + pS_r_boot + pS_circ + pS_threshold + pS_goslim + pS_fry_lead +
  plot_layout(design = diag_layout) +
  plot_annotation(tag_levels = "A", theme = theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.margin = margin(4, 4, 4, 4, "mm")
  ))
ggsave(file.path(RPT_SUPP_PNG, "SUPP_F04_diagnostics.png"), diagnostics,
  width = 360, height = 400, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(RPT_SUPP_PDF, "SUPP_F04_diagnostics.pdf"), diagnostics,
  width = 360, height = 400, units = "mm", device = pdf_device
)

methods <- pS_melov + pS_cmap + pS_asym +
  plot_layout(design = "AAAA\nBBBB\nCCCC") +
  plot_annotation(tag_levels = "A", theme = theme(
    plot.tag = element_text(face = "bold", size = 14),
    plot.margin = margin(4, 4, 4, 4, "mm")
  ))
ggsave(file.path(RPT_SUPP_PNG, "SUPP_F04_methods.png"), methods,
  width = 360, height = 340, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(RPT_SUPP_PDF, "SUPP_F04_methods.pdf"), methods,
  width = 360, height = 340, units = "mm", device = pdf_device
)

message("F04 supplementary composites saved")
message("F04 driver complete")
