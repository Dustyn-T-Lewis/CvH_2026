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

source(file.path(A, "panel_B_nes_scatter.R"))
n_sig_pw_B <- nrow(fgsea_sig)
rho_B <- as.numeric(nes_cor_all$estimate)
rho_lo_B <- nes_ci_all[1]
rho_hi_B <- nes_ci_all[2]
pw_rev_B <- rev_frac

source(file.path(A, "panel_D_fry.R"))
circ_r_D <- circ_r
fry_up_D <- fry_up
fry_dn_D <- fry_dn

source(file.path(A, "panel_E_rrho2.R"))
n_shared_E <- n_shared
max_rev_E <- max(max_UD, max_DU)
message("  self-rendered supplement: SUPP_F04_rrho2")

source(file.path(A, "panel_C_pattern_heatmap.R"))
n_total_C <- n_total
n_pw_C <- n_pw
message("  self-rendered supplement: SUPP_F04_pattern_heatmap")

# Capture the quadrant-ORA composite (panel A) for the main figure.
pA_comp <- composite

# Trajectory and fry are their own figures now; source panel_F so it self-renders.
source(file.path(A, "panel_F_trajectory.R"))

RPT_PDF <- "04_Figures/F04_Reversal/b_reports/main/pdf"
RPT_PNG <- "04_Figures/F04_Reversal/b_reports/main/png"
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

COMP_W <- 340
COMP_H <- 165
TTL_SZ <- FIG_TITLE_SIZE
SUB_SZ <- FIG_SUBTITLE_SIZE

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

# The reversal landscape: quadrant ORA (left, wider) beside the pathway-NES
# scatter (right). A blank top row reserves space for the drawn tag+titles.
layout <- paste(c(
  "############",
  rep("AAAAAAABBBBB", 11)
), collapse = "\n")

pA_comp <- pA_comp +
  plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0, "mm")))
pB <- pB + theme(plot.margin = margin(0, 2, 0, 2, "mm"))

fig <- wrap_elements(full = pA_comp) +
  wrap_elements(full = pB) +
  plot_layout(
    design = layout, widths = rep(1, 12),
    heights = c(12, rep(10, 11))
  )

SUB_OFFSET <- 0.03
Y_TOP <- 0.98

tag_block <- function(d, tag, ttl, sub, x, y) {
  d +
    draw_label(paste0(tag, "  ", ttl),
      x = x, y = y, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1
    ) +
    draw_label(sub,
      x = x, y = y - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic",
      hjust = 0, vjust = 1, colour = "grey30"
    )
}

composite_final <- ggdraw(fig) |>
  tag_block("A", ttl_A, sub_A, 0.005, Y_TOP) |>
  tag_block("B", ttl_B, sub_B, 0.585, Y_TOP)

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

# ── Residual volcano -> supplement (Resid = what training did NOT fix) ──
source(file.path(A, "panel_G_resid_volcano.R"))
pG_resid <- composite
RESID_PNG <- "04_Figures/F04_Reversal/b_reports/supp/png"
RESID_PDF <- "04_Figures/F04_Reversal/b_reports/supp/pdf"
dir.create(RESID_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RESID_PDF, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(RESID_PNG, "SUPP_F04_residual_volcano.png"), pG_resid,
  width = 220, height = 120, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(RESID_PDF, "SUPP_F04_residual_volcano.pdf"), pG_resid,
  width = 220, height = 120, units = "mm", device = pdf_device
)
message("F04 residual volcano saved to supplement")

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

# Each of these is already a patchwork; a bare + flattens the leftmost one into the
# top level, leaving 5 plots for 3 areas and silently dropping two. wrap_elements keeps
# each row as a single unit.
methods <- wrap_elements(pS_melov) + wrap_elements(pS_cmap) + wrap_elements(pS_asym) +
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
