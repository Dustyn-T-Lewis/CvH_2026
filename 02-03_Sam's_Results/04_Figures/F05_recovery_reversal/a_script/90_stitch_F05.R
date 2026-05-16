#!/usr/bin/env Rscript
# F05 Recovery Reversal — Master Orchestrator
# Sources supp panels first (so CSVs + PNGs exist for xlsx and composite),
# then main stitcher (composite + xlsx + cleanup), then supp stitcher (composite).
#
# Run order:
#   1. _supp_enrichment_heatmap.R  -> supp panel CSVs + PNGs (pre-generate)
#   2. 02_supp_panels.R            -> supp composite from pre-rendered PNGs
#   3. 01_main_panels.R            -> 5-panel main composite
#   4. Supplementary xlsx          -> build workbook
#
# Usage (from A_CvH_2026/):
#   Rscript "02-03_Sam's_Results/04_Figures/F05_recovery_reversal/a_script/90_stitch_F05.R"

setwd(rprojroot::find_rstudio_root_file())

source("02-03_Sam's_Results/04_Figures/shared/style.R")

BASE <- "02-03_Sam's_Results/04_Figures/F05_recovery_reversal"

message("=== F05: Pre-generating supp panel data + PNGs ===")
source("02-03_Sam's_Results/04_Figures/F05_recovery_reversal/a_script/_supp_enrichment_heatmap.R")

message("=== F05: Running supp diagnostics composite ===")
source("02-03_Sam's_Results/04_Figures/F05_recovery_reversal/a_script/02_supp_panels.R")

message("=== F05: Running main composite ===")
source("02-03_Sam's_Results/04_Figures/F05_recovery_reversal/a_script/01_main_panels.R")

# Supplementary Excel: one workbook, sheets keyed to figure panels
source("02-03_Sam's_Results/04_Figures/shared/figure_supplement_helpers.R")

message("=== F05 supplementary workbook ===")
enrichment_reversal_df    <- read.csv(file.path(BASE, "c_data", "panel_supp", "enrichment_reversal.csv"),
                                      stringsAsFactors = FALSE, check.names = FALSE)
reversal_pathway_stats_df <- read.csv(file.path(BASE, "c_data", "panel_supp", "reversal_pathway_stats.csv"),
                                      stringsAsFactors = FALSE, check.names = FALSE)
f05_specs <- list(
  list(name = "panel_A_ora_quadrant",       path = file.path(BASE, "c_data", "panel_A", "ora_quadrant.csv")),
  list(name = "panel_B_pattern_class",      path = file.path(BASE, "c_data", "panel_B_heatmap", "pattern_classification.csv")),
  list(name = "panel_B_sankey",             path = file.path(BASE, "c_data", "panel_B_heatmap", "sankey_links.csv")),
  list(name = "panel_B_bar",                path = file.path(BASE, "c_data", "panel_B_heatmap", "bar_data.csv")),
  list(name = "panel_C_fry_results",        path = file.path(BASE, "c_data", "panel_C_fry", "fry_results_all.csv")),
  list(name = "panel_C_fry_driving",        path = file.path(BASE, "c_data", "panel_C_fry", "driving_proteins.csv")),
  list(name = "panel_D_nes_scatter",        path = file.path(BASE, "c_data", "panel_D", "nes_scatter.csv")),
  list(name = "panel_E_rrho2_summary",      path = file.path(BASE, "c_data", "panel_E", "rrho2_summary.csv")),
  list(name = "panel_E_rrho2_hotspot",      path = file.path(BASE, "c_data", "panel_E", "rrho2_hotspot_genes.csv")),
  list(name = "panel_E_rrho2_ora_concord",  path = file.path(BASE, "c_data", "panel_E", "rrho2_ora_concordant.csv")),
  list(name = "panel_E_rrho2_ora_discord",  path = file.path(BASE, "c_data", "panel_E", "rrho2_ora_discordant.csv")),
  list(name = "panel_E_rrho2_discord_note", path = file.path(BASE, "c_data", "panel_E", "rrho2_ora_discordant_note.csv")),
  list(name = "SUPP_enrichment_reversal",   df = enrichment_reversal_df),
  list(name = "SUPP_reversal_pathway_stats", df = reversal_pathway_stats_df),
  list(name = "SUPP_ora_dedup",             path = file.path(BASE, "c_data", "panel_supp", "SUPP_ora_dedup_sensitivity.csv")),
  list(name = "SUPP_r_bootstrap",           path = file.path(BASE, "c_data", "panel_supp", "SUPP_r_bootstrap.csv")),
  list(name = "SUPP_reversal_threshold",    path = file.path(BASE, "c_data", "panel_supp", "SUPP_reversal_threshold.csv")),
  list(name = "SUPP_goslim_bars",           path = file.path(BASE, "c_data", "panel_supp", "SUPP_goslim_distribution.csv")),
  list(name = "SUPP_fry_leading",           path = file.path(BASE, "c_data", "panel_supp", "SUPP_fry_leading_edge.csv")),
  list(name = "SUPP_fry_circularity",       path = file.path(BASE, "c_data", "panel_supp", "SUPP_fry_circularity.csv"))
)
build_workbook(
  file.path(BASE, "c_data", "F05_supplementary.xlsx"),
  title = "F05 — Figure 5 source data",
  description = "Recovery-reversal diagnostics: quadrant ORA, pathway NES scatter, per-protein pattern classification, fry rotation test, RRHO2.",
  overview_df = data.frame(
    Sheet = c(
      "panel_A_ora_quadrant",
      "panel_B_pattern_class", "panel_B_sankey", "panel_B_bar",
      "panel_C_fry_results", "panel_C_fry_driving",
      "panel_D_nes_scatter",
      "panel_E_rrho2_summary",
      "panel_E_rrho2_hotspot", "panel_E_rrho2_ora_concord",
      "panel_E_rrho2_ora_discord", "panel_E_rrho2_discord_note",
      "SUPP_enrichment_reversal", "SUPP_reversal_pathway_stats",
      "SUPP_ora_dedup", "SUPP_r_bootstrap", "SUPP_reversal_threshold",
      "SUPP_goslim_bars", "SUPP_fry_leading", "SUPP_fry_circularity"),
    Description = c(
      "Panel A: ORA by reversal-quadrant scatter (Reversed Up/Down, Exacerbated Up/Down)",
      "Panel B: per-protein reversal pattern classification (heatmap source)",
      "Panel B: pathway-protein sankey links (heatmap source)",
      "Panel B: per-pattern bar chart counts",
      "Panel C: fry rotation test for reversal contrasts",
      "Panel C: reversal driving proteins (concordant sign-flip)",
      "Panel D: NES scatter (Cancer_vs_Healthy vs Training_CR) per pathway with Spearman + Fisher Z",
      "Panel E: RRHO2 Cancer_vs_Healthy vs Training_CR quadrant summary (max -log10p per quadrant)",
      "Panel E: RRHO2 hotspot genes per quadrant",
      "Panel E: ORA on RRHO2 reversed quadrant genes",
      "Panel E: ORA on RRHO2 exacerbated quadrant genes",
      "Panel E: Notes on exacerbated quadrant ORA (if applicable)",
      "SUPP: pathway-level reversal enrichment",
      "SUPP: per-pathway reversal stats (NES CvH, NES CR, sign agreement)",
      "SUPP: ORA dedup sensitivity across Jaccard cutoffs (reversal quadrants)",
      "SUPP: Pearson r bootstrap (1000 reps, 95% CI)",
      "SUPP: Reversal classification threshold sensitivity",
      "SUPP: GO Slim category distribution by reversal quadrant",
      "SUPP: Top fry driving proteins by |t-stat| in Training CR",
      "SUPP: Circularity diagnostic (protein-permuted null)"),
    stringsAsFactors = FALSE),
  sheet_specs = f05_specs
)
cleanup_after_workbook(f05_specs,
  extra_subdirs = c(file.path(BASE, "c_data", "panel_A"),
                     file.path(BASE, "c_data", "panel_B_heatmap"),
                     file.path(BASE, "c_data", "panel_C_fry"),
                     file.path(BASE, "c_data", "panel_D"),
                     file.path(BASE, "c_data", "panel_E"),
                     file.path(BASE, "c_data", "panel_supp"),
                     file.path(BASE, "c_data", "reversal_tests")))

# ── Verify outputs ────────────────────────────────────────────────────────────
RPT_PDF <- file.path(BASE, "b_reports", "main",  "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main",  "png")
SUP_PDF <- file.path(BASE, "b_reports", "supp",  "pdf")
SUP_PNG <- file.path(BASE, "b_reports", "supp",  "png")

main_pdf <- file.path(RPT_PDF, "MAIN_F05_composite.pdf")
main_png <- file.path(RPT_PNG, "MAIN_F05_composite.png")
supp_pdf <- file.path(SUP_PDF, "SUPP_F05_diagnostics.pdf")
supp_png <- file.path(SUP_PNG, "SUPP_F05_diagnostics.png")
xlsx_out <- file.path(BASE, "c_data", "F05_supplementary.xlsx")

outputs <- c(main_pdf, main_png, supp_pdf, supp_png, xlsx_out)
ok      <- file.exists(outputs)
for (i in seq_along(outputs)) {
  sz <- if (ok[i]) sprintf("%.1f KB", file.size(outputs[i]) / 1024) else "MISSING"
  message(sprintf("  %s [%s]", basename(outputs[i]), sz))
}

main_panels <- list.files(file.path(BASE, "b_reports", "main", "png", "panels"),
                           pattern = "\\.png$", full.names = FALSE)
supp_panels <- list.files(file.path(BASE, "b_reports", "supp", "png", "panels"),
                           pattern = "\\.png$", full.names = FALSE)
message(sprintf("  Main panel PNGs: %d", length(main_panels)))
message(sprintf("  Supp panel PNGs: %d", length(supp_panels)))

if (all(ok)) {
  message("\n=== F05 COMPLETE ===")
} else {
  warning("Some outputs are MISSING — check errors above")
}
