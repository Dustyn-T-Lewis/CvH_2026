#!/usr/bin/env Rscript
# F04 — Training Concordance: Master Orchestrator
# Sources supp panels first (so CSVs + PNGs exist for xlsx and composite),
# then main panels (composite + xlsx + cleanup).

setwd(rprojroot::find_rstudio_root_file())

source("02-03_Sam's_Results/04_Figures/shared/figure_supplement_helpers.R")

DAT <- "02-03_Sam's_Results/04_Figures/F04_training_concordance/c_data"

# Supp panels first (CSVs needed for xlsx)
source("02-03_Sam's_Results/04_Figures/F04_training_concordance/a_script/02_supp_panels.R")

# Main panels + composite
source("02-03_Sam's_Results/04_Figures/F04_training_concordance/a_script/01_main_panels.R")

# Build supplementary xlsx
message("=== F04 supplementary workbook ===")
enrichment_blunting_df <- read.csv(file.path(DAT, "panel_supp", "enrichment_blunting.csv"),
                                   stringsAsFactors = FALSE, check.names = FALSE)
f04_specs <- list(
  list(name = "panel_A_ora_quadrant",       path = file.path(DAT, "panel_A", "ora_quadrant.csv")),
  list(name = "panel_B_pattern_class",      path = file.path(DAT, "panel_B_heatmap", "pattern_classification.csv")),
  list(name = "panel_B_sankey",             path = file.path(DAT, "panel_B_heatmap", "sankey_links.csv")),
  list(name = "panel_B_bar",                path = file.path(DAT, "panel_B_heatmap", "bar_data.csv")),
  list(name = "panel_C_fry_results",        path = file.path(DAT, "panel_C_fry", "fry_results_all.csv")),
  list(name = "panel_C_fry_driving",        path = file.path(DAT, "panel_C_fry", "driving_proteins.csv")),
  list(name = "panel_D_nes_scatter",        path = file.path(DAT, "panel_D", "nes_scatter.csv")),
  list(name = "panel_E_rrho2_summary",      path = file.path(DAT, "panel_E", "rrho2_summary.csv")),
  list(name = "panel_E_rrho2_hotspot",      path = file.path(DAT, "panel_E", "rrho2_hotspot_genes.csv")),
  list(name = "panel_E_rrho2_ora_concord",  path = file.path(DAT, "panel_E", "rrho2_ora_concordant.csv")),
  list(name = "panel_E_rrho2_ora_discord",  path = file.path(DAT, "panel_E", "rrho2_ora_discordant.csv")),
  list(name = "SUPP_enrichment_blunting",   df = enrichment_blunting_df),
  list(name = "SUPP_ora_dedup",             path = file.path(DAT, "panel_supp", "SUPP_ora_dedup_sensitivity.csv")),
  list(name = "SUPP_rho_bootstrap",         path = file.path(DAT, "panel_supp", "SUPP_rho_bootstrap.csv")),
  list(name = "SUPP_threshold_sens",        path = file.path(DAT, "panel_supp", "SUPP_threshold_sensitivity.csv")),
  list(name = "SUPP_goslim_bars",           path = file.path(DAT, "panel_supp", "SUPP_goslim_distribution.csv")),
  list(name = "SUPP_fry_leading",           path = file.path(DAT, "panel_supp", "SUPP_fry_leading_edge.csv"))
)
build_workbook(
  file.path(DAT, "F04_supplementary.xlsx"),
  title = "F04 — Figure 4 source data",
  description = "Training-concordance diagnostics (CRE vs PLA): quadrant ORA, pathway NES scatter, per-protein pattern classification, fry rotation test, RRHO2.",
  overview_df = data.frame(
    Sheet = c(
      "panel_A_ora_quadrant",
      "panel_B_pattern_class", "panel_B_sankey", "panel_B_bar",
      "panel_C_fry_results", "panel_C_fry_driving",
      "panel_D_nes_scatter",
      "panel_E_rrho2_summary",
      "panel_E_rrho2_hotspot", "panel_E_rrho2_ora_concord", "panel_E_rrho2_ora_discord",
      "SUPP_enrichment_blunting",
      "SUPP_ora_dedup", "SUPP_rho_bootstrap", "SUPP_threshold_sens",
      "SUPP_goslim_bars", "SUPP_fry_leading"),
    Description = c(
      "Panel A: ORA by training-concordance scatter quadrant",
      "Panel B: per-protein concordance pattern classification (heatmap source)",
      "Panel B: pathway-protein sankey links (heatmap source)",
      "Panel B: per-pattern bar chart counts",
      "Panel C: fry rotation test p-values + directions per pathway",
      "Panel C: concordant driving proteins (gene, t-stats, logFC)",
      "Panel D: NES scatter (Training_CRE vs Training_PLA) per pathway, Spearman + Fisher Z",
      "Panel E: RRHO2 quadrant summary (max -log10p per quadrant)",
      "Panel E: RRHO2 hotspot genes per quadrant",
      "Panel E: ORA on RRHO2 concordant quadrant genes",
      "Panel E: ORA on RRHO2 discordant quadrant genes",
      "SUPP: blunting enrichment (Pi-score-weighted ORA on Training_PLA DEPs)",
      "SUPP: ORA dedup sensitivity across Jaccard cutoffs",
      "SUPP: Spearman rho bootstrap (1000 reps, 95% CI)",
      "SUPP: Protein counts per quadrant across significance thresholds",
      "SUPP: GO Slim category distribution by concordance quadrant",
      "SUPP: Top 20 fry driving proteins by |t-stat| in Training PLA"),
    stringsAsFactors = FALSE),
  sheet_specs = f04_specs
)
cleanup_after_workbook(f04_specs,
  extra_subdirs = c(file.path(DAT, "panel_A"),
                     file.path(DAT, "panel_B_heatmap"),
                     file.path(DAT, "panel_C_fry"),
                     file.path(DAT, "panel_D"),
                     file.path(DAT, "panel_E"),
                     file.path(DAT, "panel_supp"),
                     file.path(DAT, "supp")))

# Final cleanup: remove any leftover CSVs
remaining <- list.files(DAT, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
if (length(remaining)) {
  file.remove(remaining)
  message(sprintf("  final cleanup: removed %d leftover CSV(s)", length(remaining)))
}

message("=== F04 complete ===")
