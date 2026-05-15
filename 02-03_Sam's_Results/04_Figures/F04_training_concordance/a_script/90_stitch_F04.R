#!/usr/bin/env Rscript
# F04 Training Concordance — Master Orchestrator
#
# Run order:
#   1. 02_supp_panels.R  -> supplementary composite + CSVs
#   2. 01_main_panels.R  -> 5-panel main composite + biological summary
#
# Usage (from A_Proteomics_Analysis/):
#   Rscript A_CvH_2026/02-03_Sam\'s_Results/04_Figures/F04_training_concordance/a_script/90_stitch_F04.R

setwd(rprojroot::find_rstudio_root_file())

source("02-03_Sam's_Results/04_Figures/shared/style.R")

BASE    <- "02-03_Sam's_Results/04_Figures/F04_training_concordance"
RPT_PDF <- file.path(BASE, "b_reports", "main",  "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main",  "png")
SUP_PDF <- file.path(BASE, "b_reports", "supp",  "pdf")
SUP_PNG <- file.path(BASE, "b_reports", "supp",  "png")

message("=== F04 Training Concordance: supplementary panels ===")
source("02-03_Sam's_Results/04_Figures/F04_training_concordance/a_script/02_supp_panels.R")

message("=== F04 Training Concordance: main composite ===")
source("02-03_Sam's_Results/04_Figures/F04_training_concordance/a_script/01_main_panels.R")

# ── Verify outputs ────────────────────────────────────────────────────────────
main_pdf  <- file.path(RPT_PDF, "MAIN_F04_composite.pdf")
main_png  <- file.path(RPT_PNG, "MAIN_F04_composite.png")
supp_pdf  <- file.path(SUP_PDF, "SUPP_F04_training_concordance_diagnostics.pdf")
supp_png  <- file.path(SUP_PNG, "SUPP_F04_training_concordance_diagnostics.png")

outputs <- c(main_pdf, main_png, supp_pdf, supp_png)
ok      <- file.exists(outputs)
for (i in seq_along(outputs)) {
  sz <- if (ok[i]) sprintf("%.1f KB", file.size(outputs[i]) / 1024) else "MISSING"
  message(sprintf("  %s [%s]", basename(outputs[i]), sz))
}

# Panel-level PNGs
main_panels <- list.files(file.path(BASE, "b_reports", "main", "png", "panels"),
                           pattern = "\\.png$", full.names = FALSE)
supp_panels <- list.files(file.path(BASE, "b_reports", "supp", "png", "panels"),
                           pattern = "\\.png$", full.names = FALSE)
message(sprintf("  Main panel PNGs: %d", length(main_panels)))
message(sprintf("  Supp panel PNGs: %d", length(supp_panels)))

# c_data CSVs
data_csvs <- list.files(file.path(BASE, "c_data"), pattern = "\\.csv$",
                         recursive = TRUE, full.names = FALSE)
message(sprintf("  c_data CSVs: %d", length(data_csvs)))

if (all(ok)) {
  message("\n=== F04 COMPLETE ===")
} else {
  warning("Some outputs are MISSING — check errors above")
}
