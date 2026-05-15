#!/usr/bin/env Rscript
# F01 — Phenotype Figure: Master Orchestrator (Sam CvH)
# Runs supp panels first (CSVs written), then main panels + composites.

setwd(rprojroot::find_rstudio_root_file())

source("04_Figures/shared/figure_supplement_helpers.R")

DAT      <- "02-03_Sam's_Results/04_Figures/F01_phenotype/c_data"
DAT_SUPP <- file.path(DAT, "supp")

# ── 1. Supp panels (writes per-panel PNGs + audit CSVs) ──────────────────────
source("02-03_Sam's_Results/04_Figures/F01_phenotype/a_script/02_supp_panels.R")

# ── 2. Main panels + both composites ─────────────────────────────────────────
source("02-03_Sam's_Results/04_Figures/F01_phenotype/a_script/01_main_panels.R")

# ── 3. Build supplementary workbook ──────────────────────────────────────────
f01_specs <- list(
  list(name = "panel_A_leg_ext",    path = file.path(DAT,      "panel_A_leg_ext.csv")),
  list(name = "panel_B_dxa_lbm",    path = file.path(DAT,      "panel_B_dxa_lbm.csv")),
  list(name = "panel_C_alm",        path = file.path(DAT,      "panel_C_alm.csv")),
  list(name = "SUPP_panel_A_chest", path = file.path(DAT_SUPP, "panel_A_chest_press.csv")),
  list(name = "SUPP_panel_B_grip",  path = file.path(DAT_SUPP, "panel_B_grip.csv")),
  list(name = "SUPP_panel_C_sts",   path = file.path(DAT_SUPP, "panel_C_sts.csv")))

build_workbook(
  file.path(DAT, "F01_supplementary.xlsx"),
  title       = "F01 — Phenotype source data (Sam CvH)",
  description = "Main panels A–C and supplementary panels A–C. CRE vs PLA within SURV cohort (N≤15).",
  overview_df = data.frame(
    Sheet       = sapply(f01_specs, `[[`, "name"),
    Description = c("Leg Extension 1RM", "DXA Lean Body Mass",
                    "Appendicular Lean Mass",
                    "Chest Press 1RM", "Grip Strength",
                    "Sit-to-Stand Max Power")),
  sheet_specs = f01_specs)

cleanup_after_workbook(f01_specs, extra_subdirs = DAT_SUPP)

# ── 4. Smoke-test: verify expected outputs ────────────────────────────────────
BASE <- "02-03_Sam's_Results/04_Figures/F01_phenotype"
expected <- c(
  file.path(BASE, "b_reports", "main", "png", "panels", "MAIN_panel_A_leg_ext.png"),
  file.path(BASE, "b_reports", "main", "png", "panels", "MAIN_panel_B_dxa_lbm.png"),
  file.path(BASE, "b_reports", "main", "png", "panels", "MAIN_panel_C_alm.png"),
  file.path(BASE, "b_reports", "supp", "png", "panels", "SUPP_panel_A_chest_press.png"),
  file.path(BASE, "b_reports", "supp", "png", "panels", "SUPP_panel_B_grip.png"),
  file.path(BASE, "b_reports", "supp", "png", "panels", "SUPP_panel_C_sts.png"),
  file.path(BASE, "b_reports", "main", "pdf", "MAIN_F01_composite.pdf"),
  file.path(BASE, "b_reports", "main", "png", "MAIN_F01_composite.png"),
  file.path(BASE, "b_reports", "supp", "pdf", "SUPP_F01_composite.pdf"),
  file.path(BASE, "b_reports", "supp", "png", "SUPP_F01_composite.png"),
  file.path(DAT,  "F01_supplementary.xlsx")
)

missing <- expected[!file.exists(expected)]
if (length(missing) == 0) {
  message("F01 smoke-test PASSED — all ", length(expected), " expected outputs present")
} else {
  warning("F01 smoke-test FAILED — missing files:\n",
          paste("  -", missing, collapse = "\n"))
}

message("F01 complete")
