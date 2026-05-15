#!/usr/bin/env Rscript
# F02_DEP_overview — Master Orchestrator
# Runs supp panels first, then main panels.
#
# Run from A_CvH_2026/ (or any subdirectory — rprojroot resolves to A_CvH_2026/).

setwd(rprojroot::find_rstudio_root_file())

message("=== F02_DEP_overview: running supp panels ===")
source("02-03_Sam's_Results/04_Figures/F02_DEP_overview/a_script/02_supp_panels.R")

message("=== F02_DEP_overview: running main panels ===")
source("02-03_Sam's_Results/04_Figures/F02_DEP_overview/a_script/01_main_panels.R")

message("=== F02_DEP_overview complete ===")
