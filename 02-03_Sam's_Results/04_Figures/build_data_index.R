# Sam's-Results data index.
# Catalogues every artifact a figure script in this subdir might need.
# Source from a panel script:  source("02-03_Sam's_Results/04_Figures/build_data_index.R")
#
# Provides the named list `sam_idx` with three groups:
#   sam_idx$sam        - Sam's own pipeline outputs
#   sam_idx$ours_on_his - our DEP rerun on Sam's normalized data
#   sam_idx$comparison  - 3-way comparison artifacts

setwd(rprojroot::find_rstudio_root_file())

sam_idx <- list(
  sam = list(
    dalist_rds   = "02-03_Sam's_Results/00_input/01_normalized_DAList_SURV_stringent_muscle.RDS",
    normalized   = "02-03_Sam's_Results/00_input/01_normalized_data_SURV_stringent_muscle.csv",
    limma_xlsx   = "02-03_Sam's_Results/00_input/limma_results_all_muscle_str.xlsx"
  ),
  ours_on_his = list(
    dalist_CRvH  = "02-03_Sam's_Results/03_DEP/c_data/01_limma_DAList_CRvH.rds",
    dalist_CR    = "02-03_Sam's_Results/03_DEP/c_data/01_limma_DAList_CR.rds",
    summary_CRvH = "02-03_Sam's_Results/03_DEP/c_data/02_DA_summary_CRvH.csv",
    summary_CR   = "02-03_Sam's_Results/03_DEP/c_data/02_DA_summary_CR.csv",
    combined_CRvH = "02-03_Sam's_Results/03_DEP/c_data/03_combined_results_CRvH.csv",
    combined_CR   = "02-03_Sam's_Results/03_DEP/c_data/03_combined_results_CR.csv",
    per_contrast  = "02-03_Sam's_Results/03_DEP/c_data/04_per_contrast_results",
    xlsx_CRvH    = "02-03_Sam's_Results/03_DEP/c_data/05_results_CRvH.xlsx",
    xlsx_CR      = "02-03_Sam's_Results/03_DEP/c_data/05_results_CR.xlsx"
  ),
  comparison = list(
    xlsx      = "02-03_Sam's_Results/03_DEP/c_data/comparison_3way.xlsx",
    summary   = "02-03_Sam's_Results/03_DEP/c_data/comparison_3way_summary.csv",
    narrative = "02-03_Sam's_Results/SAM_VS_CVH_DIFF.md"
  )
)

# Sanity check: every path resolves
.missing <- unlist(sam_idx)[!file.exists(unlist(sam_idx))]
if (length(.missing) > 0) {
  warning("sam_idx: missing files:\n  ", paste(.missing, collapse = "\n  "))
}
invisible(sam_idx)
