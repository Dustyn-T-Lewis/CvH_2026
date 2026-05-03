# CvH Imputation — Supplementary Excel workbook
# Input: CSV files from c_data/benchmark/ | Output: c_data/benchmark/10_imputation_supp.xlsx

library(openxlsx)
library(readr)
library(tibble)

setwd(rprojroot::find_rstudio_root_file())
add_sheet <- function(wb, name, df, title = NULL, notes = NULL) {
  openxlsx::addWorksheet(wb, name)
  start_row <- 1L
  if (!is.null(title)) {
    openxlsx::writeData(wb, name, title, startRow = 1)
    openxlsx::addStyle(wb, name,
      openxlsx::createStyle(textDecoration = "bold", fontSize = 12),
      rows = 1, cols = 1)
    start_row <- start_row + 1L
  }
  if (!is.null(notes)) {
    for (i in seq_along(notes)) {
      openxlsx::writeData(wb, name, notes[i], startRow = start_row)
      openxlsx::addStyle(wb, name,
        openxlsx::createStyle(fontSize = 10, fontColour = "#555555", wrapText = TRUE),
        rows = start_row, cols = 1)
      start_row <- start_row + 1L
    }
    start_row <- start_row + 1L
  }
  hs <- openxlsx::createStyle(textDecoration = "bold", border = "Bottom",
                               fgFill = "#DCE6F1")
  openxlsx::writeData(wb, name, df, startRow = start_row, headerStyle = hs)
  openxlsx::freezePane(wb, name, firstActiveRow = start_row + 1L,
                        firstActiveCol = 2)
  openxlsx::setColWidths(wb, name, cols = seq_len(ncol(df)), widths = "auto")
}

DATA_DIR <- "02_Imputation/c_data/benchmark"

bench_sum  <- read_csv(file.path(DATA_DIR, "03_benchmark_summary.csv"), show_col_types = FALSE)
bench_df   <- read_csv(file.path(DATA_DIR, "04_benchmark_raw_iterations.csv"), show_col_types = FALSE)
ext_sum    <- read_csv(file.path(DATA_DIR, "05_benchmark_extended.csv"), show_col_types = FALSE)
bin_sum    <- read_csv(file.path(DATA_DIR, "06_benchmark_per_intensity.csv"), show_col_types = FALSE)
miss_class <- read_csv(file.path(DATA_DIR, "..", "02_mar_mnar_classification.csv"), show_col_types = FALSE)
mnar_audit <- read_csv(file.path(DATA_DIR, "..", "08_mnar_imputation_audit.csv"), show_col_types = FALSE)
loocv_df   <- read_csv(file.path(DATA_DIR, "11_per_sample_loocv.csv"), show_col_types = FALSE)
info_lines <- readLines(file.path(DATA_DIR, "..", "09_imputation_summary.txt"))
info_df    <- tibble(Parameter = sub(" = .*", "", info_lines),
                     Value     = sub(".* = ", "", info_lines))

wb <- createWorkbook()
readme <- tibble(
  Sheet = c("Benchmark", "Extended", "Per_Intensity", "Classification",
            "MNAR_Audit", "LOOCV", "Iterations", "Summary"),
  Description = c(
    "Method ranking by NRMSE and PSS (12 methods x 20 iterations)",
    "Top 5 methods: NRMSE + PSS + per-intensity-tertile breakdown",
    "Per-intensity bin NRMSE for all methods (low/mid/high abundance)",
    "Per-protein Complete/MAR/MNAR classification with reliability flag",
    "MNAR pre/post means, shift, Cohen's d",
    "Per-sample leave-one-out cross-validation NRMSE",
    "Raw per-iteration NRMSE and PSS for all methods",
    "Pipeline summary statistics"))
add_sheet(wb, "README", readme)
add_sheet(wb, "Benchmark", bench_sum)
add_sheet(wb, "Extended", ext_sum)
add_sheet(wb, "Per_Intensity", bin_sum)
add_sheet(wb, "Classification", miss_class)
add_sheet(wb, "MNAR_Audit", mnar_audit)
add_sheet(wb, "LOOCV", loocv_df)
add_sheet(wb, "Iterations", bench_df)
add_sheet(wb, "Summary", info_df)

saveWorkbook(wb, file.path(DATA_DIR, "10_imputation_supp.xlsx"), overwrite = TRUE)
cat(sprintf("Done: %s/10_imputation_supp.xlsx\n", DATA_DIR))
