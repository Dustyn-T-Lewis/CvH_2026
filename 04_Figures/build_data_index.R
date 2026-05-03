# Build master supplementary data index from per-figure data dictionaries
# Output: 04_Figures/supplementary_data_index.xlsx (one sheet per figure/stream)
setwd(rprojroot::find_rstudio_root_file())

library(readr)
library(openxlsx)

wb <- createWorkbook()

# F02 QC — CRvH and CR streams
for (stream in c("CRvH", "CR")) {
  fig_stream <- paste0("F02_", stream)
  dict_path <- file.path("04_Figures", "F02", stream, "c_data", "00_data_dictionary.csv")
  if (!file.exists(dict_path)) {
    message("Skipping ", fig_stream, " -- no data dictionary found")
    next
  }
  df <- read_csv(dict_path, show_col_types = FALSE)
  addWorksheet(wb, fig_stream)
  writeData(wb, fig_stream, df)
  setColWidths(wb, fig_stream, cols = 1:ncol(df), widths = c(40, 10, 60, 30))
}

# F03-F05 — CRvH and CR streams (F05 only CRvH)
for (fig in c("F03", "F04", "F05")) {
  streams <- if (fig == "F05") "CRvH" else c("CRvH", "CR")
  for (stream in streams) {
    fig_stream <- paste0(fig, "_", stream)
    dict_path <- file.path("04_Figures", fig, stream, "c_data", "00_data_dictionary.csv")
    if (!file.exists(dict_path)) {
      message("Skipping ", fig_stream, " -- no data dictionary found")
      next
    }
    df <- read_csv(dict_path, show_col_types = FALSE)
    addWorksheet(wb, fig_stream)
    writeData(wb, fig_stream, df)
    setColWidths(wb, fig_stream, cols = 1:ncol(df), widths = c(40, 10, 60, 30))
  }
}

out_path <- "04_Figures/supplementary_data_index.xlsx"
saveWorkbook(wb, out_path, overwrite = TRUE)
message("Wrote ", out_path)
