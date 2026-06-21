# Build master supplementary data index from per-figure data dictionaries
# Output: 04_Figures/supplementary_data_index.xlsx (one sheet per figure/stream)
setwd(rprojroot::find_rstudio_root_file())

library(readr)
library(openxlsx)

wb <- createWorkbook()

# Single 3-group stream (CRvH) per figure — supplements/interaction dropped.
for (fig in c("F02", "F03", "F04")) {
  fig_stream <- paste0(fig, "_CRvH")
  dict_path <- file.path("04_Figures", fig, "CRvH", "c_data", "00_data_dictionary.csv")
  if (!file.exists(dict_path)) {
    message("Skipping ", fig_stream, " -- no data dictionary found")
    next
  }
  df <- read_csv(dict_path, show_col_types = FALSE)
  addWorksheet(wb, fig_stream)
  writeData(wb, fig_stream, df)
  setColWidths(wb, fig_stream, cols = 1:ncol(df), widths = c(40, 10, 60, 30))
}

out_path <- "04_Figures/supplementary_data_index.xlsx"
saveWorkbook(wb, out_path, overwrite = TRUE)
message("Wrote ", out_path)
