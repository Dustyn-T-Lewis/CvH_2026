# Excel supplement helpers used by every figure's 90_stitch_F0x.R.
# Pattern: panel scripts write clean per-panel CSVs; the stitcher reads them
# via safe_read() and writes one panel-labeled workbook to c_data/F0X_supplementary.xlsx.
# After the workbook is saved, cleanup_after_workbook() removes the consumed
# CSVs/subdirs so each c_data/ ends up holding just the Excel.

pacman::p_load(openxlsx, readr, readxl)

add_sheet <- function(wb, name, data) {
  addWorksheet(wb, name)
  writeData(wb, name, data)
  hs <- createStyle(textDecoration = "bold")
  addStyle(wb, name, hs, rows = 1, cols = seq_len(ncol(data)), gridExpand = TRUE)
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(data)), widths = "auto")
  cat(sprintf("    + %s: %d x %d\n", name, nrow(data), ncol(data)))
}

safe_read <- function(path) {
  if (file.exists(path)) {
    as.data.frame(read_csv(path, show_col_types = FALSE))
  } else {
    cat(sprintf("    SKIP (not found): %s\n", path))
    NULL
  }
}

build_workbook <- function(out_file, title, description, overview_df, sheet_specs) {
  wb <- createWorkbook()
  addWorksheet(wb, "Overview")

  writeData(wb, "Overview", title, startRow = 1, startCol = 1, colNames = FALSE)
  mergeCells(wb, "Overview", cols = 1:2, rows = 1)
  addStyle(wb, "Overview",
    createStyle(fontSize = 14, textDecoration = "bold"),
    rows = 1, cols = 1
  )

  writeData(wb, "Overview", description, startRow = 2, startCol = 1, colNames = FALSE)
  mergeCells(wb, "Overview", cols = 1:2, rows = 2)
  addStyle(wb, "Overview",
    createStyle(
      fontSize = 11, textDecoration = "italic",
      fontColour = "#555555",
      wrapText = TRUE, valign = "top"
    ),
    rows = 2, cols = 1
  )
  setRowHeights(wb, "Overview", rows = 2, heights = 40)

  writeData(wb, "Overview", overview_df, startRow = 4, startCol = 1, colNames = TRUE)
  addStyle(wb, "Overview",
    createStyle(
      fontSize = 11, fontColour = "#FFFFFF",
      textDecoration = "bold", fgFill = "#2F4F4F",
      valign = "center"
    ),
    rows = 4, cols = 1:2, gridExpand = TRUE
  )
  n_rows <- nrow(overview_df)
  if (n_rows > 0) {
    addStyle(wb, "Overview",
      createStyle(wrapText = TRUE, valign = "top"),
      rows = 5:(4 + n_rows), cols = 1:2, gridExpand = TRUE
    )
  }
  setRowHeights(wb, "Overview", rows = 4, heights = 20)
  setColWidths(wb, "Overview", cols = 1:2, widths = c(32, 95))
  freezePane(wb, "Overview", firstActiveRow = 5)

  for (spec in sheet_specs) {
    df <- if (!is.null(spec$df)) spec$df else safe_read(spec$path)
    if (!is.null(df)) add_sheet(wb, spec$name, df)
  }
  saveWorkbook(wb, out_file, overwrite = TRUE)
  cat(sprintf("  Saved: %s (%.0f KB)\n\n", out_file, file.size(out_file) / 1e3))
}

cleanup_after_workbook <- function(sheet_specs,
                                   extra_subdirs = character(),
                                   extra_files = character(),
                                   preserve_patterns = c(
                                     "^00_input/",
                                     "^01_Filtering/",
                                     "^02_Normalization/",
                                     "^03_DEP/",
                                     "^04_Figures/shared/"
                                   )) {
  is_preserved <- function(path) {
    path <- paste0(sub("/$", "", path), "/")
    any(vapply(preserve_patterns, function(p) grepl(p, path), logical(1)))
  }
  removed <- 0L
  preserved <- 0L
  for (sp in sheet_specs) {
    if (!is.null(sp$path) && file.exists(sp$path)) {
      if (is_preserved(sp$path)) {
        preserved <- preserved + 1L
      } else {
        unlink(sp$path)
        removed <- removed + 1L
      }
    }
  }
  for (sd in extra_subdirs) {
    if (dir.exists(sd) && !is_preserved(sd)) {
      unlink(sd, recursive = TRUE)
      removed <- removed + 1L
    }
  }
  for (f in extra_files) {
    if (file.exists(f) && !is_preserved(f)) {
      unlink(f)
      removed <- removed + 1L
    }
  }
  cat(sprintf(
    "  cleanup: removed %d intermediate(s); preserved %d upstream/shared path(s)\n",
    removed, preserved
  ))
}
