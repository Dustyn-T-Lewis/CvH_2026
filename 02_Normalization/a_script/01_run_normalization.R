#!/usr/bin/env Rscript
# =============================================================================
# 01_run_normalization.R  --  CvH Stage 02: Normalization
#
# Reads the filtered (un-normalized) DAList from Stage 01, writes pre-norm QC +
# norm-method comparison, applies cycloess (Bolstad 2003) via proteoDA, writes
# post-norm QC, and exports the CANONICAL normalized matrix used by Stage 04 DEP.
# =============================================================================

suppressPackageStartupMessages({
  library(proteoDA); library(here); library(readr); library(dplyr); library(tibble)
})
set.seed(42)

in_rds    <- here("01_Filtering", "c_data", "01_DAList_filtered.rds")
report_dir <- here("02_Normalization", "b_reports")
data_dir   <- here("02_Normalization", "c_data")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir,   recursive = TRUE, showWarnings = FALSE)

dal <- readRDS(in_rds)
cat(sprintf("Loaded filtered DAList: %d proteins x %d samples\n",
            nrow(dal$data), ncol(dal$data)))

# --- Pre-norm QC + normalization-method comparison ---------------------------
write_norm_report(dal, grouping_column = "group_time",
                  output_dir = report_dir, filename = "01_norm_comparison.pdf", overwrite = TRUE)
write_qc_report(dal, color_column = "group_time",
                output_dir = report_dir, filename = "02_qc_pre.pdf", overwrite = TRUE)
saveRDS(dal, file.path(data_dir, "01_DAList_prenorm.rds"))

# --- Cycloess normalization --------------------------------------------------
dal <- normalize_data(dal, norm_method = "cycloess")
cat(sprintf("Normalized (cycloess): %d proteins x %d samples\n",
            nrow(dal$data), ncol(dal$data)))

write_qc_report(dal, color_column = "group_time",
                output_dir = report_dir, filename = "03_qc_post.pdf", overwrite = TRUE)

# --- Export canonical normalized artifacts -----------------------------------
write_csv(bind_cols(as_tibble(dal$annotation) |> select(uniprot_id, protein, gene, description),
                    as_tibble(dal$data)),
          file.path(data_dir, "02_normalized.csv"))
saveRDS(dal, file.path(data_dir, "03_DAList_normalized.rds"))

if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
writeLines(capture.output(sessionInfo()), file.path(data_dir, "sessionInfo.txt"))
cat(sprintf("Done -> %s/\n", data_dir))
