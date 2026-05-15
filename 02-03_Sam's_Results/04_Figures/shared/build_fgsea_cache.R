#!/usr/bin/env Rscript
# Build fGSEA enrichment cache for Sam's 5 CvH DEP contrasts.
# One RDS per contrast written to fgsea_cache/.
# Cache is treated as a frozen artifact: existing files are skipped.
# Delete individual RDS files (or the whole fgsea_cache/ dir) to force regeneration.
#
# Usage:
#   Rscript 04_Figures/shared/build_fgsea_cache.R
# Run from: A_CvH_2026/02-03_Sam's_Results/

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

SCRIPT_DIR   <- "04_Figures/shared"
CACHE_DIR    <- file.path(SCRIPT_DIR, "fgsea_cache")
DEP_DIR      <- "03_DEP/c_data/04_per_contrast_results"

CONTRASTS <- c(
  "Cancer_vs_Healthy",
  "Baseline_Supplement",
  "Training_CR",
  "Training_CRE",
  "Training_PLA",
  "Supplement_Interaction"
)

dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)

source(file.path(SCRIPT_DIR, "pathway_utils.R"))

# ------------------------------------------------------------------
# Check which contrasts still need building
# ------------------------------------------------------------------
todo <- character(0)
for (ctr in CONTRASTS) {
  rds_path <- file.path(CACHE_DIR, paste0(ctr, "_fgsea.rds"))
  if (file.exists(rds_path)) {
    message(sprintf("Cache present for %s — skipping (delete to regenerate)", ctr))
  } else {
    todo <- c(todo, ctr)
  }
}

if (length(todo) == 0) {
  message("All contrast caches present — nothing to do.")
  quit(save = "no", status = 0)
}

message(sprintf("Building fGSEA cache for %d contrast(s): %s",
                length(todo), paste(todo, collapse = ", ")))

# ------------------------------------------------------------------
# Build pathway collection (shared across all contrasts)
# ------------------------------------------------------------------
message("\nLoading pathway collection...")
pw_list <- build_pathway_collection(
  min_size         = 10,
  max_size         = 500,
  include_goslim   = TRUE,
  exclude_variants = FALSE
)
message(sprintf("Pathway collection ready: %d sets", length(pw_list)))

# ------------------------------------------------------------------
# Run fGSEA per contrast
# ------------------------------------------------------------------
set.seed(42)
t_start_all <- proc.time()

for (ctr in todo) {
  csv_path <- file.path(DEP_DIR, paste0(ctr, ".csv"))
  rds_path <- file.path(CACHE_DIR, paste0(ctr, "_fgsea.rds"))

  if (!file.exists(csv_path)) {
    warning(sprintf("Input CSV not found: %s — skipping", csv_path))
    next
  }

  message(sprintf("\n=== %s ===", ctr))
  t_start <- proc.time()

  dep <- read_csv(csv_path, show_col_types = FALSE)

  # Drop rows with missing gene or t-stat
  n_before <- nrow(dep)
  dep <- dep[!is.na(dep$gene) & dep$gene != "" & !is.na(dep$t), ]
  if (nrow(dep) < n_before) {
    message(sprintf("  Dropped %d rows with NA gene or t-stat (kept %d)",
                    n_before - nrow(dep), nrow(dep)))
  }

  # Build named rank vector from t-statistic
  ranks <- setNames(dep$t, dep$gene)
  # Resolve duplicate gene names: keep the one with highest |t|
  if (anyDuplicated(names(ranks))) {
    dup_genes <- names(ranks)[duplicated(names(ranks))]
    message(sprintf("  Deduplicating %d gene name duplicates (max |t|)",
                    length(unique(dup_genes))))
    dep_dd <- dep[order(abs(dep$t), decreasing = TRUE), ]
    dep_dd <- dep_dd[!duplicated(dep_dd$gene), ]
    ranks  <- setNames(dep_dd$t, dep_dd$gene)
  }
  ranks <- sort(ranks, decreasing = TRUE)

  message(sprintf("  Ranks: %d genes, t in [%.3f, %.3f]",
                  length(ranks), min(ranks), max(ranks)))

  res <- run_fgsea_deduplicated(
    ranks          = ranks,
    pathways       = pw_list,
    jaccard_cutoff = 0.5,
    nperm          = 10000,
    min_size       = 15,
    max_size       = 500
  )
  res$contrast <- ctr

  # Serialize leadingEdge list column as semicolon-delimited string
  res <- res |>
    mutate(leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";"))

  # Standardise column order
  keep_cols <- c("pathway", "pval", "padj", "log2err", "ES", "NES",
                 "size", "leadingEdge", "database", "contrast")
  res <- res[, intersect(keep_cols, names(res))]

  saveRDS(res, rds_path)

  elapsed <- (proc.time() - t_start)[["elapsed"]]
  n_sig   <- sum(!is.na(res$padj) & res$padj < 0.05, na.rm = TRUE)
  message(sprintf("  Done: %d pathways tested, %d sig (padj<0.05) | %.1f s",
                  nrow(res), n_sig, elapsed))
}

elapsed_all <- (proc.time() - t_start_all)[["elapsed"]]
message(sprintf("\nTotal fGSEA cache build time: %.1f s (%.1f min)",
                elapsed_all, elapsed_all / 60))

# ------------------------------------------------------------------
# Verify outputs
# ------------------------------------------------------------------
message("\n--- Output verification ---")
for (ctr in CONTRASTS) {
  rds_path <- file.path(CACHE_DIR, paste0(ctr, "_fgsea.rds"))
  if (file.exists(rds_path)) {
    sz  <- file.info(rds_path)$size
    dat <- readRDS(rds_path)
    n_sig <- sum(!is.na(dat$padj) & dat$padj < 0.05, na.rm = TRUE)
    message(sprintf("  %s: %d rows, %d sig | %.1f KB",
                    ctr, nrow(dat), n_sig, sz / 1024))
  } else {
    message(sprintf("  %s: MISSING", ctr))
  }
}
