#!/usr/bin/env Rscript
# Build the fGSEA cache from the primary (non-imputed) Stage 03 t-statistics,
# one ranking per contrast over the shared database panel (Hallmark + KEGG +
# Reactome + GO:BP + GO Slim), EnrichmentMap Jaccard-deduplicated. Source from
# any figure that needs it; skipped if present. Delete the CSV to regenerate.

setwd(here::here())

CACHE_PATH <- "04_Figures/shared/fgsea_CRvH.csv"
DEP_CSV <- "03_DEP/a_non_imputed/c_data/combined_results_pi.csv"

stopifnot("Stage 03 combined_results_pi.csv missing" = file.exists(DEP_CSV))

if (file.exists(CACHE_PATH)) {
  message("fGSEA cache present — delete to regenerate: ", CACHE_PATH)
} else {
  source("04_Figures/shared/pathway_utils.R")
  pacman::p_load(dplyr, readr, purrr)

  dep <- read_csv(DEP_CSV, show_col_types = FALSE)
  contrasts <- c(
    "CRvH_Baseline", "CR_Training", "Resid",
    "Baseline_Supplement", "Training_CRE", "Training_PLA",
    "Supplement_Interaction"
  )
  pw <- build_pathway_collection(min_size = 10, max_size = 500, include_goslim = TRUE)

  set.seed(42)
  cache <- map_dfr(contrasts, function(ctr) {
    ranks <- dep |>
      filter(contrast == ctr, !is.na(t), !is.na(gene)) |>
      group_by(gene) |>
      slice_max(abs(t), n = 1, with_ties = FALSE) |>
      ungroup()
    ranks <- sort(setNames(ranks$t, ranks$gene), decreasing = TRUE)
    res <- run_fgsea_deduplicated(ranks, pw,
      jaccard_cutoff = 0.5,
      nperm = 10000, min_size = 15, max_size = 500
    )
    res$contrast <- ctr
    res
  }) |>
    mutate(leadingEdge = vapply(leadingEdge, paste, character(1), collapse = ";")) |>
    select(pathway, pval, padj, log2err, ES, NES, size, leadingEdge, database, contrast)

  write_csv(cache, CACHE_PATH)
  message(sprintf(
    "Wrote %d rows across %d contrasts to %s",
    nrow(cache), length(contrasts), CACHE_PATH
  ))
}
