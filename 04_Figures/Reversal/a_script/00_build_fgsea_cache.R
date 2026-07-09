#!/usr/bin/env Rscript
# 00_build_fgsea_cache.R  --  regenerate the pathway-NES cache from the NEW DEP.
#
# Panel B (NES scatter) consumes 04_Figures/shared/fgsea_CRvH.csv, a long table
# of fgsea results for the disease axis (Cancer_vs_Healthy = CRvH_Baseline) and
# the training axis (Training_CR = CR_Training). The cache was built on the old
# DEP; this rebuilds it from 03_DEP/a_non_imputed using the same schema
# (pathway, pval, padj, log2err, ES, NES, size, leadingEdge, database, contrast)
# so the panel needs no changes. Ranks = signed limma t-statistic per gene.

setwd(here::here())
pacman::p_load(dplyr, tidyr, readr, tibble, fgsea)
source("04_Figures/shared/pathway_utils.R")
source("04_Figures/Reversal/a_script/reversal.R")

set.seed(42)
wide <- load_reversal_table("03_DEP/a_non_imputed/c_data/combined_results_pi.csv")

# axis -> gene-level signed t-stat (one value/gene = max |t|; drop NA/blank genes)
make_ranks <- function(gene, tstat) {
  d <- tibble(gene = gene, t = tstat) |>
    filter(!is.na(gene), gene != "", !is.na(t)) |>
    group_by(gene) |> slice_max(abs(t), n = 1, with_ties = FALSE) |> ungroup()
  setNames(d$t, d$gene)
}
ranks <- list(
  Cancer_vs_Healthy = make_ranks(wide$gene, wide$t_D),   # disease axis (D)
  Training_CR       = make_ranks(wide$gene, wide$t_T)    # training axis (T)
)

pw <- build_pathway_collection(min_size = 10, max_size = 500,
                               include_goslim = TRUE, exclude_variants = TRUE)

long <- bind_rows(lapply(names(ranks), function(ctr) {
  res <- fgsea::fgseaMultilevel(pathways = pw, stats = ranks[[ctr]],
                                minSize = 10, maxSize = 500, eps = 0)
  res <- as.data.frame(res)
  res$leadingEdge <- vapply(res$leadingEdge, paste, character(1), collapse = ";")
  res$database <- classify_database(res$pathway)
  res$contrast <- ctr
  res[, c("pathway", "pval", "padj", "log2err", "ES", "NES",
          "size", "leadingEdge", "database", "contrast")]
}))

write_csv(long, "04_Figures/shared/fgsea_CRvH.csv")
cat(sprintf("fgsea cache rebuilt: %d rows | %s | dbs: %s\n",
            nrow(long), paste(names(ranks), collapse = " + "),
            paste(sort(unique(long$database)), collapse = ", ")))
print(long |> count(contrast, database) |> tidyr::pivot_wider(names_from = contrast, values_from = n))
