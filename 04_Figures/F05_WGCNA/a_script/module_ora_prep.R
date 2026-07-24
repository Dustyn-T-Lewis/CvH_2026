# Per-module over-representation analysis feeding the F05 module cards.
# Background is every quantified protein assigned to a module; each module is
# tested against the full multi-DB collection and deduplicated by Jaccard.
# Cached to c_data/module_ora.csv so the panels stay cheap to re-render.

setwd(here::here())
source("04_Figures/F05_WGCNA/a_script/style.R")
source("04_Figures/shared/pathway_utils.R")

pacman::p_load(readr, dplyr, purrr, tibble)

DAT <- "04_Figures/F05_WGCNA/c_data"

module_df <- read_csv(file.path(DAT, "wgcna_module_assignments.csv"), show_col_types = FALSE)
mod_bio <- read_csv(file.path(DAT, "mod_bio_labels.csv"), show_col_types = FALSE)

universe <- unique(na.omit(module_df$gene))
pw <- build_pathway_collection(min_size = 15, max_size = 500, include_goslim = FALSE)

ora <- map_dfr(mod_bio$module_color, function(mod) {
  genes <- unique(na.omit(module_df$gene[module_df$module_color == mod]))
  message(sprintf("ORA %s (%d genes)", mod, length(genes)))
  res <- run_ora_deduplicated(
    genes = genes, universe = universe, pathways = pw,
    jaccard_cutoff = 0.5, min_size = 10, max_size = 500, padj_cutoff = 1
  )
  if (is.null(res) || nrow(res) == 0) {
    return(tibble())
  }
  res |>
    mutate(module_color = mod) |>
    select(module_color, pathway, database, pval, padj, overlap, size, odds_ratio)
})

write_csv(ora, file.path(DAT, "module_ora.csv"))
message(sprintf(
  "Wrote %d ORA rows across %d modules",
  nrow(ora), dplyr::n_distinct(ora$module_color)
))
