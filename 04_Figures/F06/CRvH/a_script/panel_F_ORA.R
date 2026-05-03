# F06 CRvH Panel F -- ORA on fry Driving Proteins (Concordance)
# Runs over-representation analysis on proteins that drive the fry concordance
# signal (set members with t in the expected concordant direction).
# ---------------------------------------------------------------------------
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F06/a_script/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
})

RPT <- "04_Figures/F06/CRvH/b_reports"
DAT <- "04_Figures/F06/CRvH/c_data"
dir.create(file.path(DAT, "panel_F_fry"), recursive = TRUE, showWarnings = FALSE)

# -- Load driving proteins from fry barcode panel ------------------------------

driving_path <- file.path(DAT, "panel_F_fry", "driving_proteins.csv")
if (!file.exists(driving_path)) {
  message("panel_F_ORA: driving_proteins.csv not found -- run panel_F_fry.R first")
  quit(save = "no", status = 0)
}

driving_df <- read_csv(driving_path, show_col_types = FALSE)

if (nrow(driving_df) < 5) {
  message("panel_F_ORA: fewer than 5 driving proteins -- skipping ORA")
  write_csv(tibble(note = "Fewer than 5 driving proteins"),
            file.path(DAT, "panel_F_fry", "driving_ora.csv"))
  quit(save = "no", status = 0)
}

# -- ORA -----------------------------------------------------------------------

dep_df   <- read_csv("03_DEP/c_data/03_combined_results_CRvH.csv", show_col_types = FALSE)
universe <- dep_df$gene

pw_collection <- build_pathway_collection(min_size = 10, max_size = 500,
                                           include_goslim = FALSE,
                                           exclude_variants = TRUE)

driving_genes <- driving_df$gene
ora_res <- tryCatch(
  run_ora_deduplicated(
    genes = driving_genes,
    universe = universe,
    pathways = pw_collection,
    jaccard_cutoff = 0.5,
    min_size = 10, max_size = 500, padj_cutoff = 0.05
  ),
  error = function(e) { message("ORA error: ", e$message); tibble() }
)

if (nrow(ora_res) == 0) {
  message("panel_F_ORA: no significant pathways at padj < 0.05")
  write_csv(tibble(note = "No pathways enriched at padj < 0.05"),
            file.path(DAT, "panel_F_fry", "driving_ora.csv"))

  message("F06 CRvH Panel F (ORA) done -- no significant pathways")
  quit(save = "no", status = 0)
}

# -- Format and export ---------------------------------------------------------

ora_res <- ora_res %>%
  mutate(
    pathway_label = clean_pathway_name(pathway),
    database = classify_database(pathway),
    neg_log10_padj = -log10(padj),
    geneID = sapply(overlapGenes, paste, collapse = "/")
  ) %>%
  arrange(desc(neg_log10_padj))

write_csv(ora_res %>% select(pathway, pathway_label, database, overlap, size,
                              padj, neg_log10_padj, geneID),
          file.path(DAT, "panel_F_fry", "driving_ora.csv"))

message("F06 CRvH Panel F (ORA) done")
