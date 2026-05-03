# F04 Panel A: Volcano Ring for Cancer_vs_Healthy
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/volcano_ring.R")

suppressPackageStartupMessages({
  library(tidyverse)
})

RPT <- "04_Figures/F04/b_reports"
DAT <- "04_Figures/F04/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_A"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# -- Load data --
dep_df <- read_csv("03_DEP/c_data/03_combined_results_CRvH.csv", show_col_types = FALSE)

fgsea_df <- read_csv("04_Figures/F02/c_data/06_panel_H_CRvH_fgsea_results.csv",
                      show_col_types = FALSE)

# -- Select ring terms --
top_terms <- select_ring_terms(fgsea_df, "Cancer_vs_Healthy",
                                databases = c("Hallmark", "GO Slim"))

ring_data <- build_ring_with_gaps(top_terms, "Cancer_vs_Healthy", fgsea_df,
                                   databases = c("Hallmark", "GO Slim"))

# -- Build volcano ring --
p_vr <- make_volcano_ring(
  de_df     = dep_df,
  go_df     = fgsea_df,
  contrast  = "Cancer_vs_Healthy",
  contrast_title    = "Cancer vs Healthy (Baseline)",
  contrast_subtitle = "CR_T1 \u2212 H_T1 | 2,582 proteins \u00b7 limma + dupCor \u00b7 missForest",
  databases = c("Hallmark", "GO Slim"),
  ring_data_override = ring_data
)

ggsave(file.path(RPT, "panel_A_volcano_CvH.pdf"), p_vr,
       width = 260, height = 240, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_A_volcano_CvH.png"), p_vr,
       width = 260, height = 240, units = "mm", dpi = 300)

# -- Export ring data --
if (!is.null(attr(p_vr, "ring_data"))) {
  attr(p_vr, "ring_data") %>%
    mutate(clean_label = as.character(clean_label)) %>%
    select(pathway, NES, padj, size, database, clean_label) %>%
    write_csv(file.path(DAT, "panel_A", "ring_terms.csv"))
}

cat("F04 Panel A done\n")
