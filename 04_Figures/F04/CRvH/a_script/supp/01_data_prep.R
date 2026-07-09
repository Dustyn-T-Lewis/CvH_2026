# Supplementary Enrichment Gallery -- Data Preparation (F04 CRvH: Concordance)
# Loads pre-computed enrichment CSV, pivots to wide, classifies patterns.
# Saves prep_concordance.rds for downstream viz scripts.
setwd(here::here())
source("04_Figures/F04/a_script/style.R")

pacman::p_load(tidyverse)

DAT <- "04_Figures/F04/CRvH/c_data/supp"
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

# Cancer Recovery concordance (CvH vs Training_CR)
conc_long <- read_csv("04_Figures/F04/CRvH/c_data/panel_supp/enrichment_concordance.csv",
                       show_col_types = FALSE)

conc_wide <- conc_long %>%
  select(pathway, pathway_label, database, contrast, NES, padj, sig, size) %>%
  pivot_wider(
    id_cols     = c(pathway, pathway_label, database),
    names_from  = contrast,
    values_from = c(NES, padj, sig, size)
  ) %>%
  rename(
    NES_CvH  = NES_Cancer_vs_Healthy,
    NES_TR   = NES_Training_CR,
    padj_CvH = padj_Cancer_vs_Healthy,
    padj_TR  = padj_Training_CR,
    sig_CvH  = sig_Cancer_vs_Healthy,
    sig_TR   = sig_Training_CR
  ) %>%
  mutate(
    sig_CvH = replace_na(sig_CvH, FALSE),
    sig_TR  = replace_na(sig_TR,  FALSE),
    pattern = case_when(
      sig_CvH & sig_TR & sign(NES_CvH) == sign(NES_TR) ~ "Concordant",
      sig_CvH & sig_TR & sign(NES_CvH) != sign(NES_TR) ~ "Discordant",
      sig_CvH & !sig_TR                                 ~ "Cancer-specific",
      !sig_CvH & sig_TR                                 ~ "Training-specific",
      TRUE                                               ~ "Other"
    ),
    pattern = factor(pattern, levels = c("Concordant", "Discordant",
                                          "Cancer-specific", "Training-specific",
                                          "Other")),
    bio_theme = classify_pathway_func(pathway),
    bio_theme = factor(bio_theme, levels = CONSOLIDATED_PATHWAY_ORDER),
    set_size  = coalesce(size_Cancer_vs_Healthy, size_Training_CR)
  )

cat(sprintf("F04 CRvH concordance: %d pathways\n", nrow(conc_wide)))
cat("Pattern counts:\n")
print(table(conc_wide$pattern))

saveRDS(conc_wide, file.path(DAT, "prep_concordance.rds"))

cat("\nData prep complete. RDS file saved to", DAT, "\n")
