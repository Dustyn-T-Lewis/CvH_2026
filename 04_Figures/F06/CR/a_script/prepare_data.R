# F06/CR: Supplement-Specific Reversal — Data Preparation
# Loads Cancer_vs_Healthy (CRvH model) + Training_CRE/PLA (CR model)
# to test whether each supplement arm reverses cancer effects.
#
# Reversal question:
#   Does Training_CRE reverse Cancer_vs_Healthy? (CRE-specific reversal)
#   Does Training_PLA reverse Cancer_vs_Healthy? (PLA-specific reversal)
#   Are they different? (supplement-moderated reversal)
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F06/a_script/style.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

RPT <- "04_Figures/F06/CR/b_reports"
DAT <- "04_Figures/F06/CR/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

# --- Load DEP results from BOTH models ---
dep_crvh <- read_csv("03_DEP/c_data/03_combined_results_CRvH.csv", show_col_types = FALSE)
dep_cr   <- read_csv("03_DEP/c_data/03_combined_results_CR.csv",   show_col_types = FALSE)

# Merge on gene — keep Cancer_vs_Healthy from CRvH + Training contrasts from CR
rev_df <- dep_crvh %>%
  select(gene, matches("Cancer_vs_Healthy")) %>%
  inner_join(
    dep_cr %>% select(gene, matches("Training_CRE|Training_PLA|Supplement_Interaction")),
    by = "gene"
  )

# --- Reversal classification ---
# CRE reversal: Cancer_vs_Healthy × Training_CRE (anti-concordance)
rev_df <- rev_df %>%
  mutate(
    # CRE reversal quadrants
    cre_quad = case_when(
      logFC_Cancer_vs_Healthy > 0 & logFC_Training_CRE < 0 ~ "Reversed (Cancer Up)",
      logFC_Cancer_vs_Healthy < 0 & logFC_Training_CRE > 0 ~ "Reversed (Cancer Down)",
      logFC_Cancer_vs_Healthy > 0 & logFC_Training_CRE > 0 ~ "Exacerbated Up",
      logFC_Cancer_vs_Healthy < 0 & logFC_Training_CRE < 0 ~ "Exacerbated Down",
      TRUE ~ "NS"
    ),
    # PLA reversal quadrants
    pla_quad = case_when(
      logFC_Cancer_vs_Healthy > 0 & logFC_Training_PLA < 0 ~ "Reversed (Cancer Up)",
      logFC_Cancer_vs_Healthy < 0 & logFC_Training_PLA > 0 ~ "Reversed (Cancer Down)",
      logFC_Cancer_vs_Healthy > 0 & logFC_Training_PLA > 0 ~ "Exacerbated Up",
      logFC_Cancer_vs_Healthy < 0 & logFC_Training_PLA < 0 ~ "Exacerbated Down",
      TRUE ~ "NS"
    )
  )

# --- Reversal statistics ---
# CRE reversal correlation
cre_r <- cor(rev_df$logFC_Cancer_vs_Healthy, rev_df$logFC_Training_CRE,
             use = "complete.obs", method = "spearman")
pla_r <- cor(rev_df$logFC_Cancer_vs_Healthy, rev_df$logFC_Training_PLA,
             use = "complete.obs", method = "spearman")

# Reversal fraction (among Cancer-sig proteins)
cancer_sig <- rev_df %>% filter(pi_score_Cancer_vs_Healthy < 0.05)
n_cancer <- nrow(cancer_sig)
cre_reversed <- sum(grepl("Reversed", cancer_sig$cre_quad))
pla_reversed <- sum(grepl("Reversed", cancer_sig$pla_quad))

cat(sprintf("F06/CR Reversal Summary:\n"))
cat(sprintf("  CRE: rho = %.3f, reversed %d/%d (%.1f%%)\n",
            cre_r, cre_reversed, n_cancer, 100 * cre_reversed / n_cancer))
cat(sprintf("  PLA: rho = %.3f, reversed %d/%d (%.1f%%)\n",
            pla_r, pla_reversed, n_cancer, 100 * pla_reversed / n_cancer))

# Save reversal data for downstream panels
write_csv(rev_df, file.path(DAT, "reversal_merged.csv"))

message("F06/CR prepare_data done")
