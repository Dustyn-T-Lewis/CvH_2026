# reversal_inputs.R  --  back-compat adapter for the Reversal panels.
#
# The panels were written against the old wide-format DEP
# (03_DEP/c_data/03_combined_results_CRvH.csv) and a removed 02_Imputation stage.
# This adapter rebuilds the SAME objects (dep_df with old column names, dal,
# imputation_df) from the new proteoDA-native pipeline so the panel bodies need
# no edits. Axis mapping: Cancer_vs_Healthy = CRvH_Baseline (D),
# Training_CR = CR_Training (T). Source after setwd(root).

source("04_Figures/Reversal/a_script/reversal.R")   # load_reversal_table + dplyr/tidyr/readr/tibble

# --- dep_df: old wide column names from the new long DEP ---------------------
.rev_wide <- load_reversal_table("03_DEP/a_non_imputed/c_data/combined_results_pi.csv")
dep_df <- .rev_wide |>
  dplyr::transmute(
    uniprot_id, gene, protein, description,
    logFC_Cancer_vs_Healthy   = logFC_D, logFC_Training_CR   = logFC_T,
    t_Cancer_vs_Healthy       = t_D,     t_Training_CR       = t_T,
    pi_score_Cancer_vs_Healthy = pi_score_D, pi_score_Training_CR = pi_score_T,
    P.Value_Cancer_vs_Healthy = P.Value_D, P.Value_Training_CR = P.Value_T,
    logFC_Resid = logFC_R, t_Resid = t_R, pi_score_Resid = pi_score_R,
    P.Value_Resid = P.Value_R
  )

# --- dal: imp4p-imputed DAList (complete matrix for fry/rotation panels) ------
dal <- readRDS("02_Normalization/imputation/c_data/DAList_imputed_imp4p.rds")

# --- imputation_df: gene -> imputed flag (had >= 1 missing pre-imputation) ----
# Replaces the old MAR/MNAR "Complete" classification: a protein is "imputed"
# in figure terms if it carried any NA in the normalized (pre-imputation) matrix.
.norm <- readRDS("02_Normalization/c_data/DAList_normalized.rds")
imputation_df <- tibble::tibble(
  uniprot_id = rownames(.norm$data),
  imputed    = rowSums(is.na(.norm$data)) > 0
) |>
  dplyr::left_join(dplyr::select(.rev_wide, uniprot_id, gene), by = "uniprot_id") |>
  dplyr::filter(!is.na(gene)) |>
  dplyr::group_by(gene) |>
  dplyr::summarise(imputed = any(imputed), .groups = "drop")

rm(.rev_wide, .norm)
