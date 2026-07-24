# Native data load for the F04 reversal panels. Builds dep_df (per-protein wide,
# with native contrast suffixes), dal (imp4p DAList for the rotation panels) and
# imputation_df (gene -> imputed flag) directly from the proteoDA-native
# pipeline. Source after setwd(root).

source("04_Figures/F04_Reversal/a_script/reversal.R")

.rev_wide <- load_reversal_table("03_DEP/a_non_imputed/c_data/combined_results_pi.csv")
dep_df <- .rev_wide |>
  dplyr::transmute(
    uniprot_id, gene, protein, description,
    logFC_CRvH_Baseline = logFC_D, logFC_CR_Training = logFC_T,
    t_CRvH_Baseline = t_D, t_CR_Training = t_T,
    pi_score_CRvH_Baseline = pi_score_D, pi_score_CR_Training = pi_score_T,
    P.Value_CRvH_Baseline = P.Value_D, P.Value_CR_Training = P.Value_T,
    logFC_Resid = logFC_R, t_Resid = t_R, pi_score_Resid = pi_score_R,
    P.Value_Resid = P.Value_R
  )

# adj.P.Val is not carried by load_reversal_table; join it per axis.
.adj <- readr::read_csv(
  "03_DEP/a_non_imputed/c_data/combined_results_pi.csv",
  show_col_types = FALSE
) |>
  dplyr::filter(contrast %in% c("CRvH_Baseline", "CR_Training", "Resid")) |>
  dplyr::select(uniprot_id, contrast, adj.P.Val) |>
  tidyr::pivot_wider(
    names_from = contrast, values_from = adj.P.Val,
    names_prefix = "adj.P.Val_"
  )
dep_df <- dplyr::left_join(dep_df, .adj, by = "uniprot_id")
rm(.adj)

dal <- readRDS("02_Normalization/imputation/c_data/DAList_imputed_imp4p.rds")

.norm <- readRDS("02_Normalization/c_data/DAList_normalized.rds")
imputation_df <- tibble::tibble(
  uniprot_id = rownames(.norm$data),
  imputed = rowSums(is.na(.norm$data)) > 0
) |>
  dplyr::left_join(dplyr::select(.rev_wide, uniprot_id, gene), by = "uniprot_id") |>
  dplyr::filter(!is.na(gene)) |>
  dplyr::group_by(gene) |>
  dplyr::summarise(imputed = any(imputed), .groups = "drop")

rm(.rev_wide, .norm)
