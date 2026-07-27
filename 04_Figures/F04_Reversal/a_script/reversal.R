# reversal.R -- the three reversal axes, pivoted onto one row per protein.
#
# Sourced by f04_data.R, which assembles the inputs every F04 panel reads.
# The reversal statistics themselves live in the panels that compute them;
# docs/decisions.md carries the shared-baseline design note and method lineage.
#
#   D = "CRvH_Baseline"  disease axis   (CR_pre  - H_pre)
#   T = "CR_Training"    training axis  (CR_post - CR_pre)
#   R = "Resid"          residual axis  (CR_post - H_pre = D + T)

pacman::p_load(dplyr, tidyr, tibble, purrr, readr)

REVERSAL_CONTRASTS <- c(D = "CRvH_Baseline", T = "CR_Training", R = "Resid")

#' Pivot the long combined_results_pi.csv into a per-protein wide frame.
#' Keeps logFC, t, P.Value, pi_score and the Pi-significance flag for each axis.
load_reversal_table <- function(combined_pi_path,
                                contrasts = REVERSAL_CONTRASTS) {
  long <- readr::read_csv(combined_pi_path, show_col_types = FALSE)
  stopifnot(all(c(
    "uniprot_id", "logFC", "t", "P.Value", "pi_score",
    "sig_pi", "contrast", "gene"
  ) %in% names(long)))
  long <- long |> filter(contrast %in% contrasts)

  # dplyr:: qualified -- AnnotationDbi (pulled in by RRHO2/GO.db) masks select()
  ann <- long |>
    dplyr::select(any_of(c("uniprot_id", "gene", "protein", "description"))) |>
    distinct(uniprot_id, .keep_all = TRUE)

  vals <- long |>
    mutate(axis = names(contrasts)[match(contrast, contrasts)]) |>
    dplyr::select(uniprot_id, axis, logFC, t, P.Value, pi_score, sig_pi) |>
    pivot_wider(
      names_from = axis,
      values_from = c(logFC, t, P.Value, pi_score, sig_pi),
      names_sep = "_"
    )

  ann |> inner_join(vals, by = "uniprot_id")
}
