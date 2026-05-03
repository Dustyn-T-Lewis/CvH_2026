# methods/17_mindet.R
# MinDet — deterministic minimum imputation
# Lazar et al. 2016
# Pure MNAR; replaces NAs with per-sample minimum detected value

impute_MinDet <- function(mat, meta, is_mnar, ...) {
  MsCoreUtils::impute_matrix(mat, method = "MinDet")
}
