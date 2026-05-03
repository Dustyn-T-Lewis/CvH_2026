# methods/19_zero.R
# Zero imputation — replace NAs with 0
# Baseline/negative control

impute_zero <- function(mat, meta, is_mnar, ...) {
  MsCoreUtils::impute_matrix(mat, method = "zero")
}
