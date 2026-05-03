# methods/16_non_imputed.R
# Non-imputed — reference method (identity)
# Returns matrix with NAs intact; limma handles per-row
# Karpievitch et al. 2012

impute_Non_imputed <- function(mat, meta, is_mnar, ...) {
  mat  # return as-is
}
