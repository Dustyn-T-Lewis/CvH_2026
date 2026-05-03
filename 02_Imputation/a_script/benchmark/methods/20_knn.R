# methods/20_knn.R
# KNN standalone — k-nearest neighbors imputation (k=10)
# Troyanskaya et al. 2001
# Pure MAR; no QRILC component (cf. KNN_QRILC hybrid)

impute_knn_standalone <- function(mat, meta, is_mnar, ...) {
  MsCoreUtils::impute_matrix(mat, method = "knn")
}
