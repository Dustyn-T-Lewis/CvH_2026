# methods/23_svd.R
# SVD imputation — singular value decomposition
# Troyanskaya et al. 2001
# Pure MAR; low-rank via pcaMethods

impute_SVD <- function(mat, meta, is_mnar, ...) {
  requireNamespace("pcaMethods", quietly = TRUE)
  pcaMethods::pca(mat, method = "svdImpute", nPcs = 2, verbose = FALSE)@completeObs
}
