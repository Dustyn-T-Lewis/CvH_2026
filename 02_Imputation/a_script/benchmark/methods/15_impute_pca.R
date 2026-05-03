# methods/15_impute_pca.R
# imputePCA — regularized iterative PCA
# Josse & Husson 2016
# Pure MAR (low-rank)

impute_imputePCA <- function(mat, meta, is_mnar, ...) {
  requireNamespace("missMDA", quietly = TRUE)
  # missMDA expects samples x variables (transposed)
  result <- missMDA::imputePCA(t(mat), ncp = 2, method = "Regularized")
  imp <- t(result$completeObs)
  dimnames(imp) <- dimnames(mat)
  imp
}
