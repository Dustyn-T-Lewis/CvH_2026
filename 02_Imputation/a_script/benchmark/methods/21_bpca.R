# methods/21_bpca.R
# BPCA standalone — Bayesian PCA imputation
# Oba et al. 2003
# Pure MAR; no QRILC component (cf. BPCA_QRILC hybrid)

impute_bpca_standalone <- function(mat, meta, is_mnar, ...) {
  MsCoreUtils::impute_matrix(mat, method = "bpca")
}
