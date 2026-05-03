# methods/07_imp4p_rf.R
# imp4p RF — condition-aware random forest imputation
# Giai Gianetto et al. 2020
# Pure MAR (proteomics-specific)

impute_imp4p_RF <- function(mat, meta, is_mnar, ...) {
  requireNamespace("imp4p", quietly = TRUE)
  conditions <- factor(meta$Group_Time)
  imp <- imp4p::impute.RF(tab = mat, conditions = conditions)
  # imp4p may leave residual NAs — QRILC fallback
  if (any(is.na(imp))) {
    imp[is.na(imp)] <- MsCoreUtils::impute_matrix(mat, method = "QRILC")[is.na(imp)]
  }
  imp
}
