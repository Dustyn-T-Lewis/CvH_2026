# methods/02_knn_qrilc.R
# KNN+QRILC hybrid imputation
# Troyanskaya et al. 2001; Lazar et al. 2016
# Assumption: MAR proteins -> KNN (k=10), MNAR proteins -> QRILC

impute_KNN_QRILC <- function(mat, meta, is_mnar, ...) {
  imp <- mat
  has_na <- rowSums(is.na(mat)) > 0
  mnar_flag <- is_mnar[rownames(mat)]
  mnar_flag[is.na(mnar_flag)] <- FALSE

  mar_idx <- which(has_na & !mnar_flag)
  mnar_idx <- which(has_na & mnar_flag)

  if (length(mar_idx) > 0) {
    imp[mar_idx, ] <- MsCoreUtils::impute_matrix(mat[mar_idx, , drop = FALSE], method = "knn")
  }
  if (length(mnar_idx) > 0) {
    imp[mnar_idx, ] <- MsCoreUtils::impute_matrix(mat[mnar_idx, , drop = FALSE], method = "QRILC")
  }
  # Fallback for any remaining NAs
  if (any(is.na(imp))) {
    fallback <- MsCoreUtils::impute_matrix(mat, method = "MinProb")
    imp[is.na(imp)] <- fallback[is.na(imp)]
  }
  imp
}
