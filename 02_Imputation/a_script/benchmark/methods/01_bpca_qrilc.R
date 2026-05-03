# methods/01_bpca_qrilc.R
# BPCA+QRILC hybrid imputation
# Oba et al. 2003; Lazar et al. 2016
# Assumption: MAR proteins -> BPCA, MNAR proteins -> QRILC

impute_BPCA_QRILC <- function(mat, meta, is_mnar, ...) {
  imp <- mat
  has_na <- rowSums(is.na(mat)) > 0
  # is_mnar may not cover all rows (complete proteins = NA in is_mnar)
  # Default: incomplete proteins not in is_mnar are treated as MAR
  mnar_flag <- is_mnar[rownames(mat)]
  mnar_flag[is.na(mnar_flag)] <- FALSE

  mar_idx <- which(has_na & !mnar_flag)
  mnar_idx <- which(has_na & mnar_flag)

  if (length(mar_idx) > 0) {
    imp[mar_idx, ] <- MsCoreUtils::impute_matrix(mat[mar_idx, , drop = FALSE], method = "bpca")
  }
  if (length(mnar_idx) > 0) {
    imp[mnar_idx, ] <- MsCoreUtils::impute_matrix(mat[mnar_idx, , drop = FALSE], method = "QRILC")
  }
  # Fallback for any remaining NAs (e.g., BPCA failure on high-missingness rows)
  if (any(is.na(imp))) {
    still_na <- which(is.na(imp), arr.ind = TRUE)
    fallback <- MsCoreUtils::impute_matrix(mat, method = "MinProb")
    imp[is.na(imp)] <- fallback[is.na(imp)]
  }
  imp
}
