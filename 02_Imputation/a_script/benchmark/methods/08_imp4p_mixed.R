# methods/08_imp4p_mixed.R
# imp4p mixed — hybrid: routes MAR->MLE, MNAR->IGCDA by MCAR probability
# Giai Gianetto et al. 2020
# Hybrid; classifier-dependent

impute_imp4p_mixed <- function(mat, meta, is_mnar, ...) {
  requireNamespace("imp4p", quietly = TRUE)
  conditions <- factor(meta$Group_Time)
  # Build per-protein MCAR probability from is_mnar:
  # prob.MCAR = 1 for MAR proteins, 0 for MNAR (controls routing)
  prob <- matrix(ifelse(is_mnar[rownames(mat)], 0, 1),
                 nrow = nrow(mat), ncol = ncol(mat))
  imp <- imp4p::impute.mix(tab = mat, conditions = conditions,
                           prob.MCAR = prob, threshold = 0.5,
                           methodMCAR = "mle", methodMNAR = "igcda",
                           progress.bar = FALSE)
  # Fallback for residual NAs
  if (any(is.na(imp))) {
    imp[is.na(imp)] <- MsCoreUtils::impute_matrix(mat, method = "QRILC")[is.na(imp)]
  }
  imp
}
