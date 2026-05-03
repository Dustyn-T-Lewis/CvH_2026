# methods/10_dreamai.R
# DreamAI — ensemble imputation
# Ma et al. 2021
# Averages KNN + MissForest + other available methods

impute_DreamAI <- function(mat, meta, is_mnar, ...) {
  # DreamAI ensemble: average of multiple imputation methods
  # If DreamAI package is available, use it; otherwise replicate ensemble logic
  if (requireNamespace("DreamAI", quietly = TRUE)) {
    result <- tryCatch({
      DreamAI::DreamAI(
        data = as.data.frame(mat),
        k = 10, maxiter_MF = 10, ntree = 100,
        maxnodes = NULL, maxiter_ADMIN = 30,
        tol = 10^(-2), gamma_ADMIN = 0,
        gamma = 50, CV = FALSE, fillmethod = "row_mean",
        maxiter_RegImpute = 10, conv_nrmse = 1e-6,
        iter_SpectroFM = 40, method = c("KNN", "MissForest", "RegImpute"),
        out = "Ensemble"
      )$Ensemble
    }, error = function(e) NULL)
    if (!is.null(result)) {
      result <- as.matrix(result)
      dimnames(result) <- dimnames(mat)
      return(result)
    }
  }

  # Fallback: manual ensemble of KNN + missForest
  imp_knn <- MsCoreUtils::impute_matrix(mat, method = "knn")
  imp_mf <- t(missForest::missForest(t(mat), verbose = FALSE, maxiter = 10)$ximp)
  imp <- (imp_knn + imp_mf) / 2
  dimnames(imp) <- dimnames(mat)
  imp
}
