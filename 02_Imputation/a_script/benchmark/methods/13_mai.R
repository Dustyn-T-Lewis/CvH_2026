# methods/13_mai.R
# MAI — mechanism-aware imputation
# Dekermanjian et al. 2022
# Auto-classifies per protein, routes to best model

impute_MAI <- function(mat, meta, is_mnar, ...) {
  requireNamespace("MAI", quietly = TRUE)
  result <- MAI::MAI(mat, MCAR_algorithm = "BPCA", MNAR_algorithm = "nsKNN",
                     verbose = FALSE)
  # MAI returns a list with Imputed_data
  imp <- result$Imputed_data
  if (is.matrix(imp)) {
    dimnames(imp) <- dimnames(mat)
    return(imp)
  }
  # If SummarizedExperiment, extract assay
  if (methods::is(imp, "SummarizedExperiment")) {
    imp <- SummarizedExperiment::assay(imp)
  }
  imp <- as.matrix(imp)
  dimnames(imp) <- dimnames(mat)
  imp
}
