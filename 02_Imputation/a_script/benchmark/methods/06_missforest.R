# methods/06_missforest.R
# missForest — random forest iterative imputation
# Stekhoven & Buhlmann 2012
# Pure MAR (nonparametric)

impute_missForest <- function(mat, meta, is_mnar, ...) {
  requireNamespace("missForest", quietly = TRUE)
  # missForest expects samples x features (transposed)
  result <- missForest::missForest(t(mat), verbose = FALSE, maxiter = 10)
  t(result$ximp)
}
