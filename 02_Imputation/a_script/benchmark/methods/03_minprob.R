# methods/03_minprob.R
# MinProb — pure MNAR imputation from left tail (q=0.01)
# Lazar et al. 2016

impute_MinProb <- function(mat, meta, is_mnar, ...) {
  MsCoreUtils::impute_matrix(mat, method = "MinProb")
}
