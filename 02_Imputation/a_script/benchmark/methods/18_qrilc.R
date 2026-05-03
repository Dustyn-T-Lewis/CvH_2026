# methods/18_qrilc.R
# QRILC — quantile regression imputation of left-censored data
# Lazar et al. 2016
# Pure MNAR; draws from fitted left tail per sample

impute_QRILC <- function(mat, meta, is_mnar, ...) {
  MsCoreUtils::impute_matrix(mat, method = "QRILC")
}
