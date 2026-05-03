# methods/22_rf.R
# RF (MsCoreUtils) — random forest via MsCoreUtils wrapper
# Stekhoven & Buhlmann 2012
# Pure MAR; cf. imp4p_RF (condition-aware) and missForest (direct)

impute_RF_MsCoreUtils <- function(mat, meta, is_mnar, ...) {
  MsCoreUtils::impute_matrix(mat, method = "RF")
}
