# methods/14_half_minimum.R
# Half-minimum — naive MNAR baseline
# Replace NAs with 0.5 * per-protein minimum observed value
# Negative control

impute_Half_minimum <- function(mat, meta, is_mnar, ...) {
  imp <- mat
  for (i in which(rowSums(is.na(mat)) > 0)) {
    obs_vals <- mat[i, !is.na(mat[i, ])]
    if (length(obs_vals) == 0) {
      imp[i, is.na(mat[i, ])] <- min(mat, na.rm = TRUE) * 0.5
    } else {
      imp[i, is.na(mat[i, ])] <- min(obs_vals) * 0.5
    }
  }
  imp
}
