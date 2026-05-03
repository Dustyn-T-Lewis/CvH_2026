# methods/04_mice.R
# MICE/PMM — predictive mean matching
# van Buuren & Groothuis-Oudshoorn 2011
# Pure MAR assumption
# NOTE: Full MICE on 2138 proteins is infeasible.
# Strategy: impute each protein with missing values using top-10 correlated
# complete proteins as predictors, via mice with pmm.

impute_MICE <- function(mat, meta, is_mnar, ...) {
  requireNamespace("mice", quietly = TRUE)
  imp <- mat
  has_na <- which(rowSums(is.na(mat)) > 0 & rowSums(is.na(mat)) < ncol(mat))
  complete <- which(rowSums(is.na(mat)) == 0)

  if (length(complete) < 5) {
    # Fallback: not enough complete proteins for correlation
    return(MsCoreUtils::impute_matrix(mat, method = "knn"))
  }

  # Precompute correlation of incomplete proteins with complete ones
  # Use pairwise.complete.obs for correlation
  for (i in has_na) {
    obs_cols <- which(!is.na(mat[i, ]))
    na_cols <- which(is.na(mat[i, ]))

    # Find top 10 most correlated complete proteins
    cors <- cor(mat[i, obs_cols], t(mat[complete, obs_cols, drop = FALSE]),
                use = "pairwise.complete.obs")
    cors[is.na(cors)] <- 0
    top10 <- complete[order(abs(cors), decreasing = TRUE)[1:min(10, length(complete))]]

    # Build small data frame: target + predictors, samples as rows
    small_df <- as.data.frame(t(mat[c(i, top10), , drop = FALSE]))
    colnames(small_df) <- c("target", paste0("p", seq_along(top10)))

    mice_fit <- tryCatch(
      mice::mice(small_df, m = 1, maxit = 3, method = "pmm",
                 printFlag = FALSE, seed = 42),
      error = function(e) NULL
    )
    if (!is.null(mice_fit)) {
      completed <- mice::complete(mice_fit, 1)
      imp[i, na_cols] <- completed$target[na_cols]
    } else {
      # Fallback: row mean
      imp[i, na_cols] <- mean(mat[i, ], na.rm = TRUE)
    }
  }
  imp
}
