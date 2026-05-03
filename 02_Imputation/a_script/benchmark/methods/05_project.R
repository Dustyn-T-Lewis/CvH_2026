# methods/05_project.R
# ProJect — correlation-weighted skew-normal prediction
# Webb-Robertson et al. 2015 (PMID 25987570)
# Model-free; no MAR/MNAR assumption
#
# Algorithm: For each protein with missing values:
#   1. Fit a skew-normal to its observed values (via sn::selm)
#   2. Find top-K correlated proteins that are complete across missing samples
#   3. For each missing sample, predict via correlation-weighted average of
#      donor proteins, then adjust to match the skew-normal marginal

impute_ProJect <- function(mat, meta, is_mnar, ...) {
  requireNamespace("sn", quietly = TRUE)
  imp <- mat
  n_prot <- nrow(mat)
  n_samp <- ncol(mat)
  K <- 20  # number of donor proteins

  # Pre-compute pairwise correlations (observed values only)
  # Use complete observations between each pair
  has_na_rows <- which(rowSums(is.na(mat)) > 0 & rowSums(is.na(mat)) < n_samp)
  complete_rows <- which(rowSums(is.na(mat)) == 0)

  if (length(complete_rows) < 5) {
    # Not enough complete proteins: fall back to row-mean imputation
    for (i in has_na_rows) {
      imp[i, is.na(mat[i, ])] <- mean(mat[i, ], na.rm = TRUE)
    }
    return(imp)
  }

  # Complete protein matrix for donor lookup
  complete_mat <- mat[complete_rows, , drop = FALSE]

  for (i in has_na_rows) {
    na_cols  <- which(is.na(mat[i, ]))
    obs_cols <- which(!is.na(mat[i, ]))
    obs_vals <- mat[i, obs_cols]

    if (length(obs_vals) < 4) {
      imp[i, na_cols] <- mean(obs_vals)
      next
    }

    # Step 1: Fit skew-normal to observed values
    sn_fit <- tryCatch({
      fit <- sn::selm(obs_vals ~ 1, family = "SN")
      cp <- sn::coef(fit, param.type = "CP")
      # CP = (mean, s.d., gamma1) in "centred parameterisation"
      list(mu = cp[1], sigma = cp[2], gamma1 = cp[3])
    }, error = function(e) {
      # Fall back to normal if sn fit fails
      list(mu = mean(obs_vals), sigma = sd(obs_vals), gamma1 = 0)
    })

    # Step 2: Find top-K correlated complete proteins
    cors <- cor(mat[i, obs_cols], t(complete_mat[, obs_cols, drop = FALSE]),
                use = "pairwise.complete.obs")
    cors[is.na(cors)] <- 0
    cors <- as.vector(cors)

    k_use <- min(K, length(complete_rows))
    top_order <- order(abs(cors), decreasing = TRUE)[1:k_use]
    top_idx <- complete_rows[top_order]
    top_cors <- cors[top_order]

    # Weights: absolute correlation, softmax-style
    w <- abs(top_cors)
    if (sum(w) < 1e-10) {
      # No correlated donors: draw from fitted marginal
      for (j in na_cols) {
        imp[i, j] <- sn_fit$mu
      }
      next
    }
    w <- w / sum(w)

    # Step 3: For each missing sample, predict via weighted donor average
    for (j in na_cols) {
      donor_vals <- mat[top_idx, j]
      raw_pred <- sum(donor_vals * w)

      # Step 4: Standardize donor prediction to match target protein's
      # skew-normal marginal. Compute donor ensemble mean/sd across all samples.
      donor_means <- rowMeans(complete_mat[top_order, , drop = FALSE])
      donor_ensemble_mean <- sum(donor_means * w)
      donor_sds <- apply(complete_mat[top_order, , drop = FALSE], 1, sd)
      donor_ensemble_sd <- sum(donor_sds * w)

      if (donor_ensemble_sd > 1e-10 && sn_fit$sigma > 1e-10) {
        # Z-score in donor space, then project to target space
        z <- (raw_pred - donor_ensemble_mean) / donor_ensemble_sd
        imp[i, j] <- sn_fit$mu + z * sn_fit$sigma
      } else {
        imp[i, j] <- raw_pred
      }
    }
  }

  # Handle any completely missing proteins (all samples NA)
  all_na <- which(rowSums(is.na(mat)) == n_samp)
  for (i in all_na) {
    imp[i, ] <- median(mat, na.rm = TRUE)
  }

  imp
}
