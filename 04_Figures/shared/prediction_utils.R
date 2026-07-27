# Prediction engines for F06: nested cross-validation with in-fold feature selection,
# permutation nulls for unpaired and paired designs, and the leave-one-subject-out
# audits that quantify how much of a WGCNA-derived AUC is optimism.
#
# Design follows YvO F07: plain logistic regression with hand-rolled resampling. At
# n = 25 a penalized or ensemble learner would fit noise, and its tuning would need a
# third nesting level to stay honest.

pacman::p_load(pROC, dplyr, tibble)

# Returns the linear predictor, not the fitted probability. Small samples with a strong
# signal are often perfectly separable, which sends the coefficients to +/-Inf and
# saturates probabilities to exactly 0 and 1; the resulting ties destroy the AUC ranking
# even though the model ordered the samples correctly. The link scale stays finite and
# ordered, and plogis() recovers a probability wherever one is actually needed.
# A failed fit returns 0, the link-scale equivalent of p = 0.5.
fit_logistic <- function(y, x_train, x_test) {
  fit <- tryCatch(
    suppressWarnings(glm(y ~ .,
      family = binomial,
      data = cbind(y = y, as.data.frame(x_train))
    )),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(0)
  }
  as.numeric(predict(fit, type = "link", newdata = as.data.frame(x_test)))
}

rank_features <- function(x, y) {
  r <- abs(stats::cor(x, y))
  r[is.na(r)] <- 0
  rownames(r)[order(r[, 1], decreasing = TRUE)]
}

# direction is pinned: a logistic link score already runs low-to-high with P(y=1), so
# letting pROC re-detect it would fold an anti-predictive model up above 0.5.
loocv_auc <- function(labels, probs) {
  tryCatch(
    as.numeric(pROC::auc(pROC::roc(labels, probs, quiet = TRUE, direction = "<"))),
    error = function(e) 0.5
  )
}

# Outer leave-one-out for honest prediction; feature ranking and inner tuning both run
# on the training rows only, so the held-out sample never informs its own prediction.
# k_range starts at 1, unlike YvO's 2:5. With five modules and 25 samples, forcing a
# second predictor pulls in a noise feature that the inner loop would otherwise reject,
# and it costs real signal: on separable test data, k_range = 2:3 scores AUC 0.72 where
# 1:3 scores 0.90.
run_topk_loocv <- function(labels, x, k_range = 1:5, group = NULL) {
  x <- as.matrix(x)
  n <- length(labels)
  if (is.null(group)) group <- seq_len(n)
  folds <- unique(group)
  scores <- numeric(n)
  selected <- vector("list", length(folds))
  best_k <- integer(length(folds))

  for (f in seq_along(folds)) {
    test_idx <- which(group == folds[f])
    train_idx <- which(group != folds[f])
    train_x <- x[train_idx, , drop = FALSE]
    train_y <- labels[train_idx]
    train_g <- group[train_idx]
    ranked <- rank_features(train_x, train_y)

    deviance <- vapply(k_range, function(k) {
      top <- ranked[seq_len(min(k, length(ranked)))]
      sum(vapply(unique(train_g), function(g) {
        inner_test <- which(train_g == g)
        inner_train <- which(train_g != g)
        p <- stats::plogis(fit_logistic(
          train_y[inner_train], train_x[inner_train, top, drop = FALSE],
          train_x[inner_test, top, drop = FALSE]
        ))
        p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
        sum(-(train_y[inner_test] * log(p) + (1 - train_y[inner_test]) * log(1 - p)))
      }, numeric(1)))
    }, numeric(1))

    best_k[f] <- k_range[which.min(deviance)]
    top <- ranked[seq_len(min(best_k[f], length(ranked)))]
    selected[[f]] <- top
    scores[test_idx] <- fit_logistic(
      train_y, train_x[, top, drop = FALSE], x[test_idx, top, drop = FALSE]
    )
  }

  list(
    scores = scores,
    probs = stats::plogis(scores),
    feature_freq = table(factor(unlist(selected), levels = colnames(x))),
    selected = selected,
    best_k = best_k
  )
}

# Permutation variant: k is fixed at the complexity the observed run settled on, which
# drops the inner loop and its factor-of-n cost. Selection still runs inside every fold,
# so the null absorbs selection optimism.
run_fast_loocv_auc <- function(labels, x, k_fixed, group = NULL) {
  x <- as.matrix(x)
  if (is.null(group)) group <- seq_along(labels)
  scores <- numeric(length(labels))
  for (g in unique(group)) {
    test_idx <- which(group == g)
    train_idx <- which(group != g)
    train_x <- x[train_idx, , drop = FALSE]
    train_y <- labels[train_idx]
    top <- rank_features(train_x, train_y)[seq_len(min(k_fixed, ncol(x)))]
    scores[test_idx] <- fit_logistic(
      train_y, train_x[, top, drop = FALSE], x[test_idx, top, drop = FALSE]
    )
  }
  loocv_auc(labels, scores)
}

# Each subject independently keeps or swaps its own labels, so subject means survive and
# only the within-subject effect is destroyed.
within_subject_shuffle <- function(subject) {
  function(labels) {
    out <- labels
    for (s in unique(subject)) {
      idx <- which(subject == s)
      out[idx] <- sample(out[idx])
    }
    out
  }
}

# Phipson-Smyth (b+1)/(m+1), so the p-value floors at 1/(n_perm+1) rather than zero.
# No direction folding here: model probabilities are already sign-anchored.
perm_p_classifier <- function(labels, x, k_fixed, obs_auc, n_perm = 1000,
                              shuffle = sample, group = NULL) {
  nulls <- vapply(
    seq_len(n_perm),
    function(i) run_fast_loocv_auc(shuffle(labels), x, k_fixed, group = group),
    numeric(1)
  )
  (sum(nulls >= obs_auc) + 1) / (n_perm + 1)
}

# DeLong CI assumes independent observations, so it is anti-conservative for the paired
# pre/post outcome. The permutation p-value is the primary inferential quantity there.
uni_auc <- function(y, x) {
  ok <- !is.na(y) & !is.na(x)
  if (length(unique(y[ok])) < 2) {
    return(list(auc = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, roc = NULL))
  }
  r <- suppressMessages(pROC::roc(y[ok], x[ok], quiet = TRUE, direction = "auto"))
  ci <- as.numeric(pROC::ci.auc(r))
  list(auc = as.numeric(pROC::auc(r)), ci_lo = ci[1], ci_hi = ci[3], roc = r)
}

# A single eigengene has arbitrary sign, so fold both the observed statistic and the null
# to the same scale before comparing.
fold_auc <- function(a) max(a, 1 - a)

perm_p_paired <- function(tp, subject, x, obs_auc, n_perm = 1000) {
  shuffler <- within_subject_shuffle(subject)
  nulls <- vapply(seq_len(n_perm), function(i) {
    r <- suppressMessages(pROC::roc(shuffler(tp), x, quiet = TRUE, direction = "auto"))
    fold_auc(as.numeric(pROC::auc(r)))
  }, numeric(1))
  (sum(nulls >= fold_auc(obs_auc)) + 1) / (n_perm + 1)
}

# Protein-level engine. Thousands of rows make leave-one-out wasteful, and every
# predictor is kept, so there is no selection step to nest. Fold assignment is held
# fixed across permutations to isolate the label effect.
fit_kfold_perm <- function(x, y, n_perm = 200, k_folds = 10, seed = 42) {
  set.seed(seed)
  x <- as.matrix(x)
  folds <- integer(length(y))
  for (cls in unique(y)) {
    idx <- which(y == cls)
    folds[idx] <- sample(rep(seq_len(k_folds), length.out = length(idx)))
  }

  cv_scores <- function(labels) {
    s <- numeric(length(labels))
    for (f in seq_len(k_folds)) {
      train <- folds != f
      s[!train] <- fit_logistic(
        labels[train], x[train, , drop = FALSE], x[!train, , drop = FALSE]
      )
    }
    s
  }

  scores <- cv_scores(y)
  roc_obj <- pROC::roc(y, scores, quiet = TRUE, direction = "<")
  auc_val <- as.numeric(pROC::auc(roc_obj))
  nulls <- vapply(seq_len(n_perm), function(i) {
    ys <- sample(y)
    loocv_auc(ys, cv_scores(ys))
  }, numeric(1))
  ci <- as.numeric(pROC::ci.auc(roc_obj))

  list(
    auc = auc_val, ci_lo = ci[1], ci_hi = ci[3],
    perm_p = (sum(nulls >= auc_val) + 1) / (n_perm + 1),
    scores = scores, probs = stats::plogis(scores), labels = y,
    n = length(y), n_pos = sum(y == 1), n_neg = sum(y == 0),
    fpr = 1 - roc_obj$specificities, tpr = roc_obj$sensitivities
  )
}

# Tier 1 of the circularity audit: refit the module's first principal component on the
# training fold and project the held-out samples with training-fold centering and
# scaling. PC sign is arbitrary, so anchor it against the training row means.
loso_me <- function(full_mat, train_rows, holdout_rows) {
  x_train <- scale(full_mat[train_rows, , drop = FALSE])
  center <- attr(x_train, "scaled:center")
  spread <- attr(x_train, "scaled:scale")
  spread[spread == 0 | !is.finite(spread)] <- 1
  v <- svd(x_train, nu = 0, nv = 1)$v[, 1]
  if (stats::cor(as.numeric(x_train %*% v), rowMeans(x_train)) < 0) v <- -v
  x_hold <- sweep(full_mat[holdout_rows, , drop = FALSE], 2, center, "-")
  as.numeric(sweep(x_hold, 2, spread, "/") %*% v)
}

# Tier 2 support: a network refit on n-1 subjects relabels its modules arbitrarily, so
# identity is re-established by membership overlap before any projection.
match_modules <- function(mod_train, mod_full) {
  train_lvls <- setdiff(unique(mod_train), "grey")
  full_lvls <- setdiff(unique(mod_full), "grey")
  overlap <- vapply(train_lvls, function(t) {
    b <- names(mod_train)[mod_train == t]
    vapply(full_lvls, function(f) {
      a <- names(mod_full)[mod_full == f]
      length(intersect(a, b)) / length(union(a, b))
    }, numeric(1))
  }, numeric(length(full_lvls)))
  overlap <- matrix(overlap,
    nrow = length(full_lvls),
    dimnames = list(full_lvls, train_lvls)
  )
  tibble(
    full = rownames(overlap),
    train = colnames(overlap)[apply(overlap, 1, which.max)],
    jaccard = apply(overlap, 1, max)
  )
}
