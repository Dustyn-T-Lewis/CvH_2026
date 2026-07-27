library(testthat)
source(here::here("04_Figures/shared/prediction_utils.R"))

make_separable <- function(n = 20, p = 6, effect = 3, seed = 1) {
  set.seed(seed)
  y <- rep(c(0, 1), length.out = n)
  x <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("f", seq_len(p))))
  x[, 1] <- x[, 1] + effect * y
  list(y = y, x = x)
}

make_noise <- function(n = 20, p = 6, seed = 2) {
  set.seed(seed)
  list(
    y = rep(c(0, 1), length.out = n),
    x = matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("f", seq_len(p))))
  )
}

test_that("run_topk_loocv recovers a planted signal", {
  d <- make_separable()
  res <- run_topk_loocv(d$y, d$x, k_range = 1:3)
  expect_length(res$scores, length(d$y))
  expect_true(all(res$probs >= 0 & res$probs <= 1))
  expect_gt(loocv_auc(d$y, res$scores), 0.85)
  expect_equal(names(which.max(res$feature_freq)), "f1")
})

# Perfectly separable data drives glm coefficients to +/-Inf, so the fitted probabilities
# collapse onto exactly 0 and 1 and tie. Ranking on the link scale keeps them ordered.
test_that("separable data does not lose AUC to probability saturation", {
  d <- make_separable(effect = 4)
  res <- run_topk_loocv(d$y, d$x, k_range = 2:3)
  expect_true(any(res$probs %in% c(0, 1)))
  expect_false(anyDuplicated(res$scores) > 0)
  expect_gt(loocv_auc(d$y, res$scores), loocv_auc(d$y, round(res$probs)) - 1e-9)
  expect_gt(loocv_auc(d$y, res$scores), 0.9)
})

test_that("shuffled labels land at chance and are not significant", {
  d <- make_noise()
  res <- run_topk_loocv(d$y, d$x, k_range = 2:3)
  auc_obs <- loocv_auc(d$y, res$scores)
  expect_lt(abs(auc_obs - 0.5), 0.25)
  p <- perm_p_classifier(d$y, d$x, k_fixed = 2, obs_auc = auc_obs, n_perm = 60)
  expect_gt(p, 0.05)
})

test_that("in-fold selection never sees the held-out label", {
  d <- make_separable(n = 16, p = 5)
  base <- run_topk_loocv(d$y, d$x, k_range = 2:3)
  flipped_y <- d$y
  flipped_y[3] <- 1 - flipped_y[3]
  flipped <- run_topk_loocv(flipped_y, d$x, k_range = 2:3)
  expect_equal(base$selected[[3]], flipped$selected[[3]])
})

# A paired design leaks if only one of a subject's two samples is held out: the partner
# sample stays in training and carries the subject's level. Grouped hold-out must strip
# that advantage, so a subject-identity-only signal should score at chance.
test_that("grouped hold-out removes the paired-partner leak", {
  set.seed(7)
  n_subj <- 12
  subject <- rep(paste0("S", seq_len(n_subj)), each = 2)
  y <- rep(c(0, 1), times = n_subj)
  subject_level <- rep(rnorm(n_subj, sd = 3), each = 2)
  x <- cbind(f1 = subject_level + rnorm(2 * n_subj, sd = 0.1))

  ungrouped <- run_topk_loocv(y, x, k_range = 1)
  grouped <- run_topk_loocv(y, x, k_range = 1, group = subject)

  expect_lt(abs(loocv_auc(y, grouped$scores) - 0.5), 0.2)
  expect_length(grouped$best_k, n_subj)
  expect_length(ungrouped$best_k, 2 * n_subj)
})

test_that("within_subject_shuffle permutes only inside a subject", {
  subject <- rep(paste0("S", 1:5), each = 2)
  labels <- rep(c(0, 1), times = 5)
  shuffler <- within_subject_shuffle(subject)
  set.seed(3)
  out <- shuffler(labels)
  expect_equal(
    tapply(out, subject, sum),
    tapply(labels, subject, sum)
  )
  expect_setequal(out, labels)
})

test_that("perm_p_paired is direction-agnostic and bounded", {
  set.seed(4)
  subject <- rep(paste0("S", 1:8), 2)
  tp <- c(rep(0, 8), rep(1, 8))
  x <- c(rnorm(8), rnorm(8, 1.5))
  obs <- uni_auc(tp, x)$auc
  p <- perm_p_paired(tp, subject, x, obs, n_perm = 200)
  expect_gte(p, 1 / 201)
  expect_lte(p, 1)
  expect_lt(p, 0.05)
})

test_that("fit_kfold_perm separates signal from noise", {
  set.seed(5)
  n <- 200
  y <- rep(c(0, 1), length.out = n)
  x <- cbind(
    signal = rnorm(n) + 2 * y,
    noise = rnorm(n)
  )
  res <- fit_kfold_perm(x, y, n_perm = 40)
  expect_gt(res$auc, 0.8)
  expect_lt(res$perm_p, 0.05)
  expect_length(res$probs, n)
})

test_that("loso_me projects the held-out sample using training-fold scaling", {
  set.seed(6)
  mat <- matrix(rnorm(20 * 8), 20, 8)
  proj <- loso_me(mat, train_rows = 2:20, holdout_rows = 1)
  expect_length(proj, 1)
  expect_true(is.finite(proj))
})

test_that("loocv_auc reports an anti-predictive score below chance", {
  y <- rep(c(0, 1), each = 15)
  anti <- ifelse(y == 1, -1, 1) + rnorm(30, sd = 0.1)
  expect_lt(loocv_auc(y, anti), 0.1)
  expect_gt(loocv_auc(y, -anti), 0.9)
})

test_that("match_modules scores identical partitions at Jaccard 1", {
  mods <- setNames(rep(c("blue", "brown"), each = 10), paste0("p", 1:20))
  out <- match_modules(mods, mods)
  expect_setequal(out$full, c("blue", "brown"))
  expect_true(all(out$jaccard == 1))
})
