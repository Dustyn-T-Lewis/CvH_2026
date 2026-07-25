# F06 supplement: how much of the module-based AUC survives when the modules are not
# allowed to see the held-out subject.
#
# Tier 1 keeps the full-sample module membership and refits only the eigengene's first
# principal component per fold. Tier 2 refits the whole network per fold and re-matches
# modules by membership overlap, which is the stricter test of module-definition
# circularity. Neither tier gates a result: both report the in-sample AUC alongside the
# out-of-fold one, and the drop between them is the quantity of interest. On 24-sample
# refits some folds will fail to reproduce a module at all, and that failure rate is
# itself reported rather than hidden.

setwd(here::here())
pacman::p_load(dplyr, tibble, readr, purrr, WGCNA, pROC)
source("04_Figures/shared/prediction_utils.R")

DAT <- "04_Figures/F06_Prediction/c_data"
F05 <- "04_Figures/F05_WGCNA/c_data"
JACCARD_FLOOR <- 0.4
SOFT_POWER <- 12L
set.seed(42)

bundle <- readRDS(file.path(DAT, "feature_bundle.rds"))
datExpr <- readRDS(file.path(F05, "datExpr.rds"))
module_colors <- readRDS(file.path(F05, "module_colors.rds"))
names(module_colors) <- colnames(datExpr)
modules <- setdiff(unique(module_colors), "grey")

summary_in <- read_csv(file.path(DAT, "classifier_summary.csv"), show_col_types = FALSE) |>
  filter(space == "module")

tier1 <- map_dfr(names(bundle$outcomes), function(outcome_name) {
  out <- bundle$outcomes[[outcome_name]]
  ids <- out$ids
  group <- if (out$paired) out$subject else ids
  map_dfr(modules, function(mod) {
    prot <- names(module_colors)[module_colors == mod]
    mat <- datExpr[ids, prot, drop = FALSE]
    scores <- numeric(length(ids))
    for (g in unique(group)) {
      hold <- which(group == g)
      scores[hold] <- loso_me(mat, setdiff(seq_along(ids), hold), hold)
    }
    obs <- uni_auc(out$y, datExpr[ids, prot, drop = FALSE] |> rowMeans())
    loso <- uni_auc(out$y, scores)
    tibble(
      tier = "tier1_pc_reprojection", outcome = outcome_name, module = mod,
      auc_insample = fold_auc(obs$auc), auc_loso = fold_auc(loso$auc),
      drop = fold_auc(obs$auc) - fold_auc(loso$auc), n = length(ids)
    )
  })
})
write_csv(tier1, file.path(DAT, "circularity_tier1.csv"))

# Tier 2. Each fold rebuilds the network from scratch at the F05 parameters, so module
# labels are arbitrary and must be re-matched by overlap before anything is projected.
refit_fold <- function(hold_rows, expr) {
  train <- setdiff(seq_len(nrow(expr)), hold_rows)
  net <- blockwiseModules(
    expr[train, , drop = FALSE],
    power = SOFT_POWER, networkType = "signed", TOMType = "signed",
    corType = "bicor", maxPOutliers = 0.05,
    minModuleSize = 30, mergeCutHeight = 0.25,
    numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE, verbose = 0
  )
  mod_train <- labels2colors(net$colors)
  names(mod_train) <- colnames(expr)
  list(mod_train = mod_train, match = match_modules(mod_train, module_colors), train = train)
}

cor <- WGCNA::cor
tier2_raw <- map_dfr(names(bundle$outcomes), function(outcome_name) {
  out <- bundle$outcomes[[outcome_name]]
  ids <- out$ids
  expr <- datExpr[ids, , drop = FALSE]
  group <- if (out$paired) out$subject else ids
  folds <- unique(group)
  scores <- matrix(NA_real_, length(ids), length(modules), dimnames = list(NULL, modules))
  jac <- matrix(NA_real_, length(folds), length(modules), dimnames = list(NULL, modules))

  for (f in seq_along(folds)) {
    hold <- which(group == folds[f])
    fit <- refit_fold(hold, expr)
    for (m in modules) {
      row <- fit$match |> filter(full == m)
      if (nrow(row) == 0) next
      jac[f, m] <- row$jaccard
      if (row$jaccard < JACCARD_FLOOR) next
      prot <- names(fit$mod_train)[fit$mod_train == row$train]
      if (length(prot) < 2) next
      scores[hold, m] <- loso_me(expr[, prot, drop = FALSE], fit$train, hold)
    }
  }

  map_dfr(modules, function(m) {
    ok <- !is.na(scores[, m])
    insample <- fold_auc(uni_auc(out$y, rowMeans(
      datExpr[ids, names(module_colors)[module_colors == m], drop = FALSE]
    ))$auc)
    loso <- if (sum(ok) > 4 && length(unique(out$y[ok])) == 2) {
      fold_auc(uni_auc(out$y[ok], scores[ok, m])$auc)
    } else {
      NA_real_
    }
    tibble(
      tier = "tier2_network_refit", outcome = outcome_name, module = m,
      auc_insample = insample, auc_loso = loso, drop = insample - loso,
      n = sum(ok), mean_jaccard = mean(jac[, m], na.rm = TRUE),
      min_jaccard = suppressWarnings(min(jac[, m], na.rm = TRUE)),
      n_folds = length(folds), n_failed = sum(jac[, m] < JACCARD_FLOOR | is.na(jac[, m]))
    )
  })
})
cor <- stats::cor
write_csv(tier2_raw, file.path(DAT, "circularity_tier2.csv"))

message(sprintf(
  "circularity ladder: tier1 median drop %.3f | tier2 median drop %.3f | tier2 mean Jaccard %.2f",
  stats::median(tier1$drop, na.rm = TRUE),
  stats::median(tier2_raw$drop, na.rm = TRUE),
  mean(tier2_raw$mean_jaccard, na.rm = TRUE)
))
