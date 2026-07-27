# F06 stage 01: run every prediction engine.
#
# Six subject-level classifiers (three feature spaces x two outcomes) plus the
# protein-level reversal arm. Sample-level power is thin (25 baseline samples, 12 paired
# subjects), so each AUC carries a permutation null and the figure is framed as
# exploratory. The protein-level arm is the well-powered one at n = 297.

setwd(here::here())
pacman::p_load(dplyr, tidyr, tibble, readr, purrr, pROC)
source("04_Figures/shared/prediction_utils.R")

DAT <- "04_Figures/F06_Prediction/c_data"
F05 <- "04_Figures/F05_WGCNA/c_data"
set.seed(42)

bundle <- readRDS(file.path(DAT, "feature_bundle.rds"))
features <- bundle$features
outcomes <- bundle$outcomes

# More candidate features means a wider k to search, but every space starts at 1: with
# five modules a single-eigengene classifier is often the honest optimum.
K_RANGE <- list(protein = c(1, 2, 5, 10, 20), module = 1:5, hallmark = c(1, 2, 5, 10))
N_PERM_HEADLINE <- 1000
N_PERM_OTHER <- 200

# Paired samples hold out by subject, not by row. Leaving one sample in place leaves its
# partner in training, and the model then scores the subject's level rather than the
# training effect.
run_one <- function(space, outcome_name) {
  out <- outcomes[[outcome_name]]
  x <- features[[space]][out$ids, , drop = FALSE]
  group <- if (out$paired) out$subject else NULL
  fit <- run_topk_loocv(out$y, x, k_range = K_RANGE[[space]], group = group)
  roc_obj <- pROC::roc(out$y, fit$scores, quiet = TRUE, direction = "<")
  ci <- as.numeric(pROC::ci.auc(roc_obj))
  list(
    space = space, outcome = outcome_name, label = out$label,
    auc = as.numeric(pROC::auc(roc_obj)), ci_lo = ci[1], ci_hi = ci[3],
    n = length(out$y), n_pos = sum(out$y), n_neg = sum(out$y == 0),
    paired = out$paired, median_k = as.integer(stats::median(fit$best_k)),
    fit = fit, roc = roc_obj, y = out$y, x = x, subject = out$subject
  )
}

grid <- tidyr::expand_grid(space = names(features), outcome = names(outcomes))
runs <- purrr::pmap(grid, function(space, outcome) run_one(space, outcome))
names(runs) <- paste(grid$space, grid$outcome, sep = "_")

# The winner earns the deep permutation budget; the rest get enough to separate a real
# signal from noise without paying for precision nobody reads.
best <- names(runs)[which.max(vapply(runs, function(r) r$auc, numeric(1)))]

runs <- imap(runs, function(r, nm) {
  shuffle <- if (r$paired) within_subject_shuffle(r$subject) else sample
  n_perm <- if (nm == best) N_PERM_HEADLINE else N_PERM_OTHER
  r$perm_p <- perm_p_classifier(
    r$y, r$x,
    k_fixed = r$median_k, obs_auc = r$auc, n_perm = n_perm, shuffle = shuffle,
    group = if (r$paired) r$subject else NULL
  )
  r$n_perm <- n_perm
  message(sprintf(
    "  %-18s AUC %.3f [%.2f-%.2f] perm p=%.4f (k=%d, %d perms)",
    nm, r$auc, r$ci_lo, r$ci_hi, r$perm_p, r$median_k, n_perm
  ))
  r
})

classifier_summary <- map_dfr(runs, function(r) {
  tibble(
    space = r$space, outcome = r$outcome, label = r$label,
    auc = r$auc, ci_lo = r$ci_lo, ci_hi = r$ci_hi,
    perm_p = r$perm_p, n_perm = r$n_perm, median_k = r$median_k,
    n = r$n, n_pos = r$n_pos, n_neg = r$n_neg, paired = r$paired
  )
}) |>
  mutate(q_bh = p.adjust(perm_p, "BH"))
write_csv(classifier_summary, file.path(DAT, "classifier_summary.csv"))

roc_curves <- map_dfr(runs, function(r) {
  tibble(
    space = r$space, outcome = r$outcome,
    fpr = 1 - r$roc$specificities, tpr = r$roc$sensitivities
  )
})
write_csv(roc_curves, file.path(DAT, "roc_curves.csv"))

feature_stability <- map_dfr(runs, function(r) {
  freq <- r$fit$feature_freq
  tibble(
    space = r$space, outcome = r$outcome,
    feature = names(freq), n_selected = as.integer(freq),
    pct_selected = 100 * as.integer(freq) / r$n
  ) |>
    filter(n_selected > 0) |>
    arrange(desc(n_selected))
})
write_csv(feature_stability, file.path(DAT, "feature_stability.csv"))

# Protein-level arm. phi = -logFC_Training / logFC_Baseline defines the outcome, so any
# logFC term would hand the model the class boundary; only module identity and abundance
# enter. This asks whether co-expression membership explains which disease-signature
# proteins training normalizes.
reversal <- read_csv(
  "04_Figures/F04_Reversal/c_data/panel_F/trajectory_clusters.csv",
  show_col_types = FALSE
)
modules <- read_csv(file.path(F05, "wgcna_module_assignments.csv"), show_col_types = FALSE)
imp_mat <- readRDS(file.path(F05, "imp_mat.rds"))

protein_df <- reversal |>
  select(uniprot_id, reversal_class) |>
  inner_join(select(modules, uniprot_id, module_color), by = "uniprot_id") |>
  mutate(
    reversed = as.integer(reversal_class == "Normalized"),
    mean_abund = rowMeans(imp_mat[uniprot_id, , drop = FALSE]),
    sd_abund = apply(imp_mat[uniprot_id, , drop = FALSE], 1, stats::sd),
    module_color = factor(module_color)
  ) |>
  filter(!is.na(module_color))

cross_tab <- protein_df |>
  count(module_color, reversal_class) |>
  pivot_wider(names_from = reversal_class, values_from = n, values_fill = 0)
write_csv(cross_tab, file.path(DAT, "reversal_crosstab.csv"))

chisq <- suppressWarnings(stats::chisq.test(
  table(protein_df$module_color, protein_df$reversed)
))

x_protein_arm <- cbind(
  model.matrix(~module_color, data = protein_df)[, -1, drop = FALSE],
  mean_abund = protein_df$mean_abund,
  sd_abund = protein_df$sd_abund
)
rev_fit <- fit_kfold_perm(x_protein_arm, protein_df$reversed, n_perm = N_PERM_OTHER)

write_csv(
  tibble(
    analysis = "reversal_from_module",
    auc = rev_fit$auc, ci_lo = rev_fit$ci_lo, ci_hi = rev_fit$ci_hi,
    perm_p = rev_fit$perm_p, n = rev_fit$n,
    n_normalized = rev_fit$n_pos, n_other = rev_fit$n_neg,
    chisq_stat = as.numeric(chisq$statistic), chisq_p = chisq$p.value
  ),
  file.path(DAT, "reversal_arm_summary.csv")
)
write_csv(
  tibble(fpr = rev_fit$fpr, tpr = rev_fit$tpr),
  file.path(DAT, "reversal_arm_roc.csv")
)

message(sprintf(
  "  reversal arm: AUC %.3f [%.2f-%.2f] perm p=%.4f | chi-sq p=%.3g | n=%d (%d normalized)",
  rev_fit$auc, rev_fit$ci_lo, rev_fit$ci_hi, rev_fit$perm_p,
  chisq$p.value, rev_fit$n, rev_fit$n_pos
))
