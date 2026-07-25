# F06 stage 02: render the prediction figure and its supplements.
#
# The comparison across feature spaces is the result, not any single AUC. Panels report
# cross-validated performance with permutation p-values; the univariate module tiles are
# in-sample by construction and are labelled as such.

setwd(here::here())
source("04_Figures/F06_Prediction/a_script/style.R")
source("04_Figures/shared/prediction_utils.R")
source("04_Figures/shared/figure_supplement_helpers.R")
pacman::p_load(ggplot2, dplyr, tidyr, tibble, readr, purrr, patchwork, ggtext)

DAT <- "04_Figures/F06_Prediction/c_data"
F05 <- "04_Figures/F05_WGCNA/c_data"
RPT_PNG <- "04_Figures/F06_Prediction/b_reports/main/png"
RPT_PDF <- "04_Figures/F06_Prediction/b_reports/main/pdf"
SUPP_PNG <- "04_Figures/F06_Prediction/b_reports/supp/png"
SUPP_PDF <- "04_Figures/F06_Prediction/b_reports/supp/pdf"
for (d in c(RPT_PNG, RPT_PDF, SUPP_PNG, SUPP_PDF)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
pdf_device <- get_pdf_device()
set.seed(42)

SPACE_LABS <- c(protein = "Proteins", module = "WGCNA modules", hallmark = "Hallmark singscore")
SPACE_COLORS <- c(protein = "#8073AC", module = "#1B7837", hallmark = "#D6604D")
OUTCOME_LABS <- c(baseline = "CR vs Ctl (baseline)", training = "CR pre vs post")

clf <- read_csv(file.path(DAT, "classifier_summary.csv"), show_col_types = FALSE) |>
  mutate(
    space = factor(space, levels = names(SPACE_LABS)),
    outcome = factor(outcome, levels = names(OUTCOME_LABS)),
    sig = if_else(perm_p < 0.05, "*", "")
  )
roc_df <- read_csv(file.path(DAT, "roc_curves.csv"), show_col_types = FALSE) |>
  mutate(
    space = factor(space, levels = names(SPACE_LABS)),
    outcome = factor(outcome, levels = names(OUTCOME_LABS))
  )
tier1 <- read_csv(file.path(DAT, "circularity_tier1.csv"), show_col_types = FALSE)
tier2 <- read_csv(file.path(DAT, "circularity_tier2.csv"), show_col_types = FALSE)
rev_sum <- read_csv(file.path(DAT, "reversal_arm_summary.csv"), show_col_types = FALSE)
rev_roc <- read_csv(file.path(DAT, "reversal_arm_roc.csv"), show_col_types = FALSE)
cross_tab <- read_csv(file.path(DAT, "reversal_crosstab.csv"), show_col_types = FALSE)
stability <- read_csv(file.path(DAT, "feature_stability.csv"), show_col_types = FALSE)

panel_auc <- function() {
  ggplot(clf, aes(auc, space, fill = space)) +
    geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey55", linewidth = 0.4) +
    geom_col(width = 0.62, colour = "grey25", linewidth = 0.3) +
    geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi),
      orientation = "y", width = 0.18, linewidth = 0.4, colour = "grey20"
    ) +
    geom_text(aes(x = ci_hi, label = sprintf("%.2f  p=%.2f%s", auc, perm_p, sig)),
      hjust = -0.08, size = 2.5, colour = "grey15"
    ) +
    facet_wrap(~outcome, ncol = 1, labeller = labeller(outcome = OUTCOME_LABS)) +
    scale_fill_manual(values = SPACE_COLORS, guide = "none") +
    scale_y_discrete(labels = SPACE_LABS) +
    scale_x_continuous(limits = c(0, 1.35), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
    labs(
      title = "Cross-validated classification", x = "AUC (LOOCV)", y = NULL,
      subtitle = "bars = AUC, whiskers = DeLong CI, p = permutation null"
    ) +
    FIG_THEME +
    theme(
      strip.text = element_text(face = "bold", size = 7),
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(size = 7, face = "bold")
    )
}

panel_roc <- function() {
  ggplot(roc_df, aes(fpr, tpr, colour = space)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey65", linewidth = 0.3) +
    geom_step(linewidth = 0.8) +
    facet_wrap(~outcome, ncol = 2, labeller = labeller(outcome = OUTCOME_LABS)) +
    scale_colour_manual(values = SPACE_COLORS, labels = SPACE_LABS, name = NULL) +
    coord_fixed() +
    labs(
      title = "ROC by feature space", x = "False positive rate", y = "True positive rate",
      subtitle = "out-of-fold predictions"
    ) +
    FIG_THEME +
    theme(
      strip.text = element_text(face = "bold", size = 7),
      legend.position = "bottom", legend.text = element_text(size = 6)
    )
}

# Univariate module tiles. These are in-sample scores: the eigengene is used directly as
# the ranking, with no model fitted, so they describe association rather than predict.
module_tiles <- tier1 |>
  select(outcome, module, auc_insample, auc_loso) |>
  left_join(select(tier2, outcome, module, mean_jaccard, n_failed, n_folds),
    by = c("outcome", "module")
  ) |>
  mutate(
    outcome = factor(outcome, levels = names(OUTCOME_LABS)),
    stable = n_failed == 0,
    lab = sprintf("%.2f", auc_insample)
  )

panel_modules <- function() {
  ggplot(module_tiles, aes(outcome, module, fill = auc_insample)) +
    geom_tile(colour = "grey80", linewidth = 0.4, width = 0.94, height = 0.9) +
    geom_tile(
      data = filter(module_tiles, stable), fill = NA, colour = "black",
      linewidth = 0.9, width = 0.94, height = 0.9
    ) +
    geom_text(aes(label = lab), size = 2.4, colour = if_else(module_tiles$auc_insample > 0.8, "white", "grey15")) +
    scale_fill_gradient2(
      low = "white", mid = "#C7E9C0", high = "#1B7837", midpoint = 0.5,
      limits = c(0.4, 1), oob = scales::squish, name = "AUC"
    ) +
    scale_x_discrete(labels = c(baseline = "CR vs Ctl", training = "pre vs post")) +
    labs(
      title = "Single-module association", x = NULL, y = NULL,
      subtitle = "in-sample AUC; box = module survives every leave-one-out refit"
    ) +
    FIG_THEME +
    theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 6.5))
}

panel_reversal <- function() {
  comp <- cross_tab |>
    pivot_longer(-module_color, names_to = "class", values_to = "n") |>
    group_by(module_color) |>
    mutate(frac = n / sum(n)) |>
    ungroup() |>
    mutate(class = factor(class, levels = c("Normalized", "Persistent", "Exacerbated")))
  ggplot(comp, aes(frac, module_color, fill = class)) +
    geom_col(width = 0.76, colour = "grey30", linewidth = 0.25) +
    scale_fill_manual(
      values = c(Normalized = "#1B7837", Persistent = "#BDBDBD", Exacerbated = "#D6604D"),
      name = NULL
    ) +
    scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
    labs(
      title = "Does module membership explain reversal?",
      subtitle = sprintf(
        "%d disease-signature proteins; grey = unassigned\n10-fold AUC %.2f [%.2f-%.2f], perm p=%.3f | chi-sq p=%.2f",
        rev_sum$n, rev_sum$auc, rev_sum$ci_lo, rev_sum$ci_hi, rev_sum$perm_p, rev_sum$chisq_p
      ),
      x = "Share of module's disease-signature proteins", y = NULL
    ) +
    FIG_THEME +
    theme(legend.position = "bottom", legend.text = element_text(size = 6))
}

main_fig <- (panel_auc() | panel_roc()) / (panel_modules() | panel_reversal()) +
  plot_layout(heights = c(1, 0.95)) +
  plot_annotation(
    title = "Can the muscle proteome classify cancer-survivor status and training response?",
    subtitle = paste(
      "Three feature spaces differing in how far their definition sits from these data.",
      "Exploratory: 25 baseline samples, 12 paired subjects; permutation nulls, not confirmatory claims."
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, colour = "grey10"),
      plot.subtitle = element_text(face = "italic", size = 8, colour = "grey40")
    )
  )

ggsave(file.path(RPT_PNG, "MAIN_F06_prediction.png"), main_fig,
  width = 280, height = 210, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(RPT_PDF, "MAIN_F06_prediction.pdf"), main_fig,
  width = 280, height = 210, units = "mm", device = pdf_device
)
message("F06 main figure saved")

supp_circ <- tier2 |>
  filter(!is.na(mean_jaccard)) |>
  mutate(outcome = recode(outcome, !!!OUTCOME_LABS)) |>
  ggplot(aes(mean_jaccard, module, fill = n_failed / n_folds)) +
  geom_col(width = 0.7, colour = "grey30", linewidth = 0.3) +
  geom_vline(xintercept = 0.4, linetype = "dashed", colour = "#D6604D") +
  geom_text(aes(label = sprintf("%d/%d folds failed", n_failed, n_folds)),
    hjust = -0.06, size = 2.2, colour = "grey25"
  ) +
  facet_wrap(~outcome) +
  scale_fill_gradient(low = "#C7E9C0", high = "#D6604D", name = "fold failure rate") +
  scale_x_continuous(limits = c(0, 1.35), expand = c(0, 0)) +
  labs(
    title = "Module reproducibility under leave-one-subject-out network refits",
    subtitle = "Dashed line = Jaccard admissibility floor (0.4). Small modules do not survive subject removal.",
    x = "Mean Jaccard vs full-sample module", y = NULL
  ) +
  FIG_THEME +
  theme(strip.text = element_text(face = "bold", size = 7))

ggsave(file.path(SUPP_PNG, "SUPP_F06_circularity.png"), supp_circ,
  width = 250, height = 110, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(SUPP_PDF, "SUPP_F06_circularity.pdf"), supp_circ,
  width = 250, height = 110, units = "mm", device = pdf_device
)

overview <- tibble::tribble(
  ~Sheet, ~Contents,
  "classifier_summary", "Cross-validated AUC, DeLong CI, permutation p and BH q per feature space and outcome.",
  "roc_curves", "Out-of-fold ROC coordinates for every classifier.",
  "feature_stability", "How often each feature was selected across outer folds.",
  "circularity_tier1", "Eigengene PC re-projection per fold: in-sample vs out-of-fold AUC.",
  "circularity_tier2", "Full WGCNA refit per fold with Jaccard module matching and fold failure counts.",
  "reversal_arm", "Protein-level test of whether module membership explains F04 reversal class.",
  "reversal_crosstab", "Module by reversal-class counts underlying the protein-level arm."
)
build_workbook(
  file.path(DAT, "F06_supplementary.xlsx"),
  "F06 Prediction — supplementary tables",
  paste(
    "Sample classification from three feature spaces (proteins, WGCNA eigengenes, Hallmark singscore).",
    "Exploratory: n=25 baseline, 12 paired subjects. All p-values are permutation-based."
  ),
  overview,
  list(
    list(name = "classifier_summary", df = clf),
    list(name = "roc_curves", df = roc_df),
    list(name = "feature_stability", df = stability),
    list(name = "circularity_tier1", df = tier1),
    list(name = "circularity_tier2", df = tier2),
    list(name = "reversal_arm", df = rev_sum),
    list(name = "reversal_crosstab", df = cross_tab)
  )
)

message("F06 render complete")
