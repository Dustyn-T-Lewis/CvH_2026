# F04 CRvH Supplementary: Directional Asymmetry in Cancer Recovery Concordance
# Tests whether concordance rate and magnitude differ by direction.
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F04/a_script/style.R")

library(tidyverse)
library(patchwork)

RPT <- "04_Figures/F04/CRvH/b_reports/supp"
DAT <- "04_Figures/F04/CRvH/c_data/supp"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

# --- Load data ---------------------------------------------------------------
res <- readr::read_csv("03_DEP/a_non_imputed/c_data/combined_results_pi.csv",
                       show_col_types = FALSE) |>
  dplyr::mutate(contrast = dplyr::recode(contrast,
                                         CRvH_Baseline = "Cancer_vs_Healthy",
                                         CR_Training   = "Training_CR")) |>
  tidyr::pivot_wider(id_cols = c(uniprot_id, gene, protein, description),
                     names_from = contrast,
                     values_from = c(logFC, t, P.Value, adj.P.Val, pi_score, sig_pi),
                     names_glue = "{.value}_{contrast}")

df <- res %>%
  select(gene, logFC_Cancer_vs_Healthy, logFC_Training_CR,
         t_Cancer_vs_Healthy, t_Training_CR) %>%
  filter(!is.na(logFC_Cancer_vs_Healthy), !is.na(logFC_Training_CR)) %>%
  mutate(
    cvh_dir    = ifelse(logFC_Cancer_vs_Healthy > 0, "CvH Up", "CvH Down"),
    concordant = sign(logFC_Cancer_vs_Healthy) == sign(logFC_Training_CR)
  )

# Concordance score: logFC_TR / logFC_CvH (positive = concordant)
df_cs <- df %>%
  filter(abs(logFC_Cancer_vs_Healthy) >= 0.01) %>%
  mutate(CS = pmin(pmax(logFC_Training_CR / logFC_Cancer_vs_Healthy, -5), 5))

# --- Test 1: Proportion test -------------------------------------------------
tab <- df %>%
  group_by(cvh_dir) %>%
  summarise(n_conc = sum(concordant), n_total = n(), .groups = "drop")

prop_res <- prop.test(tab$n_conc, tab$n_total)

# --- Test 2: Wilcoxon on |logFC_TR| by direction ----------------------------
conc_only <- df %>% filter(concordant)
wilcox_res <- wilcox.test(
  abs(logFC_Training_CR) ~ cvh_dir, data = conc_only,
  conf.int = TRUE
)

# --- Test 3: KS test on CS distributions ------------------------------------
cs_up   <- df_cs$CS[df_cs$cvh_dir == "CvH Up"]
cs_down <- df_cs$CS[df_cs$cvh_dir == "CvH Down"]
ks_res  <- ks.test(cs_up, cs_down)

# --- Test 4: Permutation test (10K) -----------------------------------------
set.seed(42)
obs_diff <- abs(diff(tab$n_conc / tab$n_total))
n_perm   <- 10000L
perm_diffs <- replicate(n_perm, {
  shuf <- sample(df$cvh_dir)
  n_c_a <- sum(df$concordant[shuf == "CvH Down"])
  n_t_a <- sum(shuf == "CvH Down")
  n_c_b <- sum(df$concordant[shuf == "CvH Up"])
  n_t_b <- sum(shuf == "CvH Up")
  abs(n_c_a / n_t_a - n_c_b / n_t_b)
})
perm_p <- mean(perm_diffs >= obs_diff)

# --- Collect test results ----------------------------------------------------
tests <- tibble(
  test = c("Proportion test", "Wilcoxon (|logFC_TR|)",
           "KS test (CS)", "Permutation (10K)"),
  statistic = c(prop_res$statistic, wilcox_res$statistic,
                ks_res$statistic, obs_diff),
  p_value = c(prop_res$p.value, wilcox_res$p.value,
              ks_res$p.value, perm_p),
  effect_size = c(
    diff(tab$n_conc / tab$n_total),
    wilcox_res$estimate,
    ks_res$statistic,
    obs_diff
  )
)
write.csv(tests, file.path(DAT, "h_directional_asymmetry_tests.csv"),
          row.names = FALSE)

# --- Test 5: Load existing ORA for side-by-side display ----------------------
ora_file <- "04_Figures/F04/CRvH/c_data/panel_E/rrho2_ora_concordant.csv"
if (file.exists(ora_file)) {
  ora <- read.csv(ora_file) %>%
    mutate(direction = quadrant)
  write.csv(ora, file.path(DAT, "h_directional_asymmetry_ora.csv"),
            row.names = FALSE)
} else {
  ora <- tibble()
  message("ORA file not found -- run panel_E.R first for ORA data")
}

# --- Panel (a): CS density by CvH direction ----------------------------------
p_a <- ggplot(df_cs, aes(x = CS, fill = cvh_dir, color = cvh_dir)) +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#2E7D32",
             linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50",
             linewidth = 0.4) +
  annotate("text", x = 1, y = Inf, label = "Perfect\nconcordance",
           hjust = -0.1, vjust = 1.3, size = 2.8, color = "#2E7D32",
           fontface = "italic") +
  annotate("text", x = 0.5, y = 0, vjust = -0.5,
           label = sprintf("KS %s | Perm %s",
                           fmt_p(ks_res$p.value), fmt_p(perm_p)),
           size = 3, fontface = "bold") +
  scale_fill_manual(values = c("CvH Up" = "#E05A4E", "CvH Down" = "#5DA5DA"),
                    name = "Cancer vs Healthy direction") +
  scale_color_manual(values = c("CvH Up" = "#E05A4E", "CvH Down" = "#5DA5DA"),
                     guide = "none") +
  labs(title = "Concordance Score Distribution",
       subtitle = "CS = logFC(Training CR) / logFC(Cancer vs Healthy), winsorized [-5, 5]",
       x = "Concordance Score", y = "Density") +
  FIG_THEME +
  theme(legend.position = "bottom") +
  labs(tag = "a")

# --- Panel (b): Concordance fractions with CIs -------------------------------
frac_df <- tab %>%
  mutate(
    frac = n_conc / n_total,
    ci_lo = mapply(function(x, n) binom.test(x, n)$conf.int[1], n_conc, n_total),
    ci_hi = mapply(function(x, n) binom.test(x, n)$conf.int[2], n_conc, n_total)
  )

p_b <- ggplot(frac_df, aes(x = cvh_dir, y = frac, fill = cvh_dir)) +
  geom_col(width = 0.6, alpha = 0.8) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.2, linewidth = 0.5) +
  geom_text(aes(label = sprintf("%d/%d\n(%.1f%%)", n_conc, n_total, 100 * frac)),
            vjust = -0.3, size = 3) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  annotate("text", x = 1.5, y = max(frac_df$ci_hi) + 0.06,
           label = sprintf("Prop test: %s", fmt_p(prop_res$p.value)),
           size = 3, fontface = "bold") +
  scale_fill_manual(values = c("CvH Up" = "#E05A4E", "CvH Down" = "#5DA5DA"),
                    guide = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Concordance Rate by Cancer Direction",
       subtitle = "Fraction of proteins with sign(logFC_TR) == sign(logFC_CvH)",
       x = NULL, y = "Concordance fraction") +
  FIG_THEME +
  labs(tag = "b")

# --- Panel (c): Direction-stratified ORA bars --------------------------------
if (nrow(ora) > 0 && "direction" %in% names(ora)) {
  ora_top <- ora %>%
    filter(nchar(direction) > 0) %>%
    group_by(direction) %>%
    slice_min(padj, n = 8, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      pathway_label = clean_pathway_name(pathway),
      neg_log_padj  = -log10(padj)
    )

  ora_colors <- c(
    "Concordant Up"   = "#E57373",
    "Concordant Down" = "#64B5F6"
  )

  p_c <- ggplot(ora_top, aes(x = neg_log_padj,
                             y = reorder_within(pathway_label, neg_log_padj, direction),
                             fill = direction)) +
    geom_col(alpha = 0.85, width = 0.7) +
    facet_wrap(~direction, scales = "free_y", ncol = 1) +
    scale_y_reordered() +
    scale_fill_manual(values = ora_colors, guide = "none") +
    labs(title = "Direction-Stratified Pathway Enrichment",
         subtitle = "Top 8 ORA pathways per concordance quadrant (padj < 0.05)",
         x = expression(-log[10](p[adj])), y = NULL) +
    FIG_THEME +
    theme(axis.text.y = element_text(size = 7),
          strip.text  = element_text(size = 9)) +
    labs(tag = "c")
} else {
  p_c <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "ORA data not available\n(run panel_E.R first)",
             size = 4, color = "grey50", fontface = "italic") +
    theme_void() +
    labs(tag = "c")
}

# --- Composite ----------------------------------------------------------------
composite <- (p_a | p_b) / p_c +
  plot_layout(heights = c(1, 1.4))

ggsave(file.path(RPT, "supp_directional_tests_composite_SUPP.pdf"), composite,
       width = 220, height = 280, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "supp_directional_tests_composite_SUPP.png"), composite,
       width = 220, height = 280, units = "mm", dpi = 300)

cat("F04 CRvH directional asymmetry analysis complete.\n")
cat(sprintf("  Proportion test: %s\n", fmt_p(prop_res$p.value)))
cat(sprintf("  KS test: %s\n", fmt_p(ks_res$p.value)))
cat(sprintf("  Permutation: %s\n", fmt_p(perm_p)))
