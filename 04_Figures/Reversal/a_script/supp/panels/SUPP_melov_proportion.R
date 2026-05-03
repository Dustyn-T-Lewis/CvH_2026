# SUPP: Melov-Style Reversal Proportion Test
# Melov et al. 2007, PLoS ONE (PMID 17520024): % of disease-DEPs reversed
# Swindell 2009, BMC Genomics (PMID 19968875): chi-square/Fisher variant
#
# Three complementary tests:
#   1. Binomial test vs 50% (H0: reversal no better than coin flip)
#   2. Fisher exact 2x2 (up/down × reversed/not) for directional independence
#   3. Permutation test (10,000×): shuffle cancer-DEP labels, recompute fraction
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
library(tidyverse)
library(patchwork)

RPT_PNG <- "04_Figures/Reversal/b_reports/supp/png/panels"
RPT_PDF <- "04_Figures/Reversal/b_reports/supp/pdf/panels"
DAT     <- "04_Figures/Reversal/c_data/panel_supp"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT,     recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

# ── Data ─────────────────────────────────────────────────────────────────────
dep_df <- read_csv("03_DEP/c_data/03_combined_results_CRvH.csv",
                   show_col_types = FALSE) %>%
  filter(!is.na(logFC_Cancer_vs_Healthy), !is.na(logFC_Training_CR))

all_df <- dep_df %>%
  mutate(cancer_sig = pi_score_Cancer_vs_Healthy < 0.05,
         reversed   = sign(logFC_Cancer_vs_Healthy) != sign(logFC_Training_CR))

cancer_dep <- all_df %>% filter(cancer_sig)
n_cancer   <- nrow(cancer_dep)
n_reversed <- sum(cancer_dep$reversed)
pct_rev    <- 100 * n_reversed / n_cancer

# ── Test 1: Binomial test vs 50% ────────────────────────────────────────────
binom_res <- binom.test(n_reversed, n_cancer, p = 0.5, alternative = "greater")
message(sprintf("  Melov: %d/%d (%.1f%%) reversed, binom p = %s",
                n_reversed, n_cancer, pct_rev,
                format.pval(binom_res$p.value, digits = 3)))

# ── Test 2: Permutation (10,000×) ───────────────────────────────────────────
# Null: randomly select n_cancer proteins from the full proteome and compute
# reversal fraction. Tests whether cancer-DEPs reverse MORE than random proteins.
set.seed(42)
B <- 10000
n_all <- nrow(all_df)
null_frac <- replicate(B, {
  idx <- sample.int(n_all, n_cancer)
  sum(all_df$reversed[idx]) / n_cancer
})

observed_frac <- n_reversed / n_cancer
perm_p <- (sum(null_frac >= observed_frac) + 1) / (B + 1)
message(sprintf("  Permutation: observed = %.3f, p = %.4f (n_perm = %d)",
                observed_frac, perm_p, B))

# ── Test 3: Fisher exact 2×2 (within cancer-DEPs, up vs down direction) ────
tab <- table(
  cancer_dep$logFC_Cancer_vs_Healthy > 0,
  cancer_dep$reversed
)
dimnames(tab) <- list(cancer_dir = c("Cancer Down", "Cancer Up"),
                      reversed = c("Not reversed", "Reversed"))
fisher_res <- fisher.test(tab)

# ── Export CSV ───────────────────────────────────────────────────────────────
melov_summary <- tibble(
  test = c("binomial_vs_50pct", "permutation_vs_random", "fisher_up_vs_down"),
  n = c(n_cancer, n_cancer, n_cancer),
  n_reversed = c(n_reversed, n_reversed, NA_integer_),
  pct_reversed = c(pct_rev, pct_rev, NA_real_),
  statistic = c(binom_res$statistic, observed_frac, fisher_res$estimate),
  p_value = c(binom_res$p.value, perm_p, fisher_res$p.value),
  ci_lower = c(binom_res$conf.int[1], quantile(null_frac, 0.025), NA_real_),
  ci_upper = c(binom_res$conf.int[2], quantile(null_frac, 0.975), NA_real_),
  note = c("one-sided, greater",
           sprintf("n_perm=%d; null=random protein selection", B),
           "Fisher exact, two-sided")
)
write_csv(melov_summary, file.path(DAT, "SUPP_melov_proportion.csv"))
write_csv(tibble(replicate = seq_len(B), null_frac = null_frac),
          file.path(DAT, "SUPP_melov_null_dist.csv"))

# ── Visualization ────────────────────────────────────────────────────────────
# Left: stacked bar (reversal by cancer direction)
bar_df <- cancer_dep %>%
  mutate(cancer_dir = ifelse(logFC_Cancer_vs_Healthy > 0, "Cancer Up", "Cancer Down"),
         status = ifelse(reversed, "Reversed", "Not reversed")) %>%
  count(cancer_dir, status) %>%
  group_by(cancer_dir) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

p_bar <- ggplot(bar_df, aes(x = cancer_dir, y = pct, fill = status)) +
  geom_col(width = 0.55, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)),
            position = position_stack(vjust = 0.5),
            size = 3.2, fontface = "bold", color = "white") +
  scale_fill_manual(values = c("Reversed" = "#2563EB", "Not reversed" = "#94A3B8"),
                    name = NULL) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey40",
             linewidth = 0.3) +
  annotate("text", x = 1.5, y = 50, label = "50% (chance)",
           vjust = -0.5, size = 2.8, fontface = "italic", color = "grey40") +
  labs(x = NULL, y = "Percentage") +
  scale_y_continuous(limits = c(0, 108), breaks = seq(0, 100, 25)) +
  FIG_THEME +
  theme(legend.position = "bottom",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 7))

# Right: permutation null distribution
p_null <- ggplot(tibble(x = null_frac), aes(x = x * 100)) +
  geom_histogram(bins = 50, fill = "grey70", color = "grey40", linewidth = 0.2) +
  geom_vline(xintercept = pct_rev, color = "#2563EB", linewidth = 0.8) +
  annotate("text", x = pct_rev, y = Inf,
           label = sprintf("Observed = %.1f%%\nperm p %s", pct_rev, fmt_p(perm_p)),
           hjust = -0.1, vjust = 1.5, size = 3, fontface = "bold", color = "#2563EB") +
  labs(x = "% Reversed (null: random proteins)", y = "Count") +
  FIG_THEME

# Combine
pS_melov <- p_bar + p_null +
  plot_layout(widths = c(1, 1.3)) +
  plot_annotation(
    title = "Melov-Style Reversal Proportion",
    subtitle = sprintf("Cancer DEPs (\u03a0 < 0.05, n = %d) | %.1f%% reversed | binom p %s | perm p %s",
                        n_cancer, pct_rev,
                        fmt_p(binom_res$p.value), fmt_p(perm_p)),
    theme = theme(plot.title = element_text(size = FIG_TITLE_SIZE, face = "bold"),
                  plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, color = "grey30")))

ggsave(file.path(RPT_PNG, "SUPP_melov_proportion.png"), pS_melov,
       width = 240, height = 100, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_melov_proportion.pdf"), pS_melov,
       width = 240, height = 100, units = "mm", device = pdf_device)

message("Done: SUPP_melov_proportion")
