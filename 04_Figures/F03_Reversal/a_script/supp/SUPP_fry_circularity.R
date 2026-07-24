# SUPP: Shared-Baseline Circularity Diagnostics
# Smyth & Altman 2013, BMC Bioinform (PMID 23705896): shared-reference designs
# Wu & Smyth 2012, NAR (PMID 22638577): inter-gene correlation in set tests
# Kim et al. 2024, Mol Sys Biol (PMID 39349762): pooled-screen shared-reference
#
# Problem: CRvH_Baseline and CR_Training share CR_T1 as baseline, creating
# a structural negative correlation between fold-change vectors. The gene-label
# shuffle (old test) produces a null centered at r=0, which is the wrong null
# because the true structural expectation is r<0 even without biological reversal.
#
# This panel runs THREE complementary diagnostics:
#   1. Background correlation: r(all proteins) as the structural baseline
#   2. Subset specificity test (10,000×): randomly draw n_dep proteins from
#      the full proteome and compute r — tests whether cancer-DEPs reverse
#      MORE than random subsets of the same size
#   3. Reversal fraction specificity (10,000×): same null, but for %reversed
#      (complements the correlation metric with a count metric)
setwd(here::here())
source("04_Figures/shared/style.R")
pacman::p_load(tidyverse, patchwork)

RPT_PNG <- "04_Figures/F03_Reversal/b_reports/supp/png/panels"
RPT_PDF <- "04_Figures/F03_Reversal/b_reports/supp/pdf/panels"
DAT     <- "04_Figures/F03_Reversal/c_data/panel_supp"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT,     recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

# ── Data ─────────────────────────────────────────────────────────────────────
source("04_Figures/F03_Reversal/a_script/f03_data.R")
dep_df <- dep_df %>%
  transmute(gene,
            logFC_CvH = logFC_CRvH_Baseline,
            logFC_TR  = logFC_CR_Training,
            pi_CvH    = pi_score_CRvH_Baseline) %>%
  filter(!is.na(logFC_CvH), !is.na(logFC_TR))

n_all <- nrow(dep_df)
cancer_dep <- dep_df %>% filter(pi_CvH < 0.05)
n_dep <- nrow(cancer_dep)

# ── Test 1: Background correlation (structural baseline) ────────────────────
r_all <- cor(dep_df$logFC_CvH, dep_df$logFC_TR, method = "pearson")
r_dep <- cor(cancer_dep$logFC_CvH, cancer_dep$logFC_TR, method = "pearson")

dep_df <- dep_df %>%
  mutate(reversed = sign(logFC_CvH) != sign(logFC_TR))
cancer_dep <- cancer_dep %>%
  mutate(reversed = sign(logFC_CvH) != sign(logFC_TR))

frac_all <- mean(dep_df$reversed)
frac_dep <- mean(cancer_dep$reversed)

message(sprintf("  Background r (all %d proteins): %.3f", n_all, r_all))
message(sprintf("  Cancer-DEP r (%d proteins):     %.3f", n_dep, r_dep))
message(sprintf("  Background %% reversed: %.1f%% | Cancer-DEP: %.1f%%",
                100 * frac_all, 100 * frac_dep))

# ── Test 2: Subset specificity permutation (correlation) ────────────────────
# Null: randomly draw n_dep proteins from the full proteome, compute r.
# This preserves the shared-baseline structure (all fold changes are real)
# and tests whether cancer-DEPs are MORE anti-correlated than random subsets.
set.seed(42)
B <- 10000
null_r <- replicate(B, {
  idx <- sample.int(n_all, n_dep)
  cor(dep_df$logFC_CvH[idx], dep_df$logFC_TR[idx], method = "pearson")
})

# One-sided: is cancer-DEP r more negative than random subsets?
p_r <- (sum(null_r <= r_dep) + 1) / (B + 1)
message(sprintf("  Subset r test: observed = %.3f, perm p = %.4f", r_dep, p_r))

# ── Test 3: Subset specificity permutation (reversal fraction) ──────────────
null_frac <- replicate(B, {
  idx <- sample.int(n_all, n_dep)
  mean(dep_df$reversed[idx])
})

p_frac <- (sum(null_frac >= frac_dep) + 1) / (B + 1)
message(sprintf("  Subset frac test: observed = %.1f%%, perm p = %.4f",
                100 * frac_dep, p_frac))

# ── Export ───────────────────────────────────────────────────────────────────
circ_summary <- tibble(
  metric = c("pearson_r", "pct_reversed"),
  background_all = c(r_all, 100 * frac_all),
  cancer_dep = c(r_dep, 100 * frac_dep),
  n_all = n_all, n_dep = n_dep,
  perm_p = c(p_r, p_frac), n_perm = B,
  note = c("one-sided: r_dep < null_r",
           "one-sided: frac_dep > null_frac")
)
write_csv(circ_summary, file.path(DAT, "SUPP_fry_circularity.csv"))
write_csv(tibble(replicate = seq_len(B), null_r = null_r, null_frac = null_frac),
          file.path(DAT, "SUPP_fry_circularity_null.csv"))

# ── Visualization ────────────────────────────────────────────────────────────
# Panel 1: Subset correlation specificity
p_r_null <- ggplot(tibble(x = null_r), aes(x = x)) +
  geom_histogram(bins = 50, fill = "grey70", colour = "grey40", linewidth = 0.2) +
  # Background r (all proteins) — the structural baseline
  geom_vline(xintercept = r_all, colour = "grey30",
             linewidth = 0.7, linetype = "dashed") +
  # Observed r (cancer-DEPs)
  geom_vline(xintercept = r_dep, colour = DIR_COLORS["Up"],
             linewidth = 0.8, linetype = "solid") +
  annotate("text", x = r_dep, y = Inf, vjust = 1.5, hjust = -0.1,
           label = sprintf("Cancer-DEP r = %.3f", r_dep),
           colour = DIR_COLORS["Up"], fontface = "bold", size = 3) +
  annotate("text", x = r_all, y = Inf, vjust = 3.0, hjust = -0.1,
           label = sprintf("Background r = %.3f", r_all),
           colour = "grey30", fontface = "bold", size = 3) +
  annotate("text", x = max(null_r), y = Inf, vjust = 1.5, hjust = 1,
           label = sprintf("perm p %s", fmt_p(p_r)),
           colour = "grey30", fontface = "bold.italic", size = 3.5) +
  labs(title    = "Subset Correlation Specificity",
       subtitle = sprintf("Null: r of %d random proteins (n = %d)", n_dep, n_all),
       x = "Pearson r (random subset)", y = "Count") +
  FIG_THEME

# Panel 2: Subset reversal fraction specificity
p_frac_null <- ggplot(tibble(x = 100 * null_frac), aes(x = x)) +
  geom_histogram(bins = 50, fill = "grey70", colour = "grey40", linewidth = 0.2) +
  geom_vline(xintercept = 100 * frac_all, colour = "grey30",
             linewidth = 0.7, linetype = "dashed") +
  geom_vline(xintercept = 100 * frac_dep, colour = "#2563EB",
             linewidth = 0.8, linetype = "solid") +
  annotate("text", x = 100 * frac_dep, y = Inf, vjust = 1.5, hjust = 1.1,
           label = sprintf("Cancer-DEP = %.1f%%", 100 * frac_dep),
           colour = "#2563EB", fontface = "bold", size = 3) +
  annotate("text", x = 100 * frac_all, y = Inf, vjust = 3.0, hjust = -0.1,
           label = sprintf("Background = %.1f%%", 100 * frac_all),
           colour = "grey30", fontface = "bold", size = 3) +
  annotate("text", x = min(100 * null_frac) + 1, y = Inf, vjust = 1.5,
           hjust = 0,
           label = sprintf("perm p %s", fmt_p(p_frac)),
           colour = "grey30", fontface = "bold.italic", size = 3.5) +
  labs(title    = "Subset Reversal Fraction Specificity",
       subtitle = sprintf("Null: %% reversed in %d random proteins", n_dep),
       x = "% Reversed (random subset)", y = "Count") +
  FIG_THEME

# Combine
pS_circ <- p_r_null + p_frac_null +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "Shared-Baseline Circularity Diagnostics",
    subtitle = sprintf(
      "Cancer-DEPs (n = %d) vs random subsets of same size from %d proteins | %s permutations",
      n_dep, n_all, format(B, big.mark = ",")),
    theme = theme(plot.title = element_text(size = FIG_TITLE_SIZE, face = "bold"),
                  plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, color = "grey30")))

ggsave(file.path(RPT_PNG, "SUPP_fry_circularity.png"), pS_circ,
       width = 280, height = 110, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_fry_circularity.pdf"), pS_circ,
       width = 280, height = 110, units = "mm", device = pdf_device)

message(sprintf("Done: SUPP_fry_circularity (r perm p %s, frac perm p %s)",
                fmt_p(p_r), fmt_p(p_frac)))
