# Pilot: Three alternative reversal quantification approaches
# 1. Melov-style reversal proportion (Melov et al. 2007 PLoS ONE)
# 2. CMap connectivity score (Subramanian et al. 2017 Cell)
# 3. Recovery score per protein (continuous metric)
setwd(here::here())
source("04_Figures/shared/style.R")
library(tidyverse)
library(patchwork)

dep_df <- read_csv("03_DEP/c_data/03_combined_results_CRvH.csv",
                    show_col_types = FALSE)
pdf_device <- get_pdf_device()
PILOT_DIR <- "04_Figures/Reversal/b_reports/pilot"
dir.create(PILOT_DIR, recursive = TRUE, showWarnings = FALSE)

# ═══════════════════════════════════════════════════════════════════════════════
# Define cancer signature (Pi < 0.05)
# ═══════════════════════════════════════════════════════════════════════════════
sig_df <- dep_df %>%
  filter(!is.na(logFC_Cancer_vs_Healthy), !is.na(logFC_Training_CR)) %>%
  mutate(
    cancer_sig   = pi_score_Cancer_vs_Healthy < 0.05,
    cancer_dir   = ifelse(logFC_Cancer_vs_Healthy > 0, "Cancer Up", "Cancer Down"),
    training_dir = ifelse(logFC_Training_CR > 0, "Training Up", "Training Down"),
    reversed     = sign(logFC_Cancer_vs_Healthy) != sign(logFC_Training_CR)
  )

cancer_dep <- sig_df %>% filter(cancer_sig)
n_cancer <- nrow(cancer_dep)
n_up     <- sum(cancer_dep$cancer_dir == "Cancer Up")
n_dn     <- sum(cancer_dep$cancer_dir == "Cancer Down")

cat(sprintf("Cancer DEPs (Pi<0.05): %d (up=%d, down=%d)\n", n_cancer, n_up, n_dn))

# ═══════════════════════════════════════════════════════════════════════════════
# 1. MELOV APPROACH: Simple reversal proportion
#    "What fraction of cancer-changed proteins show opposite-sign training change?"
#    Melov et al. 2007: reported % of age-affected transcripts reversed by exercise
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n═══ 1. MELOV REVERSAL PROPORTION ═══\n")

melov_df <- cancer_dep %>%
  group_by(cancer_dir) %>%
  summarise(
    n_total    = n(),
    n_reversed = sum(reversed),
    n_same     = sum(!reversed),
    pct_rev    = 100 * n_reversed / n_total,
    .groups    = "drop"
  )
print(melov_df)

# Binomial test: is reversal proportion > 50%?
binom_up <- binom.test(
  sum(cancer_dep$reversed[cancer_dep$cancer_dir == "Cancer Up"]),
  sum(cancer_dep$cancer_dir == "Cancer Up"),
  p = 0.5, alternative = "greater")
binom_dn <- binom.test(
  sum(cancer_dep$reversed[cancer_dep$cancer_dir == "Cancer Down"]),
  sum(cancer_dep$cancer_dir == "Cancer Down"),
  p = 0.5, alternative = "greater")
binom_all <- binom.test(sum(cancer_dep$reversed), n_cancer,
                         p = 0.5, alternative = "greater")

cat(sprintf("  All:         %.1f%% reversed (binom p=%s)\n",
            100 * sum(cancer_dep$reversed) / n_cancer, format.pval(binom_all$p.value, digits=3)))
cat(sprintf("  Cancer Up:   %.1f%% reversed (binom p=%s)\n",
            melov_df$pct_rev[melov_df$cancer_dir == "Cancer Up"],
            format.pval(binom_up$p.value, digits=3)))
cat(sprintf("  Cancer Down: %.1f%% reversed (binom p=%s)\n",
            melov_df$pct_rev[melov_df$cancer_dir == "Cancer Down"],
            format.pval(binom_dn$p.value, digits=3)))

# Visualization: stacked bar with reversal proportions
melov_long <- cancer_dep %>%
  mutate(status = ifelse(reversed, "Reversed", "Not reversed")) %>%
  count(cancer_dir, status) %>%
  group_by(cancer_dir) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup()

p_melov <- ggplot(melov_long, aes(x = cancer_dir, y = pct, fill = status)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)),
            position = position_stack(vjust = 0.5),
            size = 3.5, fontface = "bold", color = "white") +
  scale_fill_manual(values = c("Reversed" = "#2563EB", "Not reversed" = "#DC2626"),
                    name = NULL) +
  # Add binomial p-values
  annotate("text", x = 1, y = 102,
           label = sprintf("p %s", format.pval(binom_up$p.value, digits=2)),
           size = 3, fontface = "bold.italic", color = "grey30") +
  annotate("text", x = 2, y = 102,
           label = sprintf("p %s", format.pval(binom_dn$p.value, digits=2)),
           size = 3, fontface = "bold.italic", color = "grey30") +
  labs(title = "Melov-Style Reversal Proportion",
       subtitle = sprintf("Cancer DEPs (Pi < 0.05, n = %d) | Binomial test vs 50%%",
                           n_cancer),
       x = NULL, y = "Percentage") +
  scale_y_continuous(limits = c(0, 110), breaks = seq(0, 100, 25)) +
  FIG_THEME +
  theme(legend.position = "bottom")

# ═══════════════════════════════════════════════════════════════════════════════
# 2. CONNECTIVITY MAP APPROACH: Weighted KS enrichment score
#    Subramanian et al. 2017 Cell: query signature vs reference profile
#    Query = cancer up/down sets; Reference = training-ranked proteins
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n═══ 2. CMAP CONNECTIVITY SCORE ═══\n")

# Rank all proteins by training t-statistic (descending)
all_ranked <- sig_df %>%
  arrange(desc(t_Training_CR)) %>%
  mutate(rank = row_number())
n_all <- nrow(all_ranked)

# Weighted KS enrichment score (GSEA-style)
compute_es <- function(ranks_in_set, n_total, weights = NULL) {
  n_set <- length(ranks_in_set)
  if (n_set == 0) return(list(es = 0, running = rep(0, n_total)))

  # Sort set positions
  hit_pos <- sort(ranks_in_set)

  # Running enrichment score
  if (is.null(weights)) weights <- rep(1, n_set)
  weights <- abs(weights)

  running <- numeric(n_total)
  hit_indicator <- rep(FALSE, n_total)
  hit_indicator[hit_pos] <- TRUE

  # Weighted hit score
  hit_weights <- numeric(n_total)
  hit_weights[hit_pos] <- weights[order(ranks_in_set)]
  hit_cum <- cumsum(hit_weights * hit_indicator) / sum(weights)

  # Miss score
  miss_cum <- cumsum(!hit_indicator) / (n_total - n_set)

  running <- hit_cum - miss_cum

  # ES = max deviation from 0
  es <- if (abs(max(running)) > abs(min(running))) max(running) else min(running)
  list(es = es, running = running)
}

# Cancer-up set: ranks in the training-ranked list
up_ranks <- all_ranked$rank[all_ranked$cancer_sig & all_ranked$cancer_dir == "Cancer Up"]
dn_ranks <- all_ranked$rank[all_ranked$cancer_sig & all_ranked$cancer_dir == "Cancer Down"]

# Weight by |t_Cancer_vs_Healthy| at each protein's rank position
up_weights <- abs(all_ranked$t_Cancer_vs_Healthy[up_ranks])
dn_weights <- abs(all_ranked$t_Cancer_vs_Healthy[dn_ranks])

es_up <- compute_es(up_ranks, n_all, up_weights)
es_dn <- compute_es(dn_ranks, n_all, dn_weights)

# Connectivity score (Subramanian 2017):
# If ES_up and ES_down have same sign → incoherent → connectivity = 0
# If opposite signs → connectivity = (ES_up - ES_down) / 2
if (sign(es_up$es) == sign(es_dn$es)) {
  connectivity <- 0
  coherent <- FALSE
} else {
  connectivity <- (es_up$es - es_dn$es) / 2
  coherent <- TRUE
}

cat(sprintf("  ES (cancer-up):   %.3f  (expect negative = reversal)\n", es_up$es))
cat(sprintf("  ES (cancer-down): %.3f  (expect positive = reversal)\n", es_dn$es))
cat(sprintf("  Connectivity:     %.3f  (negative = reversal, coherent = %s)\n",
            connectivity, coherent))

# Permutation test for connectivity score
set.seed(42)
n_perm <- 1000
null_conn <- numeric(n_perm)
for (i in seq_len(n_perm)) {
  perm_ranks_up <- sample(n_all, length(up_ranks))
  perm_ranks_dn <- sample(n_all, length(dn_ranks))
  perm_es_up <- compute_es(perm_ranks_up, n_all)$es
  perm_es_dn <- compute_es(perm_ranks_dn, n_all)$es
  if (sign(perm_es_up) == sign(perm_es_dn)) {
    null_conn[i] <- 0
  } else {
    null_conn[i] <- (perm_es_up - perm_es_dn) / 2
  }
}
perm_p <- mean(null_conn <= connectivity)
cat(sprintf("  Permutation p-value: %.4f (n_perm = %d)\n", perm_p, n_perm))

# Visualization: running enrichment curves
run_df <- tibble(
  rank = rep(1:n_all, 2),
  ES   = c(es_up$running, es_dn$running),
  Set  = rep(c("Cancer Up", "Cancer Down"), each = n_all)
)

marks_up <- tibble(rank = up_ranks, Set = "Cancer Up")
marks_dn <- tibble(rank = dn_ranks, Set = "Cancer Down")

p_cmap_curves <- ggplot(run_df, aes(x = rank, y = ES, color = Set)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_rug(data = marks_up, aes(x = rank), sides = "b",
           color = "#E57373", alpha = 0.3, length = unit(2, "mm"), inherit.aes = FALSE) +
  geom_rug(data = marks_dn, aes(x = rank), sides = "t",
           color = "#64B5F6", alpha = 0.3, length = unit(2, "mm"), inherit.aes = FALSE) +
  scale_color_manual(values = c("Cancer Up" = "#E57373", "Cancer Down" = "#64B5F6")) +
  annotate("text", x = n_all * 0.98, y = max(run_df$ES) * 0.9,
           label = sprintf("ES(up) = %.3f\nES(dn) = %.3f\nConnectivity = %.3f\nperm p = %.3f",
                            es_up$es, es_dn$es, connectivity, perm_p),
           hjust = 1, size = 3, fontface = "bold", color = "grey25") +
  labs(title = "CMap Connectivity Score",
       subtitle = sprintf("Cancer signature (%d up, %d dn) on Training CR rank | %s",
                           length(up_ranks), length(dn_ranks),
                           ifelse(connectivity < 0, "REVERSAL", "EXACERBATION")),
       x = sprintf("Rank by Training CR t-stat (n = %d)", n_all),
       y = "Enrichment Score", color = NULL) +
  FIG_THEME +
  theme(legend.position = c(0.15, 0.85))

# Null distribution
p_cmap_null <- ggplot(tibble(null = null_conn), aes(x = null)) +
  geom_histogram(bins = 50, fill = "grey70", color = "grey50", linewidth = 0.2) +
  geom_vline(xintercept = connectivity, color = "#2563EB", linewidth = 1) +
  annotate("text", x = connectivity, y = Inf,
           label = sprintf("Observed = %.3f\np = %.3f", connectivity, perm_p),
           hjust = -0.1, vjust = 1.5, size = 3, fontface = "bold", color = "#2563EB") +
  labs(title = "Null Distribution",
       subtitle = "1000 permutations",
       x = "Connectivity Score", y = "Count") +
  FIG_THEME

# ═══════════════════════════════════════════════════════════════════════════════
# 3. RECOVERY SCORE: Continuous per-protein reversal metric
#    recovery = -logFC_TR / logFC_CvH
#    > 0 = reversal direction, = 1 = exact reversal, > 1 = overcorrection
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n═══ 3. RECOVERY SCORE ═══\n")

cancer_dep <- cancer_dep %>%
  mutate(
    recovery = -logFC_Training_CR / logFC_Cancer_vs_Healthy,
    # Cap at [-3, 3] for visualization
    recovery_cap = pmax(-3, pmin(3, recovery))
  )

# Summary stats
med_all   <- median(cancer_dep$recovery)
med_up    <- median(cancer_dep$recovery[cancer_dep$cancer_dir == "Cancer Up"])
med_dn    <- median(cancer_dep$recovery[cancer_dep$cancer_dir == "Cancer Down"])

# Wilcoxon test: is median recovery > 0?
wilcox_all <- wilcox.test(cancer_dep$recovery, mu = 0, alternative = "greater")
wilcox_up  <- wilcox.test(cancer_dep$recovery[cancer_dep$cancer_dir == "Cancer Up"],
                           mu = 0, alternative = "greater")
wilcox_dn  <- wilcox.test(cancer_dep$recovery[cancer_dep$cancer_dir == "Cancer Down"],
                           mu = 0, alternative = "greater")

cat(sprintf("  All DEPs:    median recovery = %.3f (Wilcoxon p = %s)\n",
            med_all, format.pval(wilcox_all$p.value, digits=3)))
cat(sprintf("  Cancer Up:   median recovery = %.3f (Wilcoxon p = %s)\n",
            med_up, format.pval(wilcox_up$p.value, digits=3)))
cat(sprintf("  Cancer Down: median recovery = %.3f (Wilcoxon p = %s)\n",
            med_dn, format.pval(wilcox_dn$p.value, digits=3)))

# Categorize recovery
cancer_dep <- cancer_dep %>%
  mutate(
    recovery_cat = case_when(
      recovery >= 0.75 ~ "Strong reversal (>75%)",
      recovery >= 0.25 ~ "Partial reversal (25-75%)",
      recovery >= 0    ~ "Weak reversal (0-25%)",
      recovery >= -0.25 ~ "Mild exacerbation",
      TRUE              ~ "Strong exacerbation"
    ) %>% factor(levels = c("Strong reversal (>75%)", "Partial reversal (25-75%)",
                             "Weak reversal (0-25%)", "Mild exacerbation",
                             "Strong exacerbation"))
  )

cat_counts <- cancer_dep %>% count(cancer_dir, recovery_cat) %>%
  group_by(cancer_dir) %>% mutate(pct = 100 * n / sum(n)) %>% ungroup()
cat("\nRecovery categories:\n")
print(cat_counts %>% dplyr::select(cancer_dir, recovery_cat, n, pct))

# Visualization: density + violin
p_recovery_density <- ggplot(cancer_dep, aes(x = recovery_cap, fill = cancer_dir)) +
  geom_density(alpha = 0.5, color = "grey30", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "#2563EB", linewidth = 0.6) +
  geom_vline(xintercept = med_all, color = "black", linewidth = 0.8) +
  annotate("text", x = 1.05, y = Inf, label = "Perfect\nreversal",
           hjust = 0, vjust = 1.5, size = 2.5, color = "#2563EB", fontface = "italic") +
  annotate("text", x = med_all, y = Inf,
           label = sprintf("Median = %.2f\np %s", med_all,
                            format.pval(wilcox_all$p.value, digits=2)),
           hjust = -0.1, vjust = 1.5, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("Cancer Up" = "#E57373", "Cancer Down" = "#64B5F6"),
                    name = NULL) +
  labs(title = "Recovery Score Distribution",
       subtitle = sprintf("recovery = -logFC(Training) / logFC(Cancer) | n = %d DEPs", n_cancer),
       x = "Recovery Score (0 = no change, 1 = full reversal)",
       y = "Density") +
  scale_x_continuous(breaks = seq(-3, 3, 1),
                     labels = c("-3", "-2", "-1", "0", "1", "2", "3")) +
  coord_cartesian(xlim = c(-3, 3)) +
  FIG_THEME +
  theme(legend.position = c(0.85, 0.85))

# Violin split by direction
p_recovery_violin <- ggplot(cancer_dep, aes(x = cancer_dir, y = recovery_cap,
                                             fill = cancer_dir)) +
  geom_violin(alpha = 0.4, color = "grey50", draw_quantiles = c(0.25, 0.5, 0.75)) +
  geom_jitter(aes(color = cancer_dir), width = 0.15, size = 0.5, alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = 1, linetype = "dotted", color = "#2563EB") +
  scale_fill_manual(values = c("Cancer Up" = "#E57373", "Cancer Down" = "#64B5F6")) +
  scale_color_manual(values = c("Cancer Up" = "#E57373", "Cancer Down" = "#64B5F6")) +
  # Annotate medians
  annotate("text", x = 0.6, y = med_up,
           label = sprintf("%.2f", med_up), size = 3, fontface = "bold") +
  annotate("text", x = 2.4, y = med_dn,
           label = sprintf("%.2f", med_dn), size = 3, fontface = "bold") +
  # p-values
  annotate("text", x = 1, y = 3.2,
           label = sprintf("p %s", format.pval(wilcox_up$p.value, digits=2)),
           size = 2.8, fontface = "bold.italic", color = "grey30") +
  annotate("text", x = 2, y = 3.2,
           label = sprintf("p %s", format.pval(wilcox_dn$p.value, digits=2)),
           size = 2.8, fontface = "bold.italic", color = "grey30") +
  labs(title = "Recovery by Cancer Direction",
       subtitle = "Wilcoxon test: median > 0?",
       x = NULL, y = "Recovery Score") +
  coord_cartesian(ylim = c(-3, 3.5)) +
  FIG_THEME +
  theme(legend.position = "none")

# Recovery category bar chart
cat_colors <- c("Strong reversal (>75%)" = "#1565C0",
                "Partial reversal (25-75%)" = "#64B5F6",
                "Weak reversal (0-25%)" = "#BBDEFB",
                "Mild exacerbation" = "#FFCDD2",
                "Strong exacerbation" = "#E53935")

p_recovery_cats <- ggplot(cat_counts, aes(x = cancer_dir, y = pct, fill = recovery_cat)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.2) +
  geom_text(aes(label = ifelse(pct > 5, sprintf("%.0f%%", pct), "")),
            position = position_stack(vjust = 0.5),
            size = 2.8, fontface = "bold", color = "white") +
  scale_fill_manual(values = cat_colors, name = "Recovery category") +
  labs(title = "Recovery Categories",
       subtitle = "% of cancer DEPs by reversal extent",
       x = NULL, y = "Percentage") +
  FIG_THEME +
  theme(legend.position = "right",
        legend.text = element_text(size = 7))

# ═══════════════════════════════════════════════════════════════════════════════
# COMPOSITE: All three approaches
# ═══════════════════════════════════════════════════════════════════════════════
pilot_composite <- (p_melov | p_cmap_curves | p_cmap_null) /
                   (p_recovery_density | p_recovery_violin | p_recovery_cats) +
  plot_annotation(
    title = "Reversal Quantification: Three Approaches",
    subtitle = sprintf("Cancer DEPs (Pi < 0.05, n = %d) | Training CR response",
                        n_cancer),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "grey30"),
      plot.tag = element_text(face = "bold", size = 13)))

ggsave(file.path(PILOT_DIR, "pilot_reversal_methods.png"), pilot_composite,
       width = 420, height = 280, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(PILOT_DIR, "pilot_reversal_methods.pdf"), pilot_composite,
       width = 420, height = 280, units = "mm", device = pdf_device)

cat(sprintf("\nPilot saved to: %s\n", PILOT_DIR))
