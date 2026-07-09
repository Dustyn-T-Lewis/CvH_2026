# SUPP: Directional Asymmetry in Signature Reversal
# NOVEL FINDING — no prior publication quantifies this at the proteome level.
#
# Tests whether cancer-DOWN proteins reverse at a different rate than
# cancer-UP proteins after exercise training. This asymmetry has biological
# implications: proteins suppressed by cachexia/disuse may recover more
# readily via exercise-driven anabolic signaling than proteins elevated
# by persistent oncogenic/inflammatory programs.
#
# Mechanistic basis:
#   Murgia et al. 2023, JCSM (PMID 36517414): protein-level recovery asymmetry
#   Tyagi & Pedrioli 2015, NAR (PMID 25870413): codon-bias differential recovery
#   Mahmassani et al. 2021, J Gerontol A (PMID 33705535): translational asymmetry
#
# Statistical tests:
#   1. Two-proportion z-test: reversal rate (DOWN) vs reversal rate (UP)
#   2. Chi-square test: 2×2 (cancer direction × reversal status)
#   3. Bootstrap CI for the asymmetry difference
#   4. Threshold sensitivity: does asymmetry persist across Pi / FDR thresholds?
setwd(here::here())
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
  mutate(reversed = sign(logFC_Cancer_vs_Healthy) != sign(logFC_Training_CR),
         cancer_dir = ifelse(logFC_Cancer_vs_Healthy > 0, "Cancer Up", "Cancer Down"))

# ── Primary analysis at Pi < 0.05 ───────────────────────────────────────────
cancer_dep <- all_df %>% filter(pi_score_Cancer_vs_Healthy < 0.05)

n_up   <- sum(cancer_dep$cancer_dir == "Cancer Up")
n_dn   <- sum(cancer_dep$cancer_dir == "Cancer Down")
rev_up <- sum(cancer_dep$reversed[cancer_dep$cancer_dir == "Cancer Up"])
rev_dn <- sum(cancer_dep$reversed[cancer_dep$cancer_dir == "Cancer Down"])
pct_up <- 100 * rev_up / n_up
pct_dn <- 100 * rev_dn / n_dn
delta  <- pct_dn - pct_up

message(sprintf("  Cancer Up:   %d/%d (%.1f%%) reversed", rev_up, n_up, pct_up))
message(sprintf("  Cancer Down: %d/%d (%.1f%%) reversed", rev_dn, n_dn, pct_dn))
message(sprintf("  Asymmetry:   %.1f percentage points (Down > Up)", delta))

# Test 1: Two-proportion z-test
prop_test <- prop.test(c(rev_dn, rev_up), c(n_dn, n_up), alternative = "greater")
message(sprintf("  Two-proportion z-test: p = %s", format.pval(prop_test$p.value, digits = 3)))

# Test 2: Chi-square 2×2
tab <- table(
  direction = cancer_dep$cancer_dir,
  reversed = cancer_dep$reversed
)
chi_test <- chisq.test(tab)

# Test 3: Bootstrap CI for difference in reversal rates
set.seed(42)
B <- 10000
boot_delta <- replicate(B, {
  idx <- sample.int(nrow(cancer_dep), replace = TRUE)
  bd <- cancer_dep[idx, ]
  up_b <- bd %>% filter(cancer_dir == "Cancer Up")
  dn_b <- bd %>% filter(cancer_dir == "Cancer Down")
  if (nrow(up_b) == 0 || nrow(dn_b) == 0) return(NA_real_)
  100 * mean(dn_b$reversed) - 100 * mean(up_b$reversed)
})
boot_delta <- boot_delta[!is.na(boot_delta)]
boot_ci <- quantile(boot_delta, probs = c(0.025, 0.975))

# ── Threshold sensitivity ────────────────────────────────────────────────────
thresholds <- list(
  "Pi < 0.01"    = all_df %>% filter(pi_score_Cancer_vs_Healthy < 0.01),
  "Pi < 0.05"    = all_df %>% filter(pi_score_Cancer_vs_Healthy < 0.05),
  "FDR < 0.05"   = all_df %>% filter(adj.P.Val_Cancer_vs_Healthy < 0.05),
  "FDR < 0.10"   = all_df %>% filter(adj.P.Val_Cancer_vs_Healthy < 0.10),
  "p < 0.05"     = all_df %>% filter(P.Value_Cancer_vs_Healthy < 0.05)
)

thresh_df <- map_dfr(names(thresholds), function(thr) {
  d <- thresholds[[thr]]
  if (nrow(d) < 10) return(tibble(threshold = thr))
  n_u <- sum(d$cancer_dir == "Cancer Up")
  n_d <- sum(d$cancer_dir == "Cancer Down")
  r_u <- sum(d$reversed[d$cancer_dir == "Cancer Up"])
  r_d <- sum(d$reversed[d$cancer_dir == "Cancer Down"])
  pct_u <- 100 * r_u / max(n_u, 1)
  pct_d <- 100 * r_d / max(n_d, 1)
  pt <- if (n_u >= 5 && n_d >= 5) {
    prop.test(c(r_d, r_u), c(n_d, n_u), alternative = "greater")$p.value
  } else NA_real_
  tibble(threshold = thr, n_total = n_u + n_d,
         pct_up = pct_u, pct_down = pct_d,
         delta = pct_d - pct_u, p_value = pt)
})

# ── Export ───────────────────────────────────────────────────────────────────
asymmetry_summary <- tibble(
  cancer_up_n = n_up, cancer_up_reversed = rev_up, cancer_up_pct = round(pct_up, 1),
  cancer_dn_n = n_dn, cancer_dn_reversed = rev_dn, cancer_dn_pct = round(pct_dn, 1),
  delta_pct = round(delta, 1),
  prop_test_p = prop_test$p.value,
  chi_sq_p = chi_test$p.value,
  boot_ci_lower = round(boot_ci[1], 1),
  boot_ci_upper = round(boot_ci[2], 1)
)
write_csv(asymmetry_summary, file.path(DAT, "SUPP_directional_asymmetry.csv"))
write_csv(thresh_df, file.path(DAT, "SUPP_asymmetry_threshold_sensitivity.csv"))

# ── Visualization ────────────────────────────────────────────────────────────
# Panel 1: Side-by-side bars showing reversal rate by cancer direction
bar_df <- tibble(
  direction = factor(c("Cancer Up", "Cancer Down"),
                     levels = c("Cancer Up", "Cancer Down")),
  pct = c(pct_up, pct_dn),
  n_rev = c(rev_up, rev_dn),
  n_tot = c(n_up, n_dn)
)

p_bars <- ggplot(bar_df, aes(x = direction, y = pct, fill = direction)) +
  geom_col(width = 0.55, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)", pct, n_rev, n_tot)),
            vjust = -0.3, size = 3, fontface = "bold", color = "grey15") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey50",
             linewidth = 0.3) +
  scale_fill_manual(values = c("Cancer Up" = "#E57373", "Cancer Down" = "#64B5F6")) +
  # Significance bracket
  annotate("segment", x = 1, xend = 2, y = max(pct_up, pct_dn) + 8,
           yend = max(pct_up, pct_dn) + 8, linewidth = 0.4) +
  annotate("text", x = 1.5, y = max(pct_up, pct_dn) + 10,
           label = sprintf("\u0394 = %.1f pp\np %s", delta, fmt_p(prop_test$p.value)),
           size = 3, fontface = "bold.italic", color = "grey25") +
  labs(x = NULL, y = "% Reversed") +
  scale_y_continuous(limits = c(0, max(pct_up, pct_dn) + 18),
                     breaks = seq(0, 100, 25)) +
  FIG_THEME +
  theme(legend.position = "none")

# Panel 2: Bootstrap distribution of asymmetry delta
p_boot <- ggplot(tibble(d = boot_delta), aes(x = d)) +
  geom_histogram(bins = 50, fill = "grey70", color = "grey40", linewidth = 0.2) +
  annotate("rect", xmin = boot_ci[1], xmax = boot_ci[2],
           ymin = -Inf, ymax = Inf, fill = "#2563EB", alpha = 0.12) +
  geom_vline(xintercept = delta, color = "#2563EB", linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50",
             linewidth = 0.3) +
  annotate("text", x = delta, y = Inf,
           label = sprintf("\u0394 = %.1f pp\n95%% CI [%.1f, %.1f]",
                            delta, boot_ci[1], boot_ci[2]),
           hjust = -0.1, vjust = 1.5, size = 3, fontface = "bold", color = "#2563EB") +
  labs(x = "Asymmetry \u0394 (% Down reversed \u2212 % Up reversed)",
       y = "Count") +
  FIG_THEME

# Panel 3: Threshold sensitivity
thresh_long <- thresh_df %>%
  filter(!is.na(pct_up)) %>%
  pivot_longer(cols = c(pct_up, pct_down), names_to = "dir",
               values_to = "pct") %>%
  mutate(dir = ifelse(dir == "pct_up", "Cancer Up", "Cancer Down"),
         threshold = factor(threshold, levels = names(thresholds)))

p_thresh <- ggplot(thresh_long, aes(x = threshold, y = pct, fill = dir)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.2) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey50",
             linewidth = 0.3) +
  geom_text(data = thresh_df %>% filter(!is.na(p_value)),
            aes(x = threshold, y = pmax(pct_up, pct_down) + 4,
                label = ifelse(p_value < 0.05, sprintf("p %s *", fmt_p(p_value)),
                                sprintf("p %s", fmt_p(p_value)))),
            inherit.aes = FALSE, size = 2.5, fontface = "bold.italic",
            color = "grey30") +
  scale_fill_manual(values = c("Cancer Up" = "#E57373", "Cancer Down" = "#64B5F6"),
                    name = NULL) +
  labs(x = "Cancer Signature Threshold", y = "% Reversed") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  FIG_THEME +
  theme(axis.text.x = element_text(size = 7, angle = 30, hjust = 1),
        legend.position = "bottom",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 7))

# Combine 3 panels
pS_asym <- p_bars + p_boot + p_thresh +
  plot_layout(widths = c(0.8, 1.2, 1.2)) +
  plot_annotation(
    title = "Directional Asymmetry in Signature Reversal",
    subtitle = sprintf("Cancer-Down proteins reverse %.1f%% vs Cancer-Up %.1f%% (\u0394 = %.1f pp, p %s) | 95%% bootstrap CI [%.1f, %.1f]",
                        pct_dn, pct_up, delta, fmt_p(prop_test$p.value),
                        boot_ci[1], boot_ci[2]),
    theme = theme(plot.title = element_text(size = FIG_TITLE_SIZE, face = "bold"),
                  plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, color = "grey30")))

ggsave(file.path(RPT_PNG, "SUPP_directional_asymmetry.png"), pS_asym,
       width = 320, height = 110, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_directional_asymmetry.pdf"), pS_asym,
       width = 320, height = 110, units = "mm", device = pdf_device)

message("Done: SUPP_directional_asymmetry")
