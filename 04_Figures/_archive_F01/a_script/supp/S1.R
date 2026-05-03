# Figure 1 — Supplementary S1: Pi-score Distributions
# 5 contrasts from two models. Histogram + ranked scatter.
# Outputs: s1 (patchwork composite), supp_S1_pi_score_distributions.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(patchwork)
})

RPT <- "04_Figures/F01/b_reports"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(RPT, "supp"), showWarnings = FALSE)

pdf_device <- get_pdf_device()

# Read both DEP result files
dep_crvh <- read_csv("03_DEP/c_data/03_combined_results_CRvH.csv", show_col_types = FALSE)
dep_cr   <- read_csv("03_DEP/c_data/03_combined_results_CR.csv",   show_col_types = FALSE)

CONTRASTS <- c("Cancer_vs_Healthy", "Training_CR",
               "Baseline_Supplement", "Training_CRE",
               "Training_PLA", "Supplement_Interaction")

# Extract pi_score columns from each file
pi_crvh <- dep_crvh |>
  dplyr::select(gene, starts_with("pi_score_")) |>
  pivot_longer(starts_with("pi_score_"), names_to = "contrast", values_to = "pi_score") |>
  mutate(contrast = str_remove(contrast, "pi_score_"))

pi_cr <- dep_cr |>
  dplyr::select(gene, starts_with("pi_score_")) |>
  pivot_longer(starts_with("pi_score_"), names_to = "contrast", values_to = "pi_score") |>
  mutate(contrast = str_remove(contrast, "pi_score_"))

pi_long <- bind_rows(pi_crvh, pi_cr) |>
  filter(!is.na(pi_score), contrast %in% CONTRASTS)
pi_long$contrast <- factor(pi_long$contrast, levels = CONTRASTS)

p_hist <- ggplot(pi_long, aes(x = pi_score)) +
  geom_histogram(bins = 50, fill = "grey60", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.4) +
  facet_wrap(~ contrast, scales = "free_y", ncol = 3,
             labeller = labeller(contrast = CTR_AXIS)) +
  labs(x = expression(bold(Pi*"-score")), y = "Count") +
  FIG_THEME

pi_ranked <- pi_long |>
  group_by(contrast) |> arrange(pi_score) |>
  mutate(rank = row_number()) |> ungroup()

n_sig <- pi_long |>
  group_by(contrast) |>
  summarise(n = sum(pi_score < 0.05), .groups = "drop")

p_rank <- ggplot(pi_ranked, aes(x = rank, y = pi_score)) +
  geom_point(size = 0.2, alpha = 0.4, color = "grey40") +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.4) +
  geom_text(data = n_sig, aes(label = sprintf("n = %d", n)),
            x = Inf, y = 0.10, hjust = 1.1, vjust = 0, size = scale_text(BASE_STAT - 1, 250), color = "red") +
  facet_wrap(~ contrast, scales = "free_x", ncol = 3,
             labeller = labeller(contrast = CTR_AXIS)) +
  labs(x = "Protein rank", y = expression(bold(Pi*"-score"))) +
  FIG_THEME

s1 <- p_hist / p_rank +
  plot_annotation(
    title = expression(bold("S1  ") * Pi * bold("-score distributions")),
    theme = theme(plot.title = element_text(size = 10)))

ggsave(file.path(RPT, "supp", "S1_pi_score_distributions.pdf"), s1,
       width = 250, height = 240, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "supp", "S1_pi_score_distributions.png"), s1,
       width = 250, height = 240, units = "mm", dpi = 300)

cat("Supplementary S1 done.\n")
