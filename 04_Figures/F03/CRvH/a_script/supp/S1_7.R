# F03/CRvH — Supplementary S1.7: Pi-score Distributions
# Contrasts: Cancer_vs_Healthy, Training_CR
# Outputs: s17 (patchwork composite)

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F03/a_script/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(patchwork)
})

DEP_FILE <- "03_DEP/c_data/03_combined_results_CRvH.csv"
RPT      <- "04_Figures/F03/CRvH/b_reports/supp"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

dep_df     <- read_csv(DEP_FILE, show_col_types = FALSE)
pdf_device <- get_pdf_device()

pi_long <- dep_df |>
  dplyr::select(gene, starts_with("pi_score_")) |>
  pivot_longer(starts_with("pi_score_"), names_to = "contrast",
               values_to = "pi_score") |>
  mutate(contrast = str_remove(contrast, "pi_score_")) |>
  filter(!is.na(pi_score))

p_hist <- ggplot(pi_long, aes(x = pi_score)) +
  geom_histogram(bins = 50, fill = "grey60", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red",
             linewidth = 0.4) +
  facet_wrap(~ contrast, scales = "free_y", ncol = 1,
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
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red",
             linewidth = 0.4) +
  geom_text(data = n_sig, aes(label = sprintf("n = %d", n)),
            x = Inf, y = 0.10, hjust = 1.1, vjust = 0, size = 2.5,
            color = "red") +
  facet_wrap(~ contrast, scales = "free_x", ncol = 1,
             labeller = labeller(contrast = CTR_AXIS)) +
  labs(x = "Protein rank", y = expression(bold(Pi*"-score"))) +
  FIG_THEME

s17 <- p_hist / p_rank +
  plot_annotation(
    title = expression(bold("S1.7  ") * Pi * bold("-score distributions (CRvH)")),
    theme = theme(plot.title = element_text(size = 10))
  )

ggsave(file.path(RPT, "S1_7_pi_score_distributions_SUPP.pdf"), s17,
       width = 180, height = 200, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "S1_7_pi_score_distributions_SUPP.png"), s17,
       width = 180, height = 200, units = "mm", dpi = 300)

cat("F03/CRvH Supp S1.7 done.\n")
