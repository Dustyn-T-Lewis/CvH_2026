#!/usr/bin/env Rscript
# Sourced by 02_supp_panels.R — expects style.R already loaded.
# F05 Supplementary Panel E: GO Slim Category Breakdown for Reversal
# Supports main Panel B — horizontal stacked bar showing GO Slim categories
# by reversal quadrant: Reversed Up, Reversed Down, Non-reversed.

source("02-03_Sam's_Results/04_Figures/shared/go_slim_categories.R")

suppressPackageStartupMessages({
  library(tidyverse)
})

BASE    <- "02-03_Sam's_Results/04_Figures/F05_recovery_reversal"
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
DAT     <- file.path(BASE, "c_data", "panel_supp")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT,     recursive = TRUE, showWarnings = FALSE)
pdf_device <- grDevices::pdf  # force base pdf (cairo DLL fails)

# Data
dep <- read_csv("02-03_Sam's_Results/03_DEP/c_data/F05_combined_CvHvTCR.csv",
                show_col_types = FALSE)

fc_df <- dep |>
  transmute(gene,
            logFC_CvH = logFC_Cancer_vs_Healthy,
            logFC_TCR = logFC_Training_CR,
            pi_CvH    = pi_score_Cancer_vs_Healthy,
            pi_TCR    = pi_score_Training_CR) |>
  filter(!is.na(logFC_CvH), !is.na(logFC_TCR)) |>
  filter(pi_CvH < 0.05 | pi_TCR < 0.05) |>
  mutate(quadrant = case_when(
    logFC_CvH > 0 & logFC_TCR < 0 ~ "Reversed Up",
    logFC_CvH < 0 & logFC_TCR > 0 ~ "Reversed Down",
    TRUE                            ~ "Non-reversed"
  ))

all_genes <- dep$gene[!is.na(dep$logFC_Cancer_vs_Healthy)]
fg_genes  <- fc_df$gene

# GO Slim assignment
slim_result <- assign_go_slim_consolidated(fg_genes, all_genes)

slim_merged <- fc_df |>
  left_join(slim_result, by = "gene") |>
  filter(!is.na(consolidated))

# Summary counts
slim_counts <- slim_merged |>
  count(consolidated, quadrant, name = "count") |>
  mutate(consolidated = factor(consolidated, levels = rev(CONSOLIDATED_PATHWAY_ORDER)))

write_csv(slim_counts, file.path(DAT, "SUPP_goslim_distribution.csv"))

# Plot
quad_colors <- c("Reversed Up" = "#E57373", "Reversed Down" = "#64B5F6",
                 "Non-reversed" = "grey60")
quad_levels <- c("Reversed Up", "Reversed Down", "Non-reversed")
slim_counts <- slim_counts |>
  mutate(quadrant = factor(quadrant, levels = quad_levels))

pS_goslim <- ggplot(slim_counts,
                    aes(x = count, y = consolidated, fill = quadrant)) +
  geom_col(position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = quad_colors, name = "Quadrant") +
  labs(title    = "GO Slim Category Distribution (Reversal)",
       subtitle = sprintf("DEPs (Pi < 0.05) | %d proteins | %d categories",
                          nrow(slim_merged),
                          n_distinct(slim_merged$consolidated)),
       x = "Protein count", y = NULL) +
  FIG_THEME +
  theme(legend.position    = "right",
        axis.text.y        = element_text(size = FIG_AXIS_TEXT),
        panel.grid.major.y = element_blank())

PW <- 89; PH <- 85
ggsave(file.path(RPT_PNG, "SUPP_goslim_bars.png"), pS_goslim,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_goslim_bars.pdf"), pS_goslim,
       width = PW, height = PH, units = "mm", device = pdf_device)

message("F05 SUPP Panel E (GO Slim bars) saved")

# Expose for composite
pS_goslim_title    <- "GO Slim Category Distribution (Reversal)"
pS_goslim_subtitle <- sprintf("DEPs (Pi < 0.05) | %d proteins | %d categories",
                               nrow(slim_merged),
                               n_distinct(slim_merged$consolidated))
pS_goslim          <- strip_for_composite(pS_goslim)
