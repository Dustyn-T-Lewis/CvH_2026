# Supplementary Enrichment Gallery -- NES Quadrant Scatter (F04 CRvH: Concordance)
# NES vs NES scatter with quadrant backgrounds and pattern coloring.
setwd(here::here())
source("04_Figures/F04/a_script/style.R")

pacman::p_load(tidyverse, ggrepel)

RPT <- "04_Figures/F04/CRvH/b_reports/supp"
DAT <- "04_Figures/F04/CRvH/c_data/supp"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

PATTERN_COLORS <- c(
  Concordant           = "#2E7D32",
  Discordant           = "#FF8F00",
  "Cancer-specific"    = "#4CAF50",
  "Training-specific"  = "#9C27B0",
  Other                = "grey60"
)

# Cancer_vs_Healthy vs Training_CR
conc <- readRDS(file.path(DAT, "prep_concordance.rds"))

conc <- conc %>%
  mutate(
    label = case_when(
      pattern == "Concordant" & rank(-abs(NES_CvH)) <= 5        ~ pathway_label,
      pattern == "Cancer-specific" & rank(-abs(NES_CvH)) <= 4   ~ pathway_label,
      pattern == "Training-specific" & rank(-abs(NES_TR)) <= 3  ~ pathway_label,
      pattern == "Discordant" & rank(-abs(NES_CvH)) <= 3        ~ pathway_label,
      TRUE ~ NA_character_
    )
  )

r_val  <- cor(conc$NES_CvH, conc$NES_TR, use = "complete.obs")
n_conc <- sum(sign(conc$NES_CvH) == sign(conc$NES_TR), na.rm = TRUE)
pct_conc <- round(100 * n_conc / nrow(conc), 1)

p_conc <- ggplot(conc, aes(NES_CvH, NES_TR)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "#E8F5E9", alpha = 0.4) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#E8F5E9", alpha = 0.4) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#E3F2FD", alpha = 0.4) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "#E3F2FD", alpha = 0.4) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_point(aes(color = pattern, size = set_size), alpha = 0.7) +
  geom_text_repel(aes(label = label), size = 2.2, max.overlaps = 30,
                  segment.size = 0.2, segment.color = "grey50",
                  box.padding = 0.6, min.segment.length = 0.2,
                  force = 2, force_pull = 0.5, seed = 42) +
  scale_color_manual(values = PATTERN_COLORS, name = "Pattern") +
  scale_size_continuous(range = c(1, 5), name = "Set size", guide = "none") +
  labs(
    title    = "Pathway-Level Response: Cancer vs Healthy vs Training CR",
    subtitle = sprintf("r = %.2f, %d/%d concordant (%.1f%%)",
                        r_val, n_conc, nrow(conc), pct_conc),
    x = "NES \u2014 Cancer vs Healthy",
    y = "NES \u2014 Training (CR)"
  ) +
  FIG_THEME +
  theme(legend.position = "bottom",
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3)) +
  coord_fixed()

ggsave(file.path(RPT, "b_nes_scatter_SUPP.pdf"), p_conc,
       width = 200, height = 200, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "b_nes_scatter_SUPP.png"), p_conc,
       width = 200, height = 200, units = "mm", dpi = 300)

cat("NES scatter plot saved.\n")
