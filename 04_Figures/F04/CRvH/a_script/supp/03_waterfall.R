# Supplementary Enrichment Gallery -- Concordance Ratio Lollipop (F04 CRvH)
# All pathways included; alpha encodes denominator significance.
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F04/a_script/style.R")

library(tidyverse)

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

# Concordance ratio (NES_TR / NES_CvH) for all pathways
conc <- readRDS(file.path(DAT, "prep_concordance.rds"))

conc_lollipop <- conc %>%
  filter(abs(NES_CvH) > 0.01) %>%
  mutate(
    ratio = pmin(pmax(NES_TR / NES_CvH, -1.5), 1.5),
    sig_denom = sig_CvH,
    pathway_label = clean_pathway_name(pathway, max_chars = 50),
    pathway_label = fct_reorder(pathway_label, ratio)
  )

n_pw <- nrow(conc_lollipop)
fig_h <- max(120, 5 * n_pw + 30)

p_conc <- ggplot(conc_lollipop, aes(x = ratio, y = pathway_label)) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = Inf,
           fill = "#FFCDD2", alpha = 0.08) +
  annotate("rect", xmin = 0, xmax = 1, ymin = -Inf, ymax = Inf,
           fill = "#FFE0B2", alpha = 0.08) +
  annotate("rect", xmin = 1, xmax = Inf, ymin = -Inf, ymax = Inf,
           fill = "#C8E6C9", alpha = 0.08) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#2E7D32",
             linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#E05A4E",
             linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = ratio, yend = pathway_label,
                   color = pattern, alpha = ifelse(sig_denom, 1, 0.4)),
               linewidth = 0.5) +
  geom_point(aes(color = pattern, alpha = ifelse(sig_denom, 1, 0.4)), size = 2) +
  annotate("text", x = 1, y = Inf, label = "Concordant", hjust = -0.1,
           vjust = 1.5, size = 3.2, color = "#2E7D32", fontface = "italic") +
  annotate("text", x = 0, y = Inf, label = "No training effect", hjust = 1.1,
           vjust = 1.5, size = 3.2, color = "#E05A4E", fontface = "italic") +
  scale_color_manual(values = PATTERN_COLORS, name = "Pattern") +
  scale_alpha_identity() +
  labs(
    title    = "Concordance Ratio: Training CR / Cancer vs Healthy",
    subtitle = sprintf("NES(TR) / NES(CvH) for all %d pathways (ratio capped at \u00b11.5)", n_pw),
    x = "NES ratio (Training CR / Cancer vs Healthy)",
    y = NULL
  ) +
  FIG_THEME +
  theme(
    axis.text.y = element_text(size = 6.5),
    legend.position = "bottom",
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.3)
  )

ggsave(file.path(RPT, "f_waterfall_SUPP.pdf"), p_conc,
       width = 200, height = fig_h, units = "mm", device = pdf_device,
       limitsize = FALSE)
ggsave(file.path(RPT, "f_waterfall_SUPP.png"), p_conc,
       width = 200, height = fig_h, units = "mm", dpi = 300,
       limitsize = FALSE)

cat("Waterfall lollipop plot saved.\n")
