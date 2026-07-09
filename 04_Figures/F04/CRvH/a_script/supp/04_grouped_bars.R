# Supplementary Enrichment Gallery -- Diverging Bars by Biological Theme (F04 CRvH)
# All significant pathways included (no top-N selection).
setwd(here::here())
source("04_Figures/F04/a_script/style.R")

pacman::p_load(tidyverse)

RPT <- "04_Figures/F04/CRvH/b_reports/supp"
DAT <- "04_Figures/F04/CRvH/c_data/supp"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

conc <- readRDS(file.path(DAT, "prep_concordance.rds"))

# All sig pathways (no top-N selection)
conc_sel <- conc %>%
  filter(sig_CvH | sig_TR)

# Reshape to long for paired bars
conc_bars <- conc_sel %>%
  select(pathway, pathway_label, bio_theme, NES_CvH, NES_TR, padj_CvH, padj_TR,
         sig_CvH, sig_TR) %>%
  pivot_longer(
    cols      = c(NES_CvH, NES_TR),
    names_to  = "contrast",
    values_to = "NES"
  ) %>%
  mutate(
    contrast = recode(contrast,
                      NES_CvH = "Cancer vs Healthy",
                      NES_TR  = "Training (CR)"),
    sig = ifelse(contrast == "Cancer vs Healthy", sig_CvH, sig_TR),
    alpha_val = ifelse(sig, 1, 0.4),
    pathway_label = clean_pathway_name(pathway, max_chars = 40)
  )

n_rows <- n_distinct(conc_bars$pathway)
fig_h <- max(150, 6 * n_rows + 40)

p_conc <- ggplot(conc_bars,
                  aes(x = NES,
                      y = reorder_within(pathway_label, NES, bio_theme),
                      fill = contrast, alpha = alpha_val)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  facet_grid(bio_theme ~ ., scales = "free_y", space = "free_y") +
  scale_y_reordered() +
  scale_fill_manual(values = c("Cancer vs Healthy" = unname(CONTRAST_COLORS["Cancer_vs_Healthy"]),
                                "Training (CR)"    = unname(CONTRAST_COLORS["Training_CR"])),
                    name = "Contrast") +
  scale_alpha_identity() +
  labs(
    title    = "Concordance Response by Biological Theme",
    subtitle = sprintf("All %d significant pathways; opacity reflects padj < 0.05", n_rows),
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  FIG_THEME +
  theme(
    axis.text.y    = element_text(size = 7),
    strip.text.y   = element_text(angle = 0, size = 7, hjust = 0),
    legend.position = "bottom",
    panel.spacing  = unit(2, "mm"),
    panel.border   = element_rect(color = "grey70", fill = NA, linewidth = 0.3)
  )

ggsave(file.path(RPT, "g_grouped_bars_SUPP.pdf"), p_conc,
       width = 200, height = fig_h, units = "mm", device = pdf_device,
       limitsize = FALSE)
ggsave(file.path(RPT, "g_grouped_bars_SUPP.png"), p_conc,
       width = 200, height = fig_h, units = "mm", dpi = 300,
       limitsize = FALSE)

cat("Grouped bar plot saved.\n")
