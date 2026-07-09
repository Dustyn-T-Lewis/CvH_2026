# Supplementary Enrichment Gallery -- Pathway Fate Alluvial (F04 CRvH: Concordance)
# Left = CvH direction, Right = TR direction. Flows colored by pattern.
setwd(here::here())
source("04_Figures/F04/a_script/style.R")

pacman::p_load(tidyverse, ggalluvial)

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

dir_label <- function(nes, sig) {
  case_when(
    !sig           ~ "NS",
    nes > 0        ~ "Up",
    nes < 0        ~ "Down",
    TRUE           ~ "NS"
  )
}

conc <- readRDS(file.path(DAT, "prep_concordance.rds"))

conc_sankey <- conc %>%
  mutate(
    CvH_dir = dir_label(NES_CvH, sig_CvH) %>% factor(levels = c("Up", "Down", "NS")),
    TR_dir  = dir_label(NES_TR, sig_TR)    %>% factor(levels = c("Up", "Down", "NS"))
  ) %>%
  count(CvH_dir, TR_dir, pattern) %>%
  filter(n > 0)

p_conc <- ggplot(conc_sankey,
                  aes(axis1 = CvH_dir, axis2 = TR_dir, y = n)) +
  geom_alluvium(aes(fill = pattern), width = 1/6, alpha = 0.7,
                curve_type = "sigmoid") +
  geom_stratum(width = 1/6, fill = "grey90", color = "grey40") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.5) +
  scale_x_discrete(limits = c("Cancer vs Healthy", "Training (CR)"),
                   expand = c(0.15, 0.05)) +
  scale_fill_manual(values = PATTERN_COLORS, name = "Pattern") +
  labs(
    title    = "Pathway Fate: Cancer vs Healthy \u2192 Training CR",
    subtitle = "Concordance = flows preserving direction; NS = pathway not significant",
    y = "Number of pathways"
  ) +
  FIG_THEME +
  theme(
    legend.position = "bottom",
    axis.text.y     = element_text(size = 9),
    panel.grid      = element_blank(),
    panel.border    = element_rect(color = "grey70", fill = NA, linewidth = 0.3)
  )

ggsave(file.path(RPT, "c_sankey_SUPP.pdf"), p_conc,
       width = 180, height = 140, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "c_sankey_SUPP.png"), p_conc,
       width = 180, height = 140, units = "mm", dpi = 300)

cat("Sankey alluvial plot saved.\n")
