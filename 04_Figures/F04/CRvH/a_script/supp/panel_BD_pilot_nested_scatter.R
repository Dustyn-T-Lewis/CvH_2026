# F04 CRvH Supplementary: Nested Scatter (Protein Concordance + NES Quadrant Insets)
# Main scatter: logFC Cancer_vs_Healthy vs Training_CR (protein level)
# Insets: NES scatter per quadrant showing pathway-level concordance
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F04/a_script/style.R")

library(tidyverse)

RPT <- "04_Figures/F04/CRvH/b_reports/supp"
DAT <- "04_Figures/F04/CRvH/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

# --- 1. Load data
dep <- readr::read_csv("03_DEP/a_non_imputed/c_data/combined_results_pi.csv",
                       show_col_types = FALSE) |>
  dplyr::mutate(contrast = dplyr::recode(contrast,
                                         CRvH_Baseline = "Cancer_vs_Healthy",
                                         CR_Training   = "Training_CR")) |>
  tidyr::pivot_wider(id_cols = c(uniprot_id, gene, protein, description),
                     names_from = contrast,
                     values_from = c(logFC, t, P.Value, adj.P.Val, pi_score, sig_pi),
                     names_glue = "{.value}_{contrast}")
nes <- read_csv(file.path(DAT, "panel_D", "nes_scatter.csv"), show_col_types = FALSE)

# --- 2. Build main scatter data (protein-level)
scatter_df <- dep %>%
  transmute(
    gene,
    x = logFC_Cancer_vs_Healthy,
    y = logFC_Training_CR,
    pi_CvH = pi_score_Cancer_vs_Healthy,
    pi_TR  = pi_score_Training_CR
  ) %>%
  filter(!is.na(x), !is.na(y)) %>%
  mutate(sig_class = classify_proteins_f4(pi_CvH, pi_TR))

# --- 3. Build NES quadrant data (pathway-level)
nes_wide <- nes %>%
  filter(!is.na(NES_Cancer_vs_Healthy), !is.na(NES_Training_CR)) %>%
  mutate(
    sig_CvH = padj_Cancer_vs_Healthy < 0.05,
    sig_TR  = padj_Training_CR < 0.05,
    nes_sig = case_when(
      sig_CvH & sig_TR ~ "Sig Both",
      sig_CvH          ~ "Sig Cancer only",
      sig_TR           ~ "Sig Training only",
      TRUE             ~ "NS"
    ),
    quadrant = case_when(
      NES_Cancer_vs_Healthy > 0 & NES_Training_CR > 0 ~ "CU",
      NES_Cancer_vs_Healthy < 0 & NES_Training_CR < 0 ~ "CD",
      NES_Cancer_vs_Healthy < 0 & NES_Training_CR > 0 ~ "DCdTu",
      TRUE ~ "DCuTd"
    )
  )

# NES inset colors
NES_COLORS <- c(
  "Sig Both"          = "#2E7D32",
  "Sig Cancer only"   = "#4CAF50",
  "Sig Training only" = "#9C27B0",
  "NS"                = "grey70"
)

# --- 4. Main scatter plot
p_main <- ggplot(scatter_df, aes(x, y, color = sig_class)) +
  geom_hline(yintercept = 0, color = "grey80") +
  geom_vline(xintercept = 0, color = "grey80") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 0.8, alpha = 0.5) +
  scale_color_manual(values = SIG_COLORS_F4, name = "Significance") +
  labs(
    x = expression(log[2]~FC~"(Cancer vs Healthy)"),
    y = expression(log[2]~FC~"(Training CR)"),
    title = "Protein Concordance with Pathway-Level Insets",
    subtitle = "NES scatter insets show pathway concordance per quadrant"
  ) +
  FIG_THEME +
  theme(
    legend.position   = c(0.15, 0.85),
    legend.background = element_rect(fill = alpha("white", 0.9), color = NA)
  )

# --- 5. Build mini NES scatter for each quadrant
make_inset <- function(quad_data) {
  if (nrow(quad_data) == 0) return(NULL)
  ggplot(quad_data, aes(NES_Cancer_vs_Healthy, NES_Training_CR, color = nes_sig)) +
    geom_hline(yintercept = 0, color = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, color = "grey85", linewidth = 0.3) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_color_manual(values = NES_COLORS, guide = "none") +
    annotate("text", x = Inf, y = -Inf,
             label = paste0("n=", nrow(quad_data)),
             hjust = 1.2, vjust = -0.5, size = 2.5, fontface = "bold",
             color = "grey40") +
    theme_minimal(base_size = 6) +
    theme(
      panel.border     = element_rect(color = "grey50", fill = NA, linewidth = 0.5),
      panel.background = element_rect(fill = alpha("white", 0.85), color = NA),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      axis.title       = element_blank(),
      axis.text        = element_text(size = 5, color = "grey40"),
      plot.margin      = margin(1, 1, 1, 1)
    )
}

inset_CU    <- make_inset(nes_wide %>% filter(quadrant == "CU"))
inset_CD    <- make_inset(nes_wide %>% filter(quadrant == "CD"))
inset_DCdTu <- make_inset(nes_wide %>% filter(quadrant == "DCdTu"))
inset_DCuTd <- make_inset(nes_wide %>% filter(quadrant == "DCuTd"))

# --- 6. Position insets at extreme corners (away from origin/density)
xrng <- range(scatter_df$x, na.rm = TRUE)
yrng <- range(scatter_df$y, na.rm = TRUE)

# Inset span: ~35% of each axis range, anchored at the outer edge
x_span <- diff(xrng) * 0.35
y_span <- diff(yrng) * 0.35

# Top-right corner (Concordant Up)
if (!is.null(inset_CU))
  p_main <- p_main + annotation_custom(
    ggplotGrob(inset_CU),
    xmin = xrng[2] - x_span, xmax = xrng[2],
    ymin = yrng[2] - y_span, ymax = yrng[2]
  )

# Bottom-left corner (Concordant Down)
if (!is.null(inset_CD))
  p_main <- p_main + annotation_custom(
    ggplotGrob(inset_CD),
    xmin = xrng[1], xmax = xrng[1] + x_span,
    ymin = yrng[1], ymax = yrng[1] + y_span
  )

# Top-left corner (Discordant C Down / T Up)
if (!is.null(inset_DCdTu))
  p_main <- p_main + annotation_custom(
    ggplotGrob(inset_DCdTu),
    xmin = xrng[1], xmax = xrng[1] + x_span,
    ymin = yrng[2] - y_span, ymax = yrng[2]
  )

# Bottom-right corner (Discordant C Up / T Down)
if (!is.null(inset_DCuTd))
  p_main <- p_main + annotation_custom(
    ggplotGrob(inset_DCuTd),
    xmin = xrng[2] - x_span, xmax = xrng[2],
    ymin = yrng[1], ymax = yrng[1] + y_span
  )

# --- 7. Ensure points render above insets
p_final <- p_main +
  # Re-draw points on top of annotation_custom grobs
  geom_point(data = scatter_df, aes(x, y, color = sig_class),
             size = 0.8, alpha = 0.5, show.legend = FALSE) +
  coord_cartesian(clip = "off")

# --- 8. Save
ggsave(file.path(RPT, "panel_BD_pilot_nested_scatter_SUPP.pdf"), p_final,
       width = 220, height = 220, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_BD_pilot_nested_scatter_SUPP.png"), p_final,
       width = 220, height = 220, units = "mm", dpi = 300)

message("Done: ", file.path(RPT, "panel_BD_pilot_nested_scatter_SUPP.{pdf,png}"))
