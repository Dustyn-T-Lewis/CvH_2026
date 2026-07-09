# Reversal Panel B: Pathway NES Scatter
# Cancer_vs_Healthy (x) vs Training_CR (y) at pathway level
# fGSEA cache is LONG format — pivot wider before scatter
setwd(here::here())
source("04_Figures/shared/style.R")
library(tidyverse)
library(ggrepel)

RPT_PNG <- "04_Figures/Reversal/b_reports/main/png/panels"
RPT_PDF <- "04_Figures/Reversal/b_reports/main/pdf/panels"
DAT     <- "04_Figures/Reversal/c_data"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_B"), recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

COMP_RED   <- unname(DIR_COLORS["Up"])
COMP_BLUE  <- unname(DIR_COLORS["Down"])
PG_W <- 146

# --- Load fGSEA cache (LONG format) ---
fgsea_cache <- "04_Figures/shared/fgsea_CRvH.csv"
stopifnot("fGSEA cache missing" = file.exists(fgsea_cache))
fgsea_long <- read_csv(fgsea_cache, show_col_types = FALSE)

# --- Pivot to wide ---
fgsea_wide <- fgsea_long %>%
  dplyr::select(pathway, NES, padj, size, database, contrast) %>%
  pivot_wider(id_cols = c(pathway, database),
              names_from = contrast,
              values_from = c(NES, padj, size),
              names_glue = "{.value}_{contrast}") %>%
  filter(!is.na(NES_Cancer_vs_Healthy), !is.na(NES_Training_CR))

# --- Filter to Hallmark + GO Slim ---
fgsea_hg <- fgsea_wide %>%
  filter(database %in% c("Hallmark", "GO Slim"))

fgsea_hg <- fgsea_hg %>%
  mutate(
    set_size = coalesce(size_Cancer_vs_Healthy, size_Training_CR),
    sig_1 = !is.na(padj_Cancer_vs_Healthy) & padj_Cancer_vs_Healthy < 0.05,
    sig_2 = !is.na(padj_Training_CR) & padj_Training_CR < 0.05,
    significance = case_when(
      sig_1 & sig_2 ~ "Sig Both",
      sig_1         ~ "Sig Cancer only",
      sig_2         ~ "Sig Training only",
      TRUE          ~ "NS"
    ) %>% factor(levels = names(SIG_COLORS_F4)),
    pathway_label = clean_pathway_name(pathway),
    db_shape = ifelse(database == "Hallmark", 24, 21)
  )

fgsea_sig <- fgsea_hg %>% filter(significance != "NS")

message(sprintf("  %d total pathways (Hallmark: %d, GO Slim: %d) | %d significant",
                nrow(fgsea_hg),
                sum(fgsea_hg$database == "Hallmark"),
                sum(fgsea_hg$database == "GO Slim"),
                nrow(fgsea_sig)))

# --- Spearman correlation ---
nes_cor_all <- cor.test(fgsea_hg$NES_Cancer_vs_Healthy,
                         fgsea_hg$NES_Training_CR, method = "spearman")
nes_ci_all  <- fisher_z_ci(nes_cor_all$estimate, nrow(fgsea_hg))
nes_cor_sig <- if (nrow(fgsea_sig) >= 3) {
  cor.test(fgsea_sig$NES_Cancer_vs_Healthy,
           fgsea_sig$NES_Training_CR, method = "spearman")
} else NULL

nes_lim <- max(abs(c(fgsea_hg$NES_Cancer_vs_Healthy,
                      fgsea_hg$NES_Training_CR))) * 1.15

# --- Quadrant counts ---
n_q1 <- sum(fgsea_sig$NES_Cancer_vs_Healthy > 0 & fgsea_sig$NES_Training_CR > 0)
n_q2 <- sum(fgsea_sig$NES_Cancer_vs_Healthy < 0 & fgsea_sig$NES_Training_CR > 0)
n_q3 <- sum(fgsea_sig$NES_Cancer_vs_Healthy < 0 & fgsea_sig$NES_Training_CR < 0)
n_q4 <- sum(fgsea_sig$NES_Cancer_vs_Healthy > 0 & fgsea_sig$NES_Training_CR < 0)
n_reversed <- n_q2 + n_q4  # off-diagonal
n_total_sig <- nrow(fgsea_sig)
rev_frac <- if (n_total_sig > 0) n_reversed / n_total_sig else 0

message(sprintf("  NES Spearman (all): rho = %.3f [%.3f, %.3f]",
                nes_cor_all$estimate, nes_ci_all[1], nes_ci_all[2]))

# --- Sizes ---
txt_pw   <- 2.2
txt_quad <- 2.8

# --- Label data: all significant pathways ---
label_pw <- fgsea_sig %>%
  mutate(nes_mag = abs(NES_Cancer_vs_Healthy) + abs(NES_Training_CR)) %>%
  arrange(desc(nes_mag)) %>%
  dplyr::select(-nes_mag) %>%
  mutate(
    label_fill     = SIG_LABEL_FILL_F4[as.character(significance)],
    label_text_col = SIG_LABEL_TEXT_F4[as.character(significance)]
  )

# --- Split for layered plotting ---
ns_df  <- fgsea_hg %>% filter(significance == "NS")
sig_df <- fgsea_hg %>% filter(significance != "NS")

# --- Subtitle ---
rho_sig_str <- if (!is.null(nes_cor_sig)) {
  sprintf(", \u03c1(sig) = %.2f", nes_cor_sig$estimate)
} else ""
subtitle_str <- sprintf(
  "GO Slim + Hallmark | %d pathways (%d sig.) | \u03c1 = %.2f [%.2f, %.2f], %s%s\n%.0f%% reversed pathways",
  nrow(fgsea_hg), n_total_sig,
  nes_cor_all$estimate, nes_ci_all[1], nes_ci_all[2],
  ifelse(nes_cor_all$p.value < 0.001, "p < 0.001",
         sprintf("p = %.3f", nes_cor_all$p.value)),
  rho_sig_str, rev_frac * 100)

# --- Build plot ---
pB <- ggplot(mapping = aes(x = NES_Cancer_vs_Healthy, y = NES_Training_CR)) +
  # Quadrant backgrounds: blue = reversed (off-diagonal), red = exacerbated
  annotate("rect", xmin = 0,    xmax = Inf,  ymin = -Inf, ymax = 0,
           fill = "#DCEEFF", alpha = 0.20, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = -Inf, xmax = 0,    ymin = 0,    ymax = Inf,
           fill = "#DCEEFF", alpha = 0.20, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 0,    xmax = Inf,  ymin = 0,    ymax = Inf,
           fill = "#FFE0E0", alpha = 0.20, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = -Inf, xmax = 0,    ymin = -Inf, ymax = 0,
           fill = "#FFE0E0", alpha = 0.20, color = "grey70", linewidth = 0.2) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  # NS points
  geom_point(data = ns_df, aes(shape = database),
             size = 1.0, fill = "grey70", color = "grey55", alpha = 0.40, stroke = 0.2) +
  # Sig points
  geom_point(data = sig_df, aes(fill = significance, size = set_size, shape = database),
             color = ifelse(sig_df$database == "Hallmark", "black", "grey65"),
             alpha = 0.80, stroke = 0.4) +
  scale_fill_manual(values = SIG_COLORS_F4, name = "Significance") +
  scale_shape_manual(values = c("Hallmark" = 24, "GO Slim" = 21), name = "Database") +
  scale_size_continuous(range = c(1.5, 5), name = "Set size",
                        breaks = c(20, 50, 100, 200)) +
  # Pathway labels
geom_label_repel(data = label_pw, aes(label = pathway_label),
                   fill = label_pw$label_fill, color = label_pw$label_text_col,
                   size = txt_pw, fontface = "bold",
                   max.overlaps = 50,
                   segment.size = 0.5, segment.color = "grey20",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 1.0, point.padding = 0.45,
                   force = 35, force_pull = 0.15,
                   label.padding = unit(1, "pt"),
                   label.r = unit(0.5, "pt"),
                   linewidth = 0.10, seed = 42) +
  # Quadrant labels
  annotate("label", x = nes_lim, y = nes_lim,
           label = sprintf("Exacerbated Up  n = %d", n_q1),
           hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
           color = COMP_RED, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -nes_lim, y = nes_lim,
           label = sprintf("Reversed (C\u2193 T\u2191)  n = %d", n_q2),
           hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
           color = COMP_BLUE, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -nes_lim, y = -nes_lim,
           label = sprintf("Exacerbated Dn  n = %d", n_q3),
           hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
           color = COMP_RED, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = nes_lim, y = -nes_lim,
           label = sprintf("Reversed (C\u2191 T\u2193)  n = %d", n_q4),
           hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
           color = COMP_BLUE, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt")) +
  scale_x_continuous(expand = expansion(0, 0)) +
  scale_y_continuous(expand = expansion(0, 0)) +
  coord_fixed(ratio = 1, xlim = c(-nes_lim, nes_lim),
              ylim = c(-nes_lim, nes_lim)) +
  labs(title = "Pathway NES: Cancer Recovery Reversal",
       subtitle = subtitle_str,
       x = "NES (Cancer vs Healthy)",
       y = "NES (Training CR)") +
  FIG_THEME +
  theme(
    axis.text         = element_text(size = FIG_AXIS_TEXT, face = "bold", color = "grey30"),
    axis.title        = element_text(size = FIG_AXIS_TEXT, face = "bold"),
    legend.position   = "bottom",
    legend.title      = element_text(size = FIG_LEGEND_TITLE, face = "bold", color = "grey25"),
    legend.text       = element_text(size = FIG_LEGEND_TEXT, color = "grey20"),
    legend.key.size   = unit(3, "mm"),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box        = "horizontal",
    legend.box.just   = "center",
    legend.spacing.x  = unit(3, "mm"),
    legend.box.margin = margin(t = -2),
    plot.margin       = margin(0, 0, 0, 0)
  ) +
  guides(fill  = "none",
         shape = guide_legend(nrow = 1, order = 1,
                               keyheight = unit(4, "mm"),
                               keywidth  = unit(4, "mm"),
                               override.aes = list(size = 3, fill = "grey50")),
         size  = guide_legend(nrow = 1, order = 2,
                               keyheight = unit(4, "mm"),
                               keywidth  = unit(4, "mm")))

ggsave(file.path(RPT_PNG, "MAIN_panel_B_nes_scatter.png"), pB,
       width = PG_W, height = PG_W, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "MAIN_panel_B_nes_scatter.pdf"), pB,
       width = PG_W, height = PG_W, units = "mm", device = pdf_device)

# --- Export ---
export_df <- fgsea_hg %>%
  transmute(
    pathway, pathway_label, database,
    NES_Cancer_vs_Healthy = round(NES_Cancer_vs_Healthy, 3),
    NES_Training_CR = round(NES_Training_CR, 3),
    padj_Cancer_vs_Healthy = signif(padj_Cancer_vs_Healthy, 4),
    padj_Training_CR = signif(padj_Training_CR, 4),
    significance = as.character(significance),
    set_size
  ) %>%
  arrange(significance, desc(abs(NES_Cancer_vs_Healthy) + abs(NES_Training_CR)))
write_csv(export_df, file.path(DAT, "panel_B", "nes_scatter.csv"))

# --- Export for composite ---
pB_title    <- "Pathway NES: Cancer Recovery Reversal"
pB_subtitle <- subtitle_str
pB_legend   <- NULL
pB          <- pB + labs(title = NULL, subtitle = NULL, tag = NULL)

# Backward-compatible alias
pw_rev_frac <- rev_frac

message("Reversal Panel B NES scatter done")
