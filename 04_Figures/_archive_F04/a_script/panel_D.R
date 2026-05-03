# F04 Panel D: fGSEA NES Scatter (Cancer vs Healthy × Training CR) — Reversal
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

PG_W <- 200
RPT <- "04_Figures/F04/b_reports"
DAT <- "04_Figures/F04/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_D"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

fgsea_all <- read_csv("04_Figures/F02/c_data/06_panel_H_CRvH_fgsea_results.csv",
                       show_col_types = FALSE)

fgsea_hg <- fgsea_all %>%
  filter(contrast %in% c("Cancer_vs_Healthy", "Training_CR"))

fgsea_wide <- fgsea_hg %>%
  dplyr::select(pathway, contrast, NES, padj, size, database) %>%
  pivot_wider(id_cols = c(pathway, database), names_from = contrast,
              values_from = c(NES, padj, size)) %>%
  filter(!is.na(NES_Cancer_vs_Healthy), !is.na(NES_Training_CR)) %>%
  mutate(set_size = coalesce(size_Cancer_vs_Healthy, size_Training_CR))

fgsea_wide <- fgsea_wide %>%
  mutate(
    sig_CvH = !is.na(padj_Cancer_vs_Healthy) & padj_Cancer_vs_Healthy < 0.05,
    sig_TR  = !is.na(padj_Training_CR)       & padj_Training_CR < 0.05,
    significance = case_when(
      sig_CvH & sig_TR ~ "Sig Both",
      sig_CvH          ~ "Sig Cancer only",
      sig_TR           ~ "Sig Training only",
      TRUE             ~ "NS"
    ) %>% factor(levels = names(SIG_COLORS_F4)),
    pathway_label = clean_pathway_name(pathway)
  )

fgsea_sig <- fgsea_wide %>% filter(significance != "NS")

nes_cor_all <- cor.test(fgsea_wide$NES_Cancer_vs_Healthy, fgsea_wide$NES_Training_CR,
                         method = "spearman")
nes_ci_all  <- fisher_z_ci(nes_cor_all$estimate, nrow(fgsea_wide))

nes_lim <- max(abs(c(fgsea_wide$NES_Cancer_vs_Healthy, fgsea_wide$NES_Training_CR))) * 1.15

# Reversal framing: opposite NES = reversed, same NES = exacerbated
n_rev_br <- sum(fgsea_sig$NES_Cancer_vs_Healthy > 0 & fgsea_sig$NES_Training_CR < 0)
n_rev_tl <- sum(fgsea_sig$NES_Cancer_vs_Healthy < 0 & fgsea_sig$NES_Training_CR > 0)
n_exac_tr <- sum(fgsea_sig$NES_Cancer_vs_Healthy > 0 & fgsea_sig$NES_Training_CR > 0)
n_exac_bl <- sum(fgsea_sig$NES_Cancer_vs_Healthy < 0 & fgsea_sig$NES_Training_CR < 0)

n_reversed <- n_rev_br + n_rev_tl
n_total_sig <- nrow(fgsea_sig)
pw_rev_frac <- if (n_total_sig > 0) n_reversed / n_total_sig else 0

txt_gene <- scale_text(BASE_GENE, PG_W)
txt_quad <- scale_text(BASE_QUADRANT, PG_W)

label_pw <- fgsea_sig %>%
  mutate(
    label_fill     = SIG_LABEL_FILL_F4[as.character(significance)],
    label_text_col = SIG_LABEL_TEXT_F4[as.character(significance)]
  )

ns_df  <- fgsea_wide %>% filter(significance == "NS")
sig_df <- fgsea_wide %>% filter(significance != "NS") %>%
  mutate(
    border_col = ifelse(database == "Hallmark", "black", "grey75"),
    bubble_alpha = case_when(
      significance == "Sig Both"          ~ 0.75,
      significance == "Sig Cancer only"   ~ 0.85,
      significance == "Sig Training only" ~ 0.85,
      TRUE ~ 0.60
    )
  )

subtitle_str <- sprintf(
  "%d pathways (%d sig.) | fGSEA on limma t-statistics\n\u03c1(all) = %.2f [%.2f, %.2f], p %s | %.0f%% reversed",
  nrow(fgsea_wide), n_total_sig,
  nes_cor_all$estimate, nes_ci_all[1], nes_ci_all[2],
  ifelse(nes_cor_all$p.value < 0.001, "< 0.001", sprintf("= %.3f", nes_cor_all$p.value)),
  pw_rev_frac * 100
)

pG <- ggplot(mapping = aes(x = NES_Cancer_vs_Healthy, y = NES_Training_CR)) +
  # Reversed quadrants (blue) — opposite signs
  annotate("rect", xmin = 0, xmax = Inf,  ymin = -Inf, ymax = 0,
           fill = "#DCEEFF", alpha = 0.55) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#DCEEFF", alpha = 0.55) +
  # Exacerbated quadrants (red) — same signs
  annotate("rect", xmin = 0, xmax = Inf,  ymin = 0, ymax = Inf,
           fill = "#FFE0E0", alpha = 0.55) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#FFE0E0", alpha = 0.55) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  geom_point(data = ns_df, aes(x = NES_Cancer_vs_Healthy, y = NES_Training_CR),
             size = 1.5, fill = "grey70", shape = 21, color = "grey55",
             alpha = 0.40, stroke = 0.4) +
  geom_point(data = sig_df, aes(fill = significance, size = set_size),
             shape = 21, color = sig_df$border_col,
             alpha = sig_df$bubble_alpha, stroke = 0.8) +
  scale_fill_manual(values = SIG_COLORS_F4, name = "Significance") +
  scale_size_continuous(range = c(2, 8), name = "Set size",
                        breaks = c(20, 50, 100, 200)) +
  geom_label_repel(data = label_pw, aes(label = pathway_label),
                   fill = label_pw$label_fill, color = label_pw$label_text_col,
                   size = txt_gene, fontface = "bold",
                   max.overlaps = 40,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.5, force = 3, force_pull = 0.5,
                   label.padding = unit(1.5, "pt"),
                   label.r = unit(1, "pt"),
                   label.size = 0.15, seed = 42) +
  annotate("label", x = Inf, y = -Inf,
           label = sprintf("Reversed  n = %d", n_rev_br),
           hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = Inf,
           label = sprintf("Reversed  n = %d", n_rev_tl),
           hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("Exacerbated  n = %d", n_exac_tr),
           hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = -Inf,
           label = sprintf("Exacerbated  n = %d", n_exac_bl),
           hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  scale_x_continuous(expand = expansion(0, 0)) +
  scale_y_continuous(expand = expansion(0, 0)) +
  coord_cartesian(xlim = c(-nes_lim, nes_lim), ylim = c(-nes_lim, nes_lim)) +
  labs(
    title = "Pathway-Level Reversal (fGSEA)",
    subtitle = subtitle_str,
    x = "NES (Cancer vs Healthy)",
    y = "NES (Training CR)"
  ) +
  FIG_THEME +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 8, face = "bold"),
    legend.text     = element_text(size = 7),
    legend.key.size = unit(3, "mm"),
    legend.margin   = margin(0, 0, 0, 0)
  ) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 3, alpha = 0.8)),
         size = guide_legend(nrow = 1))

ggsave(file.path(RPT, "panel_D_nes_reversal.pdf"), pG,
       width = PG_W, height = PG_W, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_D_nes_reversal.png"), pG,
       width = PG_W, height = PG_W, units = "mm", dpi = 300)

fgsea_wide %>%
  transmute(pathway, pathway_label, database,
            NES_Cancer_vs_Healthy = round(NES_Cancer_vs_Healthy, 3),
            NES_Training_CR       = round(NES_Training_CR, 3),
            padj_Cancer_vs_Healthy = signif(padj_Cancer_vs_Healthy, 4),
            padj_Training_CR       = signif(padj_Training_CR, 4),
            significance           = as.character(significance),
            set_size) %>%
  arrange(significance, desc(abs(NES_Cancer_vs_Healthy) + abs(NES_Training_CR))) %>%
  write_csv(file.path(DAT, "panel_D", "nes_reversal.csv"))

message("F04 Panel D done")
