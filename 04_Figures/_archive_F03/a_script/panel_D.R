# F03 Panel D: fGSEA NES Scatter (Training CRE vs Training PLA)
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

PG_W <- 200
RPT <- "04_Figures/F03/b_reports"
DAT <- "04_Figures/F03/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_D"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

fgsea_all <- read_csv("04_Figures/F02/c_data/06_panel_H_CR_fgsea_results.csv",
                       show_col_types = FALSE)

fgsea_hg <- fgsea_all %>%
  filter(contrast %in% c("Training_CRE", "Training_PLA"))

fgsea_wide <- fgsea_hg %>%
  dplyr::select(pathway, contrast, NES, padj, size, database) %>%
  pivot_wider(id_cols = c(pathway, database), names_from = contrast,
              values_from = c(NES, padj, size)) %>%
  filter(!is.na(NES_Training_CRE), !is.na(NES_Training_PLA)) %>%
  mutate(set_size = coalesce(size_Training_CRE, size_Training_PLA))

fgsea_wide <- fgsea_wide %>%
  mutate(
    sig_CRE = !is.na(padj_Training_CRE) & padj_Training_CRE < 0.05,
    sig_PLA = !is.na(padj_Training_PLA) & padj_Training_PLA < 0.05,
    significance = case_when(
      sig_CRE & sig_PLA ~ "Sig Both",
      sig_CRE           ~ "Sig CRE only",
      sig_PLA           ~ "Sig PLA only",
      TRUE              ~ "NS"
    ) %>% factor(levels = c("Sig Both", "Sig CRE only", "Sig PLA only", "NS")),
    pathway_label = clean_pathway_name(pathway)
  )

fgsea_sig <- fgsea_wide %>% filter(significance != "NS")

message(sprintf("  %d total pathways | %d significant", nrow(fgsea_wide), nrow(fgsea_sig)))

nes_cor_all <- cor.test(fgsea_wide$NES_Training_CRE, fgsea_wide$NES_Training_PLA,
                         method = "spearman")
nes_ci_all  <- fisher_z_ci(nes_cor_all$estimate, nrow(fgsea_wide))
nes_cor_sig <- if (nrow(fgsea_sig) >= 3) {
  cor.test(fgsea_sig$NES_Training_CRE, fgsea_sig$NES_Training_PLA, method = "spearman")
} else NULL

nes_lim <- max(abs(c(fgsea_wide$NES_Training_CRE, fgsea_wide$NES_Training_PLA))) * 1.15

n_conc_tr <- sum(fgsea_sig$NES_Training_CRE > 0 & fgsea_sig$NES_Training_PLA > 0)
n_conc_bl <- sum(fgsea_sig$NES_Training_CRE < 0 & fgsea_sig$NES_Training_PLA < 0)
n_disc_q2 <- sum(fgsea_sig$NES_Training_CRE < 0 & fgsea_sig$NES_Training_PLA > 0)
n_disc_q4 <- sum(fgsea_sig$NES_Training_CRE > 0 & fgsea_sig$NES_Training_PLA < 0)

n_conc_pw   <- n_conc_tr + n_conc_bl
n_total_sig <- nrow(fgsea_sig)
pw_conc_frac <- if (n_total_sig > 0) n_conc_pw / n_total_sig else 0

NES_SIG_COLORS <- c(
  "Sig Both"     = "#2E7D32",
  "Sig CRE only" = "#2166AC",
  "Sig PLA only" = "#D6604D",
  "NS"           = "grey70"
)

NES_SIG_LABEL_FILL <- c(
  "Sig Both"     = scales::alpha("#2E7D32", 0.75),
  "Sig CRE only" = scales::alpha("#2166AC", 0.75),
  "Sig PLA only" = scales::alpha("#D6604D", 0.75),
  "NS"           = scales::alpha("grey70",  0.75)
)
NES_SIG_LABEL_TEXT <- c(
  "Sig Both" = "white", "Sig CRE only" = "white",
  "Sig PLA only" = "white", "NS" = "white"
)

txt_gene <- scale_text(BASE_GENE, PG_W)
txt_quad <- scale_text(BASE_QUADRANT, PG_W)

label_pw <- fgsea_sig %>%
  mutate(
    label_fill     = NES_SIG_LABEL_FILL[as.character(significance)],
    label_text_col = NES_SIG_LABEL_TEXT[as.character(significance)]
  )

ns_df  <- fgsea_wide %>% filter(significance == "NS")
sig_df <- fgsea_wide %>% filter(significance != "NS") %>%
  mutate(
    border_col = ifelse(database == "Hallmark", "black", "grey75"),
    bubble_alpha = case_when(
      significance == "Sig Both"     ~ 0.75,
      significance == "Sig CRE only" ~ 0.85,
      significance == "Sig PLA only" ~ 0.85,
      TRUE ~ 0.60
    )
  )

rho_sig_str <- if (!is.null(nes_cor_sig)) sprintf(", rho(sig) = %.2f", nes_cor_sig$estimate) else ""
subtitle_str <- sprintf(
  "%d pathways (%d sig.) | fGSEA on limma t-statistics\n\u03c1(all) = %.2f [%.2f, %.2f], p %s%s | %.0f%% concordant",
  nrow(fgsea_wide), n_total_sig,
  nes_cor_all$estimate, nes_ci_all[1], nes_ci_all[2],
  ifelse(nes_cor_all$p.value < 0.001, "< 0.001", sprintf("= %.3f", nes_cor_all$p.value)),
  rho_sig_str, pw_conc_frac * 100
)

pG <- ggplot(mapping = aes(x = NES_Training_CRE, y = NES_Training_PLA)) +
  annotate("rect", xmin = 0, xmax = Inf,  ymin = 0, ymax = Inf,
           fill = "#FFE0E0", alpha = 0.55) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#FFE0E0", alpha = 0.55) +
  annotate("rect", xmin = 0, xmax = Inf,  ymin = -Inf, ymax = 0,
           fill = "#DCEEFF", alpha = 0.55) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#DCEEFF", alpha = 0.55) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  geom_point(data = ns_df, aes(x = NES_Training_CRE, y = NES_Training_PLA),
             size = 1.5, fill = "grey70", shape = 21, color = "grey55",
             alpha = 0.40, stroke = 0.4) +
  geom_point(data = sig_df, aes(fill = significance, size = set_size),
             shape = 21, color = sig_df$border_col,
             alpha = sig_df$bubble_alpha, stroke = 0.8) +
  scale_fill_manual(values = NES_SIG_COLORS, name = "Significance") +
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
  annotate("label", x = Inf, y = Inf,
           label = sprintf("Concordant Up  n = %d", n_conc_tr),
           hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = -Inf,
           label = sprintf("Concordant Down  n = %d", n_conc_bl),
           hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = Inf,
           label = sprintf("Discordant  n = %d", n_disc_q2),
           hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = Inf, y = -Inf,
           label = sprintf("Discordant  n = %d", n_disc_q4),
           hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  scale_x_continuous(expand = expansion(0, 0)) +
  scale_y_continuous(expand = expansion(0, 0)) +
  coord_cartesian(xlim = c(-nes_lim, nes_lim), ylim = c(-nes_lim, nes_lim)) +
  labs(
    title = "Pathway-Level Concordance (fGSEA)",
    subtitle = subtitle_str,
    x = "NES (Training CRE)",
    y = "NES (Training PLA)"
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

ggsave(file.path(RPT, "panel_D_nes_scatter.pdf"), pG,
       width = PG_W, height = PG_W, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_D_nes_scatter.png"), pG,
       width = PG_W, height = PG_W, units = "mm", dpi = 300)

fgsea_wide %>%
  transmute(
    pathway, pathway_label, database,
    NES_Training_CRE = round(NES_Training_CRE, 3),
    NES_Training_PLA = round(NES_Training_PLA, 3),
    padj_Training_CRE = signif(padj_Training_CRE, 4),
    padj_Training_PLA = signif(padj_Training_PLA, 4),
    significance      = as.character(significance),
    set_size
  ) %>%
  arrange(significance, desc(abs(NES_Training_CRE) + abs(NES_Training_PLA))) %>%
  write_csv(file.path(DAT, "panel_D", "nes_scatter.csv"))

message("F03 Panel D done")
