# F04 CRvH Supplementary: Volcano + Flanking fGSEA Pathway Bars
# Minimal volcanos with top fGSEA pathways as bars growing outward
# Down-regulated pathways (blue bars) grow right->left from volcano left edge
# Up-regulated pathways (red bars) grow left->right from volcano right edge
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F04/a_script/style.R")

library(tidyverse)
library(patchwork)

RPT <- "04_Figures/F04/CRvH/b_reports/supp"
DAT <- "04_Figures/F04/CRvH/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

# --- Load data
dep <- readr::read_csv("03_DEP/a_non_imputed/c_data/combined_results_pi.csv",
                       show_col_types = FALSE) |>
  dplyr::mutate(contrast = dplyr::recode(contrast,
                                         CRvH_Baseline = "Cancer_vs_Healthy",
                                         CR_Training   = "Training_CR")) |>
  tidyr::pivot_wider(id_cols = c(uniprot_id, gene, protein, description),
                     names_from = contrast,
                     values_from = c(logFC, t, P.Value, adj.P.Val, pi_score, sig_pi),
                     names_glue = "{.value}_{contrast}")
fgsea_all <- read_csv(file.path(DAT, "shared", "fgsea_tstat_CRvH.csv"),
                       show_col_types = FALSE)

CONTRASTS  <- c("Cancer_vs_Healthy", "Training_CR")
CTR_LABELS <- c("Cancer vs Healthy", "Training (CR)")
N_PATHWAYS <- 7  # top pathways per direction

PW  <- 180  # panel width mm
txt <- scale_text(3.8, PW)

# --- Build one row: left bars | volcano | right bars
build_volcano_row <- function(ctr, ctr_label, show_xlab = FALSE) {

  # --- Volcano data
  lfc_col  <- paste0("logFC_", ctr)
  pval_col <- paste0("P.Value_", ctr)
  pi_col   <- paste0("pi_score_", ctr)

  vdf <- dep %>%
    transmute(gene,
              logFC     = .data[[lfc_col]],
              pval      = .data[[pval_col]],
              pi        = .data[[pi_col]],
              neg_log10p = -log10(pval),
              sig       = !is.na(pi) & pi < 0.05,
              direction = case_when(
                !sig       ~ "NS",
                logFC > 0  ~ "Up",
                TRUE       ~ "Down"
              )) %>%
    filter(!is.na(logFC), !is.na(pval))

  # --- Minimal volcano
  p_volc <- ggplot(vdf, aes(logFC, neg_log10p, color = direction)) +
    geom_point(size = 0.4, alpha = 0.5) +
    scale_color_manual(values = DIR_COLORS, guide = "none") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey60", linewidth = 0.3) +
    labs(x     = if (show_xlab) expression(log[2]~FC) else NULL,
         y     = if (show_xlab) expression(-log[10]~p) else NULL,
         title = ctr_label) +
    FIG_THEME +
    theme(
      plot.title    = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.text     = if (show_xlab)
                        element_text(size = 8, color = "grey15")
                      else element_blank(),
      axis.title    = if (show_xlab)
                        element_text(size = 9, face = "bold")
                      else element_blank(),
      axis.ticks    = if (show_xlab) element_line() else element_blank(),
      panel.border  = element_rect(color = "grey70", linewidth = 0.4, fill = NA),
      plot.margin   = margin(2, 0, 2, 0)
    )

  # --- fGSEA bars
  fgsea_ctr <- fgsea_all %>%
    filter(contrast == ctr, padj < 0.05) %>%
    mutate(direction     = ifelse(NES > 0, "Up", "Down"),
           pathway_label = str_wrap(clean_pathway_name(pathway), 25),
           abs_NES       = abs(NES))

  up_df <- fgsea_ctr %>%
    filter(direction == "Up") %>%
    arrange(desc(abs_NES)) %>%
    slice_head(n = N_PATHWAYS) %>%
    mutate(y = row_number())

  dn_df <- fgsea_ctr %>%
    filter(direction == "Down") %>%
    arrange(desc(abs_NES)) %>%
    slice_head(n = N_PATHWAYS) %>%
    mutate(y = row_number())

  # Shared axis scale
  max_nes <- max(c(up_df$abs_NES, dn_df$abs_NES, 2), na.rm = TRUE) * 1.1
  n_rows  <- max(nrow(up_df), nrow(dn_df), 1)

  # Stars helper
  sig_star <- function(p) {
    ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
  }

  # --- Right bars (Up, red)
  if (nrow(up_df) > 0) {
    p_right <- ggplot(up_df, aes(y = y)) +
      geom_rect(aes(xmin = 0, xmax = abs_NES,
                     ymin = y - 0.4, ymax = y + 0.4),
                fill = DIR_COLORS["Up"], color = NA) +
      geom_text(aes(x = 0.1, y = y, label = pathway_label),
                hjust = 0, size = txt * 0.7, color = "white", fontface = "bold",
                lineheight = 0.85) +
      geom_text(aes(x = abs_NES + 0.05, y = y, label = sig_star(padj)),
                hjust = 0, size = txt * 0.8, fontface = "bold", color = "grey30") +
      scale_x_continuous(limits = c(0, max_nes), expand = c(0, 0)) +
      scale_y_reverse(limits = c(n_rows + 0.5, 0.5)) +
      theme_void() +
      theme(plot.margin = margin(2, 4, 2, 0))
  } else {
    p_right <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No sig.\npathways",
               size = txt, color = "grey50", fontface = "italic") +
      theme_void() +
      theme(plot.margin = margin(2, 4, 2, 0))
  }

  # --- Left bars (Down, blue) -- grow right to left
  if (nrow(dn_df) > 0) {
    p_left <- ggplot(dn_df, aes(y = y)) +
      geom_rect(aes(xmin = max_nes - abs_NES, xmax = max_nes,
                     ymin = y - 0.4, ymax = y + 0.4),
                fill = DIR_COLORS["Down"], color = NA) +
      geom_text(aes(x = max_nes - 0.1, y = y, label = pathway_label),
                hjust = 1, size = txt * 0.7, color = "white", fontface = "bold",
                lineheight = 0.85) +
      geom_text(aes(x = max_nes - abs_NES - 0.05, y = y, label = sig_star(padj)),
                hjust = 1, size = txt * 0.8, fontface = "bold", color = "grey30") +
      scale_x_continuous(limits = c(0, max_nes), expand = c(0, 0)) +
      scale_y_reverse(limits = c(n_rows + 0.5, 0.5)) +
      theme_void() +
      theme(plot.margin = margin(2, 0, 2, 4))
  } else {
    p_left <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No sig.\npathways",
               size = txt, color = "grey50", fontface = "italic") +
      theme_void() +
      theme(plot.margin = margin(2, 0, 2, 4))
  }

  # Combine: left bars | volcano | right bars
  p_left + p_volc + p_right +
    plot_layout(widths = c(1.3, 1, 1.3), nrow = 1)
}

# --- Build 2 rows (CvH and Training_CR -- no Interaction in CRvH model)
rows <- lapply(seq_along(CONTRASTS), function(i) {
  build_volcano_row(CONTRASTS[i], CTR_LABELS[i],
                    show_xlab = (i == length(CONTRASTS)))
})

# Stack vertically
pilot <- wrap_plots(rows, ncol = 1) +
  plot_annotation(
    title    = "Volcano + fGSEA Pathway Bars (Cancer Recovery Concordance)",
    subtitle = "Top 7 pathways per direction (padj < 0.05) | bar length = |NES|",
    theme = theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(face = "bold.italic", size = 9,
                                   color = "grey30")
    )
  )

ggsave(file.path(RPT, "panel_A_pilot_volcano_fgsea_SUPP.pdf"), pilot,
       width = PW, height = 300, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "panel_A_pilot_volcano_fgsea_SUPP.png"), pilot,
       width = PW, height = 300, units = "mm", dpi = 300,
       limitsize = FALSE)

message("F04 CRvH Panel A pilot (volcano + fGSEA bars) done")
