# F05 CRvH Panel B ORA: Scatter + Flanking ORA Bars Composite
# Threshold-free ORA on all proteins per concordance quadrant.
# Concordance framing: Cancer_vs_Healthy (x) vs Training_CR (y).
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F05/a_script/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
  library(ggrepel)
  library(patchwork)
})

RPT <- "04_Figures/F05/CRvH/b_reports"
DAT <- "04_Figures/F05/CRvH/c_data"
dir.create(file.path(DAT, "panel_B"), recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

COMP_RED   <- unname(DIR_COLORS["Up"])     # #D6604D
COMP_BLUE  <- unname(DIR_COLORS["Down"])   # #4393C3
COMP_RED_LIGHT  <- scales::alpha(COMP_RED, 0.30)
COMP_BLUE_LIGHT <- scales::alpha(COMP_BLUE, 0.30)
N_SHOW    <- 5

# -- Data --
dep_df <- read_csv("03_DEP/c_data/03_combined_results_CRvH.csv", show_col_types = FALSE)
imputation_df <- read_csv("02_Imputation/c_data/02_mar_mnar_classification.csv",
                           show_col_types = FALSE) %>%
  transmute(gene, imputed = classification != "Complete")

scatter_df <- dep_df %>%
  transmute(gene,
            logFC_CvH = logFC_Cancer_vs_Healthy, logFC_TR = logFC_Training_CR,
            pi_CvH = pi_score_Cancer_vs_Healthy, pi_TR = pi_score_Training_CR) %>%
  filter(!is.na(logFC_CvH), !is.na(logFC_TR)) %>%
  left_join(imputation_df, by = "gene") %>%
  mutate(
    imputed      = replace_na(imputed, FALSE),
    sig_class    = classify_proteins_f4(pi_CvH, pi_TR),
    is_sig       = sig_class != "NS",
    quadrant     = case_when(
      logFC_CvH > 0 & logFC_TR > 0 ~ "Concordant Up",
      logFC_CvH < 0 & logFC_TR < 0 ~ "Concordant Down",
      logFC_CvH > 0 & logFC_TR < 0 ~ "Discordant (Cancer Up / Training Down)",
      TRUE                          ~ "Discordant (Cancer Down / Training Up)"),
    border_col   = "grey75",
    point_size   = ifelse(sig_class == "NS", 1.8, 2.3),
    point_stroke = ifelse(sig_class == "NS", 0.6, 0.9),
    bubble_alpha = case_when(
      sig_class == "NS"       ~ 0.30,
      sig_class == "Sig Both" ~ 0.75,
      TRUE                    ~ 0.85))

universe <- scatter_df$gene
message(sprintf("  Total proteins: %d | Significant: %d",
                nrow(scatter_df), sum(scatter_df$is_sig)))

# -- ORA --
pw_collection <- build_pathway_collection(min_size = 15, max_size = 500,
                                           include_goslim = FALSE,
                                           exclude_variants = TRUE)

run_set_ora <- function(genes, set_name) {
  if (length(genes) < 5) return(tibble())
  res <- tryCatch(
    run_ora_deduplicated(genes = genes, universe = universe,
                          pathways = pw_collection, jaccard_cutoff = 0.5,
                          min_size = 15, max_size = 500, padj_cutoff = 1),
    error = function(e) { message("  ORA error: ", e$message); tibble() })
  if (nrow(res) == 0) return(tibble())
  res %>%
    mutate(set = set_name, pathway_label = clean_pathway_name(pathway),
           neg_log10_padj = -log10(padj), significant = padj < 0.05) %>%
    arrange(desc(neg_log10_padj)) %>%
    slice_head(n = N_SHOW)
}

message("\n--- Quadrant ORA (threshold-free) ---")
ora_q1 <- run_set_ora(scatter_df$gene[scatter_df$quadrant == "Concordant Up"],
                       "Concordant Up")
ora_q2 <- run_set_ora(scatter_df$gene[scatter_df$quadrant == "Discordant (Cancer Down / Training Up)"],
                       "Discordant (Cancer Down / Training Up)")
ora_q3 <- run_set_ora(scatter_df$gene[scatter_df$quadrant == "Concordant Down"],
                       "Concordant Down")
ora_q4 <- run_set_ora(scatter_df$gene[scatter_df$quadrant == "Discordant (Cancer Up / Training Down)"],
                       "Discordant (Cancer Up / Training Down)")

for (qn in c("Concordant Up", "Concordant Down",
             "Discordant (Cancer Down / Training Up)",
             "Discordant (Cancer Up / Training Down)"))
  message(sprintf("  %s: %d genes", qn,
                  sum(scatter_df$quadrant == qn)))

all_quad_ora <- bind_rows(ora_q1, ora_q2, ora_q3, ora_q4)
if (nrow(all_quad_ora) > 0)
  write_csv(all_quad_ora, file.path(DAT, "panel_B", "ora_quadrant.csv"))

# -- Scatter panel (composite version) --
xlim_range <- c(-3, 3)
ylim_range <- c(-2, 2)

ns_df  <- filter(scatter_df, sig_class == "NS")
sig_df <- filter(scatter_df, sig_class != "NS")

q_df <- scatter_df %>%
  mutate(q = case_when(
    logFC_CvH > 0 & logFC_TR > 0 ~ "Q1",
    logFC_CvH < 0 & logFC_TR < 0 ~ "Q3",
    logFC_CvH > 0 & logFC_TR < 0 ~ "Q4",
    TRUE ~ "Q2"))
q_counts <- q_df %>% count(q) %>% deframe()
q_sig    <- q_df %>% filter(sig_class != "NS") %>% count(q) %>% deframe()
for (qq in c("Q1","Q2","Q3","Q4")) if (is.na(q_sig[qq])) q_sig[qq] <- 0

label_df <- sig_df %>%
  group_by(sig_class) %>%
  arrange(desc(abs(logFC_CvH) + abs(logFC_TR))) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(label_fill     = SIG_LABEL_FILL_F4[as.character(sig_class)],
         label_text_col = SIG_LABEL_TEXT_F4[as.character(sig_class)])

txt_gene <- scale_text(BASE_GENE, 200) * 1.25
txt_quad <- scale_text(BASE_QUADRANT, 200) * 1.3

p_scatter <- ggplot(mapping = aes(x = logFC_CvH, y = logFC_TR)) +
  annotate("rect", xmin = 0, xmax = Inf,  ymin = 0, ymax = Inf,
           fill = "#FFE0E0", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#FFE0E0", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 0, xmax = Inf,  ymin = -Inf, ymax = 0,
           fill = "#DCEEFF", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#DCEEFF", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  geom_point(data = ns_df, color = "grey80", fill = "grey85", shape = 21,
             size = 1.0, alpha = 0.3, stroke = 0.2) +
  geom_point(data = sig_df, aes(fill = sig_class), shape = 21,
             size = sig_df$point_size, color = sig_df$border_col,
             alpha = sig_df$bubble_alpha, stroke = sig_df$point_stroke) +
  scale_fill_manual(values = SIG_COLORS_F4, name = "Significance") +
  geom_label_repel(data = label_df, aes(label = gene),
                   fill = label_df$label_fill, color = label_df$label_text_col,
                   size = txt_gene, fontface = "italic", max.overlaps = 40,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.6, force = 3, force_pull = 0.5,
                   label.padding = unit(2, "pt"), label.r = unit(1, "pt"),
                   label.size = 0.15, seed = 42,
                   xlim = xlim_range * 0.9, ylim = ylim_range * 0.9) +
  annotate("label", x = xlim_range[2] / 2, y = ylim_range[2] / 2,
           label = sprintf("Concordant Up\nn = %s/%s", q_sig["Q1"], q_counts["Q1"]),
           hjust = 0.5, vjust = 0.5, size = txt_quad, fontface = "bold",
           color = COMP_RED, fill = alpha("white", 0.9),
           label.padding = unit(3, "pt"), lineheight = 0.85) +
  annotate("label", x = xlim_range[1] / 2, y = ylim_range[1] / 2,
           label = sprintf("Concordant Down\nn = %s/%s", q_sig["Q3"], q_counts["Q3"]),
           hjust = 0.5, vjust = 0.5, size = txt_quad, fontface = "bold",
           color = COMP_RED, fill = alpha("white", 0.9),
           label.padding = unit(3, "pt"), lineheight = 0.85) +
  annotate("label", x = xlim_range[1] / 2, y = ylim_range[2] / 2,
           label = sprintf("Discordant (C\u2193 T\u2191)\nn = %s/%s", q_sig["Q2"], q_counts["Q2"]),
           hjust = 0.5, vjust = 0.5, size = txt_quad, fontface = "bold",
           color = COMP_BLUE, fill = alpha("white", 0.9),
           label.padding = unit(3, "pt"), lineheight = 0.85) +
  annotate("label", x = xlim_range[2] / 2, y = ylim_range[1] / 2,
           label = sprintf("Discordant (C\u2191 T\u2193)\nn = %s/%s", q_sig["Q4"], q_counts["Q4"]),
           hjust = 0.5, vjust = 0.5, size = txt_quad, fontface = "bold",
           color = COMP_BLUE, fill = alpha("white", 0.9),
           label.padding = unit(3, "pt"), lineheight = 0.85) +
  annotate("label", x = xlim_range[2] * 0.35, y = 0,
           label = expression(log[2]*FC ~ "(Cancer vs Healthy)"),
           hjust = 0.5, vjust = 0.5, size = 3.5, color = "grey30", fontface = "bold",
           fill = alpha("white", 0.8), label.size = 0, label.padding = unit(2, "pt")) +
  annotate("label", x = 0, y = ylim_range[2] * 0.75,
           label = expression(log[2]*FC ~ "(Training CR)"),
           hjust = 0.5, vjust = 0.5, size = 3.5, color = "grey30", fontface = "bold",
           fill = alpha("white", 0.8), label.size = 0, label.padding = unit(2, "pt"),
           angle = 90) +
  scale_x_continuous(breaks = seq(-3, 3, 1)) +
  scale_y_continuous(breaks = seq(-2, 2, 1)) +
  coord_cartesian(xlim = xlim_range, ylim = ylim_range, expand = FALSE) +
  labs(x = NULL, y = NULL) +
  FIG_THEME +
  theme(plot.title    = element_blank(),
        plot.subtitle = element_blank(),
        axis.text.x  = element_text(size = 8, color = "grey40"),
        axis.text.y  = element_text(size = 8, color = "grey40"),
        axis.ticks   = element_line(color = "grey40", linewidth = 0.3),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin  = margin(2, 0, 2, 0, "mm"),
        legend.position    = "bottom",
        legend.title       = element_text(size = 10, face = "bold"),
        legend.text        = element_text(size = 9),
        legend.key.size    = unit(4, "mm"),
        legend.margin      = margin(2, 0, 0, 0),
        legend.box.margin  = margin(0, 0, 5, 0)) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 4, alpha = 0.8)))

# -- Half-bar builder --
make_half_bars <- function(df, fill_color, side, ylim,
                            display_labels = character(0)) {
  bar_h  <- 0.24
  n_bars <- if (is.null(df) || nrow(df) == 0) 0L else min(nrow(df), 5L)

  if (n_bars == 0)
    return(ggplot() + theme_void() +
             scale_y_continuous(limits = ylim, expand = c(0, 0)))

  if (ylim[1] >= 0) {
    y_pos <- rev(seq(0.3, 1.7, length.out = 5))[seq_len(n_bars)]
  } else {
    y_pos <- seq(-0.3, -1.7, length.out = 5)[seq_len(n_bars)]
  }

  bars <- df %>%
    arrange(desc(neg_log10_padj)) %>%
    slice_head(n = 5) %>%
    mutate(
      y            = y_pos,
      bar_fill     = ifelse(significant, scales::alpha(fill_color, 0.85),
                            scales::alpha(fill_color, 0.30)),
      display_name = ifelse(pathway_label %in% names(display_labels),
                            display_labels[pathway_label],
                            pathway_label),
      display_name = ifelse(!grepl("\n", display_name),
                            stringr::str_wrap(display_name, width = 22),
                            display_name),
      label_size   = ifelse(significant, 4.5, 3.4),
      star         = sig_stars(padj))

  x_max  <- max(bars$neg_log10_padj)
  is_upper <- ylim[1] >= 0

  p <- ggplot(bars, aes(y = y)) +
    geom_rect(aes(xmin = 0, xmax = neg_log10_padj,
                  ymin = y - bar_h / 2, ymax = y + bar_h / 2),
              fill = bars$bar_fill, color = "black", linewidth = 0.3) +
    geom_text(data = bars %>% filter(significant),
              aes(x = neg_log10_padj / 2, y = y, label = display_name),
              hjust = 0.5, size = 4.5, fontface = "bold",
              color = "white", lineheight = 0.85) +
    geom_text(data = bars %>% filter(!significant),
              aes(x = neg_log10_padj / 2, y = y, label = display_name),
              hjust = 0.5, size = 3.4, fontface = "bold",
              color = "grey30", lineheight = 0.85) +
    geom_text(aes(x = neg_log10_padj + x_max * 0.05, label = star),
              hjust = 0, size = 3.5, fontface = "bold", color = "black") +
    annotate("segment", x = 0, xend = x_max, y = -Inf, yend = -Inf,
             color = "grey40", linewidth = 0.3) +
    scale_fill_identity() +
    labs(x = if (!is_upper) expression(-log[10](p[adj])) else NULL,
         y = NULL) +
    theme_minimal(base_size = 8) +
    theme(panel.grid   = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x  = if (is_upper) element_blank()
                          else element_text(size = 7, face = "bold"),
          axis.title.x = if (is_upper) element_blank()
                          else element_text(size = 8, face = "bold"),
          axis.line.x  = element_blank(),
          axis.ticks.x = if (is_upper) element_blank()
                          else element_line(color = "grey40", linewidth = 0.3),
          plot.margin  = if (side == "left") margin(2, 0, 2, 3, "mm")
                         else margin(2, 3, 2, 0, "mm"))

  if (side == "left") {
    p + scale_x_reverse(limits = c(x_max * 1.25, 0),
                         expand = expansion(mult = c(0.10, 0))) +
        scale_y_continuous(limits = ylim, expand = c(0, 0))
  } else {
    p + scale_x_continuous(limits = c(0, x_max * 1.25),
                            expand = expansion(mult = c(0, 0.10))) +
        scale_y_continuous(limits = ylim, expand = c(0, 0))
  }
}

# Left:  upper = Q2 (Discordant C Down/T Up, BLUE), lower = Q3 (Concordant Down, RED)
# Right: upper = Q1 (Concordant Up, RED),            lower = Q4 (Discordant C Up/T Down, BLUE)
p_ul <- make_half_bars(ora_q2, COMP_BLUE_LIGHT, "left",  c(0, 2))
p_ll <- make_half_bars(ora_q3, COMP_RED_LIGHT,  "left",  c(-2, 0))
p_ur <- make_half_bars(ora_q1, COMP_RED_LIGHT,  "right", c(0, 2))
p_lr <- make_half_bars(ora_q4, COMP_BLUE_LIGHT, "right", c(-2, 0))

# -- Composite --
design <- c(
  area(1, 1),          # upper-left bars
  area(1, 2, 2, 2),   # scatter (spans both rows)
  area(1, 3),          # upper-right bars
  area(2, 1),          # lower-left bars
  area(2, 3)           # lower-right bars
)
n_total   <- nrow(scatter_df)
n_sig     <- sum(scatter_df$is_sig)
n_enrich  <- if (nrow(all_quad_ora) > 0) sum(all_quad_ora$significant) else 0L
r_pearson <- cor(scatter_df$logFC_CvH, scatter_df$logFC_TR, use = "complete.obs")

composite <- p_ul + p_scatter + p_ur + p_ll + p_lr +
  plot_layout(design = design, widths = c(1.4, 2, 1.4)) +
  plot_annotation(
    title    = "CRvH Concordance: Quadrant ORA",
    subtitle = sprintf("Threshold-free quadrant ORA (hypergeometric) | N = %d proteins | %d DEPs (Pi < 0.05) | %d enriched terms (FDR < 0.05) | r = %.2f | Hallmark + Reactome + KEGG + GO:BP",
                        n_total, n_sig, n_enrich, r_pearson),
    theme    = theme(plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
                     plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey30")))

COMP_W <- 450; COMP_H <- 220
ggsave(file.path(RPT, "panel_B_ORA_composite_MAIN.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_B_ORA_composite_MAIN.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message("\nF05 CRvH Panel B ORA composite done")
