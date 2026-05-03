# F05/CR Panel B ORA: Scatter + Flanking ORA Bars Composite
# Threshold-free ORA on all proteins per quadrant, displayed as bars flanking
# the concordance scatter. Red = same-direction (concordant), Blue = opposing.
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F05/a_script/style.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
  library(ggrepel)
  library(patchwork)
})

RPT <- "04_Figures/F05/CR/b_reports"
DAT <- "04_Figures/F05/CR/c_data"
dir.create(file.path(DAT, "panel_B"), recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

COMP_RED   <- unname(DIR_COLORS["Up"])     # #D6604D
COMP_BLUE  <- unname(DIR_COLORS["Down"])   # #4393C3
COMP_RED_LIGHT  <- scales::alpha(COMP_RED, 0.30)
COMP_BLUE_LIGHT <- scales::alpha(COMP_BLUE, 0.30)
N_SHOW    <- 5

# -- Data ----------------------------------------------------------------------
dep_df <- read_csv("03_DEP/c_data/03_combined_results_CR.csv", show_col_types = FALSE)
imputation_df <- read_csv("02_Imputation/c_data/02_mar_mnar_classification.csv",
                           show_col_types = FALSE) %>%
  transmute(gene, imputed = classification != "Complete")

scatter_df <- dep_df %>%
  transmute(gene,
            logFC_CRE = logFC_Training_CRE, logFC_PLA = logFC_Training_PLA,
            pi_CRE = pi_score_Training_CRE, pi_PLA = pi_score_Training_PLA,
            pi_Int = pi_score_Supplement_Interaction) %>%
  filter(!is.na(logFC_CRE), !is.na(logFC_PLA)) %>%
  left_join(imputation_df, by = "gene") %>%
  mutate(
    imputed      = replace_na(imputed, FALSE),
    sig_class    = classify_proteins_f3(pi_CRE, pi_PLA, pi_Int),
    is_sig       = sig_class != "NS",
    quadrant     = case_when(
      logFC_CRE > 0 & logFC_PLA > 0 ~ "Concordant Up",
      logFC_CRE < 0 & logFC_PLA < 0 ~ "Concordant Down",
      logFC_CRE > 0 & logFC_PLA < 0 ~ "Discordant (CRE Up / PLA Down)",
      TRUE                           ~ "Discordant (CRE Down / PLA Up)"),
    border_col   = ifelse(imputed, "grey50", "grey75"),
    point_size   = ifelse(sig_class == "NS", 1.8, 2.3),
    point_stroke = ifelse(sig_class == "NS", 0.6, 0.9),
    bubble_alpha = case_when(
      sig_class == "NS"          ~ 0.30,
      sig_class == "Interaction" ~ 0.55,
      sig_class == "Sig Both"    ~ 0.75,
      TRUE                       ~ 0.85))

universe <- scatter_df$gene
message(sprintf("  Total proteins: %d | Significant: %d",
                nrow(scatter_df), sum(scatter_df$is_sig)))

# -- ORA -----------------------------------------------------------------------
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
ora_q2 <- run_set_ora(scatter_df$gene[scatter_df$quadrant == "Discordant (CRE Down / PLA Up)"],
                       "Discordant (CRE Down / PLA Up)")
ora_q3 <- run_set_ora(scatter_df$gene[scatter_df$quadrant == "Concordant Down"],
                       "Concordant Down")
ora_q4 <- run_set_ora(scatter_df$gene[scatter_df$quadrant == "Discordant (CRE Up / PLA Down)"],
                       "Discordant (CRE Up / PLA Down)")

for (qn in c("Concordant Up", "Concordant Down",
             "Discordant (CRE Down / PLA Up)", "Discordant (CRE Up / PLA Down)"))
  message(sprintf("  %s: %d genes", qn,
                  sum(scatter_df$quadrant == qn)))

all_quad_ora <- bind_rows(ora_q1, ora_q2, ora_q3, ora_q4)
if (nrow(all_quad_ora) > 0)
  write_csv(all_quad_ora, file.path(DAT, "panel_B", "ora_quadrant.csv"))

# -- Scatter panel (composite version -- no title/subtitle) --------------------
x_range <- range(scatter_df$logFC_CRE, na.rm = TRUE)
y_range <- range(scatter_df$logFC_PLA, na.rm = TRUE)
x_pad   <- diff(x_range) * 0.10
y_pad   <- diff(y_range) * 0.10
xlim_range <- c(x_range[1] - x_pad, x_range[2] + x_pad)
ylim_range <- c(y_range[1] - y_pad, y_range[2] + y_pad)

ns_df  <- filter(scatter_df, sig_class == "NS")
sig_df <- filter(scatter_df, sig_class != "NS")

q_df <- scatter_df %>%
  mutate(q = case_when(
    logFC_CRE > 0 & logFC_PLA > 0 ~ "Q1",
    logFC_CRE < 0 & logFC_PLA < 0 ~ "Q3",
    logFC_CRE > 0 & logFC_PLA < 0 ~ "Q4",
    TRUE ~ "Q2"))
q_counts <- q_df %>% count(q) %>% deframe()
q_sig    <- q_df %>% filter(sig_class != "NS") %>% count(q) %>% deframe()
for (qq in c("Q1","Q2","Q3","Q4")) if (is.na(q_sig[qq])) q_sig[qq] <- 0

label_df <- sig_df %>%
  group_by(sig_class) %>%
  arrange(desc(abs(logFC_CRE) + abs(logFC_PLA))) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(label_fill     = SIG_LABEL_FILL_F3[as.character(sig_class)],
         label_text_col = SIG_LABEL_TEXT_F3[as.character(sig_class)])

txt_gene <- scale_text(BASE_GENE, 200) * 1.25
txt_quad <- scale_text(BASE_QUADRANT, 200) * 1.3

p_scatter <- ggplot(mapping = aes(x = logFC_CRE, y = logFC_PLA)) +
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
  scale_fill_manual(values = SIG_COLORS_F3, name = "Significance") +
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
           label = sprintf("Discordant (CRE\u2193 PLA\u2191)\nn = %s/%s", q_sig["Q2"], q_counts["Q2"]),
           hjust = 0.5, vjust = 0.5, size = txt_quad, fontface = "bold",
           color = COMP_BLUE, fill = alpha("white", 0.9),
           label.padding = unit(3, "pt"), lineheight = 0.85) +
  annotate("label", x = xlim_range[2] / 2, y = ylim_range[1] / 2,
           label = sprintf("Discordant (CRE\u2191 PLA\u2193)\nn = %s/%s", q_sig["Q4"], q_counts["Q4"]),
           hjust = 0.5, vjust = 0.5, size = txt_quad, fontface = "bold",
           color = COMP_BLUE, fill = alpha("white", 0.9),
           label.padding = unit(3, "pt"), lineheight = 0.85) +
  annotate("label", x = xlim_range[2] * 0.35, y = 0,
           label = expression(log[2]*FC ~ "(Training CRE)"),
           hjust = 0.5, vjust = 0.5, size = 3.5, color = "grey30", fontface = "bold",
           fill = alpha("white", 0.8), label.size = 0, label.padding = unit(2, "pt")) +
  annotate("label", x = 0, y = ylim_range[2] * 0.75,
           label = expression(log[2]*FC ~ "(Training PLA)"),
           hjust = 0.5, vjust = 0.5, size = 3.5, color = "grey30", fontface = "bold",
           fill = alpha("white", 0.8), label.size = 0, label.padding = unit(2, "pt"),
           angle = 90) +
  scale_x_continuous(breaks = pretty(xlim_range)) +
  scale_y_continuous(breaks = pretty(ylim_range)) +
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

# -- Half-bar builder ----------------------------------------------------------
make_half_bars <- function(df, fill_color, side, ylim) {
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
      display_name = stringr::str_wrap(pathway_label, width = 22),
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

# Compute ylim from scatter data
y_half <- max(abs(ylim_range))
p_ul <- make_half_bars(ora_q2, COMP_BLUE_LIGHT, "left",  c(0, y_half))
p_ll <- make_half_bars(ora_q3, COMP_RED_LIGHT,  "left",  c(-y_half, 0))
p_ur <- make_half_bars(ora_q1, COMP_RED_LIGHT,  "right", c(0, y_half))
p_lr <- make_half_bars(ora_q4, COMP_BLUE_LIGHT, "right", c(-y_half, 0))

# -- Composite -----------------------------------------------------------------
design <- c(
  area(1, 1),          # upper-left bars
  area(1, 2, 2, 2),   # scatter (spans both rows)
  area(1, 3),          # upper-right bars
  area(2, 1),          # lower-left bars
  area(2, 3)           # lower-right bars
)
n_total   <- nrow(scatter_df)
n_sig_c   <- sum(scatter_df$is_sig)
n_enrich  <- if (nrow(all_quad_ora) > 0) sum(all_quad_ora$significant) else 0L
r_pearson <- cor(scatter_df$logFC_CRE, scatter_df$logFC_PLA, use = "complete.obs")

composite <- p_ul + p_scatter + p_ur + p_ll + p_lr +
  plot_layout(design = design, widths = c(1.4, 2, 1.4)) +
  plot_annotation(
    title    = "Training Concordance: Quadrant ORA (CRE vs PLA)",
    subtitle = sprintf("Threshold-free quadrant ORA (hypergeometric) | N = %d proteins | %d DEPs (Pi < 0.05) | %d enriched terms (FDR < 0.05) | r = %.2f",
                        n_total, n_sig_c, n_enrich, r_pearson),
    theme    = theme(plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
                     plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey30")))

COMP_W <- 450; COMP_H <- 220
ggsave(file.path(RPT, "panel_B_ORA_composite_MAIN.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_B_ORA_composite_MAIN.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message("\nF05/CR Panel B ORA composite done")
