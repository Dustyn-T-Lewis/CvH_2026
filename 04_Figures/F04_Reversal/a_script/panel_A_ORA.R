# Reversal Panel A: Quadrant ORA Scatter + Flanking Bars
# Cancer Recovery Reversal: CRvH_Baseline (x) vs CR_Training (y)
# Threshold-free ORA on all proteins per quadrant
# Blue = reversed (off-diagonal), Red = exacerbated (diagonal)
setwd(here::here())
source("04_Figures/shared/style.R")
source("04_Figures/shared/stats.R")
source("04_Figures/shared/pathway_utils.R")
pacman::p_load(tidyverse, fgsea, ggrepel, patchwork)

RPT_PNG <- "04_Figures/F04_Reversal/b_reports/main/png/panels"
RPT_PDF <- "04_Figures/F04_Reversal/b_reports/main/pdf/panels"
DAT <- "04_Figures/F04_Reversal/c_data"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_A"), recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

COMP_RED <- unname(DIR_COLORS["Up"])
COMP_BLUE <- unname(DIR_COLORS["Down"])
N_SHOW <- 5

# Data
source("04_Figures/F04_Reversal/a_script/f04_data.R") # dep_df + imputation_df

scatter_df <- dep_df %>%
  transmute(gene,
    logFC_CvH = logFC_CRvH_Baseline,
    logFC_TR  = logFC_CR_Training,
    pi_CvH    = pi_score_CRvH_Baseline,
    pi_TR     = pi_score_CR_Training
  ) %>%
  filter(!is.na(logFC_CvH), !is.na(logFC_TR)) %>%
  left_join(imputation_df, by = "gene") %>%
  mutate(
    imputed = replace_na(imputed, FALSE),
    sig_class = classify_proteins_f4(pi_CvH, pi_TR),
    is_sig = sig_class != "NS",
    quadrant = case_when(
      logFC_CvH > 0 & logFC_TR < 0 ~ "Reversed (Cancer Up / Training Down)",
      logFC_CvH < 0 & logFC_TR > 0 ~ "Reversed (Cancer Down / Training Up)",
      logFC_CvH > 0 & logFC_TR > 0 ~ "Exacerbated Up",
      TRUE ~ "Exacerbated Down"
    )
  )

universe <- unique(scatter_df$gene)
message(sprintf(
  "  Total proteins: %d | Significant: %d",
  nrow(scatter_df), sum(scatter_df$is_sig)
))

# ORA
pw_collection <- build_pathway_collection(
  min_size = 15, max_size = 500,
  include_goslim = FALSE,
  exclude_variants = TRUE
)

run_set_ora <- function(genes, set_name) {
  if (length(genes) < 5) {
    return(tibble())
  }
  res <- tryCatch(
    run_ora_deduplicated(
      genes = genes, universe = universe,
      pathways = pw_collection, jaccard_cutoff = 0.5,
      min_size = 15, max_size = 500, padj_cutoff = 1
    ),
    error = function(e) {
      message("  ORA error: ", e$message)
      tibble()
    }
  )
  if (nrow(res) == 0) {
    return(tibble())
  }
  res %>%
    mutate(
      set = set_name, pathway_label = clean_pathway_name(pathway),
      neg_log10_padj = -log10(padj), significant = padj < 0.05
    ) %>%
    arrange(desc(neg_log10_padj)) %>%
    slice_head(n = N_SHOW)
}

message("\n--- Quadrant ORA (threshold-free) ---")
ora_tl <- run_set_ora(
  scatter_df$gene[scatter_df$quadrant == "Reversed (Cancer Down / Training Up)"],
  "Reversed (Cancer Down / Training Up)"
)
ora_tr <- run_set_ora(
  scatter_df$gene[scatter_df$quadrant == "Exacerbated Up"],
  "Exacerbated Up"
)
ora_bl <- run_set_ora(
  scatter_df$gene[scatter_df$quadrant == "Exacerbated Down"],
  "Exacerbated Down"
)
ora_br <- run_set_ora(
  scatter_df$gene[scatter_df$quadrant == "Reversed (Cancer Up / Training Down)"],
  "Reversed (Cancer Up / Training Down)"
)

all_quad_ora <- bind_rows(ora_tl, ora_tr, ora_bl, ora_br)
if (nrow(all_quad_ora) > 0) {
  write_csv(all_quad_ora, file.path(DAT, "panel_A", "ora_quadrant.csv"))
}

# Scatter panel
xlim_range <- range(scatter_df$logFC_CvH, na.rm = TRUE) * 1.15
ylim_range <- range(scatter_df$logFC_TR, na.rm = TRUE) * 1.15

ns_df <- filter(scatter_df, sig_class == "NS")
sig_df <- filter(scatter_df, sig_class != "NS")

q_df <- scatter_df %>%
  mutate(q = case_when(
    logFC_CvH > 0 & logFC_TR < 0 ~ "BR",
    logFC_CvH < 0 & logFC_TR > 0 ~ "TL",
    logFC_CvH > 0 & logFC_TR > 0 ~ "TR",
    TRUE ~ "BL"
  ))
q_counts <- q_df %>%
  count(q) %>%
  deframe()
q_sig <- q_df %>%
  filter(sig_class != "NS") %>%
  count(q) %>%
  deframe()
for (qq in c("BR", "TL", "TR", "BL")) if (is.na(q_sig[qq])) q_sig[qq] <- 0

label_df <- sig_df %>%
  group_by(sig_class) %>%
  arrange(desc(abs(logFC_CvH) + abs(logFC_TR))) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    label_fill = SIG_LABEL_FILL_F4[as.character(sig_class)],
    label_text_col = SIG_LABEL_TEXT_F4[as.character(sig_class)]
  )

txt_gene <- 2.8
txt_quad <- 2.8

# Center-axis tick positions
x_breaks <- pretty(xlim_range, n = 6)
y_breaks <- pretty(ylim_range, n = 6)
x_tick_df <- tibble(
  x = x_breaks[x_breaks != 0], y = 0,
  label = as.character(x_breaks[x_breaks != 0])
)
y_tick_df <- tibble(
  x = 0, y = y_breaks[y_breaks != 0],
  label = as.character(y_breaks[y_breaks != 0])
)

p_scatter <- ggplot(mapping = aes(x = logFC_CvH, y = logFC_TR)) +
  # Quadrant backgrounds: blue = reversed, red = exacerbated
  annotate("rect",
    xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
    fill = "#DCEEFF", alpha = 0.55, color = "grey70", linewidth = 0.2
  ) +
  annotate("rect",
    xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
    fill = "#DCEEFF", alpha = 0.55, color = "grey70", linewidth = 0.2
  ) +
  annotate("rect",
    xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
    fill = "#FFE0E0", alpha = 0.55, color = "grey70", linewidth = 0.2
  ) +
  annotate("rect",
    xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
    fill = "#FFE0E0", alpha = 0.55, color = "grey70", linewidth = 0.2
  ) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_abline(
    slope = -1, intercept = 0, linetype = "dashed",
    color = "black", linewidth = 0.3
  ) +
  # Center-axis tick labels
  geom_text(
    data = x_tick_df, aes(x = x, y = y, label = label),
    vjust = 1.5, size = 2.5, color = "grey40", fontface = "bold"
  ) +
  geom_text(
    data = y_tick_df, aes(x = x, y = y, label = label),
    hjust = -0.5, size = 2.5, color = "grey40", fontface = "bold"
  ) +
  # NS points
  geom_point(
    data = ns_df, color = "grey80", fill = "grey85", shape = 21,
    size = 0.35, alpha = 0.3, stroke = 0.10
  ) +
  # Sig points
  geom_point(
    data = sig_df, aes(fill = sig_class), shape = 21,
    size = ifelse(sig_df$sig_class == "NS", 0.6, 0.9),
    color = ifelse(sig_df$imputed, "black", "grey75"),
    alpha = case_when(
      sig_df$sig_class == "NS" ~ 0.30,
      sig_df$sig_class == "Sig Both" ~ 0.75,
      TRUE ~ 0.85
    ),
    stroke = ifelse(sig_df$sig_class == "NS", 0.4, 0.6)
  ) +
  scale_fill_manual(values = SIG_COLORS_F4, name = "Significance") +
  # Gene labels
  geom_label_repel(
    data = label_df, aes(label = gene),
    fill = label_df$label_fill, color = label_df$label_text_col,
    size = txt_gene, fontface = "italic", max.overlaps = 50,
    segment.size = 0.3, segment.color = "grey25",
    min.segment.length = 0, show.legend = FALSE,
    box.padding = 0.3, point.padding = 0.3,
    force = 6, force_pull = 0.3,
    label.padding = unit(1, "pt"), label.r = unit(0.5, "pt"),
    linewidth = 0.15, seed = 42
  ) +
  # Quadrant labels
  annotate("label",
    x = xlim_range[1], y = ylim_range[2],
    label = sprintf("Reversed (Cancer\u2193 Tr\u2191)\n%s/%s", q_sig["TL"], q_counts["TL"]),
    hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
    color = COMP_BLUE, fill = alpha("white", 0.92),
    label.padding = unit(2.5, "pt"), lineheight = 0.9
  ) +
  annotate("label",
    x = xlim_range[2], y = ylim_range[2],
    label = sprintf("Exacerbated Up\n%s/%s", q_sig["TR"], q_counts["TR"]),
    hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
    color = COMP_RED, fill = alpha("white", 0.92),
    label.padding = unit(2.5, "pt"), lineheight = 0.9
  ) +
  annotate("label",
    x = xlim_range[1], y = ylim_range[1],
    label = sprintf("%s/%s\nExacerbated Down", q_sig["BL"], q_counts["BL"]),
    hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
    color = COMP_RED, fill = alpha("white", 0.92),
    label.padding = unit(2.5, "pt"), lineheight = 0.9
  ) +
  annotate("label",
    x = xlim_range[2], y = ylim_range[1],
    label = sprintf("%s/%s\nReversed (Cancer\u2191 Tr\u2193)", q_sig["BR"], q_counts["BR"]),
    hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
    color = COMP_BLUE, fill = alpha("white", 0.92),
    label.padding = unit(2.5, "pt"), lineheight = 0.9
  ) +
  # Axis titles inside plot
  annotate("text",
    x = mean(xlim_range[xlim_range > 0]), y = 0,
    label = expression(log[2] * FC ~ "(Cancer vs Healthy)"),
    hjust = 0.5, vjust = -0.4, size = 3.0, color = "grey30", fontface = "bold"
  ) +
  annotate("text",
    x = 0, y = mean(ylim_range[ylim_range > 0]),
    label = expression(log[2] * FC ~ "(Training CR)"),
    hjust = 0.5, vjust = -0.4, size = 3.0, color = "grey30", fontface = "bold",
    angle = 90
  ) +
  coord_cartesian(xlim = xlim_range, ylim = ylim_range, expand = FALSE) +
  labs(x = NULL, y = NULL) +
  FIG_THEME +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(2, 0, 0, 0, "mm"),
    legend.position = "none"
  )

# Custom Significance key
key_lvls <- c("Sig Both", "Sig Cancer only", "Sig Training only")
key_display <- c("Sig Both", "Sig Cancer", "Sig Training")
key_df <- tibble(
  category = factor(key_lvls, levels = key_lvls),
  display  = key_display,
  fill_col = unname(SIG_COLORS_F4[key_lvls]),
  x        = c(1.25, 1.90, 2.55),
  y        = 0
)
p_key <- ggplot(key_df, aes(x = x, y = y)) +
  geom_point(aes(fill = category),
    shape = 21, size = 2.5,
    color = "grey50", stroke = 0.6, alpha = 0.85,
    show.legend = FALSE
  ) +
  geom_text(aes(label = display),
    nudge_x = 0.06, hjust = 0,
    size = 2.5, fontface = "bold", color = "grey25"
  ) +
  scale_fill_manual(values = setNames(key_df$fill_col, key_df$category)) +
  scale_x_continuous(limits = c(0.2, 4.0), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.15, 0.15), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(-24, 0, 0, 0, "mm"))

# Half-bar builder
make_half_bars <- function(df, fill_color, side, ylim) {
  bar_h <- 0.42
  n_bars <- if (is.null(df) || nrow(df) == 0) 0L else min(nrow(df), 5L)

  if (n_bars == 0) {
    return(ggplot() +
      theme_void() +
      scale_y_continuous(limits = ylim, expand = c(0, 0)))
  }

  y_pos <- if (ylim[1] >= 0) {
    rev(seq(0.3, 2.3, length.out = 5))[seq_len(n_bars)]
  } else {
    seq(-0.3, -2.5, length.out = 5)[seq_len(n_bars)]
  }

  bars <- df %>%
    arrange(desc(neg_log10_padj)) %>%
    slice_head(n = 5) %>%
    mutate(
      y = y_pos,
      bar_fill = ifelse(significant, scales::alpha(fill_color, 0.85),
        scales::alpha(fill_color, 0.30)
      ),
      name = pathway_label,
      p_lab = ifelse(padj < 0.001, "p<0.001", sprintf("p=%.3f", padj)),
      p_face = ifelse(significant, "bold", "plain")
    )

  x_max <- max(bars$neg_log10_padj)
  x_display_max <- x_max * 1.45
  is_upper <- ylim[1] >= 0
  out_hjust <- if (side == "left") 1 else 0
  brk_fn <- function(limits) {
    b <- scales::pretty_breaks(n = 3)(limits)
    b[b != 0]
  }

  # Names fit inside long bars (ggfittext reflows and shrinks, auto-contrast);
  # short bars carry the name outside. The p-value always sits past the bar tip,
  # bold when FDR < 0.05.
  bars <- bars %>%
    mutate(
      fits   = neg_log10_padj >= x_max * 0.5,
      p_x    = neg_log10_padj + x_max * 0.03,
      name_x = neg_log10_padj + x_max * 0.27
    )

  p <- ggplot(bars, aes(y = y)) +
    geom_rect(
      aes(
        xmin = 0, xmax = neg_log10_padj,
        ymin = y - bar_h / 2, ymax = y + bar_h / 2
      ),
      fill = bars$bar_fill, color = "black", linewidth = 0.3
    ) +
    ggfittext::geom_fit_text(
      data = ~ dplyr::filter(.x, fits),
      aes(
        xmin = 0, xmax = neg_log10_padj,
        ymin = y - bar_h / 2, ymax = y + bar_h / 2, label = name
      ),
      reflow = TRUE, grow = FALSE, contrast = TRUE, fontface = "bold",
      min.size = 3
    ) +
    geom_text(
      data = ~ dplyr::filter(.x, !fits),
      aes(x = name_x, y = y, label = name), hjust = out_hjust,
      size = 2.2, fontface = "bold", color = "grey15", lineheight = 0.85
    ) +
    geom_text(aes(x = p_x, y = y, label = p_lab, fontface = p_face),
      hjust = out_hjust, size = 2.0, color = "grey25"
    ) +
    labs(
      x = if (!is_upper) expression(-log[10](p[adj])) else NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(
        size = FIG_AXIS_TEXT, face = "bold",
        margin = margin(t = 0, unit = "mm")
      ),
      axis.title.x = if (!is_upper) {
        element_text(
          size = 8, face = "bold",
          margin = margin(t = 0, unit = "mm")
        )
      } else {
        element_blank()
      },
      axis.line.x = element_line(color = "grey50", linewidth = 0.3),
      axis.ticks.x = element_line(color = "grey50", linewidth = 0.3),
      plot.margin = if (is_upper && side == "left") {
        margin(4, 0, 0, 3, "mm")
      } else if (is_upper) {
        margin(4, 3, 0, 0, "mm")
      } else if (side == "left") {
        margin(2, 0, 0, 3, "mm")
      } else {
        margin(2, 3, 0, 0, "mm")
      }
    )

  if (side == "left") {
    p + scale_x_reverse(
      limits = c(x_display_max, 0),
      breaks = brk_fn,
      expand = expansion(mult = c(0, 0))
    ) +
      scale_y_continuous(limits = ylim, expand = c(0, 0)) +
      coord_cartesian(clip = "off")
  } else {
    p + scale_x_continuous(
      limits = c(0, x_display_max),
      breaks = brk_fn,
      expand = expansion(mult = c(0, 0))
    ) +
      scale_y_continuous(limits = ylim, expand = c(0, 0)) +
      coord_cartesian(clip = "off")
  }
}

# 4 half-bar panels
p_ul <- make_half_bars(
  ora_tl, scales::alpha(COMP_BLUE, 0.30), "left",
  c(0, abs(ylim_range[2]))
)
p_ll <- make_half_bars(
  ora_bl, scales::alpha(COMP_RED, 0.30), "left",
  c(-abs(ylim_range[2]), 0)
)
p_ur <- make_half_bars(
  ora_tr, scales::alpha(COMP_RED, 0.30), "right",
  c(0, abs(ylim_range[2]))
)
p_lr <- make_half_bars(
  ora_br, scales::alpha(COMP_BLUE, 0.30), "right",
  c(-abs(ylim_range[2]), 0)
)

# Composite
design <- c(
  area(1, 1), # p_ul
  area(1, 2, 2, 2), # p_scatter (rows 1-2, center)
  area(1, 3), # p_ur
  area(2, 1), # p_ll
  area(2, 3), # p_lr
  area(3, 1, 3, 3) # key spans full width
)
n_total <- nrow(scatter_df)
n_sig <- sum(scatter_df$is_sig)
n_enrich <- if (nrow(all_quad_ora) > 0) sum(all_quad_ora$significant) else 0L
r_spear <- cor(scatter_df$logFC_CvH, scatter_df$logFC_TR,
  use = "complete.obs",
  method = "spearman"
)

composite <- p_ul + p_scatter + p_ur + p_ll + p_lr + p_key +
  plot_layout(
    design = design,
    widths = c(70, 100, 70) / 240,
    heights = c(85, 85, 8) / 178
  ) +
  plot_annotation(
    title = "Cancer Recovery Reversal: Quadrant ORA",
    subtitle = sprintf(
      "Threshold-free ORA | N = %d | %d DEPs (\u03a0 < 0.05) | %d enriched (FDR < 0.05) | \u03c1 = %.2f",
      n_total, n_sig, n_enrich, r_spear
    ),
    theme = theme(
      plot.title = element_text(size = FIG_TITLE_SIZE, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, hjust = 0, color = "grey30"),
      plot.title.position = "panel"
    )
  )

COMP_W <- 200
COMP_H <- 120
ggsave(file.path(RPT_PNG, "MAIN_panel_A_ORA_composite.png"), composite,
  width = COMP_W, height = COMP_H, units = "mm", dpi = 300
)
ggsave(file.path(RPT_PDF, "MAIN_panel_A_ORA_composite.pdf"), composite,
  width = COMP_W, height = COMP_H, units = "mm", device = pdf_device
)

composite <- composite &
  labs(title = NULL, subtitle = NULL, tag = NULL) &
  theme(legend.position = "none")
composite <- composite +
  plot_annotation(
    title = NULL, subtitle = NULL,
    theme = theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )
  )
