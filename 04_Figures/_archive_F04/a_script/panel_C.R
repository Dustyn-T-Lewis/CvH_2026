# F04 Panel C: Recovery Pattern Heatmap + Sankey + Scatterpie
# Patterns: Reversed, Partially Reversed, Persistent, Exacerbated
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")
source("04_Figures/shared/go_slim_categories.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggnewscale)
})

PF_W <- 360

RPT <- "04_Figures/F04/b_reports"
DAT <- "04_Figures/F04/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_C"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# --- Load data ---
dep_df <- read_csv("03_DEP/c_data/03_combined_results_CRvH.csv", show_col_types = FALSE)
imputation_df <- read_csv("02_Imputation/c_data/02_mar_mnar_classification.csv",
                           show_col_types = FALSE) %>%
  transmute(gene, imputed = classification != "Complete")

# --- Pattern classification ---
expanded_df <- dep_df %>%
  filter(!is.na(logFC_Cancer_vs_Healthy), !is.na(logFC_Training_CR)) %>%
  mutate(
    pattern = case_when(
      pi_score_Cancer_vs_Healthy < 0.05 & pi_score_Training_CR < 0.05 &
        sign(logFC_Cancer_vs_Healthy) != sign(logFC_Training_CR) ~ "Reversed",
      pi_score_Cancer_vs_Healthy < 0.05 & pi_score_Training_CR >= 0.05 &
        sign(logFC_Cancer_vs_Healthy) != sign(logFC_Training_CR) ~ "Partially Reversed",
      pi_score_Cancer_vs_Healthy < 0.05 & pi_score_Training_CR < 0.05 &
        sign(logFC_Cancer_vs_Healthy) == sign(logFC_Training_CR) ~ "Exacerbated",
      pi_score_Cancer_vs_Healthy < 0.05 ~ "Persistent",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(pattern)) %>%
  mutate(pattern = factor(pattern, levels = PATTERN_ORDER_F4)) %>%
  left_join(imputation_df, by = "gene") %>%
  mutate(imputed = replace_na(imputed, FALSE))

# Sort key
expanded_df <- expanded_df %>%
  mutate(sort_key = case_when(
    pattern %in% c("Reversed", "Partially Reversed") ~ -logFC_Cancer_vs_Healthy,
    pattern == "Exacerbated" ~ logFC_Cancer_vs_Healthy,
    pattern == "Persistent"  ~ abs(logFC_Cancer_vs_Healthy)
  ))

n_expanded <- nrow(expanded_df)

# GO Slim classification
expanded_genes <- expanded_df$gene
all_genes      <- dep_df$gene
go_slim_result <- assign_go_slim_consolidated(expanded_genes, all_genes)

links_1to1 <- go_slim_result %>%
  rename(pathway = consolidated) %>%
  mutate(pathway = as.character(pathway)) %>%
  left_join(expanded_df %>% select(gene, pattern, sort_key,
                                    logFC_Cancer_vs_Healthy, logFC_Training_CR),
            by = "gene")

gene_order <- links_1to1 %>%
  mutate(pat_idx = match(as.character(pattern), PATTERN_ORDER_F4)) %>%
  arrange(pat_idx, desc(sort_key)) %>%
  pull(gene)

n_genes <- length(gene_order)
Y_SPAN  <- n_genes

# --- Heatmap geometry ---
n_pat_breaks <- length(PATTERN_ORDER_F4) - 1
gap_size     <- 1.2
total_gap    <- n_pat_breaks * gap_size
gene_h       <- (Y_SPAN - total_gap) / n_genes

gene_bar_df <- tibble(gene = gene_order, idx = seq_along(gene_order)) %>%
  mutate(
    ymax  = Y_SPAN - (idx - 1) * gene_h,
    ymin  = ymax - gene_h,
    y_ctr = (ymin + ymax) / 2
  ) %>%
  left_join(links_1to1 %>% select(gene, pattern) %>% distinct(), by = "gene")

cat_breaks <- gene_bar_df %>%
  mutate(pat_str = as.character(pattern)) %>%
  mutate(is_break = pat_str != lag(pat_str, default = pat_str[1]))
cum_gap <- cumsum(cat_breaks$is_break) * gap_size
gene_bar_df <- gene_bar_df %>%
  mutate(ymax = ymax - cum_gap, ymin = ymin - cum_gap,
         y_ctr = (ymin + ymax) / 2)

HM_X_CVH <- 1
HM_X_TR  <- 1.92
HM_HW    <- 0.45
ANNO_X   <- HM_X_TR + HM_HW + 0.01 + 0.12
ANNO_HW  <- 0.12
S_X_PW   <- 5.2
S_PW_HW  <- 0.12

heatmap_df <- expanded_df %>%
  filter(gene %in% gene_order) %>%
  left_join(gene_bar_df %>% select(gene, y_ctr, ymin, ymax), by = "gene") %>%
  select(gene, y_ctr, ymin, ymax, logFC_Cancer_vs_Healthy, logFC_Training_CR) %>%
  pivot_longer(cols = c(logFC_Cancer_vs_Healthy, logFC_Training_CR),
               names_to = "contrast", values_to = "logFC") %>%
  mutate(x = ifelse(contrast == "logFC_Cancer_vs_Healthy", HM_X_CVH, HM_X_TR))

fc_max <- max(abs(heatmap_df$logFC), na.rm = TRUE)

anno_rects <- gene_bar_df %>%
  mutate(fill = PATTERN_COLS_F4[as.character(pattern)])

cat_heatmap_pos <- gene_bar_df %>%
  group_by(pattern) %>%
  summarise(y_bot = min(ymin), y_top = max(ymax), .groups = "drop") %>%
  mutate(y_ctr = (y_top + y_bot) / 2, bar_h = y_top - y_bot)

group_outlines <- cat_heatmap_pos %>%
  mutate(xmin = HM_X_CVH - HM_HW, xmax = ANNO_X + ANNO_HW,
         outline_col = PATTERN_COLS_F4[as.character(pattern)])

S_Y_MIN  <- min(gene_bar_df$ymin) - gene_h * 2
S_Y_MAX  <- max(gene_bar_df$ymax) + gene_h * 2
S_Y_SPAN <- S_Y_MAX - S_Y_MIN

pat_counts <- links_1to1 %>% count(pattern) %>% deframe()

PATTERN_SUBTITLES_F4_local <- c(
  "Reversed"            = "Opposite sign, both sig",
  "Partially Reversed"  = "Opposite sign, cancer sig",
  "Persistent"          = "Cancer sig, not reversed",
  "Exacerbated"         = "Same sign, both sig"
)

inline_label_df <- cat_heatmap_pos %>%
  mutate(
    x     = HM_X_CVH - HM_HW - 0.15,
    n     = ifelse(as.character(pattern) %in% names(pat_counts),
                   pat_counts[as.character(pattern)], 0L),
    label = sprintf("%s\n%s (n=%d)",
                    as.character(pattern),
                    PATTERN_SUBTITLES_F4_local[as.character(pattern)],
                    n),
    color = PATTERN_COLS_F4[as.character(pattern)]
  )

# --- Sankey ribbons ---
cat_pw_counts <- links_1to1 %>%
  filter(pathway != "Other") %>%
  count(pattern, pathway, name = "n_proteins") %>%
  filter(n_proteins > 0)

cat_totals <- links_1to1 %>%
  filter(pathway != "Other") %>%
  count(pattern, name = "total") %>%
  arrange(match(as.character(pattern), PATTERN_ORDER_F4))

pw_1to1_counts <- links_1to1 %>% count(pathway, name = "n_1to1")
pw_order_df <- pw_1to1_counts %>% arrange(desc(n_1to1))
pw_order <- pw_order_df$pathway
if ("Other" %in% pw_order) pw_order <- c(setdiff(pw_order, "Other"), "Other")

pw_totals <- links_1to1 %>%
  filter(pathway != "Other") %>%
  count(pathway, name = "total") %>%
  mutate(pw_idx = match(pathway, pw_order)) %>%
  filter(!is.na(pw_idx)) %>%
  arrange(pw_idx)

n_cats <- nrow(cat_totals)
n_pws  <- nrow(pw_totals)

cat_bars <- cat_totals %>%
  left_join(cat_heatmap_pos, by = "pattern") %>%
  mutate(fill = PATTERN_COLS_F4[as.character(pattern)])

pw_gap_frac <- 0.03
pw_usable   <- S_Y_SPAN * (1 - pw_gap_frac * max(n_pws - 1, 0) / max(n_pws, 1))
pw_gap_size <- if (n_pws > 1) (S_Y_SPAN - pw_usable) / (n_pws - 1) else 0

pw_bars <- pw_totals %>%
  mutate(
    bar_h = pw_usable / n_pws,
    y_top = S_Y_MAX - cumsum(c(0, head(bar_h + pw_gap_size, -1))),
    y_bot = y_top - bar_h,
    y_ctr = (y_top + y_bot) / 2
  )

active_pw_order <- pw_totals$pathway
pw_bars$fill <- CONSOLIDATED_COLORS[pw_bars$pathway]
pw_bars <- pw_bars %>%
  mutate(pw_label = stringr::str_wrap(pathway, width = 22))

txt_pf <- scale_text(BASE_STAT, PF_W)

cat_cum <- setNames(rep(0, n_cats), as.character(cat_totals$pattern))
pw_cum  <- setNames(rep(0, n_pws), pw_totals$pathway)
ribbon_list <- list()
ribbon_idx  <- 0

for (ct in as.character(cat_totals$pattern)) {
  ct_contribs <- cat_pw_counts %>%
    filter(as.character(pattern) == ct) %>%
    arrange(match(pathway, active_pw_order))
  for (r in seq_len(nrow(ct_contribs))) {
    pw <- ct_contribs$pathway[r]; n <- ct_contribs$n_proteins[r]
    if (n == 0) next
    ribbon_idx <- ribbon_idx + 1
    ct_row <- cat_bars %>% filter(as.character(pattern) == ct)
    frac_ct <- n / ct_row$total
    y0_top <- ct_row$y_top - cat_cum[ct] * ct_row$bar_h
    y0_bot <- y0_top - frac_ct * ct_row$bar_h
    cat_cum[ct] <- cat_cum[ct] + frac_ct
    pw_row <- pw_bars %>% filter(pathway == pw)
    frac_pw <- n / pw_row$total
    y1_top <- pw_row$y_top - pw_cum[pw] * pw_row$bar_h
    y1_bot <- y1_top - frac_pw * pw_row$bar_h
    pw_cum[pw] <- pw_cum[pw] + frac_pw
    ribbon_list[[ribbon_idx]] <- make_sigmoid_ribbon(
      x0 = ANNO_X + ANNO_HW, x1 = S_X_PW - S_PW_HW,
      y0_top = y0_top, y0_bot = y0_bot,
      y1_top = y1_top, y1_bot = y1_bot,
      ribbon_id = paste0(ct, "->", pw)
    ) %>% mutate(pattern = ct, pathway = pw, fill_col = PATTERN_COLS_F4[ct])
  }
}
ribbons_df <- bind_rows(ribbon_list)

# --- Scatterpie ---
pie_counts <- links_1to1 %>%
  filter(pathway != "Other") %>%
  count(pathway, pattern, name = "n") %>%
  complete(pathway = active_pw_order,
           pattern = factor(PATTERN_ORDER_F4, levels = PATTERN_ORDER_F4),
           fill = list(n = 0L))

pie_totals <- pie_counts %>%
  group_by(pathway) %>%
  summarise(total = sum(n), .groups = "drop")

max_total <- max(pie_totals$total)

PIE_X_START <- S_X_PW + S_PW_HW + 0.3
PIE_X_SCALE <- 3.4 / max_total
PIE_X_END   <- PIE_X_START + 3.4
PIE_R       <- pw_bars$bar_h[1] * 0.38

header_y <- S_Y_MIN - gene_h * 0.5
f_height <- max(300, n_genes * 0.8 + 80)

x_range <- (PIE_X_END + 1.5) - (HM_X_CVH - HM_HW - 3.0)
y_range <- (S_Y_MAX + 2) - (header_y - 7)
aspect  <- (PF_W / x_range) / (f_height / y_range)
pie_r_x <- PIE_R / aspect

pie_arcs <- pie_counts %>%
  left_join(pie_totals, by = "pathway") %>%
  filter(total > 0) %>%
  left_join(pw_bars %>% select(pathway, y_ctr), by = "pathway") %>%
  group_by(pathway) %>%
  mutate(
    frac       = n / total,
    end_frac   = cumsum(frac),
    start_frac = lag(end_frac, default = 0),
    start      = start_frac * 2 * pi,
    end        = end_frac * 2 * pi,
    x0         = PIE_X_START + total * PIE_X_SCALE
  ) %>%
  ungroup() %>%
  filter(n > 0)

n_arc_pts  <- 30
pie_polys  <- list()
for (i in seq_len(nrow(pie_arcs))) {
  row   <- pie_arcs[i, ]
  theta <- seq(row$start, row$end, length.out = n_arc_pts)
  xs    <- c(row$x0, row$x0 + pie_r_x * sin(theta), row$x0)
  ys    <- c(row$y_ctr, row$y_ctr + PIE_R * cos(theta), row$y_ctr)
  pie_polys[[i]] <- tibble(
    x = xs, y = ys,
    wedge_id = paste0(row$pathway, "_", as.character(row$pattern)),
    pattern  = as.character(row$pattern),
    pathway  = row$pathway
  )
}
pie_polys_df <- bind_rows(pie_polys)

pie_labels <- pie_totals %>%
  left_join(pw_bars %>% select(pathway, y_ctr), by = "pathway") %>%
  mutate(
    label = paste0("n=", total),
    x     = PIE_X_START + total * PIE_X_SCALE + pie_r_x + 0.15
  )

count_ticks <- pretty(c(0, max_total), n = 4)
count_ticks <- count_ticks[count_ticks >= 0 & count_ticks <= max_total * 1.1]
tick_x      <- PIE_X_START + count_ticks * PIE_X_SCALE
tick_df     <- tibble(x = tick_x, label = as.character(as.integer(count_ticks)))
axis_y      <- S_Y_MIN - 2
tick_len    <- 1.5

# --- Combined panel ---
pF <- ggplot() +
  geom_rect(data = heatmap_df,
            aes(xmin = x - HM_HW, xmax = x + HM_HW,
                ymin = ymin, ymax = ymax, fill = logFC),
            color = NA, linewidth = 0) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-fc_max, fc_max),
    name = expression(log[2]~FC)
  ) +
  ggnewscale::new_scale_fill() +
  geom_rect(data = anno_rects,
            aes(xmin = ANNO_X - ANNO_HW, xmax = ANNO_X + ANNO_HW,
                ymin = ymin, ymax = ymax, fill = fill),
            color = NA, linewidth = 0) +
  geom_rect(data = group_outlines,
            aes(xmin = xmin, xmax = xmax, ymin = y_bot, ymax = y_top),
            fill = NA, color = "black", linewidth = 0.5) +
  geom_polygon(data = ribbons_df,
               aes(x = x, y = y, group = ribbon_id, fill = fill_col),
               alpha = 0.30, color = NA) +
  geom_rect(data = pw_bars,
            aes(xmin = S_X_PW - S_PW_HW, xmax = S_X_PW + S_PW_HW,
                ymin = y_bot, ymax = y_top, fill = fill),
            color = "black", linewidth = 0.3) +
  geom_text(data = pw_bars,
            aes(x = S_X_PW - S_PW_HW - 0.05, y = y_ctr, label = pw_label),
            hjust = 1, size = txt_pf, fontface = "bold",
            color = "grey20", lineheight = 0.85) +
  geom_polygon(data = pie_polys_df,
               aes(x = x, y = y, group = wedge_id,
                   fill = PATTERN_COLS_F4[pattern]),
               color = "white", linewidth = 0.3) +
  geom_text(data = pie_labels,
            aes(x = x, y = y_ctr, label = label),
            hjust = 0, size = txt_pf * 0.65, color = "grey40") +
  geom_text(data = inline_label_df,
            aes(x = x, y = y_ctr, label = label, color = color),
            hjust = 1, size = txt_pf * 0.85, fontface = "bold",
            lineheight = 0.85, show.legend = FALSE) +
  scale_color_identity() +
  scale_fill_identity() +
  annotate("text", x = HM_X_CVH, y = header_y, label = "Cancer\nvs Healthy",
           size = txt_pf, fontface = "bold", lineheight = 0.85) +
  annotate("text", x = HM_X_TR, y = header_y, label = "Training\n(CR)",
           size = txt_pf, fontface = "bold", lineheight = 0.85) +
  annotate("segment", x = PIE_X_START, xend = PIE_X_END,
           y = axis_y, yend = axis_y, color = "black", linewidth = 0.5) +
  geom_segment(data = tick_df, aes(x = x, xend = x),
               y = axis_y, yend = axis_y - tick_len,
               color = "black", linewidth = 0.4) +
  geom_text(data = tick_df, aes(x = x, label = label),
            y = axis_y - tick_len - 1.5, size = txt_pf, color = "grey30") +
  annotate("text", x = (PIE_X_START + PIE_X_END) / 2,
           y = axis_y - tick_len - 5.0,
           label = "Protein count",
           size = txt_pf, fontface = "bold", color = "grey20") +
  scale_y_continuous(breaks = NULL,
    limits = c(axis_y - tick_len - 7.5, S_Y_MAX + 2), expand = c(0, 0)) +
  scale_x_continuous(
    limits = c(HM_X_CVH - HM_HW - 3.0, PIE_X_END + 1.5), expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title = "Cancer Recovery: Functional Classification",
       subtitle = sprintf("%d cancer-signature proteins | %d GO Slim categories | 4 recovery patterns",
                          nrow(expanded_df), length(active_pw_order))) +
  FIG_THEME +
  theme(
    axis.text.y  = element_blank(), axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
    panel.grid   = element_blank(), panel.border = element_blank(),
    legend.position = "none",
    plot.margin  = margin(8, 5, 8, 1)
  )

# --- Data exports ---
expanded_df %>%
  transmute(gene, pattern = as.character(pattern), sort_key = round(sort_key, 4),
            logFC_Cancer_vs_Healthy = round(logFC_Cancer_vs_Healthy, 4),
            logFC_Training_CR       = round(logFC_Training_CR, 4),
            pi_score_CvH            = round(pi_score_Cancer_vs_Healthy, 6),
            pi_score_TR             = round(pi_score_Training_CR, 6),
            imputed) %>%
  arrange(match(pattern, PATTERN_ORDER_F4), desc(sort_key)) %>%
  write_csv(file.path(DAT, "panel_C", "pattern_classification.csv"))

links_1to1 %>%
  transmute(gene, pathway, go_slim_id = slim,
            pattern = as.character(pattern),
            logFC_Cancer_vs_Healthy = round(logFC_Cancer_vs_Healthy, 4),
            logFC_Training_CR = round(logFC_Training_CR, 4)) %>%
  arrange(match(pattern, PATTERN_ORDER_F4), pathway, gene) %>%
  write_csv(file.path(DAT, "panel_C", "sankey_links.csv"))

ggsave(file.path(RPT, "panel_C_recovery_heatmap.pdf"), pF,
       width = PF_W, height = f_height, units = "mm", device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "panel_C_recovery_heatmap.png"), pF,
       width = PF_W, height = f_height, units = "mm", dpi = 300, limitsize = FALSE)

message("F04 Panel C done")
