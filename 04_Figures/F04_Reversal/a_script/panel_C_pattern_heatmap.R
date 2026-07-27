# Reversal Panel C: Pattern Heatmap + Sankey
# Per-protein reversal classification with GO Slim pathway bars
# Source LAST — AnnotationDbi masks dplyr::select
setwd(here::here())
source("04_Figures/shared/style.R")
source("04_Figures/shared/go_slim_categories.R")

pacman::p_load(tidyverse)

RPT_PNG <- "04_Figures/F04_Reversal/b_reports/supp/png"
RPT_PDF <- "04_Figures/F04_Reversal/b_reports/supp/pdf"
DAT <- "04_Figures/F04_Reversal/c_data"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_C_heatmap"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# 1. LOAD & CLASSIFY
source("04_Figures/F04_Reversal/a_script/f04_data.R") # dep_df (old column names)

sig_df <- dep_df %>%
  filter(pi_score_CRvH_Baseline < 0.05 | pi_score_CR_Training < 0.05) %>%
  mutate(
    quadrant = case_when(
      logFC_CRvH_Baseline > 0 & logFC_CR_Training < 0 ~ "Reversed Up",
      logFC_CRvH_Baseline < 0 & logFC_CR_Training > 0 ~ "Reversed Down",
      TRUE ~ "Non-reversed"
    ),
    sig_cat = case_when(
      pi_score_CRvH_Baseline < 0.05 & pi_score_CR_Training < 0.05 ~ "Both",
      pi_score_CRvH_Baseline < 0.05 ~ "Cancer",
      pi_score_CR_Training < 0.05 ~ "Tr.(CR)",
      TRUE ~ "NS"
    )
  )

QUAD_ORDER <- c("Reversed Up", "Reversed Down", "Non-reversed")
QUAD_COLORS <- c(
  "Reversed Up" = "#2563EB", "Reversed Down" = "#DC2626",
  "Non-reversed" = "#FFB74D"
)
QUAD_BG <- c(
  "Reversed Up" = scales::alpha("#2563EB", 0.08),
  "Reversed Down" = scales::alpha("#DC2626", 0.08),
  "Non-reversed" = scales::alpha("#FFB74D", 0.08)
)
ENDPOINT_COLORS <- QUAD_COLORS
SIG_COLORS_HM <- c(
  "Both" = "#2E7D32", "Cancer" = "#4CAF50",
  "Tr.(CR)" = "#9C27B0", "NS" = "grey70"
)

# GO Slim assignment
go_result <- assign_go_slim_consolidated(sig_df$gene, dep_df$gene)
sig_df <- sig_df %>%
  left_join(go_result %>% dplyr::select(gene, consolidated), by = "gene") %>%
  mutate(pathway = ifelse(is.na(consolidated), "Other", as.character(consolidated)))

sig_df <- sig_df %>%
  mutate(quadrant = factor(quadrant, levels = QUAD_ORDER)) %>%
  arrange(quadrant, pathway, desc(logFC_CRvH_Baseline))

n_total <- nrow(sig_df)
message(sprintf(
  "  %d significant proteins across %d quadrants", n_total,
  n_distinct(sig_df$quadrant)
))

# 2. Y-COORDINATE LAYOUT
ROW_H <- 0.078

quad_counts <- sig_df %>%
  count(quadrant, .drop = FALSE) %>%
  mutate(quadrant = factor(quadrant, levels = QUAD_ORDER)) %>%
  arrange(quadrant)

y_pos <- numeric(n_total)
quad_starts <- numeric(nrow(quad_counts))
quad_ends <- numeric(nrow(quad_counts))
idx <- 1
current_y <- 0

for (q in seq_len(nrow(quad_counts))) {
  nq <- quad_counts$n[q]
  quad_starts[q] <- current_y
  if (nq > 0) {
    for (pp in seq_len(nq)) {
      y_pos[idx] <- current_y + (pp - 0.5) * ROW_H
      idx <- idx + 1
    }
  }
  quad_ends[q] <- current_y + nq * ROW_H
  current_y <- current_y + nq * ROW_H
}
total_h <- current_y
sig_df$y <- y_pos
names(quad_starts) <- QUAD_ORDER
names(quad_ends) <- QUAD_ORDER

BAR_FRAC <- 1.0
BAR_YMIN <- 0
BAR_YMAX <- total_h * BAR_FRAC

# 3. PATHWAY LAYOUT
pw_counts <- sig_df %>%
  filter(pathway != "Other") %>%
  count(pathway, name = "n_prot") %>%
  arrange(desc(n_prot)) %>%
  filter(n_prot >= 2)
n_pw <- nrow(pw_counts)

row_height <- (BAR_YMAX - BAR_YMIN) / n_pw
pw_counts$y_center <- BAR_YMIN + row_height * (seq_len(n_pw) - 0.5)
pw_counts$y_top <- BAR_YMIN + row_height * (seq_len(n_pw) - 1)
pw_counts$y_bot <- BAR_YMIN + row_height * seq_len(n_pw)
BAR_H <- row_height * 0.78

dom_quad <- sig_df %>%
  filter(pathway %in% pw_counts$pathway) %>%
  group_by(pathway, quadrant) %>%
  summarise(
    n = n(),
    lfc_sum = sum(abs(logFC_CRvH_Baseline) + abs(logFC_CR_Training)),
    .groups = "drop"
  ) %>%
  group_by(pathway) %>%
  mutate(is_max = n == max(n), n_tied = sum(is_max)) %>%
  arrange(pathway, desc(n), desc(lfc_sum)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(dom_quad = as.character(quadrant)) %>%
  dplyr::select(pathway, dom_quad)

pw_counts <- pw_counts %>% left_join(dom_quad, by = "pathway")

# 4. X-COORDINATE LAYOUT
STRIP_W <- 0.10
TILE_W <- 0.70

X_SIG <- 0.8
X_COL1 <- X_SIG + STRIP_W / 2 + TILE_W / 2 + 0.01
X_COL2 <- X_COL1 + TILE_W + 0.01
X_QUAD <- X_COL2 + TILE_W / 2 + STRIP_W / 2 + 0.01
HEAT_RIGHT <- X_QUAD + STRIP_W / 2

X_SANK_L <- HEAT_RIGHT + 0.08
X_SANK_R <- 3.2
X_BAR_L <- 3.3
BAR_SCALE <- 0.055

count_max <- max(pw_counts$n_prot)
X_BAR_MAX <- max(X_BAR_L + 40 * BAR_SCALE, X_BAR_L + count_max * BAR_SCALE)

PW_OUT <- 178
PH_OUT <- 130

# 5. STACKED BAR DATA
bar_data <- sig_df %>%
  filter(pathway %in% pw_counts$pathway) %>%
  count(pathway, quadrant, name = "n_seg") %>%
  left_join(pw_counts %>% dplyr::select(pathway, y_center, n_prot), by = "pathway") %>%
  group_by(pathway) %>%
  arrange(pathway, desc(n_seg)) %>%
  mutate(
    cum_n = cumsum(n_seg) - n_seg,
    xmin = X_BAR_L + cum_n * BAR_SCALE,
    xmax = X_BAR_L + (cum_n + n_seg) * BAR_SCALE,
    ymin = y_center - BAR_H / 2,
    ymax = y_center + BAR_H / 2
  ) %>%
  ungroup()

bg_stripes <- pw_counts %>%
  transmute(
    xmin = X_BAR_L - 0.05, xmax = X_BAR_MAX + 0.05,
    ymin = y_top, ymax = y_bot,
    fill = QUAD_BG[dom_quad]
  )

pw_labels <- pw_counts %>%
  transmute(
    x = X_BAR_L + n_prot * BAR_SCALE + 0.08, y = y_center,
    label = pathway
  )

count_ticks <- tibble(
  val = pretty(c(0, count_max), n = 4),
  x = X_BAR_L + val * BAR_SCALE,
  y_tick_top = BAR_YMAX, y_tick_bot = BAR_YMAX + ROW_H * 1.6,
  y_label = BAR_YMAX + ROW_H * 3.5
) %>% filter(val >= 0, val <= count_max)

# 6. SANKEY
flow_df <- sig_df %>%
  filter(pathway %in% pw_counts$pathway) %>%
  count(quadrant, pathway, name = "n_flow") %>%
  filter(n_flow > 0)

source_bands <- flow_df %>%
  group_by(quadrant) %>%
  mutate(
    total_q = sum(n_flow), frac = n_flow / total_q,
    q_start = quad_starts[as.character(quadrant)],
    q_end = quad_ends[as.character(quadrant)],
    q_height = q_end - q_start
  ) %>%
  arrange(quadrant, match(pathway, pw_counts$pathway)) %>%
  mutate(
    cum_frac = cumsum(frac) - frac,
    src_top = q_start + cum_frac * q_height,
    src_bot = q_start + (cum_frac + frac) * q_height
  ) %>%
  ungroup()

target_bands <- bar_data %>%
  group_by(pathway) %>%
  arrange(pathway, desc(n_seg)) %>%
  mutate(
    frac = n_seg / sum(n_seg),
    cum_frac = cumsum(frac) - frac,
    tgt_top = ymin + cum_frac * (ymax - ymin),
    tgt_bot = ymin + (cum_frac + frac) * (ymax - ymin)
  ) %>%
  ungroup() %>%
  dplyr::select(pathway, quadrant, tgt_top, tgt_bot)

ribbon_df <- source_bands %>%
  dplyr::select(quadrant, pathway, n_flow, src_top, src_bot) %>%
  left_join(target_bands, by = c("quadrant", "pathway"))

all_ribbons <- pmap_dfr(ribbon_df, function(quadrant, pathway, n_flow,
                                            src_top, src_bot, tgt_top, tgt_bot) {
  rid <- paste(quadrant, pathway, sep = "___")
  df <- make_sigmoid_ribbon(X_SANK_L, X_SANK_R, src_top, src_bot, tgt_top, tgt_bot,
    n_pts = 60, ribbon_id = rid
  )
  df$quadrant <- quadrant
  df$pathway <- pathway
  df
})

endpoint_bars <- bar_data %>%
  transmute(
    xmin = X_SANK_R - 0.04, xmax = X_SANK_R + 0.04,
    ymin, ymax, quadrant = as.character(quadrant)
  )

# 7. HEATMAP
fc_max <- max(abs(c(sig_df$logFC_CRvH_Baseline, sig_df$logFC_CR_Training)),
  na.rm = TRUE
)

lfc_to_color <- function(v, fc_max) {
  v <- pmax(-fc_max, pmin(fc_max, v))
  ifelse(v >= 0,
    scales::seq_gradient_pal("#FFFFFF", "#B2182B")(v / fc_max),
    scales::seq_gradient_pal("#2166AC", "#FFFFFF")((v + fc_max) / fc_max)
  )
}

heat_tiles <- bind_rows(
  sig_df %>% transmute(
    x = X_COL1, y, w = TILE_W, h = ROW_H,
    fill = lfc_to_color(logFC_CRvH_Baseline, fc_max)
  ),
  sig_df %>% transmute(
    x = X_COL2, y, w = TILE_W, h = ROW_H,
    fill = lfc_to_color(logFC_CR_Training, fc_max)
  )
)
sig_tiles <- sig_df %>%
  transmute(x = X_SIG, y, w = STRIP_W, h = ROW_H, fill = SIG_COLORS_HM[sig_cat])
quad_tiles <- sig_df %>%
  transmute(
    x = X_QUAD, y, w = STRIP_W, h = ROW_H,
    fill = QUAD_COLORS[as.character(quadrant)]
  )

divider_ys <- quad_ends[1:(length(QUAD_ORDER) - 1)]
divider_ys <- divider_ys[divider_ys > 0 & divider_ys < total_h]

col_headers <- tibble(
  x = c(X_COL1, X_COL2), y = total_h + ROW_H * 2.2,
  label = c("Cancer", "Tr.(CR)"),
  color = unname(CONTRAST_COLORS[c("CRvH_Baseline", "CR_Training")])
)


FONT_UNI <- 2.5
FONT_BAR <- 2.0
FONT_PW <- 2.2

# 9. RENDER
p <- ggplot() +
  geom_rect(
    data = bg_stripes,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = bg_stripes$fill, color = "grey70", linewidth = 0.2
  ) +
  geom_rect(
    data = heat_tiles,
    aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2),
    fill = heat_tiles$fill, color = NA
  ) +
  geom_rect(
    data = sig_tiles,
    aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2),
    fill = sig_tiles$fill, color = NA
  ) +
  geom_rect(
    data = quad_tiles,
    aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2),
    fill = quad_tiles$fill, color = NA
  ) +
  geom_segment(
    data = tibble(y = divider_ys),
    aes(
      x = X_SIG - STRIP_W / 2, xend = X_QUAD + STRIP_W / 2,
      y = y, yend = y
    ),
    color = "grey30", linewidth = 0.4
  ) +
  geom_text(
    data = col_headers, aes(x = x, y = y, label = label),
    size = FONT_UNI, fontface = "bold", color = col_headers$color
  ) +
  geom_polygon(
    data = all_ribbons, aes(x = x, y = y, group = ribbon_id),
    fill = QUAD_COLORS[all_ribbons$quadrant], alpha = 0.40, color = NA
  ) +
  geom_rect(
    data = endpoint_bars,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = ENDPOINT_COLORS[endpoint_bars$quadrant], color = NA
  ) +
  geom_rect(
    data = bar_data,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = QUAD_COLORS[as.character(bar_data$quadrant)],
    color = "black", linewidth = 0.3
  ) +
  geom_text(
    data = bar_data,
    aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = n_seg),
    size = FONT_BAR, fontface = "bold", color = "white"
  ) +
  geom_text(
    data = pw_labels, aes(x = x, y = y, label = label),
    size = FONT_PW, hjust = 0, fontface = "bold", color = "grey15",
    lineheight = 0.8
  ) +
  annotate("segment",
    x = X_BAR_L, xend = X_BAR_MAX,
    y = BAR_YMAX, yend = BAR_YMAX, color = "grey20", linewidth = 0.5
  ) +
  geom_segment(
    data = count_ticks,
    aes(x = x, xend = x, y = y_tick_top, yend = y_tick_bot),
    color = "grey20", linewidth = 0.3
  ) +
  geom_text(
    data = count_ticks, aes(x = x, y = y_label, label = val),
    size = FONT_UNI * 0.78 + 0.5, fontface = "bold", color = "grey20"
  ) +
  annotate("text",
    x = X_BAR_L + (count_max / 2) * BAR_SCALE,
    y = BAR_YMAX + ROW_H * 4.4,
    label = "Protein count", size = FONT_UNI, fontface = "bold",
    color = "grey20"
  ) +
  scale_y_reverse() +
  coord_cartesian(
    xlim = c(0.0, X_BAR_MAX + 2.0),
    ylim = c(BAR_YMAX + ROW_H * 7.5, -ROW_H * 0.05),
    expand = FALSE
  ) +
  labs(
    title = "Reversal Pattern Classification",
    subtitle = sprintf("%d proteins | GO Slim | %d pathways", n_total, n_pw)
  ) +
  theme_void() +
  theme(
    plot.margin = margin(6, -30, 38, -12, "mm"),
    plot.title = element_text(
      face = "bold", size = FIG_TITLE_SIZE, hjust = 0,
      margin = margin(l = 31.5, unit = "mm")
    ),
    plot.subtitle = element_text(
      face = "italic", size = FIG_SUBTITLE_SIZE,
      hjust = 0, color = "grey40",
      margin = margin(l = 31.5, unit = "mm")
    ),
    plot.title.position = "panel"
  )

# 10. SAVE
ggsave(file.path(RPT_PNG, "SUPP_F04_pattern_heatmap.png"), p,
  width = PW_OUT, height = PH_OUT, units = "mm", dpi = 300
)
ggsave(file.path(RPT_PDF, "SUPP_F04_pattern_heatmap.pdf"), p,
  width = PW_OUT, height = PH_OUT, units = "mm", device = pdf_device
)

# 11. DATA EXPORTS
sig_df %>%
  transmute(gene,
    quadrant = as.character(quadrant), sig_cat, pathway,
    logFC_CRvH_Baseline = round(logFC_CRvH_Baseline, 4),
    logFC_CR_Training = round(logFC_CR_Training, 4)
  ) %>%
  write_csv(file.path(DAT, "panel_C_heatmap", "pattern_classification.csv"))
flow_df %>% write_csv(file.path(DAT, "panel_C_heatmap", "sankey_links.csv"))
bar_data %>%
  dplyr::select(pathway, quadrant, n_seg, xmin, xmax) %>%
  write_csv(file.path(DAT, "panel_C_heatmap", "bar_data.csv"))

p <- p + labs(title = NULL, subtitle = NULL, tag = NULL) +
  coord_cartesian(
    xlim = c(0.0, X_BAR_MAX + 2.0),
    ylim = c(BAR_YMAX + ROW_H * 7.5, -ROW_H * 0.05),
    expand = FALSE
  ) +
  theme(plot.margin = margin(2, -30, 6, -12, "mm"))

message("Reversal Panel C (pattern heatmap) done")
