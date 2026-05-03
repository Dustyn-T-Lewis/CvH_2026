# F06 CRvH Panel B: CvH Concordance Pattern Heatmap
# Layout: [legends] [sig|heatmap|quad] [sankey] [stacked bars in striped bg]
# 3 groups: Concordant Up, Concordant Down, Discordant
# Classifies DEPs by Cancer_vs_Healthy and Training_CR significance
# Outputs: panel_B_pattern_heatmap_MAIN.pdf/.png
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F06/a_script/style.R")
source("04_Figures/shared/go_slim_categories.R")

library(tidyverse)

RPT <- "04_Figures/F06/CRvH/b_reports"
DAT <- "04_Figures/F06/CRvH/c_data"
dir.create(file.path(DAT, "panel_B_heatmap"), recursive = TRUE, showWarnings = FALSE)
dir.create(RPT, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# =============================================================================
# 1. LOAD & CLASSIFY -- 3 groups
# =============================================================================
dep_df <- read_csv("03_DEP/c_data/03_combined_results_CRvH.csv", show_col_types = FALSE)

sig_df <- dep_df %>%
  filter(!is.na(logFC_Cancer_vs_Healthy), !is.na(logFC_Training_CR)) %>%
  filter(pi_score_Cancer_vs_Healthy < 0.05 | pi_score_Training_CR < 0.05) %>%
  mutate(
    quadrant = case_when(
      logFC_Cancer_vs_Healthy > 0 & logFC_Training_CR > 0 ~ "Concordant Up",
      logFC_Cancer_vs_Healthy < 0 & logFC_Training_CR < 0 ~ "Concordant Down",
      TRUE ~ "Discordant"
    ),
    sig_cat = case_when(
      pi_score_Cancer_vs_Healthy < 0.05 & pi_score_Training_CR < 0.05 ~ "Both",
      pi_score_Cancer_vs_Healthy < 0.05 ~ "CvH",
      pi_score_Training_CR < 0.05       ~ "Tr.(CR)",
      TRUE ~ "NS"
    )
  )

QUAD_ORDER <- c("Concordant Up", "Concordant Down", "Discordant")
QUAD_COLORS <- c("Concordant Up" = "#D32F2F", "Concordant Down" = "#1976D2",
                 "Discordant" = "#388E3C")
QUAD_BG <- c("Concordant Up" = "#FFCDD2", "Concordant Down" = "#BBDEFB",
             "Discordant" = "#C8E6C9", "Tied" = "#EEEEEE")
ENDPOINT_COLORS <- c("Concordant Up" = "#8B0000", "Concordant Down" = "#0D47A1",
                     "Discordant" = "#1B5E20")

QUAD_LABELS <- c("Concordant Up"   = "Concordant:\nCvH\u2191 Tr.(CR)\u2191",
                 "Concordant Down" = "Concordant:\nCvH\u2193 Tr.(CR)\u2193",
                 "Discordant"      = "Discordant")

SIG_COLORS <- c("Both" = "#2E7D32", "CvH" = "#4CAF50",
                "Tr.(CR)" = "#9C27B0", "NS" = "grey70")

go_result <- assign_go_slim_consolidated(sig_df$gene, dep_df$gene)
sig_df <- sig_df %>%
  left_join(go_result %>% select(gene, consolidated), by = "gene") %>%
  mutate(pathway = ifelse(is.na(consolidated), "Other", as.character(consolidated)))

sig_df <- sig_df %>%
  mutate(quadrant = factor(quadrant, levels = QUAD_ORDER)) %>%
  arrange(quadrant, pathway, desc(logFC_Cancer_vs_Healthy))

n_total <- nrow(sig_df)
message(sprintf("  %d significant proteins across %d quadrants", n_total,
                n_distinct(sig_df$quadrant)))

# =============================================================================
# 2. Y-COORDINATE LAYOUT
# =============================================================================
ROW_H <- 0.08

quad_counts <- sig_df %>% count(quadrant, .drop = FALSE) %>%
  mutate(quadrant = factor(quadrant, levels = QUAD_ORDER)) %>% arrange(quadrant)

y_pos <- numeric(n_total)
quad_starts <- numeric(nrow(quad_counts))
quad_ends   <- numeric(nrow(quad_counts))
idx <- 1; current_y <- 0

for (q in seq_len(nrow(quad_counts))) {
  nq <- quad_counts$n[q]
  quad_starts[q] <- current_y
  if (nq > 0) for (p in seq_len(nq)) {
    y_pos[idx] <- current_y + (p - 0.5) * ROW_H
    idx <- idx + 1
  }
  quad_ends[q] <- current_y + nq * ROW_H
  current_y <- current_y + nq * ROW_H
}
total_h <- current_y
sig_df$y <- y_pos
names(quad_starts) <- QUAD_ORDER
names(quad_ends)   <- QUAD_ORDER

# =============================================================================
# 3. PATHWAY LAYOUT
# =============================================================================
pw_counts <- sig_df %>%
  filter(pathway != "Other") %>%
  count(pathway, name = "n_prot") %>%
  arrange(desc(n_prot)) %>%
  filter(n_prot >= 2)
n_pw <- nrow(pw_counts)

row_height <- total_h / n_pw
pw_counts$y_center <- row_height * (seq_len(n_pw) - 0.5)
pw_counts$y_top    <- row_height * (seq_len(n_pw) - 1)
pw_counts$y_bot    <- row_height * seq_len(n_pw)
BAR_H <- row_height * 0.60

dom_quad <- sig_df %>%
  filter(pathway %in% pw_counts$pathway) %>%
  group_by(pathway, quadrant) %>%
  summarise(n = n(),
            lfc_sum = sum(abs(logFC_Cancer_vs_Healthy) + abs(logFC_Training_CR)),
            .groups = "drop") %>%
  group_by(pathway) %>%
  mutate(is_max = n == max(n), n_tied = sum(is_max)) %>%
  arrange(pathway, desc(n), desc(lfc_sum)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(dom_quad = as.character(quadrant)) %>%
  select(pathway, dom_quad)

pw_counts <- pw_counts %>% left_join(dom_quad, by = "pathway")

# =============================================================================
# 4. X-COORDINATE LAYOUT
# =============================================================================
STRIP_W <- 0.10; TILE_W <- 0.70

X_SIG   <- 0.8
X_COL1  <- X_SIG + STRIP_W/2 + TILE_W/2 + 0.01
X_COL2  <- X_COL1 + TILE_W + 0.01
X_QUAD  <- X_COL2 + TILE_W/2 + STRIP_W/2 + 0.01
HEAT_LEFT  <- X_SIG - STRIP_W/2
HEAT_RIGHT <- X_QUAD + STRIP_W/2

X_SANK_L  <- HEAT_RIGHT + 0.08
X_SANK_R  <- 3.2
X_BAR_L   <- 3.3
BAR_SCALE <- 0.18

count_max <- max(pw_counts$n_prot)
X_BAR_MAX <- max(X_BAR_L + 17 * BAR_SCALE, X_BAR_L + count_max * BAR_SCALE)

PW_OUT <- 210; PH_OUT <- 145

# =============================================================================
# 5. STACKED BAR DATA
# =============================================================================
bar_data <- sig_df %>%
  filter(pathway %in% pw_counts$pathway) %>%
  count(pathway, quadrant, name = "n_seg") %>%
  left_join(pw_counts %>% select(pathway, y_center, n_prot), by = "pathway") %>%
  group_by(pathway) %>%
  arrange(pathway, desc(n_seg)) %>%
  mutate(
    cum_n = cumsum(n_seg) - n_seg,
    xmin = X_BAR_L + cum_n * BAR_SCALE,
    xmax = X_BAR_L + (cum_n + n_seg) * BAR_SCALE,
    ymin = y_center - BAR_H / 2,
    ymax = y_center + BAR_H / 2
  ) %>% ungroup()

bg_stripes <- pw_counts %>%
  transmute(
    xmin = X_BAR_L - 0.05, xmax = X_BAR_MAX + 0.05,
    ymin = y_top, ymax = y_bot,
    fill = QUAD_BG[dom_quad]
  )

pw_labels <- pw_counts %>%
  transmute(x = X_BAR_L + n_prot * BAR_SCALE + 0.08, y = y_center,
            label = ifelse(nchar(pathway) > 18,
                           stringr::str_wrap(pathway, width = 18), pathway))

count_ticks <- tibble(
  val = pretty(c(0, count_max), n = 4),
  x = X_BAR_L + val * BAR_SCALE,
  y_tick_top = total_h, y_tick_bot = total_h + ROW_H * 1.5,
  y_label = total_h + ROW_H * 3
) %>% filter(val >= 0, val <= count_max)

# =============================================================================
# 6. SANKEY
# =============================================================================
flow_df <- sig_df %>%
  filter(pathway %in% pw_counts$pathway) %>%
  count(quadrant, pathway, name = "n_flow") %>%
  filter(n_flow > 0)

source_bands <- flow_df %>%
  group_by(quadrant) %>%
  mutate(
    total_q = sum(n_flow), frac = n_flow / total_q,
    q_start = quad_starts[as.character(quadrant)],
    q_end   = quad_ends[as.character(quadrant)],
    q_height = q_end - q_start
  ) %>%
  arrange(quadrant, match(pathway, pw_counts$pathway)) %>%
  mutate(
    cum_frac = cumsum(frac) - frac,
    src_top = q_start + cum_frac * q_height,
    src_bot = q_start + (cum_frac + frac) * q_height
  ) %>% ungroup()

target_bands <- bar_data %>%
  group_by(pathway) %>%
  arrange(pathway, desc(n_seg)) %>%
  mutate(
    frac = n_seg / sum(n_seg),
    cum_frac = cumsum(frac) - frac,
    tgt_top = ymin + cum_frac * (ymax - ymin),
    tgt_bot = ymin + (cum_frac + frac) * (ymax - ymin)
  ) %>% ungroup() %>%
  select(pathway, quadrant, tgt_top, tgt_bot)

ribbon_df <- source_bands %>%
  select(quadrant, pathway, n_flow, src_top, src_bot) %>%
  left_join(target_bands, by = c("quadrant", "pathway"))

all_ribbons <- pmap_dfr(ribbon_df, function(quadrant, pathway, n_flow,
                                              src_top, src_bot, tgt_top, tgt_bot) {
  rid <- paste(quadrant, pathway, sep = "___")
  df <- make_sigmoid_ribbon(X_SANK_L, X_SANK_R, src_top, src_bot, tgt_top, tgt_bot,
                            n_pts = 60, ribbon_id = rid)
  df$quadrant <- quadrant; df$pathway <- pathway; df
})

endpoint_bars <- bar_data %>%
  transmute(xmin = X_SANK_R - 0.04, xmax = X_SANK_R + 0.04,
            ymin, ymax, quadrant = as.character(quadrant))

# =============================================================================
# 7. HEATMAP TILE DATA
# =============================================================================
fc_max <- max(abs(c(sig_df$logFC_Cancer_vs_Healthy, sig_df$logFC_Training_CR)), na.rm = TRUE)

lfc_to_color <- function(v, fc_max) {
  v <- pmax(-fc_max, pmin(fc_max, v))
  ifelse(v >= 0,
         scales::seq_gradient_pal("#FFFFFF", "#B2182B")(v / fc_max),
         scales::seq_gradient_pal("#2166AC", "#FFFFFF")((v + fc_max) / fc_max))
}

heat_tiles <- bind_rows(
  sig_df %>% transmute(x = X_COL1, y, w = TILE_W, h = ROW_H,
                        fill = lfc_to_color(logFC_Cancer_vs_Healthy, fc_max)),
  sig_df %>% transmute(x = X_COL2, y, w = TILE_W, h = ROW_H,
                        fill = lfc_to_color(logFC_Training_CR, fc_max))
)
sig_tiles <- sig_df %>%
  transmute(x = X_SIG, y, w = STRIP_W, h = ROW_H, fill = SIG_COLORS[sig_cat])
quad_tiles <- sig_df %>%
  transmute(x = X_QUAD, y, w = STRIP_W, h = ROW_H,
            fill = QUAD_COLORS[as.character(quadrant)])

divider_ys <- quad_ends[1:(length(QUAD_ORDER) - 1)]
divider_ys <- divider_ys[divider_ys > 0 & divider_ys < total_h]

col_headers <- tibble(x = c(X_COL1, X_COL2), y = -ROW_H * 2,
                      label = c("CvH", "Tr.(CR)"))

# =============================================================================
# 8. LEGENDS
# =============================================================================
LEG_X <- -0.5

n_g <- 50
grad_xs <- seq(HEAT_LEFT, HEAT_RIGHT, length.out = n_g)
grad_h_legend <- tibble(
  xmin = grad_xs,
  xmax = lead(grad_xs, default = max(grad_xs) + diff(grad_xs)[1]),
  fv = seq(-fc_max, fc_max, length.out = n_g),
  fill = lfc_to_color(seq(-fc_max, fc_max, length.out = n_g), fc_max)
)
GRAD_Y <- total_h + ROW_H * 2

FONT_UNI <- 2.8
FONT_BAR <- 2.6

sig_leg_cats <- c("Both", "CvH", "Tr.(CR)")
sig_leg_y0 <- total_h * 0.07
SQ_MM <- 2.5
plot_x_range <- (X_BAR_MAX + 1.0) - (-0.8)
plot_y_range <- (total_h + ROW_H * 12) - (-ROW_H * 4)
SQ_X <- SQ_MM * plot_x_range / PW_OUT
SQ_Y <- SQ_MM * plot_y_range / PH_OUT

sig_leg <- tibble(cat = sig_leg_cats,
                  sq_left = LEG_X,
                  sq_right = LEG_X + SQ_X,
                  x_center = LEG_X + SQ_X / 2,
                  y = sig_leg_y0 + seq(0, by = ROW_H * 5, length.out = length(sig_leg_cats)))

quad_leg_y0 <- total_h * 0.50
quad_leg <- tibble(
  cat = QUAD_ORDER,
  display = QUAD_LABELS[QUAD_ORDER],
  sq_left = LEG_X,
  sq_right = LEG_X + SQ_X,
  x_center = LEG_X + SQ_X / 2,
  y = quad_leg_y0 + seq(0, by = ROW_H * 8, length.out = length(QUAD_ORDER))
)

# =============================================================================
# 9. RENDER
# =============================================================================
p <- ggplot() +
  geom_rect(data = bg_stripes,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = bg_stripes$fill, color = "grey70", linewidth = 0.2) +
  geom_rect(data = heat_tiles,
            aes(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2),
            fill = heat_tiles$fill, color = NA) +
  geom_rect(data = sig_tiles,
            aes(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2),
            fill = sig_tiles$fill, color = NA) +
  geom_rect(data = quad_tiles,
            aes(xmin = x - w/2, xmax = x + w/2, ymin = y - h/2, ymax = y + h/2),
            fill = quad_tiles$fill, color = NA) +
  geom_segment(data = tibble(y = divider_ys),
               aes(x = X_SIG - STRIP_W/2, xend = X_QUAD + STRIP_W/2,
                   y = y, yend = y),
               color = "grey30", linewidth = 0.4) +
  geom_text(data = col_headers, aes(x = x, y = y, label = label),
            size = FONT_UNI, fontface = "bold", color = "grey20") +
  geom_polygon(data = all_ribbons, aes(x = x, y = y, group = ribbon_id),
               fill = QUAD_COLORS[all_ribbons$quadrant], alpha = 0.40, color = NA) +
  geom_rect(data = endpoint_bars,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = ENDPOINT_COLORS[endpoint_bars$quadrant], color = NA) +
  geom_rect(data = bar_data,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = QUAD_COLORS[as.character(bar_data$quadrant)],
            color = "black", linewidth = 0.3) +
  geom_text(data = bar_data,
            aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = n_seg),
            size = FONT_BAR, fontface = "bold", color = "white") +
  geom_text(data = pw_labels, aes(x = x, y = y, label = label),
            size = FONT_UNI, hjust = 0, fontface = "bold", color = "grey15",
            lineheight = 0.8) +
  annotate("segment", x = X_BAR_L, xend = X_BAR_MAX,
           y = total_h, yend = total_h, color = "grey20", linewidth = 0.5) +
  geom_segment(data = count_ticks,
               aes(x = x, xend = x, y = y_tick_top, yend = y_tick_bot),
               color = "grey20", linewidth = 0.3) +
  geom_text(data = count_ticks, aes(x = x, y = y_label, label = val),
            size = FONT_UNI, fontface = "bold", color = "grey20") +
  annotate("text", x = (X_BAR_L + X_BAR_MAX) / 2, y = total_h + ROW_H * 6,
           label = "Protein count", size = FONT_UNI, fontface = "bold", color = "grey20") +
  # logFC gradient horizontal
  geom_rect(data = grad_h_legend,
            aes(xmin = xmin, xmax = xmax, ymin = GRAD_Y, ymax = GRAD_Y + ROW_H * 2),
            fill = grad_h_legend$fill, color = NA) +
  annotate("text", x = HEAT_LEFT, y = GRAD_Y + ROW_H * 4,
           label = sprintf("%.1f", -fc_max), size = FONT_UNI, fontface = "bold",
           color = "grey20", hjust = 0) +
  annotate("text", x = HEAT_RIGHT, y = GRAD_Y + ROW_H * 4,
           label = sprintf("+%.1f", fc_max), size = FONT_UNI, fontface = "bold",
           color = "grey20", hjust = 1) +
  annotate("text", x = (HEAT_LEFT + HEAT_RIGHT) / 2, y = GRAD_Y + ROW_H * 4,
           label = "0", size = FONT_UNI * 0.7, fontface = "bold", color = "grey30") +
  annotate("text", x = (HEAT_LEFT + HEAT_RIGHT) / 2, y = GRAD_Y + ROW_H * 6.5,
           label = "logFC", size = FONT_UNI, fontface = "bold", color = "grey15") +
  # Sig legend
  annotate("text", x = LEG_X, y = sig_leg_y0 - ROW_H * 6,
           label = "Significance", size = FONT_UNI, fontface = "bold",
           color = "grey15", hjust = 0) +
  geom_rect(data = sig_leg,
            aes(xmin = sq_left, xmax = sq_right,
                ymin = y - SQ_Y/2, ymax = y + SQ_Y/2),
            fill = SIG_COLORS[sig_leg$cat], color = "black", linewidth = 0.3) +
  geom_text(data = sig_leg, aes(x = sq_right + 0.06, y = y, label = cat),
            size = FONT_UNI, fontface = "bold", hjust = 0, color = "grey20") +
  # Quad legend
  annotate("text", x = LEG_X, y = quad_leg_y0 - ROW_H * 6,
           label = "Quadrant", size = FONT_UNI, fontface = "bold",
           color = "grey15", hjust = 0) +
  geom_rect(data = quad_leg,
            aes(xmin = sq_left, xmax = sq_right,
                ymin = y - SQ_Y/2, ymax = y + SQ_Y/2),
            fill = QUAD_COLORS[quad_leg$cat], color = "black", linewidth = 0.3) +
  geom_text(data = quad_leg, aes(x = sq_right + 0.06, y = y, label = display),
            size = FONT_UNI, fontface = "bold", hjust = 0, color = "grey20",
            lineheight = 0.75) +
  scale_y_reverse() +
  coord_cartesian(xlim = c(-0.8, X_BAR_MAX + 1.0),
                  ylim = c(total_h + ROW_H * 12, -ROW_H * 4),
                  expand = FALSE) +
  labs(title = "CRvH Concordance Patterns",
       subtitle = sprintf("%d DEPs | GO Slim | %d pathways", n_total, n_pw)) +
  theme_void() +
  theme(plot.margin = margin(3, 3, 3, 3, "mm"),
        plot.title = element_text(face = "bold", size = 12, hjust = 0),
        plot.subtitle = element_text(face = "italic", size = 9,
                                     hjust = 0, color = "grey40"))

# =============================================================================
# 10. SAVE
# =============================================================================
ggsave(file.path(RPT, "panel_B_pattern_heatmap_MAIN.pdf"), p,
       width = PW_OUT, height = PH_OUT, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_B_pattern_heatmap_MAIN.png"), p,
       width = PW_OUT, height = PH_OUT, units = "mm", dpi = 300)

# =============================================================================
# 11. DATA EXPORTS
# =============================================================================
sig_df %>%
  transmute(gene, quadrant = as.character(quadrant), sig_cat, pathway,
            logFC_Cancer_vs_Healthy = round(logFC_Cancer_vs_Healthy, 4),
            logFC_Training_CR = round(logFC_Training_CR, 4)) %>%
  write_csv(file.path(DAT, "panel_B_heatmap", "pattern_classification.csv"))
flow_df %>% write_csv(file.path(DAT, "panel_B_heatmap", "sankey_links.csv"))
bar_data %>%
  select(pathway, quadrant, n_seg, xmin, xmax) %>%
  write_csv(file.path(DAT, "panel_B_heatmap", "bar_data.csv"))

message("F06 CRvH Panel B (pattern heatmap) done")
