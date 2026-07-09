# Supplementary Enrichment Gallery -- UpSet Intersection Plot (F04 CRvH: Concordance)
# Custom ggplot dual-panel (bars + dot matrix) matching F03 Panel B style.
setwd(here::here())
source("04_Figures/F04/a_script/style.R")

pacman::p_load(dplyr, tibble, tidyr, patchwork, ComplexHeatmap)


RPT <- "04_Figures/F04/CRvH/b_reports/supp"
DAT <- "04_Figures/F04/CRvH/c_data/supp"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()
PW <- 200

conc <- readRDS(file.path(DAT, "prep_concordance.rds"))

# Build binary matrix for ComplexHeatmap
SET_NAMES <- c("Cancer vs Healthy", "Training (CR)")
SET_COLORS <- c(
  "Cancer vs Healthy" = unname(CONTRAST_COLORS["Cancer_vs_Healthy"]),
  "Training (CR)"     = unname(CONTRAST_COLORS["Training_CR"])
)

bin_mat <- cbind(
  `Cancer vs Healthy` = as.integer(conc$sig_CvH),
  `Training (CR)`     = as.integer(conc$sig_TR)
)
rownames(bin_mat) <- conc$pathway

cm <- make_comb_mat(bin_mat, mode = "intersect")
cs <- comb_size(cm)
keep <- cs > 0 & comb_degree(cm) > 0
cm_sub <- cm[keep]

comb_names_vec    <- comb_name(cm_sub)
n_comb            <- length(comb_names_vec)
set_names_ordered <- set_name(cm_sub)
counts            <- comb_size(cm_sub)

# Order by size descending
comb_ord       <- order(-counts)
counts_ord     <- counts[comb_ord]
comb_names_ord <- comb_names_vec[comb_ord]
n_int          <- length(comb_ord)

comb_deg_ord <- vapply(comb_names_ord, function(cn) {
  sum(as.integer(strsplit(cn, "")[[1]]))
}, integer(1))

set_y_levels <- rev(SET_NAMES)

# Bar data -- single bar per intersection
bar_df <- tibble(
  x         = seq_len(n_int),
  count     = counts_ord,
  is_single = comb_deg_ord == 1,
  fill      = vapply(seq_len(n_int), function(i) {
    bits <- as.logical(as.integer(strsplit(comb_names_ord[i], "")[[1]]))
    if (sum(bits) == 1) unname(SET_COLORS[set_names_ordered[bits]]) else "grey40"
  }, character(1))
)

# Dot matrix data
dot_df <- expand_grid(x = seq_len(n_int), set = SET_NAMES)
dot_df$active <- vapply(seq_len(nrow(dot_df)), function(r) {
  bits <- as.integer(strsplit(comb_names_ord[dot_df$x[r]], "")[[1]])
  as.logical(bits[match(dot_df$set[r], set_names_ordered)])
}, logical(1))
dot_df$set  <- factor(dot_df$set, levels = set_y_levels)
dot_df$ynum <- as.numeric(dot_df$set)

# Connecting segments for multi-set intersections
seg_list <- vector("list", n_int)
for (i in seq_len(n_int)) {
  bits <- as.logical(as.integer(strsplit(comb_names_ord[i], "")[[1]]))
  ypos <- match(set_names_ordered[bits], set_y_levels)
  if (length(ypos) > 1)
    seg_list[[i]] <- tibble(x = i, ymin = min(ypos), ymax = max(ypos))
}
seg_df <- bind_rows(seg_list)

# Stripe fills for dot matrix rows
stripe_fills <- adjustcolor(unname(SET_COLORS[set_y_levels]), alpha.f = 0.20)

# Background rects for single-set bars
bar_bg_list <- lapply(seq_len(n_int), function(i) {
  bits <- as.logical(as.integer(strsplit(comb_names_ord[i], "")[[1]]))
  if (sum(bits) != 1) return(NULL)
  tibble(xmin = i - 0.5, xmax = i + 0.5,
         fill = unname(SET_COLORS[set_names_ordered[bits]]))
})
bar_bg <- bind_rows(Filter(Negate(is.null), bar_bg_list))

lbl_sz <- scale_text(BASE_COUNT, PW)

p_bars <- ggplot(bar_df, aes(x, count)) +
  {if (nrow(bar_bg) > 0)
    geom_rect(data = bar_bg,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = bar_bg$fill, alpha = 0.20,
              color = "grey70", linewidth = 0.2, inherit.aes = FALSE)} +
  geom_col(fill = bar_df$fill, color = "black", linewidth = 0.3, width = 0.6) +
  geom_text(data = bar_df |> filter(count > 0, is_single),
            aes(label = count, y = count / 2),
            size = lbl_sz - 0.8, color = "white", fontface = "bold") +
  geom_text(data = bar_df |> filter(count > 0, !is_single),
            aes(label = count, y = count + max(counts_ord) * 0.03),
            size = lbl_sz - 0.8, color = "black", fontface = "bold") +
  scale_x_continuous(expand = expansion(add = 0.3)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(y = "Intersection\nsize",
       title = "Pathway Overlap (UpSet)",
       subtitle = "Cancer Recovery Concordance: sig pathways across contrasts") +
  FIG_THEME +
  theme(axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        plot.margin = margin(2, 5.5, 0, 5.5))

n_sets <- length(SET_NAMES)
p_dots <- ggplot() +
  annotate("rect",
           xmin = 0.4, xmax = n_int + 0.6,
           ymin = seq_along(set_y_levels) - 0.40,
           ymax = seq_along(set_y_levels) + 0.40,
           fill = stripe_fills) +
  {if (nrow(seg_df) > 0)
    geom_segment(data = seg_df,
                 aes(x = x, xend = x, y = ymin, yend = ymax),
                 linewidth = 0.6, color = "grey25")} +
  geom_point(data = dot_df |> filter(!active),
             aes(x = x, y = ynum), color = "grey78", size = 1.6) +
  geom_point(data = dot_df |> filter(active),
             aes(x = x, y = ynum), color = "grey15", size = 1.6) +
  scale_x_continuous(expand = expansion(add = 0.3)) +
  scale_y_continuous(breaks = seq_along(set_y_levels), labels = set_y_levels,
                     expand = expansion(add = 0)) +
  coord_cartesian(ylim = c(0.55, n_sets + 0.45)) +
  labs(x = NULL, y = NULL) +
  FIG_THEME +
  theme(axis.text.x  = element_blank(), axis.ticks = element_blank(),
        panel.grid   = element_blank(),
        panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.3),
        axis.text.y  = element_text(size = FIG_AXIS_TEXT - 1, face = "bold",
                                     margin = margin(r = 1)),
        plot.margin  = margin(0, 5.5, 0, 5.5))

p_combined <- (p_bars / p_dots) + plot_layout(heights = c(0.80, 0.20))

ggsave(file.path(RPT, "a_upset_SUPP.pdf"), p_combined,
       width = PW, height = 140, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "a_upset_SUPP.png"), p_combined,
       width = PW, height = 140, units = "mm", dpi = 300)

cat("UpSet plot saved.\n")
