# Figure 2 — Panel C (CRvH): UpSet Plot (Dual Bar + Dot Matrix)
# Outputs: pC_crvh_bars, pC_crvh_dots (ggplot objects)
# STANDALONE: generates sig_sets_crvh, dir_map_crvh, all_genes_crvh locally (no panel_A dependency)

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
  library(patchwork)
  library(ComplexHeatmap)
  library(purrr)
})

CONTRASTS <- c("Cancer_vs_Healthy", "Training_CR")

# ============================================================================
# Generate required data locally (formerly from Panel A)
# ============================================================================
DEP_FILE <- "03_DEP/c_data/03_combined_results_CRvH.csv"
dep_df <- read_csv(DEP_FILE, show_col_types = FALSE)

SET_LABELS_crvh <- CTR_SHORT[CONTRASTS]
all_genes_crvh  <- unique(dep_df$gene[!is.na(dep_df$gene)])

sig_sets_crvh <- list()
dir_map_crvh  <- list()

for (ctr in CONTRASTS) {
  pi_vals  <- dep_df[[paste0("pi_score_", ctr)]]
  lfc_vals <- dep_df[[paste0("logFC_", ctr)]]
  is_sig   <- !is.na(pi_vals) & pi_vals < 0.05
  sig_sets_crvh[[ctr]] <- dep_df$gene[is_sig]
  dir_map_crvh[[ctr]]  <- setNames(ifelse(lfc_vals[is_sig] > 0, "Up", "Down"),
                                    dep_df$gene[is_sig])
}

pi_total_crvh  <- sum(sapply(sig_sets_crvh, length))
fdr_total_crvh <- sum(sapply(CONTRASTS, function(ctr) {
  fdr_col <- paste0("adj.P.Val_", ctr)
  if (fdr_col %in% names(dep_df)) sum(dep_df[[fdr_col]] < 0.10, na.rm = TRUE) else 0
}))

SET_DISPLAY_COLORS_crvh <- setNames(
  unname(CONTRAST_COLORS[CONTRASTS]),
  unname(SET_LABELS_crvh[CONTRASTS])
)

RPT <- "04_Figures/F02/b_reports"
DAT <- "04_Figures/F02/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()
PG_W <- 170

bin_mat <- sapply(sig_sets_crvh, function(s) as.integer(all_genes_crvh %in% s))
rownames(bin_mat) <- all_genes_crvh
colnames(bin_mat) <- SET_LABELS_crvh[colnames(bin_mat)]
label_to_contrast <- setNames(names(SET_LABELS_crvh), SET_LABELS_crvh)

cm <- make_comb_mat(bin_mat, mode = "intersect")
cs <- comb_size(cm)
keep <- cs > 0 & comb_degree(cm) > 0
cm_sub <- cm[keep]

comb_names_vec    <- comb_name(cm_sub)
n_comb            <- length(comb_names_vec)
set_names_ordered <- set_name(cm_sub)
up_counts   <- numeric(n_comb)
down_counts <- numeric(n_comb)

for (i in seq_len(n_comb)) {
  members <- extract_comb(cm_sub, comb_names_vec[i])
  if (length(members) == 0) next
  bits <- as.logical(as.integer(strsplit(comb_names_vec[i], "")[[1]]))
  active_contrasts <- label_to_contrast[set_names_ordered[bits]]
  gene_dirs <- sapply(members, function(g) {
    dirs <- sapply(active_contrasts, function(ctr) {
      if (g %in% names(dir_map_crvh[[ctr]])) dir_map_crvh[[ctr]][g] else NA
    })
    dirs <- dirs[!is.na(dirs)]
    if (all(dirs == "Up")) "Up" else if (all(dirs == "Down")) "Down" else "Mixed"
  })
  up_counts[i]   <- sum(gene_dirs == "Up")
  down_counts[i] <- sum(gene_dirs == "Down")
}

n_mixed <- sum(sapply(seq_len(n_comb), function(i) {
  members <- extract_comb(cm_sub, comb_names_vec[i])
  if (length(members) == 0) return(0L)
  bits <- as.logical(as.integer(strsplit(comb_names_vec[i], "")[[1]]))
  active_contrasts <- label_to_contrast[set_names_ordered[bits]]
  gene_dirs <- sapply(members, function(g) {
    dirs <- sapply(active_contrasts, function(ctr) {
      if (g %in% names(dir_map_crvh[[ctr]])) dir_map_crvh[[ctr]][g] else NA
    })
    dirs <- dirs[!is.na(dirs)]
    if (all(dirs == "Up")) "Up" else if (all(dirs == "Down")) "Down" else "Mixed"
  })
  sum(gene_dirs == "Mixed")
}))

display_total <- up_counts + down_counts
keep_display  <- display_total > 0
comb_ord <- which(keep_display)[order(-display_total[keep_display])]
up_ord   <- up_counts[comb_ord]
down_ord <- down_counts[comb_ord]

set_display_order <- unname(SET_LABELS_crvh[CONTRASTS])
n_int          <- length(comb_ord)
comb_names_ord <- comb_names_vec[comb_ord]
set_order_ch   <- set_name(cm_sub)
set_y_levels   <- rev(set_display_order)

comb_deg_ord <- vapply(comb_names_ord, function(cn) {
  sum(as.integer(strsplit(cn, "")[[1]]))
}, integer(1))

bar_long <- tibble(
  x         = rep(seq_len(n_int), 2),
  direction = factor(rep(c("Up", "Down"), each = n_int), levels = c("Down", "Up")),
  count     = c(up_ord, down_ord),
  is_single = rep(comb_deg_ord == 1, 2)
)

dot_df <- expand_grid(x = seq_len(n_int), set = set_display_order)
dot_df$active <- vapply(seq_len(nrow(dot_df)), function(r) {
  bits <- as.integer(strsplit(comb_names_ord[dot_df$x[r]], "")[[1]])
  as.logical(bits[match(dot_df$set[r], set_order_ch)])
}, logical(1))
dot_df$set  <- factor(dot_df$set, levels = set_y_levels)
dot_df$ynum <- as.numeric(dot_df$set)

seg_list <- vector("list", n_int)
for (i in seq_len(n_int)) {
  bits <- as.logical(as.integer(strsplit(comb_names_ord[i], "")[[1]]))
  ypos <- match(set_order_ch[bits], set_y_levels)
  if (length(ypos) > 1)
    seg_list[[i]] <- tibble(x = i, ymin = min(ypos), ymax = max(ypos))
}
seg_df <- bind_rows(seg_list)

stripe_fills <- adjustcolor(unname(SET_DISPLAY_COLORS_crvh[set_y_levels]), alpha.f = 0.20)

bar_bg_list <- lapply(seq_len(n_int), function(i) {
  bits <- as.logical(as.integer(strsplit(comb_names_ord[i], "")[[1]]))
  if (sum(bits) != 1) return(NULL)
  tibble(xmin = i - 0.5, xmax = i + 0.5,
         fill = unname(SET_DISPLAY_COLORS_crvh[set_order_ch[bits]]))
})
bar_bg <- bind_rows(compact(bar_bg_list))

lbl_sz <- scale_text(BASE_COUNT, PG_W)

pC_crvh_bars <- ggplot(bar_long, aes(x, count, fill = direction)) +
  geom_rect(data = bar_bg,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = bar_bg$fill, alpha = 0.20,
            color = "grey70", linewidth = 0.2, inherit.aes = FALSE) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.3) +
  geom_text(data = \(d) d |> filter(count > 0, is_single),
            aes(label = count, y = count / 2),
            position = position_dodge(width = 0.7), vjust = 0.5,
            size = lbl_sz - 0.3, color = "white", fontface = "bold") +
  geom_text(data = \(d) d |> filter(count > 0, !is_single),
            aes(label = count, y = count + 2),
            position = position_dodge(width = 0.7), vjust = 0,
            size = lbl_sz - 0.3, color = "black", fontface = "bold") +
  scale_fill_manual(values = c(Up = unname(DIR_COLORS["Up"]),
                                Down = unname(DIR_COLORS["Down"]))) +
  scale_x_continuous(expand = expansion(add = 0.3)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(y = "Intersection\nsize",
       title = "Contrast Overlap (CRvH)",
       subtitle = sprintf(
         "\u03A0 < 0.05 DEPs; %d \u03A0 / %d FDR total (%d mixed excl.)",
         pi_total_crvh, fdr_total_crvh, n_mixed
       ),
       tag = "C1") +
  FIG_THEME +
  theme(axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = margin(2, 5.5, 0, 5.5))

pC_crvh_dots <- ggplot() +
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
  coord_cartesian(ylim = c(0.55, 2.45)) +
  labs(x = NULL, y = NULL) +
  FIG_THEME +
  theme(axis.text.x  = element_blank(), axis.ticks = element_blank(),
        panel.grid   = element_blank(),
        panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.3),
        axis.text.y  = element_text(size = FIG_AXIS_TEXT - 1, face = "bold",
                                     margin = margin(r = 1)),
        plot.margin  = margin(0, 5.5, 0, 5.5))

# --- Fisher enrichment tests for all pairwise overlaps ---
overlap_tests <- list()
contrast_pairs <- combn(CONTRASTS, 2, simplify = FALSE)
n_bg <- length(all_genes_crvh)
for (pair in contrast_pairs) {
  a <- sig_sets_crvh[[pair[1]]]; b <- sig_sets_crvh[[pair[2]]]
  n_both <- length(intersect(a, b))
  n_a    <- length(a); n_b <- length(b)
  mat <- matrix(c(n_both, n_a - n_both, n_b - n_both,
                  n_bg - n_a - n_b + n_both), nrow = 2)
  ft <- fisher.test(mat, alternative = "greater")
  overlap_tests[[paste(pair, collapse = " & ")]] <- data.frame(
    set_A = pair[1], set_B = pair[2],
    n_A = n_a, n_B = n_b, overlap = n_both,
    expected = round(n_a * n_b / n_bg, 1),
    odds_ratio = round(ft$estimate, 2), p_value = ft$p.value
  )
}
overlap_df <- bind_rows(overlap_tests)
overlap_df$p_bh <- p.adjust(overlap_df$p_value, method = "BH")
write.csv(overlap_df, file.path(DAT, "panel_C_CRvH_overlap_enrichment.csv"), row.names = FALSE)

pC_crvh_combined <- (pC_crvh_bars / pC_crvh_dots) + plot_layout(heights = c(0.80, 0.20))

ggsave(file.path(RPT, "panel_C_CRvH_upset.pdf"), pC_crvh_combined,
       width = PG_W, height = 90, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_C_CRvH_upset.png"), pC_crvh_combined,
       width = PG_W, height = 90, units = "mm", dpi = 300)

cat("Panel C (CRvH) done.\n")
