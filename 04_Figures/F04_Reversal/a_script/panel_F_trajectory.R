# Reversal Panel F: Trajectory clustering of the disease signature
# Soft (fuzzy c-means) clustering of each disease-signature protein's standardized
# abundance profile across the ordered axis H_pre -> CR_pre -> CR_post, then a
# reversal-class label per cluster + per-cluster ORA. Pattern-first complement to
# the per-contrast scatter (A) and NES scatter (B).
#   Engine = e1071::cmeans (the fuzzy c-means Mfuzz wraps; Futschik & Carlisle 2005).
#   Descriptive only -- inference stays with fry / RRHO2 / permutation null.
setwd(here::here())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")
pacman::p_load(tidyverse, e1071, ggrepel, patchwork)

RPT_PNG <- "04_Figures/F04_Reversal/b_reports/main/png/panels"
RPT_PDF <- "04_Figures/F04_Reversal/b_reports/main/pdf/panels"
DAT <- "04_Figures/F04_Reversal/c_data/panel_F"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()
set.seed(42)

ORDER <- c("H_pre", "CR_pre", "CR_post") # healthy -> diseased -> trained
AXISLAB <- c(H_pre = "Healthy", CR_pre = "CR pre", CR_post = "CR post")
K <- 6 # soft clusters
M_FUZZ <- 1.25 # fuzzifier (Mfuzz-style)
PHI_BAND <- 0.25 # engine band (REVERSAL_PHI_BAND)
CLASS_COLORS <- c(
  Normalized = "#1B7837", Persistent = "#878787",
  Exacerbated = "#D6604D"
)

# ── Data: group-mean trajectories of the disease signature ───────────────────
source("04_Figures/F04_Reversal/a_script/f04_data.R") # dep_df, dal
stopifnot(all(ORDER %in% dal$metadata$group_time))

gmean <- sapply(ORDER, function(g) {
  rowMeans(dal$data[, dal$metadata$Col_ID[dal$metadata$group_time == g], drop = FALSE])
})
rownames(gmean) <- rownames(dal$data) # uniprot_id

# disease signature = Pi < 0.05 on CRvH_Baseline. Per-protein reversal class
# from the engine band on the actual D/T logFC (phi = -T/D), so the colour is the
# rigorous reversal call (consistent with fry / permutation null), not a cluster guess.
sig <- dep_df |>
  filter(pi_score_CRvH_Baseline < 0.05) |>
  mutate(
    phi = -logFC_CR_Training / logFC_CRvH_Baseline,
    reversal_class = case_when(
      phi >= PHI_BAND ~ "Normalized",
      phi <= -PHI_BAND ~ "Exacerbated",
      TRUE ~ "Persistent"
    )
  ) |>
  distinct(uniprot_id, gene, phi, reversal_class)
gmean <- gmean[rownames(gmean) %in% sig$uniprot_id, , drop = FALSE]

# z-score each protein's 3-point profile (drop constant rows)
z <- t(apply(gmean, 1, function(r) (r - mean(r)) / sd(r)))
z <- z[is.finite(rowSums(z)), , drop = FALSE]
message(sprintf("  Trajectory clustering: %d disease-signature proteins", nrow(z)))

# ── Fuzzy c-means soft clustering ────────────────────────────────────────────
fcm <- cmeans(z, centers = K, m = M_FUZZ, iter.max = 200)
memb_max <- apply(fcm$membership, 1, max)
assign <- tibble(
  uniprot_id = rownames(z),
  cluster = fcm$cluster, membership = memb_max
) |>
  left_join(sig, by = "uniprot_id") # per-protein phi + reversal_class
write_csv(assign, file.path(DAT, "trajectory_clusters.csv"))

# centroid shape per cluster -> dominant-pattern label for the facet strip only
cent <- as_tibble(fcm$centers, .name_repair = "minimal")
names(cent) <- ORDER
cent <- cent |>
  mutate(
    cluster = row_number(),
    d_pre = CR_pre - H_pre,
    phi_c = -(CR_post - CR_pre) / d_pre,
    pattern = case_when(
      abs(d_pre) < 0.4 ~ "training-emergent",
      phi_c >= PHI_BAND ~ "normalizing",
      phi_c <= -PHI_BAND ~ "exacerbating",
      TRUE ~ "persistent"
    )
  )

clab <- cent |>
  mutate(
    n = as.integer(table(factor(assign$cluster, levels = seq_len(K)))[cluster]),
    lab = sprintf("Cluster %d · %s (n=%d)", cluster, pattern, n)
  )
class_n <- assign |> count(reversal_class)
rev_dir <- mean(assign$phi > 0, na.rm = TRUE) # directional reverse fraction (phi > 0)

# ── Per-cluster ORA (top 5 pathways per cluster) ─────────────────────────────
universe <- unique(dep_df$gene)
pw <- build_pathway_collection(
  min_size = 15, max_size = 500,
  include_goslim = FALSE, exclude_variants = TRUE
)
cluster_ora <- map(seq_len(K), function(k) {
  g <- assign$gene[assign$cluster == k]
  if (length(g) < 5) {
    return(tibble())
  }
  r <- tryCatch(
    run_ora_deduplicated(
      genes = g, universe = universe, pathways = pw,
      jaccard_cutoff = 0.5, min_size = 15, max_size = 500, padj_cutoff = 1
    ),
    error = function(e) tibble()
  )
  if (nrow(r) == 0) {
    return(tibble())
  }
  r |>
    arrange(padj) |>
    slice_head(n = 5) |>
    mutate(pathway_label = clean_pathway_name(pathway))
})
top_path <- map_dfr(seq_len(K), function(k) {
  d <- cluster_ora[[k]]
  if (nrow(d) == 0) {
    return(tibble(cluster = k, pathway_label = NA_character_, padj = NA_real_))
  }
  d |>
    slice_head(n = 1) |>
    transmute(cluster = k, pathway_label, padj)
})
if (nrow(top_path)) write_csv(top_path, file.path(DAT, "cluster_top_pathway.csv"))

class_n <- assign |> count(reversal_class)
rev_dir <- mean(assign$phi > 0, na.rm = TRUE)

# ── Long frame for trajectory plotting ───────────────────────────────────────
long <- as_tibble(z, rownames = "uniprot_id") |>
  pivot_longer(all_of(ORDER), names_to = "cond", values_to = "zval") |>
  left_join(select(assign, uniprot_id, cluster, membership, reversal_class), by = "uniprot_id") |>
  mutate(cond = factor(cond, levels = ORDER))
cent_long <- cent |>
  select(cluster, all_of(ORDER)) |>
  pivot_longer(all_of(ORDER), names_to = "cond", values_to = "zval") |>
  mutate(cond = factor(cond, levels = ORDER))

# ── Per-cluster unit: trajectory line (left) + top-5 ORA bars (right) ─────────
CLUSTER_COLORS <- c(
  "#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377"
)
cluster_color <- function(k) CLUSTER_COLORS[((k - 1) %% length(CLUSTER_COLORS)) + 1]

traj_one <- function(k) {
  d <- filter(long, cluster == k)
  cl <- filter(cent_long, cluster == k)
  ggplot(d, aes(cond, zval, group = uniprot_id)) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
    geom_line(aes(alpha = membership),
      color = cluster_color(k), linewidth = 0.25
    ) +
    geom_line(
      data = cl, aes(cond, zval, group = cluster), inherit.aes = FALSE,
      color = "grey10", linewidth = 1.0
    ) +
    scale_alpha_continuous(range = c(0.12, 0.6), guide = "none") +
    scale_x_discrete(labels = AXISLAB) +
    labs(title = clab$lab[clab$cluster == k], x = NULL, y = "z") +
    FIG_THEME +
    theme(
      plot.title = element_text(size = 8, face = "bold"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(size = 6),
      plot.margin = margin(2, 1, 0, 0, "mm")
    )
}

ora_bars_one <- function(k) {
  df <- cluster_ora[[k]]
  if (nrow(df) == 0) {
    return(ggplot() +
      theme_void() +
      annotate("text",
        x = 0.5, y = 0.5, label = "No sig. pathways",
        size = 2.6, color = "grey60"
      ))
  }
  col <- cluster_color(k)
  bars <- df |>
    mutate(
      neg_log_padj = -log10(pmax(padj, 1e-20)),
      significant = padj < 0.05,
      bar_fill = ifelse(significant, scales::alpha(col, 0.85),
        scales::alpha(col, 0.30)
      ),
      name = pathway_label,
      p_lab = ifelse(padj < 0.001, "p<0.001", sprintf("p=%.3f", padj)),
      p_face = ifelse(significant, "bold", "plain"),
      y = rev(row_number()), bar_h = 0.85
    )
  x_max <- max(bars$neg_log_padj)
  bars <- bars |>
    mutate(
      fits = neg_log_padj >= x_max * 0.5,
      p_x = neg_log_padj + x_max * 0.03,
      name_x = neg_log_padj + x_max * 0.28
    )
  ggplot(bars, aes(y = y)) +
    geom_rect(aes(xmin = 0, xmax = neg_log_padj, ymin = y - bar_h / 2, ymax = y + bar_h / 2),
      fill = bars$bar_fill, color = "black", linewidth = 0.3
    ) +
    ggfittext::geom_fit_text(
      data = ~ filter(.x, fits),
      aes(
        xmin = 0, xmax = neg_log_padj, ymin = y - bar_h / 2, ymax = y + bar_h / 2,
        label = name
      ),
      reflow = TRUE, grow = FALSE, contrast = TRUE, fontface = "bold", min.size = 3
    ) +
    geom_text(
      data = ~ filter(.x, !fits), aes(x = name_x, y = y, label = name),
      hjust = 0, size = 2.0, fontface = "bold", color = "grey15", lineheight = 0.85
    ) +
    geom_text(aes(x = p_x, y = y, label = p_lab, fontface = p_face),
      hjust = 0, size = 1.9, color = "grey25"
    ) +
    labs(x = expression(-log[10](p[adj])), y = NULL) +
    scale_x_continuous(
      limits = c(0, x_max * 1.5), breaks = scales::pretty_breaks(3),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(limits = c(0.3, nrow(bars) + 0.7), expand = c(0, 0)) +
    theme_minimal(base_size = 8) +
    theme(
      panel.grid = element_blank(), axis.text.y = element_blank(),
      axis.ticks.y = element_blank(), axis.title.y = element_blank(),
      axis.text.x = element_text(size = 6), axis.title.x = element_text(size = 6),
      axis.line.x = element_line(color = "grey40", linewidth = 0.3),
      plot.margin = margin(2, 2, 0, 0, "mm")
    )
}

units <- map(seq_len(K), function(k) {
  (traj_one(k) | ora_bars_one(k)) + plot_layout(widths = c(1, 1.35))
})

traj_fig <- wrap_plots(units, ncol = 2) +
  plot_annotation(
    title = "Cancer Recovery Reversal: trajectory clusters + pathway ORA",
    subtitle = sprintf(
      "Fuzzy c-means k=%d | %d disease DEPs | %.0f%% reverse the disease signature | one colour per cluster, dominant pattern in each title",
      K, nrow(z), 100 * rev_dir
    ),
    theme = theme(
      plot.title = element_text(size = FIG_TITLE_SIZE, face = "bold"),
      plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, color = "grey30")
    )
  )

COMP_W <- 340
COMP_H <- 230
ggsave(file.path(RPT_PNG, "MAIN_F04_trajectory.png"), traj_fig,
  width = COMP_W, height = COMP_H, units = "mm", dpi = 300
)
ggsave(file.path(RPT_PDF, "MAIN_F04_trajectory.pdf"), traj_fig,
  width = COMP_W, height = COMP_H, units = "mm", device = pdf_device
)
message("Reversal trajectory figure (clusters + ORA) done")
