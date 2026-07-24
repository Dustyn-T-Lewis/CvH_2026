# Global-proteome overview panel builders (Mito F01 layout), shared by the main
# (pooled CR vs Healthy) and supplement (CRE/PLA/Healthy) variants. Each builder
# takes a contrast set so main and supp render from one definition.

pacman::p_load(dplyr, tidyr, ggplot2, eulerr, ggplotify, patchwork, vegan)

DEP_CSV <- "03_DEP/a_non_imputed/c_data/combined_results_pi.csv"
FGSEA_CSV <- "04_Figures/shared/fgsea_CRvH.csv"
IMP_RDS <- "02_Normalization/imputation/c_data/DAList_imputed_missforest.rds"

MAIN_CTR <- c(CRvH_Baseline = "CR vs H", CR_Training = "Training", Resid = "Residual")
SUPP_CTR <- c(
  Baseline_Supplement = "CRE vs PLA", Training_CRE = "Tr. CRE",
  Training_PLA = "Tr. PLA", Supplement_Interaction = "Interaction"
)
# keyed by display label so the fill scales match the `ctr` factor directly
MAIN_PAL <- c("CR vs H" = "#D6604D", "Training" = "#9C27B0", "Residual" = "#00897B")
SUPP_PAL <- c(
  "CRE vs PLA" = "#00897B", "Tr. CRE" = "#2166AC",
  "Tr. PLA" = "#D6604D", "Interaction" = "#FF8F00"
)

load_dep <- function(ctr_map) {
  readr::read_csv(DEP_CSV, show_col_types = FALSE) |>
    filter(contrast %in% names(ctr_map)) |>
    mutate(ctr = factor(unname(ctr_map[contrast]), levels = unname(ctr_map)))
}

# Panel B — nested-threshold DEP counts (p / FDR / Pi telescoped per contrast)
panel_dep_counts <- function(dep, pal) {
  n_total <- dplyr::n_distinct(dep$uniprot_id)
  tiers <- dep |>
    group_by(ctr) |>
    summarise(
      `p < 0.05` = sum(P.Value < 0.05, na.rm = TRUE),
      `FDR < 0.10` = sum(adj.P.Val < 0.10, na.rm = TRUE),
      `Pi < 0.05` = sum(sig_pi != 0, na.rm = TRUE),
      .groups = "drop"
    ) |>
    pivot_longer(-ctr, names_to = "tier", values_to = "n") |>
    mutate(
      tier = factor(tier, levels = c("p < 0.05", "FDR < 0.10", "Pi < 0.05")),
      pct = 100 * n / n_total,
      alpha = c("p < 0.05" = 0.20, "FDR < 0.10" = 0.5, "Pi < 0.05" = 1)[as.character(tier)]
    )
  lab <- filter(tiers, tier == "Pi < 0.05")

  ggplot(tiers, aes(ctr, pct)) +
    geom_col(aes(fill = ctr, alpha = alpha),
      position = "identity",
      width = 0.9, colour = "black", linewidth = 0.3
    ) +
    geom_text(data = lab, aes(label = n), hjust = -0.15, size = 2.4, fontface = "bold") +
    scale_alpha_identity() +
    scale_fill_manual(values = pal, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    coord_flip() +
    labs(
      title = "DEP counts", subtitle = "p / FDR / Π nested per contrast",
      x = NULL, y = "% of proteome"
    ) +
    FIG_THEME +
    theme(axis.text.y = element_text(face = "bold"), panel.grid.major.y = element_blank())
}

# Panel C — log2FC distribution per contrast, median |log2FC| annotated
panel_effect <- function(dep, pal) {
  stat <- dep |>
    group_by(ctr) |>
    summarise(lab = sprintf("med|FC| %.2f", median(abs(logFC), na.rm = TRUE)), .groups = "drop")
  ggplot(dep, aes(logFC)) +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey55") +
    geom_histogram(aes(fill = ctr),
      bins = 40, colour = "white",
      linewidth = 0.1, alpha = 0.85
    ) +
    geom_label(
      data = stat, aes(x = -Inf, y = Inf, label = lab), inherit.aes = FALSE,
      hjust = -0.05, vjust = 1.1, size = 2, fontface = "bold",
      fill = alpha("white", 0.8), label.size = 0
    ) +
    facet_wrap(~ctr, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = pal, guide = "none") +
    scale_x_continuous(limits = c(-1, 1), oob = scales::squish) +
    labs(
      title = "Effect size", subtitle = "log2FC per contrast",
      x = expression(bold(log[2] ~ FC)), y = NULL
    ) +
    FIG_THEME +
    theme(
      strip.text = element_text(size = FIG_STRIP_SIZE),
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      panel.grid = element_blank()
    )
}

# Panel D — area-proportional overlap of the Pi-significant protein sets
panel_overlap <- function(dep, pal) {
  sets <- split(dep$uniprot_id[dep$sig_pi != 0], dep$ctr[dep$sig_pi != 0])
  sets <- sets[lengths(sets) > 0]
  fit <- eulerr::euler(sets, shape = "ellipse")
  grob <- plot(fit,
    fills = list(fill = unname(pal[names(sets)]), alpha = 0.5),
    edges = list(col = unname(pal[names(sets)]), lwd = 1.2),
    labels = list(fontsize = 6, fontfamily = "Helvetica", font = 2),
    quantities = list(fontsize = 6, fontfamily = "Helvetica"),
    legend = FALSE
  )
  ggplotify::as.ggplot(grob) +
    labs(title = "DEP overlap", subtitle = "Π < 0.05 sets") +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
      plot.subtitle = element_text(face = "italic", size = FIG_SUBTITLE_SIZE, colour = "grey30")
    )
}

# Panel E — up / down split within each Pi-significant set
panel_direction <- function(dep, pal) {
  dir_df <- dep |>
    filter(sig_pi != 0) |>
    mutate(direction = ifelse(sig_pi > 0, "Up", "Down")) |>
    count(ctr, direction)
  ggplot(dir_df, aes(ctr, n, fill = direction)) +
    geom_col(
      position = position_dodge(preserve = "single", width = 0.8),
      width = 0.74, colour = "black", linewidth = 0.2
    ) +
    geom_text(aes(y = n / 2, label = n),
      position = position_dodge(preserve = "single", width = 0.8),
      size = 2, colour = "white", fontface = "bold"
    ) +
    scale_fill_manual(values = DIR_COLORS[c("Up", "Down")], name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = "Direction", subtitle = "Up / Down per contrast",
      x = NULL, y = "Proteins (Π < 0.05)"
    ) +
    FIG_THEME +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      axis.title.y = element_text(margin = margin(r = 1)),
      legend.position = c(0.82, 0.92), legend.key.size = unit(2.4, "mm"),
      legend.background = element_rect(fill = alpha("white", 0.85), colour = NA)
    )
}

# Panel A — baseline PCA (T1 only); mode "main" pools CR, "supp" splits CRE/PLA
panel_pca <- function(mode = c("main", "supp")) {
  mode <- match.arg(mode)
  dal <- readRDS(IMP_RDS)
  meta <- as.data.frame(dal$metadata)
  base <- meta$Timepoint == "T1"
  meta <- meta[base, ]
  mat <- as.matrix(dal$data)[, meta$Col_ID]

  if (mode == "main") {
    grp <- factor(ifelse(meta$Group_Time == "H_T1", "Healthy", "CR"),
      levels = c("CR", "Healthy")
    )
    pal <- c(CR = "#D6604D", Healthy = "#4DAF4A")
    shp <- c(CR = 16, Healthy = 15)
    title <- "Sample PCA"
  } else {
    grp <- factor(
      dplyr::recode(meta$Group_Time,
        CRE_T1 = "CRE", PLA_T1 = "PLA", H_T1 = "Healthy"
      ),
      levels = c("CRE", "PLA", "Healthy")
    )
    pal <- c(CRE = "#2166AC", PLA = "#D6604D", Healthy = "#4DAF4A")
    shp <- c(CRE = 16, PLA = 17, Healthy = 15)
    title <- "Sample PCA (arms)"
  }

  pca <- prcomp(t(mat), center = TRUE, scale. = TRUE)
  var_pct <- round(100 * summary(pca)$importance[2, 1:2], 1)
  df <- data.frame(pca$x[, 1:2], group = grp)
  d <- dist(scale(t(mat)))
  set.seed(42)
  perm <- adonis2(d ~ grp, permutations = 999)

  ggplot(df, aes(PC1, PC2, colour = group, shape = group)) +
    stat_ellipse(aes(fill = group),
      geom = "polygon", alpha = 0.10,
      level = 0.80, show.legend = FALSE
    ) +
    stat_ellipse(level = 0.80, linewidth = 0.35, linetype = "dashed", show.legend = FALSE) +
    geom_point(size = 2.2, alpha = 0.9) +
    scale_colour_manual(values = pal, name = NULL) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_shape_manual(values = shp, name = NULL) +
    labs(
      title = title,
      subtitle = sprintf("PERMANOVA R² = %.2f, %s", perm$R2[1], fmt_p(perm$`Pr(>F)`[1])),
      x = sprintf("PC1 (%.1f%%)", var_pct[1]), y = sprintf("PC2 (%.1f%%)", var_pct[2])
    ) +
    FIG_THEME +
    theme(
      legend.position = c(0.72, 0.14),
      legend.background = element_rect(fill = alpha("white", 0.7), colour = NA)
    )
}

# Panel F — significant pathways per contrast, diverging Up/Down, stacked by database
panel_pathway <- function(ctr_map) {
  fg <- readr::read_csv(FGSEA_CSV, show_col_types = FALSE) |>
    filter(
      contrast %in% names(ctr_map), !is.na(padj), padj < 0.05,
      database %in% DB_ORDER
    ) |>
    mutate(
      ctr = factor(unname(ctr_map[contrast]), levels = unname(ctr_map)),
      database = factor(database, levels = DB_ORDER),
      dir = ifelse(NES > 0, "Up", "Down")
    ) |>
    count(ctr, database, dir) |>
    mutate(y = ifelse(dir == "Up", n, -n))

  ggplot(fg, aes(ctr, y, fill = database)) +
    geom_col(
      width = 0.74, colour = "black", linewidth = 0.2,
      position = position_stack(reverse = TRUE)
    ) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey35") +
    scale_fill_manual(
      values = DB_COLORS, breaks = DB_ORDER, name = NULL,
      guide = guide_legend(ncol = 1)
    ) +
    scale_y_continuous(labels = abs) +
    labs(
      title = "Pathways", subtitle = "fGSEA FDR < 0.05; up / down",
      x = NULL, y = "Sig. pathways"
    ) +
    FIG_THEME +
    guides(fill = guide_legend(ncol = 1)) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
      legend.background = element_rect(
        fill = alpha("white", 0.75), colour = "grey70", linewidth = 0.3
      )
    )
}
