#!/usr/bin/env Rscript
# volcano_ring.R — Circular volcano-in-ring composite plot utility
# Standard Cartesian ggplot with ggforce::geom_arc_bar(); NO coord_polar().
# Sources style.R for palettes and sizing constants.

pacman::p_load(tidyverse, ggforce, patchwork, scales)

if (!exists("FIG_THEME")) source("04_Figures/shared/style.R")

clean_ring_label <- function(name) {
  name %>%
    str_remove("^HALLMARK_") %>%
    str_remove("^GOBP_") %>%
    str_remove("^GOCC_") %>%
    str_remove("^GOMF_") %>%
    str_remove("^GOSLIM_") %>%
    str_replace_all("_", " ") %>%
    str_to_title() %>%
    # Capitalisation fixes
    str_replace("Mtorc1", "mTORC1") %>%
    str_replace("Myc ", "MYC ") %>%
    str_replace("E2f ", "E2F ") %>%
    str_replace("Dna ", "DNA ") %>%
    str_replace("Rna ", "RNA ") %>%
    str_replace("Trna ", "tRNA ") %>%
    str_replace("Tnfa ", "TNFa ") %>%
    str_replace("Uv ", "UV ") %>%
    str_replace("G2m ", "G2M ") %>%
    str_replace("Il6 ", "IL6 ") %>%
    str_replace("Il2 ", "IL2 ") %>%
    str_replace("Kras ", "KRAS ") %>%
    str_replace("P53 ", "p53 ") %>%
    str_replace("Tgf ", "TGF ") %>%
    str_replace("Nf Kb", "NF-kB") %>%
    str_replace("Atp ", "ATP ") %>%
    str_replace("Nadh ", "NADH ") %>%
    str_replace("Ifn", "IFN") %>%
    str_replace("Pi3k", "PI3K") %>%
    str_replace("Akt", "AKT") %>%
    str_replace("Mtor", "mTOR") %>%
    str_replace("Oxidative Phosphorylation", "OXPHOS") %>%
    # Smart abbreviations for ring labels
    str_replace("Mitochondrial", "Mito.") %>%
    str_replace("Ubiquinone", "UQ") %>%
    str_replace("Organization", "Org.") %>%
    str_replace("Cytoskeleton", "Cytoskel.") %>%
    str_replace("Microtubule", "MT") %>%
    str_replace("Respiratory", "Resp.") %>%
    str_replace("Electron Transport", "ETC") %>%
    str_replace("Synthesis Coupled", "Synth.-Coupled") %>%
    str_replace("Ubiquitin Dependent", "Ub-Dep.") %>%
    str_replace("Proteasome Mediated", "Proteasome-Med.") %>%
    str_replace("Proteasomal", "Proteas.") %>%
    str_replace("Phosphorylation", "Phosph.") %>%
    str_replace("Modification", "Mod.") %>%
    str_replace("Intracellular", "Intracell.") %>%
    str_replace("Regulation Of", "Reg.") %>%
    str_replace("Signaling Pathway", "Signaling") %>%
    str_replace("Biosynthetic Process", "Biosynthesis") %>%
    str_replace("Catabolic Process", "Catabolism") %>%
    str_replace("Metabolic Process", "Metabolism") %>%
    str_replace("Based Process", "Process") %>%
    str_replace("Response To", "Resp. to") %>%
    str_replace("Extracellular Matrix", "ECM") %>%
    str_replace("Epithelial Mesenchymal Transition", "EMT") %>%
    str_wrap(width = 22)
}

prepare_ring_data <- function(go_df,
                              contrast,
                              n_terms      = 12,
                              gap_degrees  = 3,
                              start_offset = 0,
                              databases    = c("Hallmark", "GO Slim")) {

  ring <- go_df %>%
    filter(contrast == !!contrast,
           database %in% databases,
           padj < 0.05) %>%
    arrange(padj) %>%
    slice_head(n = n_terms)

  n <- nrow(ring)
  if (n == 0) {
    warning("prepare_ring_data: no significant terms for '", contrast, "'")
    return(tibble())
  }

  arc_width_deg <- (360 - n * gap_degrees) / n

  ring %>%
    mutate(
      term_idx    = row_number(),
      start_deg   = start_offset + (term_idx - 1) * (arc_width_deg + gap_degrees),
      end_deg     = start_deg + arc_width_deg,
      mid_deg     = (start_deg + end_deg) / 2,
      start_rad   = start_deg * pi / 180,
      end_rad     = end_deg   * pi / 180,
      mid_rad     = mid_deg   * pi / 180,
      arc_r1_var  = 5.6,
      clean_label = clean_ring_label(pathway),
      gene_list   = str_split(leadingEdge, ";")
    )
}

build_tick_data <- function(ring_data,
                            de_df,
                            contrast,
                            tick_r0 = 4.4,
                            tick_r1 = 4.8) {

  if (nrow(ring_data) == 0) return(tibble())

  logfc_col <- paste0("logFC_", contrast)

  gene_lfc <- de_df %>%
    dplyr::select(gene, lfc = all_of(logfc_col)) %>%
    filter(!is.na(lfc)) %>%
    distinct(gene, .keep_all = TRUE)

  pad_rad <- 0.5 * pi / 180

  map_dfr(seq_len(nrow(ring_data)), function(i) {
    row <- ring_data[i, ]
    genes_in_arc <- intersect(row$gene_list[[1]], gene_lfc$gene)
    n_genes <- length(genes_in_arc)
    if (n_genes == 0) return(tibble())

    arc_start <- row$start_rad + pad_rad
    arc_end   <- row$end_rad   - pad_rad
    if (arc_end <= arc_start) arc_end <- arc_start + pad_rad

    tick_angles <- seq(arc_start, arc_end, length.out = n_genes)

    matched <- gene_lfc %>%
      filter(gene %in% genes_in_arc) %>%
      `[`(match(genes_in_arc, .$gene), ) %>%
      filter(!is.na(gene))

    n_final <- nrow(matched)
    if (n_final == 0) return(tibble())
    tick_angles <- tick_angles[seq_len(n_final)]

    tibble(
      gene      = matched$gene,
      logFC     = matched$lfc,
      direction = ifelse(matched$lfc > 0, "Up", "Down"),
      angle_rad = tick_angles,
      x0        = tick_r0 * sin(tick_angles),
      y0        = tick_r0 * cos(tick_angles),
      x1        = tick_r1 * sin(tick_angles),
      y1        = tick_r1 * cos(tick_angles),
      term_idx  = row$term_idx,
      pathway   = row$pathway
    )
  })
}

build_volcano_layers <- function(de_df,
                                 contrast,
                                 volcano_radius = 3.5,
                                 fc_thresh      = log2(1.5),
                                 p_thresh       = 0.05,
                                 up_color       = DIR_COLORS["Up"],
                                 down_color     = DIR_COLORS["Down"],
                                 ns_color       = DIR_COLORS["NS"],
                                 point_size     = 0.6,
                                 point_alpha    = 0.5,
                                 count_label_size = 2.8) {

  logfc_col <- paste0("logFC_", contrast)
  pval_col  <- paste0("P.Value_", contrast)
  pi_col    <- paste0("pi_score_", contrast)

  vdf <- de_df %>%
    transmute(
      gene       = gene,
      logFC      = .data[[logfc_col]],
      pvalue     = .data[[pval_col]],
      pi_score   = .data[[pi_col]],
      neg_log10p = -log10(pvalue)
    ) %>%
    filter(!is.na(logFC), !is.na(pvalue), is.finite(neg_log10p)) %>%
    mutate(
      direction = case_when(
        pi_score < 0.05 & logFC > 0 ~ "Up",
        pi_score < 0.05 & logFC < 0 ~ "Down",
        TRUE                         ~ "NS"
      )
    )

  n_up   <- sum(vdf$direction == "Up")
  n_down <- sum(vdf$direction == "Down")

  x_data_max <- max(abs(vdf$logFC), na.rm = TRUE)
  y_data_max <- max(vdf$neg_log10p, na.rm = TRUE)

  margin <- 0.92
  vr <- volcano_radius * margin

  scale_x_fn <- function(val) val / x_data_max * vr
  scale_y_fn <- function(val) val / y_data_max * vr

  vdf <- vdf %>%
    mutate(
      x_plot = logFC / x_data_max * vr,
      y_plot = (neg_log10p / y_data_max) * 2 * vr - vr
    )

  vdf_ns  <- vdf %>% filter(direction == "NS")
  vdf_sig <- vdf %>% filter(direction != "NS")

  layers <- list(
    ns_points = geom_point(
      data = vdf_ns, aes(x = x_plot, y = y_plot),
      color = ns_color, size = point_size * 0.9, alpha = point_alpha * 0.6,
      inherit.aes = FALSE
    ),
    sig_points = geom_point(
      data = vdf_sig, aes(x = x_plot, y = y_plot, color = direction),
      size = point_size * 1.3, alpha = point_alpha * 1.3, stroke = 0.3,
      inherit.aes = FALSE
    ),
    color_scale = scale_color_manual(
      values = c(Up = unname(up_color), Down = unname(down_color)),
      guide  = "none"
    ),
    x_axis_line = annotate(
      "segment", x = -vr * 0.42, xend = vr * 0.42, y = -vr, yend = -vr,
      linewidth = 0.3, linetype = "dashed", color = "grey50",
      arrow = arrow(ends = "both", length = unit(1.2, "mm"), type = "closed")
    ),
    x_axis_up = annotate(
      "text", x = vr * 0.45, y = -vr,
      label = "up", size = count_label_size * 0.6,
      color = unname(up_color), fontface = "italic", hjust = 0
    ),
    x_axis_down = annotate(
      "text", x = -vr * 0.45, y = -vr,
      label = "down", size = count_label_size * 0.6,
      color = unname(down_color), fontface = "italic", hjust = 1
    ),
    x_axis_label = annotate(
      "text", x = 0, y = -vr - 0.35,
      label = "log2 FC", size = count_label_size * 0.6,
      color = "grey40", fontface = "italic"
    ),
    y_axis_line = annotate(
      "segment", x = 0, xend = 0, y = -vr, yend = vr * 0.96,
      linewidth = 0.3, linetype = "dashed", color = "grey50",
      arrow = arrow(ends = "last", length = unit(1.2, "mm"), type = "closed")
    ),
    y_axis_label = annotate(
      "text", x = 0, y = vr * 1.04,
      label = expression(-log[10]~italic(p)), size = count_label_size * 0.6,
      color = "grey40"
    ),
    n_up_box = annotate(
      "label", x = vr * 0.5, y = vr * 0.85,
      label = n_up, size = count_label_size,
      color = "black", fill = alpha(up_color, 0.9), fontface = "bold",
      label.padding = unit(2.5, "pt"), label.r = unit(2, "pt"),
      linewidth = 0.4
    ),
    n_up_text = annotate(
      "text", x = vr * 0.5, y = vr * 0.85,
      label = n_up, size = count_label_size,
      color = "white", fontface = "bold"
    ),
    n_down_box = annotate(
      "label", x = -vr * 0.5, y = vr * 0.85,
      label = n_down, size = count_label_size,
      color = "black", fill = alpha(down_color, 0.9), fontface = "bold",
      label.padding = unit(2.5, "pt"), label.r = unit(2, "pt"),
      linewidth = 0.4
    ),
    n_down_text = annotate(
      "text", x = -vr * 0.5, y = vr * 0.85,
      label = n_down, size = count_label_size,
      color = "white", fontface = "bold"
    )
  )

  attr(layers, "x_data_max") <- x_data_max
  attr(layers, "y_data_max") <- y_data_max
  attr(layers, "n_up")       <- n_up
  attr(layers, "n_down")     <- n_down
  attr(layers, "scale_x")    <- scale_x_fn
  attr(layers, "scale_y")    <- scale_y_fn

  layers
}

build_ring_layers <- function(ring_data,
                              tick_data,
                              tick_r0    = 4.4,
                              tick_r1    = 4.8,
                              arc_r0     = 4.8,
                              arc_r1     = 5.6,
                              up_color   = DIR_COLORS["Up"],
                              down_color = DIR_COLORS["Down"]) {

  if (nrow(ring_data) == 0) return(list())

  layers <- list()

  layers$tick_bg <- geom_arc_bar(
    data = ring_data,
    aes(x0 = 0, y0 = 0, r0 = tick_r0, r = tick_r1,
        start = start_rad, end = end_rad),
    fill = "grey93", color = "grey78", linewidth = 0.15,
    inherit.aes = FALSE
  )



  if (nrow(tick_data) > 0) {
    layers$ticks <- geom_segment(
      data = tick_data,
      aes(x = x0, y = y0, xend = x1, yend = y1, color = direction),
      linewidth = 0.15, alpha = 0.7,
      inherit.aes = FALSE
    )
  }

  layers$enrich_arcs <- geom_arc_bar(
    data = ring_data,
    aes(x0 = 0, y0 = 0, r0 = arc_r0, r = arc_r1_var,
        start = start_rad, end = end_rad, fill = NES),
    color = "grey40", linewidth = 0.2,
    inherit.aes = FALSE
  )

  layers$fill_scale <- scale_fill_gradientn(
    colours = c("#08306B", "#4393C3", "white", "#D6604D", "#67000D"),
    values  = scales::rescale(c(-3, -1.5, 0, 1.5, 3)),
    limits  = c(-3, 3),
    oob     = scales::squish,
    name    = "NES"
  )

  layers
}

build_label_layer <- function(ring_data,
                              label_r    = 6.8,
                              label_pad  = 0.5,
                              label_size = 2.8) {

  if (nrow(ring_data) == 0) return(list())

  n_up   <- sum(ring_data$NES > 0)
  n_down <- sum(ring_data$NES <= 0)

  up_df <- ring_data %>% filter(NES > 0) %>%
    mutate(idx = row_number(),
           legend_label = clean_label)
  down_df <- ring_data %>% filter(NES <= 0) %>%
    mutate(idx = row_number(),
           legend_label = clean_label)

  layers <- list()
  y_top <- label_r * 0.8
  fs <- label_size * 0.6

  if (n_up > 0) {
    up_df$y_pos <- seq(y_top, -y_top, length.out = n_up)
    layers$up_legend <- geom_text(
      data = up_df,
      aes(x = label_r + 0.3, y = y_pos, label = legend_label),
      hjust = 0, size = fs, color = "grey25", fontface = "bold",
      lineheight = 0.85, inherit.aes = FALSE
    )
    layers$up_header <- annotate(
      "text", x = label_r + 0.3, y = y_top + 0.6,
      label = "Upregulated", hjust = 0, size = fs * 1.1,
      fontface = "bold.italic", color = DIR_COLORS["Up"]
    )
  }

  if (n_down > 0) {
    down_df$y_pos <- seq(y_top, -y_top, length.out = n_down)
    layers$down_legend <- geom_text(
      data = down_df,
      aes(x = -(label_r + 0.3), y = y_pos, label = legend_label),
      hjust = 1, size = fs, color = "grey25", fontface = "bold",
      lineheight = 0.85, inherit.aes = FALSE
    )
    layers$down_header <- annotate(
      "text", x = -(label_r + 0.3), y = y_top + 0.6,
      label = "Downregulated", hjust = 1, size = fs * 1.1,
      fontface = "bold.italic", color = DIR_COLORS["Down"]
    )
  }

  layers
}

# min_size excludes small gene sets prone to tissue-irrelevant GO artifacts
# (Reimand et al. 2019, Nat Protocols S3.4).
# n_each = NULL shows all significant terms (GO Slim + Hallmark are small,
# non-redundant collections — typically 9-15 sig terms per contrast).
select_ring_terms <- function(go_df, contrast_name, n_each = NULL,
                              databases = c("Hallmark", "GO Slim"),
                              min_size = 15) {
  sig <- go_df %>%
    filter(contrast == contrast_name, database %in% databases,
           padj < 0.05, size >= min_size) %>%
    arrange(padj)

  up   <- sig %>% filter(NES > 0)
  down <- sig %>% filter(NES < 0)

  if (!is.null(n_each)) {
    up   <- up   %>% slice_head(n = n_each)
    down <- down %>% slice_head(n = n_each)
  }

  bind_rows(up, down)
}

center_ring_angles <- function(ring, n_up) {
  n <- nrow(ring)
  if (n < 2 || n_up < 1) return(ring)
  # Midpoint of the Up block = average of first arc start and last Up arc end
  up_mid <- (ring$start_deg[1] + ring$end_deg[min(n_up, n)]) / 2
  offset <- 90 - up_mid   # center Up block at 90° (right side)
  ring$start_deg <- ring$start_deg + offset
  ring$end_deg   <- ring$end_deg   + offset
  ring$mid_deg   <- ring$mid_deg   + offset
  ring$start_rad <- ring$start_deg * pi / 180
  ring$end_rad   <- ring$end_deg   * pi / 180
  ring$mid_rad   <- ring$mid_deg   * pi / 180
  ring
}

build_ring_with_gaps <- function(top_terms, contrast_name, go_df,
                                 n_each = NULL,
                                 databases = c("Hallmark", "GO Slim")) {
  real_rows <- go_df %>%
    filter(contrast == contrast_name, pathway %in% top_terms$pathway)
  # Save real padj before overwriting with fake ordering values
  padj_lookup <- real_rows %>%
    dplyr::select(pathway, padj) %>%
    distinct(pathway, .keep_all = TRUE)
  go_subset <- real_rows %>%
    mutate(padj = match(pathway, top_terms$pathway) * 1e-10)

  ring <- prepare_ring_data(
    go_df = go_subset, contrast = contrast_name,
    n_terms = nrow(top_terms), gap_degrees = 3, start_offset = 0,
    databases = databases
  )

  # Restore real padj for proportional arc widths
  ring$padj <- padj_lookup$padj[match(ring$pathway, padj_lookup$pathway)]

  n <- nrow(ring)
  n_up <- sum(ring$NES > 0)

  if (n >= 2) {
    gap_normal <- 3; gap_split <- 8
    gaps <- rep(gap_normal, n)
    if (n_up > 0 && n_up < n) gaps[n_up] <- gap_split
    gaps[n] <- gap_split
    arc_budget <- 360 - sum(gaps)
    arc_widths <- rep(arc_budget / n, n)
    # Variable radial height: more significant terms extend further outward
    min_height <- 0.05; max_height <- 1.6
    neg_lp <- -log10(pmax(ring$padj, .Machine$double.xmin))
    scaled <- (neg_lp - min(neg_lp)) / (max(neg_lp) - min(neg_lp))
    ring$arc_r1_var <- 4.8 + min_height +
      (max_height - min_height) * sqrt(scaled)
    cum_offset <- 0
    for (i in seq_len(n)) {
      if (i > 1) cum_offset <- cum_offset + arc_widths[i - 1] + gaps[i - 1]
      ring$start_deg[i] <- cum_offset
      ring$end_deg[i]   <- ring$start_deg[i] + arc_widths[i]
      ring$mid_deg[i]   <- (ring$start_deg[i] + ring$end_deg[i]) / 2
      ring$start_rad[i] <- ring$start_deg[i] * pi / 180
      ring$end_rad[i]   <- ring$end_deg[i]   * pi / 180
      ring$mid_rad[i]   <- ring$mid_deg[i]   * pi / 180
    }
    ring <- center_ring_angles(ring, n_up)
  }
  ring
}

make_volcano_ring <- function(de_df,
                              go_df,
                              contrast,
                              title              = NULL,
                              contrast_title     = NULL,
                              contrast_subtitle  = NULL,
                              title_size         = 12,
                              n_terms            = 12,
                              gap_degrees        = 3,
                              start_offset       = 0,
                              databases          = c("Hallmark", "GO Slim"),
                              volcano_radius     = 3.5,
                              tick_r0            = 4.4,
                              tick_r1            = 4.8,
                              arc_r0             = 4.8,
                              arc_r1             = 5.6,
                              label_r            = 6.8,
                              fc_thresh          = log2(1.5),
                              p_thresh           = 0.05,
                              up_color           = DIR_COLORS["Up"],
                              down_color         = DIR_COLORS["Down"],
                              ns_color           = DIR_COLORS["NS"],
                              point_size         = 0.6,
                              point_alpha        = 0.5,
                              label_size         = 2.8,
                              count_label_size   = 2.8,
                              ring_data_override = NULL) {

  if (!is.null(ring_data_override)) {
    ring_data <- ring_data_override
  } else {
    ring_data <- prepare_ring_data(
      go_df = go_df, contrast = contrast, n_terms = n_terms,
      gap_degrees = gap_degrees, start_offset = start_offset,
      databases = databases
    )
  }

  tick_data <- build_tick_data(
    ring_data = ring_data, de_df = de_df, contrast = contrast,
    tick_r0 = tick_r0, tick_r1 = tick_r1
  )

  volcano_layers <- build_volcano_layers(
    de_df = de_df, contrast = contrast, volcano_radius = volcano_radius,
    fc_thresh = fc_thresh, p_thresh = p_thresh,
    up_color = up_color, down_color = down_color, ns_color = ns_color,
    point_size = point_size, point_alpha = point_alpha,
    count_label_size = count_label_size
  )

  ring_layers <- build_ring_layers(
    ring_data = ring_data, tick_data = tick_data,
    tick_r0 = tick_r0, tick_r1 = tick_r1,
    arc_r0 = arc_r0, arc_r1 = arc_r1,
    up_color = up_color, down_color = down_color
  )

  label_layers <- build_label_layer(
    ring_data = ring_data, label_r = label_r, label_size = label_size
  )

  p <- ggplot() +
    ring_layers$tick_bg +
    volcano_layers +
    ring_layers$ticks +
    ring_layers$enrich_arcs +
    ring_layers$fill_scale +
    label_layers +
    coord_fixed(
      xlim = c(-(label_r + 5), label_r + 5),
      ylim = c(-(label_r + 1.5), label_r + 1.8),
      clip = "off"
    ) +
    theme_void() +
    theme(plot.margin = margin(1, 2, 2, 2, "mm"),
          legend.position = "right",
          legend.title = element_text(size = 7, face = "bold", color = "grey30"),
          legend.text  = element_text(size = 6, color = "grey40"),
          legend.key.width  = unit(2, "mm"),
          legend.key.height = unit(12, "mm"),
          legend.margin = margin(l = 0, r = 0)) +
    guides(color = "none",
           fill = guide_colorbar(direction = "vertical"))

  if (!is.null(title)) {
    p <- p + ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size = title_size, face = "bold",
                                      margin = margin(b = 1)))
  }

  # Contrast title/subtitle rendered inside the plot (above outermost ring)
  if (!is.null(contrast_title)) {
    p <- p + annotate("text", x = 0, y = label_r + 1.5,
                       label = contrast_title,
                       size = title_size / .pt, fontface = "bold", hjust = 0.5)
    if (!is.null(contrast_subtitle)) {
      p <- p + annotate("text", x = 0, y = label_r + 0.7,
                         label = contrast_subtitle,
                         size = (title_size - 3) / .pt, fontface = "italic",
                         color = "grey30", hjust = 0.5)
    }
  }

  attr(p, "ring_data")  <- ring_data
  attr(p, "tick_data")  <- tick_data
  attr(p, "n_up")       <- attr(volcano_layers, "n_up")
  attr(p, "n_down")     <- attr(volcano_layers, "n_down")

  p
}

