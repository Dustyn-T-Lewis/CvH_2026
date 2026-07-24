# Module-card columns for F05_WGCNA (Mito F03 layout). Each column is faceted by
# module on the same order with blank strips, so one row reads left to right as a
# single module: protein count, module-trait heatmap, member-response fry tiles,
# eigengene trajectory with post-hoc brackets, and top-5 ORA pathways.

pacman::p_load(ggplot2, dplyr, tidyr, stringr, patchwork, ggtext, ggfittext)

CTR_KEEP <- c("CRvH_Baseline", "CR_Training")
CTR_LABS <- c(CRvH_Baseline = "CR vs Ctl", CR_Training = "Training")

tint_colour <- function(hex, amt = 0.62) {
  v <- grDevices::col2rgb(hex) / 255
  grDevices::rgb(t(v + (1 - v) * amt))
}
fmt_fdr <- function(q) ifelse(is.na(q), "", ifelse(q < 0.01, "<.01", sprintf("%.2f", q)))

card_facets <- function(spacing = 1.0) {
  list(
    facet_wrap(~module_color, ncol = 1, drop = FALSE),
    theme(
      strip.text = element_blank(),
      panel.spacing.y = unit(spacing, "mm"),
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 9, colour = "grey15"),
      plot.subtitle = element_text(face = "italic", size = 6.6, colour = "grey40")
    )
  )
}

card_count <- function(mod_bio, modules) {
  d <- mod_bio |>
    filter(module_color %in% modules) |>
    mutate(module_color = factor(module_color, levels = modules))
  ggplot(d, aes(n_proteins, 1)) +
    geom_col(aes(fill = as.character(module_color)),
      width = 1, colour = "grey35", linewidth = 0.25, orientation = "y"
    ) +
    ggtext::geom_richtext(
      aes(x = 0, y = 1.75, label = sprintf(
        "<b>%s: %s</b> · %d", module_id, module_color, n_proteins
      )),
      hjust = 0, vjust = 1, size = 2.4, fill = NA, label.colour = NA,
      label.padding = unit(rep(0, 4), "pt")
    ) +
    scale_fill_identity() +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
    scale_y_continuous(limits = c(0.4, 1.85), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(title = "Proteins", subtitle = "count per module", x = NULL, y = NULL) +
    FIG_THEME +
    card_facets() +
    theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 6.5, colour = "grey35"),
      legend.position = "none", plot.margin = margin(2, 2, 2, 4)
    )
}

# MITO-style stacked member response: r_equiv module-trait tiles (top row) over
# fry tiles (bottom row, signed -log10 q), per contrast, with two colorbars.
card_response <- function(lmm, settests, modules) {
  fmt2 <- function(p) {
    ifelse(is.na(p), "", ifelse(p < 0.001, "<.001", sprintf("%.2f", p)))
  }
  as_card <- function(d) {
    mutate(d,
      module_color = factor(module_color, levels = modules),
      contrast = factor(contrast, levels = CTR_KEEP),
      row = factor(row, levels = c("trait", "fry"))
    )
  }
  trait <- lmm |>
    mutate(module_color = sub("^ME", "", module)) |>
    filter(module_color %in% modules, contrast %in% CTR_KEEP) |>
    transmute(module_color, contrast,
      row = "trait", value = r_equiv,
      sig = case_when(p_bh < 0.05 ~ "FDR", p_raw < 0.05 ~ "nominal", TRUE ~ "ns"),
      lab = sprintf("%.2f%s", r_equiv, ifelse(p_bh < 0.05, sig_stars(p_bh), ""))
    ) |>
    as_card()
  fry <- settests |>
    filter(module_color %in% modules, contrast %in% CTR_KEEP) |>
    transmute(module_color, contrast,
      row = "fry",
      value = if_else(direction == "Up", 1, -1) * -log10(pmax(fry_fdr, 1e-4)),
      sig = if_else(!is.na(fry_fdr) & fry_fdr < 0.05, "FDR", "ns"),
      lab = sprintf("p=%s\nq=%s", fmt2(fry_p), fmt2(fry_fdr))
    ) |>
    as_card()
  boxes <- bind_rows(trait, fry) |> filter(sig != "ns")

  ggplot(mapping = aes(contrast, row)) +
    geom_tile(
      data = trait, aes(fill = value), colour = "grey80",
      linewidth = 0.4, width = 0.96, height = 0.92
    ) +
    scale_fill_gradient2(
      low = "#1B7837", mid = "white", high = "#762A83",
      midpoint = 0, limits = c(-0.8, 0.8), oob = scales::squish,
      name = expression(r[equiv]),
      guide = guide_colorbar(
        order = 1, barwidth = unit(12, "mm"),
        barheight = unit(2, "mm"), title.vjust = 1
      )
    ) +
    ggnewscale::new_scale_fill() +
    geom_tile(
      data = fry, aes(fill = value), colour = "grey80",
      linewidth = 0.4, width = 0.96, height = 0.92
    ) +
    scale_fill_gradient2(
      low = "#4393C3", mid = "white", high = "#D6604D",
      midpoint = 0, limits = c(-2, 2), oob = scales::squish,
      name = "fry ±log10 q",
      guide = guide_colorbar(
        order = 2, barwidth = unit(12, "mm"),
        barheight = unit(2, "mm"), title.vjust = 1
      )
    ) +
    geom_tile(
      data = filter(boxes, sig == "nominal"), fill = NA, colour = "grey30",
      linetype = "dashed", linewidth = 0.5, width = 0.96, height = 0.92
    ) +
    geom_tile(
      data = filter(boxes, sig == "FDR"), fill = NA, colour = "black",
      linewidth = 0.9, width = 0.96, height = 0.92
    ) +
    geom_text(
      data = trait, aes(label = lab), size = 1.9,
      colour = if_else(abs(trait$value) > 0.5, "white", "grey15")
    ) +
    geom_text(
      data = fry, aes(label = lab), size = 1.5, lineheight = 0.9,
      colour = if_else(fry$sig == "FDR", "white", "grey30")
    ) +
    scale_x_discrete(labels = CTR_LABS) +
    scale_y_discrete(
      limits = c("fry", "trait"), labels = c(fry = "fry", trait = "r_eq")
    ) +
    labs(
      title = "Member response", subtitle = "r_equiv (top) · fry (below)",
      x = NULL, y = NULL
    ) +
    FIG_THEME +
    card_facets() +
    theme(
      axis.text.y = element_text(size = 5.5, face = "bold", colour = "grey30"),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 6.5, face = "bold", angle = 20, hjust = 1),
      legend.position = "bottom", legend.title = element_text(size = 6),
      legend.text = element_text(size = 5), legend.box = "vertical",
      legend.margin = margin(0, 0, 0, 0), legend.spacing.y = unit(0.5, "mm"),
      plot.margin = margin(2, 2, 2, 2)
    )
}

POSTHOC_SPANS <- tibble::tribble(
  ~contrast, ~x1, ~x2, ~level, ~span_label,
  "CR_Training", 2, 3, 1, "Training",
  "CRvH_Baseline", 1, 3, 2, "CR vs Ctl"
)

card_trajectory <- function(traj, lmm, modules) {
  x_lev <- c("Ctl", "CR_T1", "CR_T2")
  d <- traj |>
    filter(module_color %in% modules) |>
    mutate(
      module_color = factor(module_color, levels = modules),
      traj = factor(traj, levels = x_lev), gx = as.integer(traj),
      point_col = if_else(traj == "Ctl", unname(GROUP_COLORS["H_T1"]), as.character(module_color))
    )
  span <- d |>
    group_by(module_color) |>
    summarise(
      top = max(mean_eig + se),
      rng = pmax(max(mean_eig + se) - min(mean_eig - se), 1e-6), .groups = "drop"
    )
  brk <- lmm |>
    mutate(module_color = sub("^ME", "", module)) |>
    filter(module_color %in% modules, contrast %in% POSTHOC_SPANS$contrast, p_raw < 0.05) |>
    inner_join(POSTHOC_SPANS, by = "contrast") |>
    inner_join(span, by = "module_color") |>
    mutate(
      module_color = factor(module_color, levels = modules),
      yb = top + (0.10 + 0.30 * (level - 1)) * rng, ytick = yb - 0.05 * rng,
      ylab = yb + 0.02 * rng, xmid = (x1 + x2) / 2,
      col = if_else(p_bh < 0.05, "grey10", "grey45"),
      lab = sprintf("<b>%s</b><br>q=%.3f%s", span_label, p_bh, ifelse(p_bh < 0.05, " ✱", ""))
    )
  p <- ggplot(d, aes(gx, mean_eig)) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.2) +
    geom_errorbar(aes(ymin = mean_eig - se, ymax = mean_eig + se),
      width = 0.16, linewidth = 0.4, colour = "grey35"
    ) +
    geom_line(
      data = filter(d, traj != "Ctl"),
      aes(group = module_color, colour = as.character(module_color)), linewidth = 1.1
    ) +
    geom_point(aes(fill = point_col), shape = 21, size = 2, stroke = 0.4, colour = "grey20") +
    scale_colour_identity() +
    scale_fill_identity() +
    scale_x_continuous(
      breaks = seq_along(x_lev), labels = c("Ctl", "CR T1", "CR T2"),
      limits = c(0.55, 3.4), expand = expansion(mult = 0.03)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.48))) +
    labs(
      title = "Module eigengene",
      subtitle = "mean ± SE; brackets = post-hoc (✱ q < 0.05)", x = NULL, y = "eigengene (a.u.)"
    ) +
    FIG_THEME +
    card_facets() +
    theme(
      axis.text.x = element_text(face = "bold", size = 7),
      axis.text.y = element_text(size = 6, colour = "grey45"),
      legend.position = "none", plot.margin = margin(2, 2, 2, 2)
    )
  if (nrow(brk)) {
    p <- p +
      geom_segment(data = brk, aes(x = x1, xend = x2, y = yb, yend = yb, colour = I(col)), inherit.aes = FALSE, linewidth = 0.35) +
      geom_segment(data = brk, aes(x = x1, xend = x1, y = yb, yend = ytick, colour = I(col)), inherit.aes = FALSE, linewidth = 0.35) +
      geom_segment(data = brk, aes(x = x2, xend = x2, y = yb, yend = ytick, colour = I(col)), inherit.aes = FALSE, linewidth = 0.35) +
      ggtext::geom_richtext(
        data = brk, aes(x = xmid, y = ylab, label = lab, colour = I(col)), inherit.aes = FALSE,
        vjust = 0, size = 1.6, lineheight = 1.05, fill = NA, label.colour = NA,
        label.padding = unit(rep(0, 4), "pt")
      )
  }
  p
}

card_ora <- function(ora, modules, top_n = 5L) {
  d <- ora |>
    filter(module_color %in% modules) |>
    group_by(module_color) |>
    slice_min(padj, n = top_n, with_ties = FALSE) |>
    ungroup() |>
    mutate(
      module_color = factor(module_color, levels = modules),
      database = ifelse(database %in% names(DB_COLORS), database, "Other"),
      lp = -log10(padj), sig_fdr = padj < 0.05,
      fill_col = if_else(sig_fdr, as.character(module_color),
        vapply(as.character(module_color), tint_colour, character(1))
      ),
      name = clean_pathway_name(pathway),
      qlab = ifelse(padj < 0.01, "< 0.01", sprintf("%.2f", padj)),
      outside = padj > 0.75
    ) |>
    arrange(module_color, lp) |>
    mutate(rid = paste(module_color, dplyr::row_number()), row = factor(rid, levels = rid))
  db_of <- stats::setNames(as.character(d$database), as.character(d$row))
  db_label <- function(x) {
    sprintf("<span style='color:%s'>%s</span>", DB_COLORS[db_of[x]], db_of[x])
  }

  ggplot(d, aes(lp, row, fill = fill_col)) +
    geom_col(width = 0.86, colour = "grey35", linewidth = 0.2, orientation = "y") +
    ggfittext::geom_fit_text(
      data = ~ filter(.x, !outside), aes(xmin = 0, xmax = lp, label = name),
      reflow = TRUE, grow = FALSE, contrast = TRUE, fontface = "bold",
      size = 6, min.size = 3, padding.x = grid::unit(0.6, "mm"), padding.y = grid::unit(0.3, "mm")
    ) +
    geom_text(
      data = ~ filter(.x, outside),
      aes(x = lp, label = sprintf("%s   %s", qlab, stringr::str_trunc(name, 44))),
      hjust = -0.1, size = 1.95, colour = "grey30"
    ) +
    geom_text(
      data = ~ filter(.x, !outside),
      aes(x = lp, label = qlab, fontface = if_else(sig_fdr, "bold", "plain")),
      hjust = -0.2, size = 2.2, colour = "grey15"
    ) +
    facet_wrap(~module_color, ncol = 1, scales = "free", drop = FALSE) +
    scale_fill_identity() +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.16)),
      breaks = function(lims) Filter(function(b) b >= 0, scales::breaks_pretty(3)(lims))
    ) +
    scale_y_discrete(labels = db_label, expand = expansion(add = 0.5)) +
    coord_cartesian(clip = "off") +
    labs(
      title = sprintf("Top %d ORA pathways", top_n),
      subtitle = "over-representation (BH-FDR); bold q = FDR < 0.05",
      x = expression(-log[10] ~ FDR), y = NULL
    ) +
    FIG_THEME +
    theme(
      strip.text = element_blank(), panel.spacing.y = unit(1, "mm"), panel.grid = element_blank(),
      plot.title = element_text(face = "bold", size = 9, colour = "grey15"),
      plot.subtitle = element_text(face = "italic", size = 6.6, colour = "grey40"),
      axis.text.y = ggtext::element_markdown(size = 6, face = "bold", hjust = 1),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 6.5, colour = "grey35"),
      axis.title.x = element_text(size = 7, colour = "grey35"),
      legend.position = "none", plot.margin = margin(2, 3, 2, 1)
    )
}

assemble_card <- function(mod_bio, lmm, settests, traj, ora, modules) {
  wrap_plots(
    list(
      card_count(mod_bio, modules),
      card_response(lmm, settests, modules),
      card_trajectory(traj, lmm, modules),
      card_ora(ora, modules)
    ),
    nrow = 1, widths = c(0.55, 0.72, 1, 1.7)
  )
}
