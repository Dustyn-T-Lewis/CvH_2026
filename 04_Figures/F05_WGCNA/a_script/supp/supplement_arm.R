# F05 supplement: the module response with creatine and placebo kept apart.
#
# The main card pools both arms into one CR group, because only 6 subjects per arm are
# paired. This supplement shows what that pooling hides: the two arms' eigengene
# trajectories, and the CR-only 2x2 contrasts tested per module. At this n the arm
# comparison is badly underpowered, so it is presented as a completeness check rather
# than a test of the supplement.

render_supplement_arm <- function(dat, supp_png, supp_pdf, pdf_device) {
  traj <- read_csv(file.path(dat, "trajectory_eigengenes_supp.csv"), show_col_types = FALSE)
  settests <- read_csv(file.path(dat, "module_set_tests_supp.csv"), show_col_types = FALSE)
  nes <- read_csv(file.path(dat, "module_fgsea_nes_supp.csv"), show_col_types = FALSE)
  mod_bio <- read_csv(file.path(dat, "mod_bio_labels.csv"), show_col_types = FALSE)
  modules <- mod_bio |>
    arrange(desc(n_proteins)) |>
    pull(module_color)

  stage_levels <- c("Ctl", "pre", "post")
  arm_colors <- c(Ctl = "#4DAF4A", CRE = "#2166AC", PLA = "#D6604D")

  # The control group belongs to neither arm, so it is duplicated onto both lines to
  # give each arm a shared starting point rather than a floating point.
  ctl <- traj |> filter(arm == "Ctl")
  traj_lines <- bind_rows(
    mutate(ctl, arm_line = "CRE"),
    mutate(ctl, arm_line = "PLA"),
    traj |> filter(arm != "Ctl") |> mutate(arm_line = arm)
  ) |>
    filter(module_color %in% modules) |>
    mutate(
      module_color = factor(module_color, levels = modules),
      stage = factor(stage, levels = stage_levels),
      gx = as.integer(stage)
    )

  p_traj <- ggplot(traj_lines, aes(gx, mean_eig, colour = arm_line)) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.2) +
    geom_errorbar(aes(ymin = mean_eig - se, ymax = mean_eig + se),
      width = 0.12, linewidth = 0.35
    ) +
    geom_line(linewidth = 0.9) +
    geom_point(aes(fill = arm), shape = 21, size = 1.9, stroke = 0.4, colour = "grey20") +
    facet_wrap(~module_color, ncol = 1, scales = "free_y") +
    scale_colour_manual(values = arm_colors[c("CRE", "PLA")], name = NULL) +
    scale_fill_manual(values = arm_colors, guide = "none") +
    scale_x_continuous(
      breaks = seq_along(stage_levels), labels = c("Ctl", "pre", "post"),
      limits = c(0.8, 3.2)
    ) +
    labs(
      title = "Eigengene by supplement arm",
      subtitle = "mean ± SE; both arms share the control point",
      x = NULL, y = "eigengene (a.u.)"
    ) +
    FIG_THEME +
    theme(
      strip.text = element_text(face = "bold", size = 6.5),
      axis.text.x = element_text(face = "bold", size = 6.5),
      legend.position = "bottom", legend.text = element_text(size = 6)
    )

  ctr_levels <- c(
    "Baseline_Supplement", "Training_CRE", "Training_PLA", "Supplement_Interaction"
  )
  ctr_labs <- c(
    Baseline_Supplement = "CRE vs PLA\n(baseline)", Training_CRE = "Training\n(CRE)",
    Training_PLA = "Training\n(PLA)", Supplement_Interaction = "Interaction"
  )
  fmt2 <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<.001", sprintf("%.2f", p)))
  as_card <- function(d) {
    mutate(d,
      module_color = factor(module_color, levels = modules),
      contrast = factor(contrast, levels = ctr_levels),
      row = factor(row, levels = c("NES", "fry"))
    )
  }
  nes_d <- nes |>
    rename(module_color = module) |>
    filter(module_color %in% modules) |>
    transmute(module_color, contrast,
      row = "NES", value = NES, lab = sprintf("%.1f", NES),
      sig = !is.na(padj) & padj < 0.05
    ) |>
    as_card()
  fry_d <- settests |>
    filter(module_color %in% modules) |>
    transmute(module_color, contrast,
      row = "fry",
      value = if_else(direction == "Up", 1, -1) * -log10(pmax(fry_fdr, 1e-4)),
      lab = sprintf("p=%s\nq=%s", fmt2(fry_p), fmt2(fry_fdr)),
      sig = !is.na(fry_fdr) & fry_fdr < 0.05
    ) |>
    as_card()

  p_tiles <- ggplot(mapping = aes(contrast, row)) +
    geom_tile(data = nes_d, aes(fill = value), colour = "grey25", linewidth = 0.3) +
    scale_fill_gradientn(
      colours = c("#1B7837", "white", "#762A83"), limits = c(-3.5, 3.5),
      oob = scales::squish, name = "fGSEA NES"
    ) +
    ggnewscale::new_scale_fill() +
    geom_tile(data = fry_d, aes(fill = value), colour = "grey25", linewidth = 0.3) +
    scale_fill_gradient2(
      low = "#4393C3", mid = "white", high = "#D6604D", midpoint = 0,
      limits = c(-2, 2), oob = scales::squish, name = "fry ±log10 q"
    ) +
    geom_tile(
      data = bind_rows(nes_d, fry_d) |> filter(sig), fill = NA,
      colour = "black", linewidth = 0.9
    ) +
    geom_text(data = nes_d, aes(label = lab), size = 1.8, colour = "grey15") +
    geom_text(data = fry_d, aes(label = lab), size = 1.4, lineheight = 0.9, colour = "grey25") +
    facet_wrap(~module_color, ncol = 1) +
    scale_x_discrete(labels = ctr_labs) +
    scale_y_discrete(limits = c("fry", "NES")) +
    labs(
      title = "CR-only 2x2 contrasts", x = NULL, y = NULL,
      subtitle = "creatine vs placebo, 6 paired subjects per arm"
    ) +
    FIG_THEME +
    theme(
      strip.text = element_blank(),
      axis.text.x = element_text(size = 5.5, face = "bold", angle = 20, hjust = 1),
      axis.text.y = element_text(size = 5.5, face = "bold"),
      legend.position = "bottom", legend.title = element_text(size = 6),
      legend.text = element_text(size = 5), legend.box = "vertical"
    )

  fig <- (p_traj | p_tiles) +
    plot_layout(widths = c(1, 0.85)) +
    plot_annotation(
      title = "Module response split by supplement arm",
      subtitle = paste(
        "The main figure pools creatine and placebo into one CR group; 6 paired subjects per arm.",
        "\nNo fry contrast reaches FDR < 0.05. NES boxes mark the competitive test, which large",
        "modules pass easily and which does not on its own show an arm effect."
      ),
      theme = theme(
        plot.title = element_text(face = "bold", size = 13, colour = "grey10"),
        plot.subtitle = element_text(face = "italic", size = 8, colour = "grey40")
      )
    )

  ggsave(file.path(supp_png, "SUPP_F05_supplement_arm.png"), fig,
    width = 220, height = 42 + 32 * length(modules), units = "mm", dpi = 300, bg = "white"
  )
  ggsave(file.path(supp_pdf, "SUPP_F05_supplement_arm.pdf"), fig,
    width = 220, height = 42 + 32 * length(modules), units = "mm", device = pdf_device
  )
  fig
}
