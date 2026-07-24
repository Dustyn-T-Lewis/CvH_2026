# Construction supplement panels: scale-free fit vs soft power and the module
# dendrogram. Single-cohort network, so module preservation (an external
# validation cohort) is not computed; noted in the figure caption.

pacman::p_load(ggplot2, dplyr, WGCNA, ggplotify)

panel_scale_free <- function(sft_df, chosen_power, rsq_cut = 0.87) {
  d <- sft_df |>
    mutate(r2 = -sign(slope) * SFT.R.sq, chosen = Power == chosen_power)
  ggplot(d, aes(Power, r2)) +
    geom_hline(yintercept = rsq_cut, linetype = "dashed", colour = "grey55") +
    geom_line(colour = "grey70", linewidth = 0.3) +
    geom_text(aes(label = Power, colour = chosen), size = 2.6, fontface = "bold") +
    scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey35"), guide = "none") +
    labs(
      title = "Scale-free topology fit",
      subtitle = sprintf("soft power = %d at signed R² > %.2f", chosen_power, rsq_cut),
      x = "Soft threshold (power)", y = expression(signed ~ R^2)
    ) +
    FIG_THEME
}

panel_dendro <- function(net) {
  colors <- WGCNA::labels2colors(net$colors)
  ggplotify::as.ggplot(function() {
    WGCNA::plotDendroAndColors(
      net$dendrograms[[1]],
      colors[net$blockGenes[[1]]],
      groupLabels = "Module",
      dendroLabels = FALSE, addGuide = TRUE,
      hang = 0.03, guideHang = 0.05, main = ""
    )
  })
}
