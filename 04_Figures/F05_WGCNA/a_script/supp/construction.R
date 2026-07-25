# F05 supplements: how the network was built, and how the module-phenotype association
# holds up across the primary and sensitivity arms. Sourced by 02_clustering.R, which
# owns the paths and the cached tables.

render_construction_supp <- function(dat, supp_png, supp_pdf, pdf_device) {
  w <- readRDS(file.path(dat, "wgcna_network.rds"))
  sft_df <- readRDS(file.path(dat, "wgcna", "sft_fitIndices.rds"))
  p <- panel_scale_free(sft_df, w$chosen_power) +
    panel_dendro(w$net) +
    plot_annotation(
      title = "WGCNA network construction",
      subtitle = "Module preservation not computed: single-cohort network, no external validation set.",
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    ) +
    plot_layout(widths = c(1, 1.4))
  ggsave(file.path(supp_png, "SUPP_F05_construction.png"), p,
    width = 240, height = 110, units = "mm", dpi = 300, bg = "white"
  )
  ggsave(file.path(supp_pdf, "SUPP_F05_construction.pdf"), p,
    width = 240, height = 110, units = "mm", device = pdf_device
  )
  p
}

short_trait <- function(x) {
  recode(sub("^(pre_|post_)", "", x),
    grip_lbs = "Grip", ALM_kg = "ALM", sts_max_pwr = "STS pwr", LBM_kg = "LBM",
    chest_press_lbs = "Chest", leg_ext_lbs = "Leg ext", age = "Age"
  )
}

pheno_heat <- function(df, value_col, title, subtitle, lim) {
  d <- df |>
    filter(module != "grey") |>
    mutate(
      trait = short_trait(trait), value = .data[[value_col]],
      sig = !is.na(padj) & padj < 0.05, lab = sprintf("%.2f", value)
    )
  ggplot(d, aes(trait, module, fill = value)) +
    geom_tile(colour = "grey85", linewidth = 0.3) +
    geom_tile(data = filter(d, sig), fill = NA, colour = "black", linewidth = 0.8) +
    geom_text(aes(label = lab),
      size = 2.2,
      colour = if_else(abs(d$value) > lim * 0.6, "white", "grey15")
    ) +
    scale_fill_gradient2(
      low = "#1B7837", mid = "white", high = "#762A83", midpoint = 0,
      limits = c(-lim, lim), oob = scales::squish
    ) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL, fill = value_col) +
    FIG_THEME +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
}

render_phenotype_supp <- function(cor_pheno, cor_matched, lmm_pheno,
                                  supp_png, supp_pdf, pdf_device) {
  p <- pheno_heat(cor_pheno, "r", "Baseline (primary)", "T1 only, pre_ outcomes, Pearson", 1) +
    pheno_heat(cor_matched, "r", "Matched Pearson", "all samples, T1→pre_/T2→post_", 1) +
    pheno_heat(lmm_pheno, "beta", "Matched LMM", "std β, (1|subject)", 1) +
    plot_annotation(
      title = "Module–phenotype association: primary vs sensitivity arms",
      subtitle = "Boxes mark BH-FDR < 0.05.",
      theme = theme(plot.title = element_text(face = "bold", size = 13))
    )
  ggsave(file.path(supp_png, "SUPP_F05_phenotype_methods.png"), p,
    width = 280, height = 110, units = "mm", dpi = 300, bg = "white"
  )
  ggsave(file.path(supp_pdf, "SUPP_F05_phenotype_methods.pdf"), p,
    width = 280, height = 110, units = "mm", device = pdf_device
  )
  p
}
