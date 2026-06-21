# Figure 6 — Panel B: WGCNA Per-Module Triptych
# Layout per row: z-score heatmap (5 groups) | eigengene dynamics | ORA bars
# Top 4 key modules ordered by LMM significance
# 5 group columns: CRE_T1, CRE_T2, PLA_T1, PLA_T2, H_T1
# Eigengene: ME mean +/- SE per group_time
# ORA bars use per-module x-scales (no artificial cap)

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F06/a_script/style.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

RPT <- "04_Figures/F06/b_reports"

RPT_PDF       <- file.path(RPT, "main", "pdf")

RPT_PNG       <- file.path(RPT, "main", "png")

RPT_SUPP_PDF  <- file.path(RPT, "supp", "pdf")

RPT_SUPP_PNG  <- file.path(RPT, "supp", "png")
dir.create(RPT_PDF,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PNG, recursive = TRUE, showWarnings = FALSE)
DAT <- "04_Figures/F06/c_data"

pdf_device <- get_pdf_device()

message("Panel B: WGCNA per-module triptych...")

# --- Load data ---
z_scores  <- read_csv(file.path(DAT, "03_panel_B_heatmap_zscores.csv"), show_col_types = FALSE)
me_data   <- read_csv(file.path(DAT, "03_panel_B_eigengene_data.csv"), show_col_types = FALSE)
enrich    <- read_csv(file.path(DAT, "03_panel_B_triptych_enrichment.csv"), show_col_types = FALSE)
mod_bio   <- read_csv(file.path(DAT, "mod_bio_labels.csv"), show_col_types = FALSE)
lmm_audit <- read_csv(file.path(DAT, "wgcna_lmm_contrast_audit.csv"), show_col_types = FALSE)

if (!"display_label" %in% colnames(mod_bio)) {
  mod_bio <- mod_bio %>%
    mutate(display_label = paste0(module_id, ": ", module_color, " (n=", n_proteins, ")"))
}
mod_labels <- setNames(mod_bio$display_label, mod_bio$module_color)

# Key modules (from file or fallback)
km_file <- file.path(DAT, "key_modules.txt")
KEY_MODULES <- if (file.exists(km_file)) {
  readLines(km_file)
} else {
  head(mod_bio$module_color, 4)
}
KEY_MODULES <- KEY_MODULES[nzchar(trimws(KEY_MODULES))]

# Order by module size (largest first) — fallback ordering
mod_order <- mod_bio %>%
  filter(module_color %in% KEY_MODULES) %>%
  arrange(desc(n_proteins)) %>%
  pull(module_color)

# --- Dimensions ---
PB_W <- 300
PB_H <- 320

txt_heat  <- scale_text(BASE_GENE, PB_W) * 0.7
txt_axis  <- scale_text(BASE_STAT, PB_W) * 1.0
txt_title <- scale_text(BASE_GENE, PB_W) * 1.3
txt_bar   <- scale_text(BASE_GENE, PB_W) * 0.95
txt_sig   <- scale_text(BASE_GENE, PB_W) * 0.85

# Group ordering: 5 group_time levels
group_order  <- c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1")
group_labels <- c(CRE_T1 = "CRE T1", CRE_T2 = "CRE T2",
                  PLA_T1 = "PLA T1", PLA_T2 = "PLA T2",
                  H_T1   = "H T1")

# LMM stats for bracket annotations (key contrasts)
lmm_stats <- lmm_audit %>%
  filter(contrast %in% c("Cancer_vs_Healthy", "Training_CR",
                          "Training_CRE", "Training_PLA")) %>%
  mutate(module = gsub("^ME", "", module)) %>%
  dplyr::select(module, contrast, p_bh)

# --- Build one triptych row per module ---
build_row <- function(mod, show_xlab = FALSE) {
  label <- if (mod %in% names(mod_labels)) mod_labels[mod] else mod
  n_mod <- z_scores %>% filter(module == mod) %>% distinct(gene) %>% nrow()
  title_txt <- paste0(label, " (n=", n_mod, ")")

  # -- Heatmap: z-scores (5 group_time columns) --
  z_mod <- z_scores %>%
    filter(module == mod) %>%
    mutate(group = factor(group, levels = group_order))

  # Order genes by CRE_T1 z-score (baseline reference)
  gene_order <- z_mod %>%
    filter(group == "CRE_T1") %>%
    arrange(z) %>%
    pull(gene)
  z_mod$gene <- factor(z_mod$gene, levels = gene_order)

  p_heat <- ggplot(z_mod, aes(x = group, y = gene, fill = z)) +
    geom_tile() +
    scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                         midpoint = 0, limits = c(-2, 2), oob = scales::squish,
                         guide = "none") +
    scale_x_discrete(labels = group_labels, position = "bottom") +
    labs(title = title_txt, y = NULL, x = NULL) +
    FIG_THEME +
    theme(
      plot.title   = element_text(size = txt_title, face = "bold"),
      axis.text.x  = if (show_xlab) element_text(size = txt_axis * 0.75, angle = 45, hjust = 1)
                      else element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks   = element_blank(),
      panel.border = element_blank(),
      plot.margin  = margin(2, 1, 2, 2)
    )

  # -- Eigengene dynamics: group means +/- SE --
  me_mod <- me_data %>%
    filter(module == mod) %>%
    mutate(group_time = factor(group_time, levels = group_order))

  # Summary: mean +/- SE per group_time
  me_summary <- me_mod %>%
    group_by(group_time) %>%
    summarise(
      mean_me = mean(eigengene, na.rm = TRUE),
      se_me   = sd(eigengene, na.rm = TRUE) / sqrt(sum(!is.na(eigengene))),
      .groups = "drop"
    )

  # LMM p-values for the 3-group contrasts
  p_cvh <- lmm_stats %>% filter(module == mod, contrast == "Cancer_vs_Healthy") %>% pull(p_bh)
  p_tr  <- lmm_stats %>% filter(module == mod, contrast == "Training_CR") %>% pull(p_bh)

  fmt_sig <- function(p) {
    if (length(p) == 0 || is.na(p)) return("ns")
    if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
  }

  y_range <- range(c(me_summary$mean_me - me_summary$se_me,
                      me_summary$mean_me + me_summary$se_me), na.rm = TRUE)

  p_eigen <- ggplot(me_summary, aes(x = group_time, y = mean_me)) +
    geom_errorbar(aes(ymin = mean_me - se_me, ymax = mean_me + se_me),
                  width = 0.2, linewidth = 0.5) +
    geom_point(aes(color = group_time), size = 3) +
    geom_line(data = me_summary %>% filter(group_time %in% c("CRE_T1", "CRE_T2")),
              aes(group = 1), color = GROUP_COLORS["CRE_T1"], linewidth = 0.8) +
    geom_line(data = me_summary %>% filter(group_time %in% c("PLA_T1", "PLA_T2")),
              aes(group = 1), color = GROUP_COLORS["PLA_T1"], linewidth = 0.8) +
    scale_color_manual(values = GROUP_COLORS, guide = "none") +
    # Significance annotations (3-group contrasts: pooled training + cancer vs healthy)
    annotate("text", x = 2.5, y = y_range[2] * 1.05,
             label = paste0("Tr.(CR) ", fmt_sig(p_tr)), size = txt_sig, fontface = "bold",
             color = "grey25") +
    annotate("text", x = 5, y = y_range[2] * 1.05,
             label = paste0("CRvH ", fmt_sig(p_cvh)), size = txt_sig, fontface = "bold",
             color = "grey25") +
    labs(y = "Eigengene", x = NULL) +
    FIG_THEME +
    theme(
      axis.text.x  = if (show_xlab) element_text(size = txt_axis * 0.75, angle = 45, hjust = 1)
                      else element_blank(),
      axis.text.y  = element_text(size = txt_axis * 0.75),
      axis.title.y = element_text(size = txt_axis * 0.8),
      panel.border = element_blank(),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
      plot.margin  = margin(2, 1, 2, 1)
    )

  # -- ORA bars: top 5, per-module x-scale --
  bar_data <- enrich %>%
    filter(module == mod, p.adjust < 0.05) %>%
    arrange(p.adjust) %>%
    head(5) %>%
    mutate(
      neg_log10_p = -log10(p.adjust),
      clean_name  = clean_pathway_name(Description),
      db_fill     = DB_COLORS[database]
    ) %>%
    mutate(clean_name = factor(clean_name, levels = rev(clean_name)))

  if (nrow(bar_data) == 0) {
    p_bars <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "No sig.\nenrichment",
               size = txt_bar, color = "grey50") +
      theme_void() + theme(plot.margin = margin(2, 2, 2, 1))
  } else {
    p_bars <- ggplot(bar_data, aes(x = neg_log10_p, y = clean_name)) +
      geom_col(aes(fill = db_fill), color = "black", linewidth = 0.3, width = 0.7) +
      geom_text(aes(label = clean_name, x = 0.3), hjust = 0, size = txt_bar,
                fontface = "bold",
                color = ifelse(bar_data$db_fill %in% c("#AA336A", "#1565C0", "#00796B"),
                               "white", "grey20")) +
      scale_fill_identity() +
      scale_x_continuous(
        expand = expansion(mult = c(0, 0.05)),
        breaks = scales::breaks_pretty(n = 3),
        name = if (show_xlab) expression(-log[10](p[adj])) else NULL
      ) +
      scale_y_discrete(labels = NULL) +
      labs(y = NULL) +
      FIG_THEME +
      theme(
        axis.text.x  = if (show_xlab) element_text(size = txt_axis * 0.75) else element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        axis.line.x  = element_line(color = "black", linewidth = 0.3),
        panel.grid   = element_blank(),
        plot.margin  = margin(2, 2, 2, 1)
      )
  }

  # Combine row
  p_heat + p_eigen + p_bars + plot_layout(widths = c(3, 2, 3))
}

# --- Build all rows ---
rows <- lapply(seq_along(mod_order), function(i) {
  build_row(mod_order[i], show_xlab = (i == length(mod_order)))
})

# Z-score legend
z_legend <- ggplot(data.frame(z = seq(-2, 2, length.out = 100)),
                   aes(x = z, y = 1, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                       midpoint = 0, limits = c(-2, 2),
                       name = "Z-score", guide = guide_colorbar(
                         barwidth = unit(40, "mm"), barheight = unit(3, "mm"))) +
  theme_void() +
  theme(legend.position = "bottom", legend.text = element_text(size = txt_axis * 0.7))

# Group color legend
group_leg_df <- tibble(
  group = factor(group_order, levels = group_order),
  y = 1, x = seq_along(group_order)
)
group_legend <- ggplot(group_leg_df, aes(x = x, y = y, fill = group)) +
  geom_tile() +
  scale_fill_manual(values = GROUP_COLORS, labels = group_labels, name = "Group") +
  theme_void() +
  theme(legend.position = "bottom", legend.text = element_text(size = txt_axis * 0.7))

# --- Assemble ---
triptych <- wrap_plots(rows, ncol = 1) /
  (wrap_elements(z_legend) | wrap_elements(group_legend)) +
  plot_layout(heights = c(rep(1, length(mod_order)), 0.12))

ggsave(file.path(RPT_PDF, "panel_B_triptych_MAIN.pdf"), triptych,
       width = PB_W, height = PB_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT_PNG, "panel_B_triptych_MAIN.png"), triptych,
       width = PB_W, height = PB_H, units = "mm",
       dpi = 300, limitsize = FALSE)

message("  Panel B (WGCNA triptych) saved")
