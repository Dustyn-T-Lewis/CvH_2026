# Figure 8 — Panel A: Module-Trait Heatmap (Two-Model LMM Contrasts)
# Layout: gene counts (left) | CRvH model (2 cols) | CR model (4 cols)
# CvH has no continuous phenotype traits — LMM contrasts only
# BH correction: per-model (CRvH: 2*M tests, CR: 4*M tests)
# Two-tier display: solid border = FDR < 0.05; dot = nominal p < 0.05
# Generates: panel_A_module_trait_MAIN.pdf/.png, c_data/01_panel_A_*.csv

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F08/a_script/style.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(patchwork)
})

RPT <- "04_Figures/F08/b_reports"
DAT <- "04_Figures/F08/c_data"
SRC <- "04_Figures/F06/c_data"  # panel-ready data from 05_WGCNA
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

message("Panel A: Module-trait heatmap (CvH)...")

# --- Data loading ---
lmm_audit <- read_csv(file.path(SRC, "wgcna_lmm_contrast_audit.csv"),
                       show_col_types = FALSE)
module_df      <- read_csv(file.path(SRC, "wgcna_module_assignments.csv"), show_col_types = FALSE)
mod_bio_labels <- read_csv(file.path(SRC, "mod_bio_labels.csv"), show_col_types = FALSE)

if (!"display_label" %in% colnames(mod_bio_labels)) {
  mod_bio_labels <- mod_bio_labels %>%
    mutate(display_label = paste0(module_id, ": ", module_color, " (n=", n_proteins, ")"))
}
mod_display_vec <- setNames(mod_bio_labels$display_label, paste0("ME", mod_bio_labels$module_color))

# --- Build r_equiv and p-value matrices ---
crvh_contrasts <- c("Cancer_vs_Healthy", "Training_CR")
cr_contrasts   <- c("Baseline_Supplement", "Training_CRE", "Training_PLA", "Supplement_Interaction")
all_contrasts  <- c(crvh_contrasts, cr_contrasts)

lmm_r <- lmm_audit %>%
  dplyr::select(module, contrast, r_equiv) %>%
  pivot_wider(names_from = contrast, values_from = r_equiv) %>%
  column_to_rownames("module") %>% as.matrix()

lmm_p <- lmm_audit %>%
  dplyr::select(module, contrast, p_bh) %>%
  pivot_wider(names_from = contrast, values_from = p_bh) %>%
  column_to_rownames("module") %>% as.matrix()

lmm_p_raw <- lmm_audit %>%
  dplyr::select(module, contrast, p_raw) %>%
  pivot_wider(names_from = contrast, values_from = p_raw) %>%
  column_to_rownames("module") %>% as.matrix()

# Non-grey modules, ordered by size (largest at top)
non_grey <- mod_bio_labels %>%
  arrange(desc(n_proteins)) %>%
  pull(module_color) %>%
  paste0("ME", .)
non_grey <- intersect(non_grey, rownames(lmm_r))

# --- Dimensions ---
PA_W <- 320
PA_H <- 220

txt_cell   <- scale_text(BASE_GENE, PA_W) * 0.85
txt_count  <- scale_text(BASE_COUNT, PA_W) * 0.85
txt_axis   <- scale_text(BASE_STAT, PA_W) * 1.5
txt_brack  <- scale_text(BASE_GENE, PA_W) * 1.15

# --- Column positions with gap between models ---
gap <- 0.4
xpos <- c(
  1, 2,                                               # CRvH (2)
  2 + gap + 1, 2 + gap + 2, 2 + gap + 3, 2 + gap + 4 # CR (4)
)
trait_xpos <- setNames(xpos, all_contrasts)

# --- Contrast display labels ---
col_labels <- c(
  Cancer_vs_Healthy      = "CR vs H",
  Training_CR            = "Tr.(CR)",
  Baseline_Supplement    = "BL(CRE-PLA)",
  Training_CRE           = "Tr.(CRE)",
  Training_PLA           = "Tr.(PLA)",
  Supplement_Interaction = "CRE\u00d7PLA"
)

# --- Build tile data ---
heat_df <- expand.grid(
  module = non_grey, trait = all_contrasts,
  stringsAsFactors = FALSE
) %>%
  mutate(
    cor = mapply(function(m, c) {
      if (m %in% rownames(lmm_r) && c %in% colnames(lmm_r)) lmm_r[m, c] else NA_real_
    }, module, trait),
    pval = mapply(function(m, c) {
      if (m %in% rownames(lmm_p) && c %in% colnames(lmm_p)) lmm_p[m, c] else NA_real_
    }, module, trait),
    pval_raw = mapply(function(m, c) {
      if (m %in% rownames(lmm_p_raw) && c %in% colnames(lmm_p_raw)) lmm_p_raw[m, c] else NA_real_
    }, module, trait),
    stars = sig_stars(pval),
    label = sprintf("%.2f%s", cor, ifelse(stars == "ns", "", paste0("\n", stars))),
    trait_label = col_labels[trait],
    xpos  = trait_xpos[trait],
    model = ifelse(trait %in% crvh_contrasts, "CRvH", "CR")
  )

heat_df$module <- factor(heat_df$module, levels = rev(non_grey))

mod_color_raw <- setNames(gsub("^ME", "", non_grey), non_grey)

# --- Gene counts bar plot ---
mod_counts <- module_df %>%
  filter(module_color != "grey") %>%
  count(module_color, name = "n_proteins") %>%
  mutate(module = paste0("ME", module_color)) %>%
  filter(module %in% non_grey) %>%
  mutate(module = factor(module, levels = rev(non_grey)))

p_counts <- ggplot(mod_counts, aes(x = n_proteins, y = module)) +
  geom_col(fill = mod_counts$module_color, color = "black",
           linewidth = 0.3, width = 0.65) +
  geom_text(aes(label = n_proteins, x = n_proteins / 2),
            size = txt_count, fontface = "bold",
            color = ifelse(mod_counts$module_color %in% LIGHT_MODULES,
                           "grey30", "white")) +
  scale_x_reverse(expand = c(0, 0), limits = c(max(mod_counts$n_proteins) * 1.02, 0)) +
  scale_y_discrete(labels = function(x) {
    lbl <- mod_display_vec[x]
    ifelse(is.na(lbl), x, lbl)
  }) +
  labs(y = NULL, x = "Protein Count") +
  FIG_THEME +
  theme(axis.text.y        = element_text(size = txt_cell * 1.3, face = "bold"),
        axis.ticks.y       = element_blank(),
        axis.text.x        = element_text(size = txt_axis),
        axis.title.x       = element_text(size = txt_axis),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        panel.border       = element_blank(),
        axis.line.x        = element_line(color = "black", linewidth = 0.3),
        legend.position    = "none",
        plot.margin        = margin(2, 0, 2, 2))

# --- Model brackets ---
xmin_all <- min(xpos) - 0.55
xmax_all <- max(xpos) + 0.55

n_crvh <- sum(lmm_audit$model == "CRvH") / length(non_grey)
n_cr   <- sum(lmm_audit$model == "CR") / length(non_grey)

brackets <- tribble(
  ~label,                          ~start,  ~end,
  "CRvH Model\n(n=39, all)",     xpos[1], xpos[2],
  "CR Model\n(n\u224828, CR only)", xpos[3], xpos[6]
)
brackets <- brackets %>% mutate(mid = (start + end) / 2)

p_brackets <- ggplot() +
  geom_segment(data = brackets,
               aes(x = start - 0.4, xend = end + 0.4, y = 0.28, yend = 0.28),
               linewidth = 0.4, color = "grey30") +
  geom_segment(data = brackets,
               aes(x = start - 0.4, xend = start - 0.4, y = 0.28, yend = 0.14),
               linewidth = 0.4, color = "grey30") +
  geom_segment(data = brackets,
               aes(x = end + 0.4, xend = end + 0.4, y = 0.28, yend = 0.14),
               linewidth = 0.4, color = "grey30") +
  geom_text(data = brackets,
            aes(x = mid, y = 0.62, label = label),
            size = txt_brack, fontface = "bold", color = "grey25",
            lineheight = 0.85) +
  scale_x_continuous(limits = c(xmin_all, xmax_all), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(2, 2, 0, 0))

# --- Heatmap (two-tier significance) ---
p_heat <- ggplot(heat_df, aes(x = xpos, y = module, fill = cor)) +
  # Section background shading
  annotate("rect", xmin = xpos[1] - 0.5, xmax = xpos[2] + 0.5,
           ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Cancer_vs_Healthy"], alpha = 0.04) +
  annotate("rect", xmin = xpos[3] - 0.5, xmax = xpos[6] + 0.5,
           ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Training_CRE"], alpha = 0.04) +
  geom_tile(color = "black", linewidth = 0.3) +
  # Red border for highest r per row, blue for lowest
  geom_tile(data = heat_df %>% group_by(module) %>%
              filter(cor == max(cor, na.rm = TRUE)) %>% slice(1) %>% ungroup(),
            color = "#B2182B", linewidth = 1.3, fill = NA) +
  geom_tile(data = heat_df %>% group_by(module) %>%
              filter(cor == min(cor, na.rm = TRUE)) %>% slice(1) %>% ungroup(),
            color = "#2166AC", linewidth = 1.3, fill = NA) +
  # Solid thick border for BH FDR < 0.05
  geom_tile(data = heat_df %>% filter(pval < 0.05),
            color = "black", linewidth = 1.0, fill = NA) +
  # Dot for nominal p < 0.05 (but FDR >= 0.05)
  geom_point(data = heat_df %>% filter(pval_raw < 0.05 & pval >= 0.05),
             shape = 16, size = 0.8, color = "grey30") +
  # Text: bold white for FDR sig, plain for non-sig
  geom_text(data = heat_df %>% filter(pval >= 0.05),
            aes(label = label), size = txt_cell, color = "black", lineheight = 0.85) +
  geom_text(data = heat_df %>% filter(pval < 0.05),
            aes(label = label), size = txt_cell * 1.15,
            fontface = "bold", color = "white", lineheight = 0.85) +
  scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                       midpoint = 0, limits = c(-0.8, 0.8),
                       oob = scales::squish,
                       name = expression(r[equiv]),
                       na.value = "white",
                       guide = guide_colorbar(
                         barwidth = unit(60, "mm"),
                         barheight = unit(4, "mm"),
                         title.position = "left",
                         title.vjust = 0.75
                       )) +
  scale_x_continuous(
    breaks = xpos,
    labels = col_labels[all_contrasts],
    limits = c(xmin_all, xmax_all),
    expand = c(0, 0)
  ) +
  scale_y_discrete(labels = NULL) +
  labs(x = NULL, y = NULL) +
  FIG_THEME +
  theme(axis.text.x        = element_text(angle = 45, hjust = 1,
                                           size = txt_cell * 1.3, face = "bold"),
        axis.text.y        = element_blank(),
        axis.ticks.y       = element_blank(),
        axis.ticks.x       = element_blank(),
        panel.grid         = element_blank(),
        panel.border       = element_blank(),
        legend.position    = "bottom",
        legend.title       = element_text(size = txt_axis, face = "bold"),
        legend.text        = element_text(size = txt_cell),
        plot.margin        = margin(0, 2, 2, -3))

# --- Layout assembly ---
design <- c(
  area(1, 4, 1, 12),   # brackets (above heatmap)
  area(2, 1, 10, 3),   # gene counts
  area(2, 4, 10, 12)   # heatmap
)

pA <- wrap_elements(p_brackets) + p_counts + p_heat +
  plot_layout(design = design)

n_mods <- length(non_grey)
n_tests <- nrow(lmm_audit)

pA <- pA +
  plot_annotation(
    title    = "WGCNA Module-Trait Associations (CvH)",
    subtitle = paste0(n_mods, " modules (signed Pearson) | LMM: lmer + KR df, BH-corrected per model"),
    caption  = paste0("Values: r-equiv. | BH per-model (CRvH: ", n_mods * 2,
                      " tests; CR: ", n_mods * 4, " tests) | ",
                      "Solid border = FDR < .05 | Dot = nominal p < .05 | ",
                      "Red border = max r per row | Blue border = min r per row"),
    theme = theme(
      plot.title    = element_text(face = "bold", size = txt_axis * 2.5),
      plot.subtitle = element_text(size = txt_cell * 2, color = "grey30",
                                   face = "italic"),
      plot.caption  = element_text(size = txt_cell * 1.5, color = "grey40",
                                   hjust = 0),
      plot.margin   = margin(2, 2, 2, 2)
    )
  )

# --- Save ---
write_csv(heat_df %>% dplyr::select(module, trait, model, cor, pval, pval_raw, stars),
          file.path(DAT, "01_panel_A_heatmap_data.csv"))

ggsave(file.path(RPT, "panel_A_module_trait_MAIN.pdf"), pA,
       width = PA_W, height = PA_H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "panel_A_module_trait_MAIN.png"), pA,
       width = PA_W, height = PA_H, units = "mm",
       dpi = 300, limitsize = FALSE)

# --- Supplementary stat audit (CIs for LMM contrasts) ---
ci_list <- list()
for (i in seq_len(nrow(lmm_audit))) {
  row <- lmm_audit[i, ]
  if (!(row$module %in% non_grey)) next
  ci_est <- row$estimate + c(-1, 1) * qt(0.975, row$df) * row$SE
  ci_list <- c(ci_list, list(tibble(
    module   = row$module,
    contrast = row$contrast,
    model    = row$model,
    n        = ifelse(row$model == "CRvH", 39L, 28L),
    r_equiv  = round(row$r_equiv, 4),
    estimate = round(row$estimate, 5),
    ci_lo    = round(ci_est[1], 4),
    ci_hi    = round(ci_est[2], 4),
    p_raw    = row$p_raw,
    p_bh     = row$p_bh
  )))
}
ci_df <- bind_rows(ci_list)
write_csv(ci_df, file.path(DAT, "01_panel_A_correlation_CIs.csv"))

message("  Panel A saved")
