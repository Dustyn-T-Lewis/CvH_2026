# Figure 6 — Panel A: Module-Trait Heatmap (Two-Section LMM Contrasts)
# Layout: gene counts (flanking bar) | 3-group LMM contrasts (CR vs H, Training CR)
# No baseline/change trait sections (CvH has no phenotype columns)
# BH correction: per-model (CRvH 2*M tests, CR 4*M tests)
# Two-tier display: solid border = FDR < 0.05; dashed = nominal p < 0.05
# Generates: panel_A_module_trait_MAIN.pdf/.png, c_data/01_panel_A_*.csv

setwd(here::here())
source("04_Figures/F06/a_script/style.R")

pacman::p_load(readr, dplyr, tidyr, tibble, stringr, patchwork)

RPT <- "04_Figures/F06/b_reports"

RPT_PDF       <- file.path(RPT, "main", "pdf")

RPT_PNG       <- file.path(RPT, "main", "png")

RPT_SUPP_PDF  <- file.path(RPT, "supp", "pdf")

RPT_SUPP_PNG  <- file.path(RPT, "supp", "png")
DAT <- "04_Figures/F06/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

message("Panel A: Module-trait heatmap...")

# --- Data loading ---

# LMM contrasts (lmer eigengene ~ group_time + (1|subject), Kenward-Roger df)
lmm_audit <- read_csv(file.path(DAT, "wgcna_lmm_contrast_audit.csv"),
                       show_col_types = FALSE)

# Module metadata
module_df      <- read_csv(file.path(DAT, "wgcna_module_assignments.csv"), show_col_types = FALSE)
mod_bio_labels <- read_csv(file.path(DAT, "mod_bio_labels.csv"), show_col_types = FALSE)

if (!"display_label" %in% colnames(mod_bio_labels)) {
  mod_bio_labels <- mod_bio_labels %>%
    mutate(display_label = paste0(module_id, ": ", module_color, " (n=", n_proteins, ")"))
}
mod_display_vec <- setNames(mod_bio_labels$display_label, paste0("ME", mod_bio_labels$module_color))

# --- Reshape LMM to matrices ---
# 3-group model contrasts (supplements pooled; no CR supplement/interaction model)
crvh_contrasts <- c("Cancer_vs_Healthy", "Training_CR")

# r_equiv matrix
lmm_r <- lmm_audit %>%
  dplyr::select(module, contrast, r_equiv) %>%
  pivot_wider(names_from = contrast, values_from = r_equiv) %>%
  column_to_rownames("module") %>% as.matrix()

# BH p-value matrix
lmm_p <- lmm_audit %>%
  dplyr::select(module, contrast, p_bh) %>%
  pivot_wider(names_from = contrast, values_from = p_bh) %>%
  column_to_rownames("module") %>% as.matrix()

# Raw p-value matrix (for nominal tier)
lmm_p_raw <- lmm_audit %>%
  dplyr::select(module, contrast, p_raw) %>%
  pivot_wider(names_from = contrast, values_from = p_raw) %>%
  column_to_rownames("module") %>% as.matrix()

# Module ordering: by size (largest at top), exclude grey
mod_order <- mod_bio_labels %>%
  arrange(desc(n_proteins)) %>%
  pull(module_color) %>%
  paste0("ME", .)
mod_order <- intersect(mod_order, rownames(lmm_r))

# Gene counts per module
gene_counts <- mod_bio_labels %>%
  mutate(me_name = paste0("ME", module_color)) %>%
  dplyr::select(me_name, n_proteins) %>%
  deframe()

# --- Contrast display labels ---
contrast_labels <- c(
  Cancer_vs_Healthy = "CR vs H",
  Training_CR       = "Tr.(CR)"
)

# --- Build tile dataframe ---
all_contrasts <- crvh_contrasts

tile_df <- expand.grid(
  module   = mod_order,
  contrast = all_contrasts,
  stringsAsFactors = FALSE
) %>% as_tibble() %>%
  mutate(
    r      = mapply(function(m, c) {
      if (m %in% rownames(lmm_r) && c %in% colnames(lmm_r)) lmm_r[m, c] else NA_real_
    }, module, contrast),
    p_bh   = mapply(function(m, c) {
      if (m %in% rownames(lmm_p) && c %in% colnames(lmm_p)) lmm_p[m, c] else NA_real_
    }, module, contrast),
    p_raw  = mapply(function(m, c) {
      if (m %in% rownames(lmm_p_raw) && c %in% colnames(lmm_p_raw)) lmm_p_raw[m, c] else NA_real_
    }, module, contrast),
    model  = "CRvH",
    # Two-tier significance
    sig_tier = case_when(
      p_bh < 0.05  ~ "FDR",
      p_raw < 0.05 ~ "nominal",
      TRUE         ~ "ns"
    ),
    label  = paste0(sprintf("%.2f", r),
                    ifelse(p_bh < 0.001, "\n***",
                    ifelse(p_bh < 0.01,  "\n**",
                    ifelse(p_bh < 0.05,  "\n*", "")))),
    mod_label = mod_display_vec[module],
    ctr_label = contrast_labels[contrast]
  )

# Factor levels for ordering
tile_df$module   <- factor(tile_df$module, levels = rev(mod_order))
tile_df$contrast <- factor(tile_df$contrast, levels = all_contrasts)

# --- Gene count bar ---
count_df <- tibble(
  module = factor(mod_order, levels = rev(mod_order)),
  n = gene_counts[mod_order]
)

p_counts <- ggplot(count_df, aes(x = n, y = module)) +
  geom_col(fill = "grey70", width = 0.7) +
  geom_text(aes(label = n), hjust = -0.1, size = 2.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.3))) +
  labs(x = "Proteins", y = NULL) +
  FIG_THEME +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid   = element_blank(),
    panel.border = element_blank(),
    axis.line.x  = element_line(color = "black", linewidth = 0.3)
  )

# --- Main heatmap ---
p_heat <- ggplot(tile_df, aes(x = contrast, y = module)) +
  geom_tile(aes(fill = r), color = "white", linewidth = 0.5) +
  # FDR border: solid
  geom_tile(data = tile_df %>% filter(sig_tier == "FDR"),
            aes(fill = r), color = "black", linewidth = 0.8) +
  # Nominal border: dashed (approximate with thin dotted rect)
  geom_tile(data = tile_df %>% filter(sig_tier == "nominal"),
            aes(fill = r), color = "grey30", linewidth = 0.5, linetype = "dashed") +
  geom_text(aes(label = label), size = 2.3, lineheight = 0.85) +
  scale_fill_gradient2(
    low = "#4393C3", mid = "white", high = "#D6604D",
    midpoint = 0, limits = c(-1, 1),
    name = expression(r[equiv]),
    guide = guide_colorbar(barwidth = unit(3, "mm"), barheight = unit(40, "mm"))
  ) +
  scale_x_discrete(labels = contrast_labels, position = "top") +
  scale_y_discrete(labels = function(x) {
    lbl <- mod_display_vec[x]
    ifelse(is.na(lbl), x, lbl)
  }) +
  labs(y = NULL, x = NULL,
       title = "Module-LMM Contrast Associations",
       subtitle = "Solid border = FDR < 0.05; dashed = nominal p < 0.05") +
  coord_cartesian(clip = "off") +
  FIG_THEME +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 0, size = 8),
    axis.text.y  = element_text(size = 7.5),
    plot.margin  = margin(15, 5, 5, 5),
    legend.position = "right"
  )

# --- Composite ---
panel_A <- p_heat + p_counts + plot_layout(widths = c(6, 1))

W <- 240; H <- 200

ggsave(file.path(RPT_PDF, "panel_A_module_trait_MAIN.pdf"), panel_A,
       width = W, height = H, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT_PNG, "panel_A_module_trait_MAIN.png"), panel_A,
       width = W, height = H, units = "mm",
       dpi = 300, limitsize = FALSE)

# --- Save audit CSV ---
write_csv(tile_df %>% dplyr::select(module, contrast, model, r, p_bh, p_raw, sig_tier),
          file.path(DAT, "01_panel_A_heatmap_data.csv"))

message("  Panel A (module-trait heatmap) saved")
