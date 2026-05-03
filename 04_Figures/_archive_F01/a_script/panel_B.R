# Figure 1 — Panel B: CV Scatter Triptych
# B1/B2: Per-protein CV% T1 vs T2 (Creatine / Placebo).
# B3: DeltaCV CRE vs DeltaCV PLA.
# Outputs: pB (combined ggplot), panel_B_cv_scatter.pdf/.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(cowplot)
})

PB_SUB <- 80; PB_W <- 380; PB_H <- 120

RPT_DIR <- "04_Figures/F01/b_reports"
DAT_DIR <- "04_Figures/F01/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

norm_df <- read_csv("01_normalization/c_data/02_normalized.csv",
                    show_col_types = FALSE)

ann_cols   <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(norm_df), ann_cols)

# --- Load metadata from CSV (CvH pattern) ---
meta <- read_csv("00_input/CvH_meta.csv", show_col_types = FALSE) |>
  filter(Col_ID %in% samp_names,
         Supplement %in% c("CRE", "PLA"))

# Build CR pooled metadata (CRE+PLA at each timepoint)
meta_cr <- meta |>
  mutate(supplement = factor("CR", levels = c("CRE", "PLA", "CR")),
         timepoint  = factor(Timepoint, levels = c("T1", "T2")))

meta$supplement <- factor(meta$Supplement, levels = c("CRE", "PLA"))
meta$timepoint  <- factor(meta$Timepoint,  levels = c("T1", "T2"))

pdf_device <- get_pdf_device()

# --- CV on linear scale per Brenes 2024 ---
lin_mat <- 2^as.matrix(norm_df[, samp_names])

compute_cv <- function(mat, idx) {
  sub <- mat[, idx, drop = FALSE]
  apply(sub, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) return(NA_real_)
    sd(x) / mean(x) * 100
  })
}

# --- CV per supplement x timepoint (CRE, PLA, CR pooled) ---
scatter_list <- lapply(c("CRE", "PLA", "CR"), function(supp) {
  if (supp == "CR") {
    t1_idx <- meta_cr$Col_ID[meta_cr$timepoint == "T1"]
    t2_idx <- meta_cr$Col_ID[meta_cr$timepoint == "T2"]
  } else {
    t1_idx <- meta$Col_ID[meta$supplement == supp & meta$timepoint == "T1"]
    t2_idx <- meta$Col_ID[meta$supplement == supp & meta$timepoint == "T2"]
  }
  cv_t1  <- compute_cv(lin_mat, t1_idx)
  cv_t2  <- compute_cv(lin_mat, t2_idx)
  tibble(uniprot_id = norm_df$uniprot_id, gene = norm_df$gene,
         cv_t1 = cv_t1, cv_t2 = cv_t2, supplement = supp)
})

SUPP_LABELS_B <- c(CRE = "Creatine", PLA = "Placebo", CR = "Cancer Recovery")
scatter_df <- bind_rows(scatter_list) |>
  filter(!is.na(cv_t1), !is.na(cv_t2)) |>
  mutate(
    delta_cv = cv_t2 - cv_t1,
    max_cv   = pmax(cv_t1, cv_t2),
    supp_code = supplement,
    supplement = factor(supplement, levels = c("CR", "CRE", "PLA"),
                        labels = SUPP_LABELS_B[c("CR", "CRE", "PLA")])
  )

# --- Color caps at 98th percentile ---
max_cv_cap <- quantile(scatter_df$max_cv, 0.98, na.rm = TRUE)
scatter_df$max_cv_capped <- pmin(scatter_df$max_cv, max_cv_cap)

cv_cap <- quantile(abs(scatter_df$delta_cv), 0.98, na.rm = TRUE)
scatter_df$delta_cv_capped <- pmin(pmax(scatter_df$delta_cv, -cv_cap), cv_cap)

# --- Top 15 labels per supplement by max CV ---
top_cv_labels <- bind_rows(
  scatter_df |> filter(supplement == "Cancer Recovery") |>
    slice_max(max_cv, n = 15, with_ties = FALSE),
  scatter_df |> filter(supplement == "Creatine") |>
    slice_max(max_cv, n = 15, with_ties = FALSE),
  scatter_df |> filter(supplement == "Placebo") |>
    slice_max(max_cv, n = 15, with_ties = FALSE)
)

# --- Pearson correlations ---
r_cr  <- cor(scatter_df$cv_t1[scatter_df$supplement == "Cancer Recovery"],
             scatter_df$cv_t2[scatter_df$supplement == "Cancer Recovery"],
             use = "complete.obs")
r_cre <- cor(scatter_df$cv_t1[scatter_df$supplement == "Creatine"],
             scatter_df$cv_t2[scatter_df$supplement == "Creatine"],
             use = "complete.obs")
r_pla <- cor(scatter_df$cv_t1[scatter_df$supplement == "Placebo"],
             scatter_df$cv_t2[scatter_df$supplement == "Placebo"],
             use = "complete.obs")

r_annotations <- tibble(
  supplement = factor(c("Cancer Recovery", "Creatine", "Placebo"),
                      levels = c("Cancer Recovery", "Creatine", "Placebo")),
  label      = sprintf("r = %.2f", c(r_cr, r_cre, r_pla))
)

# --- DeltaCV wide: CRE vs PLA ---
delta_wide <- scatter_df |>
  select(uniprot_id, gene, delta_cv, supp_code) |>
  pivot_wider(id_cols = c(uniprot_id, gene),
              names_from = supp_code, values_from = delta_cv,
              names_prefix = "dcv_") |>
  filter(!is.na(dcv_CRE), !is.na(dcv_PLA)) |>
  mutate(
    dist_origin = sqrt(dcv_CRE^2 + dcv_PLA^2),
    mean_dcv    = (dcv_CRE + dcv_PLA) / 2,
    mean_dcv_capped = pmin(pmax(mean_dcv, -cv_cap), cv_cap)
  )

top_delta <- delta_wide |>
  slice_max(dist_origin, n = 15, with_ties = FALSE)

delta_abs_max <- max(abs(c(delta_wide$dcv_CRE, delta_wide$dcv_PLA)), na.rm = TRUE) * 1.15

r_delta <- cor(delta_wide$dcv_CRE, delta_wide$dcv_PLA,
               use = "complete.obs")

# --- Shared theme ---
theme_B <- FIG_THEME +
  theme(
    legend.position  = "none",
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(hjust = 0.5, size = FIG_STRIP_SIZE,
                                    face = "bold")
  )

axis_max_cv <- 300

# --- B1/B2: CV scatter by supplement ---
pB12 <- ggplot(scatter_df, aes(x = cv_t1, y = cv_t2)) +
  facet_wrap(~supplement, nrow = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey50", linewidth = 0.4) +
  geom_point(aes(color = max_cv_capped), alpha = 0.35, size = 0.9) +
  geom_label_repel(data = top_cv_labels,
                   aes(label = gene, fill = max_cv_capped),
                   color = "white", fontface = "bold",
                   size = scale_text(BASE_GENE, PB_SUB),
                   label.padding = unit(1, "pt"),
                   label.size = 0.3, max.overlaps = 20,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, seed = 42, show.legend = FALSE) +
  geom_text(data = r_annotations, aes(label = label),
            x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1,
            size = scale_text(BASE_STAT + 0.5, PB_SUB),
            color = "grey30", fontface = "bold", inherit.aes = FALSE) +
  scale_color_viridis_c(option = "inferno", direction = -1,
                        begin = 0.1, end = 0.85,
                        name = "CV%",
                        guide = guide_colorbar(barwidth = unit(2, "mm"),
                                               barheight = unit(12, "mm"),
                                               title.position = "top",
                                               title.hjust = 0.5)) +
  scale_fill_viridis_c(option = "inferno", direction = -1,
                       begin = 0.1, end = 0.85,
                       name = "CV%",
                       guide = "none") +
  coord_equal(xlim = c(0, axis_max_cv), ylim = c(0, axis_max_cv)) +
  labs(title = "Per-Protein Variability (CV%)",
       subtitle = "T1 vs T2; labels = top 15 by max CV%",
       x = expression(bold(CV * "%"[T1])),
       y = expression(bold(CV * "%"[T2])),
       tag = "B") +
  theme_B +
  theme(plot.title    = element_text(hjust = 0, size = FIG_TITLE_SIZE, face = "bold"),
        plot.subtitle = element_text(hjust = 0, size = FIG_SUBTITLE_SIZE,
                                     face = "bold.italic", color = "grey30"),
        strip.text    = element_text(face = "bold", size = FIG_STRIP_SIZE),
        legend.position  = c(0.97, 0.02),
        legend.justification = c(1, 0),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.title = element_text(face = "bold", size = FIG_LEGEND_TITLE),
        legend.text  = element_text(size = FIG_LEGEND_TEXT),
        legend.key.size = unit(3, "mm"),
        plot.margin = margin(5.5, 0, 5.5, 5.5))

# --- B3: DeltaCV CRE vs PLA ---
pB3 <- ggplot(delta_wide, aes(x = dcv_CRE, y = dcv_PLA)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(aes(color = mean_dcv_capped), alpha = 0.4, size = 0.9) +
  geom_label_repel(data = top_delta, aes(label = gene, fill = mean_dcv_capped),
                   color = "white", fontface = "bold",
                   size = scale_text(BASE_GENE, PB_SUB),
                   label.padding = unit(1, "pt"), label.size = 0.3,
                   max.overlaps = 25,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, seed = 44, show.legend = FALSE) +
  annotate("text", x = -Inf, y = Inf, label = sprintf("r = %.2f", r_delta),
           hjust = -0.05, vjust = 1.5,
           size = scale_text(BASE_STAT + 0.5, PB_SUB),
           color = "grey30", fontface = "bold") +
  scale_color_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B",
                        midpoint = 0, limits = c(-cv_cap, cv_cap),
                        name = expression(bold(Delta * "CV%")),
                        guide = guide_colorbar(barwidth = unit(2, "mm"),
                                               barheight = unit(12, "mm"),
                                               title.position = "top",
                                               title.hjust = 0.5)) +
  scale_fill_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B",
                       midpoint = 0, limits = c(-cv_cap, cv_cap),
                       guide = "none") +
  coord_equal(xlim = c(-delta_abs_max, delta_abs_max),
              ylim = c(-delta_abs_max, delta_abs_max)) +
  labs(x = expression(bold(Delta * "CV%"[Creatine])),
       y = expression(bold(Delta * "CV%"[Placebo])),
       title = "Training Response") +
  theme_B +
  theme(legend.position  = c(0.97, 0.02),
        legend.justification = c(1, 0),
        legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
        legend.title = element_text(face = "bold", size = FIG_LEGEND_TITLE),
        legend.text  = element_text(size = FIG_LEGEND_TEXT),
        legend.key.size = unit(3, "mm"),
        axis.title.y = element_text(margin = margin(r = 0, l = 0)),
        plot.margin = margin(5.5, 5.5, 5.5, 0))

# --- Audit CSVs ---
write.csv(scatter_df |> select(gene, cv_t1, cv_t2, delta_cv,
                               supplement = supplement),
          file.path(DAT_DIR, "panel_B_cv_scatter.csv"), row.names = FALSE)
write.csv(delta_wide |> select(gene, dcv_CRE, dcv_PLA,
                               mean_dcv, dist_origin),
          file.path(DAT_DIR, "panel_B3_delta_cv.csv"), row.names = FALSE)

# --- Combine: 2:1 ratio (pB12 has 2 facets, pB3 has 1) ---
pB <- cowplot::plot_grid(
  pB12, pB3,
  nrow = 1, rel_widths = c(3, 1), align = "h", axis = "tb"
)

ggsave(file.path(RPT_DIR, "panel_B_cv_scatter.pdf"), pB,
       width = PB_W, height = PB_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_B_cv_scatter.png"), pB,
       width = PB_W, height = PB_H, units = "mm", dpi = 300)

cat("Panel B done.\n")
