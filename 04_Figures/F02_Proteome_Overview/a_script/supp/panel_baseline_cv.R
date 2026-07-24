# F02 CRvH — Panel F: Baseline CV% Scatter (Biopsy Method Comparison)
# CR_T1 (microneedle, n=14) vs H_T1 (Bergstrom, n=10)
# Plasma proteins highlighted; key blood markers labeled.
# Outputs: baseline_cv.png, audit CSVs

setwd(here::here())
source("04_Figures/shared/style.R")

pacman::p_load(dplyr, tidyr, readr, ggplot2, ggrepel)

PF_W <- 160
PF_H <- 150

RPT <- "04_Figures/F02_Proteome_Overview/b_reports/supp/panels"
DAT <- "04_Figures/F02_Proteome_Overview/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

# --- Load data ---
.dal <- readRDS("02_Normalization/c_data/DAList_normalized.rds")
norm_df <- tibble::as_tibble(cbind(
  .dal$annotation[, c("uniprot_id", "protein", "gene", "description")],
  as.data.frame(.dal$data)
))
meta <- read_csv("00_input/CvH_meta.csv", show_col_types = FALSE)
hpa <- read.delim("00_input/HPA_annotations.tsv",
  check.names = FALSE, stringsAsFactors = FALSE
)

ann_cols <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(norm_df), ann_cols)
meta <- meta |> filter(Col_ID %in% samp_names)

# --- Sample groups ---
cr_t1_ids <- meta$Col_ID[meta$Group %in% c("CR_CRE", "CR_PLA") & meta$Timepoint == "T1"]
h_t1_ids <- meta$Col_ID[meta$Group == "PPS"]

# --- CV on linear scale (Brenes 2024) ---
lin_mat <- 2^as.matrix(norm_df[, samp_names])

compute_cv <- function(mat, idx) {
  sub <- mat[, idx, drop = FALSE]
  apply(sub, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) {
      return(NA_real_)
    }
    sd(x) / mean(x) * 100
  })
}

cv_cr <- compute_cv(lin_mat, cr_t1_ids)
cv_h <- compute_cv(lin_mat, h_t1_ids)

# --- HPA annotation ---
hpa_sub <- hpa[, c("Gene", "Protein class")]
names(hpa_sub) <- c("gene", "protein_class")
hpa_sub <- hpa_sub[!duplicated(hpa_sub$gene), ]

scatter_df <- tibble(
  gene  = norm_df$gene,
  cv_cr = cv_cr,
  cv_h  = cv_h
) |>
  filter(!is.na(cv_cr), !is.na(cv_h)) |>
  left_join(hpa_sub, by = "gene") |>
  mutate(
    is_plasma  = grepl("Plasma proteins", protein_class, fixed = TRUE),
    is_plasma  = replace_na(is_plasma, FALSE),
    delta_cv   = cv_cr - cv_h,
    max_cv     = pmax(cv_cr, cv_h)
  )

# --- Blood markers ---
blood_markers <- c(
  "HBB", "HBA1", "HBA2", "ALB", "TF", "HP", "SERPINA1",
  "C3", "A2M", "FGA", "FGB", "FGG", "APOA1", "HPX"
)
blood_df <- scatter_df |> filter(gene %in% blood_markers)

# --- Correlations ---
n_all <- sum(!is.na(scatter_df$cv_cr) & !is.na(scatter_df$cv_h))
r_all <- cor(scatter_df$cv_cr, scatter_df$cv_h, use = "complete.obs")
ci_all <- fisher_z_ci(r_all, n_all)

plasma_df <- scatter_df |> filter(is_plasma)
nonplasma_df <- scatter_df |> filter(!is_plasma)

n_p <- nrow(plasma_df)
r_p <- cor(plasma_df$cv_cr, plasma_df$cv_h, use = "complete.obs")
ci_p <- fisher_z_ci(r_p, n_p)

n_np <- nrow(nonplasma_df)
r_np <- cor(nonplasma_df$cv_cr, nonplasma_df$cv_h, use = "complete.obs")
ci_np <- fisher_z_ci(r_np, n_np)

# Paired Wilcoxon: is CV systematically different between groups?
w_test <- wilcox.test(scatter_df$cv_cr, scatter_df$cv_h, paired = TRUE)

# Plasma vs non-plasma delta_cv comparison
w_plasma <- wilcox.test(plasma_df$delta_cv, nonplasma_df$delta_cv)

med_cr <- median(scatter_df$cv_cr, na.rm = TRUE)
med_h <- median(scatter_df$cv_h, na.rm = TRUE)

cor_label <- paste0(
  sprintf("All: r = %.2f [%.2f, %.2f] (n=%d)", r_all, ci_all["lo"], ci_all["hi"], n_all), "\n",
  sprintf("Plasma: r = %.2f [%.2f, %.2f] (n=%d)", r_p, ci_p["lo"], ci_p["hi"], n_p), "\n",
  sprintf("Non-plasma: r = %.2f [%.2f, %.2f] (n=%d)", r_np, ci_np["lo"], ci_np["hi"], n_np)
)

# --- Plot ---
axis_max <- quantile(pmax(scatter_df$cv_cr, scatter_df$cv_h), 0.995, na.rm = TRUE)

PLASMA_COL <- "#C62828"

pF <- ggplot(scatter_df, aes(x = cv_h, y = cv_cr)) +
  geom_abline(
    slope = 1, intercept = 0, linetype = "dashed",
    color = "grey50", linewidth = 0.4
  ) +
  geom_point(
    data = scatter_df |> filter(!is_plasma),
    color = "grey60", alpha = 0.3, size = 0.8
  ) +
  geom_point(
    data = scatter_df |> filter(is_plasma),
    color = PLASMA_COL, alpha = 0.45, size = 1.0
  ) +
  geom_label_repel(
    data = blood_df, aes(label = gene),
    fill = PLASMA_COL, color = "white", fontface = "bold",
    size = scale_text(BASE_GENE, PF_W),
    label.padding = unit(1.5, "pt"), label.size = 0.3,
    max.overlaps = 20, segment.size = 0.3, segment.color = "grey40",
    min.segment.length = 0, seed = 42, show.legend = FALSE
  ) +
  annotate("label",
    x = Inf, y = -Inf, label = cor_label,
    hjust = 1.05, vjust = -0.3,
    size = scale_text(BASE_STAT, PF_W) * 0.85,
    color = "grey20", fontface = "bold",
    fill = alpha("white", 0.9), linewidth = 0,
    label.padding = unit(3, "pt")
  ) +
  coord_equal(xlim = c(0, axis_max), ylim = c(0, axis_max)) +
  scale_color_identity() +
  labs(
    title = "Baseline CV% by Biopsy Method",
    subtitle = sprintf(
      "%s proteins | CR n=%d (microneedle), H n=%d (Bergstrom) | Med CV: CR=%.0f%%, H=%.0f%% | W p %s",
      format(nrow(scatter_df), big.mark = ","),
      length(cr_t1_ids), length(h_t1_ids),
      med_cr, med_h, fmt_p(w_test$p.value)
    ),
    x = "CV% (Healthy, Bergstrom needle)",
    y = "CV% (CR Baseline, microneedle)",
    tag = "F"
  ) +
  FIG_THEME +
  theme(
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    legend.position  = "none"
  )

# --- Save ---
ggsave(file.path(RPT, "baseline_cv.png"), pF,
  width = PF_W, height = PF_H, units = "mm", dpi = 300
)

write_csv(
  scatter_df |> dplyr::select(gene, cv_cr, cv_h, delta_cv, is_plasma, protein_class),
  file.path(DAT, "audit_baseline_cv.csv")
)

write_csv(
  tibble(
    subset = c("all", "plasma", "non_plasma"),
    n = c(n_all, n_p, n_np),
    r = round(c(r_all, r_p, r_np), 4),
    ci_lo = round(c(ci_all["lo"], ci_p["lo"], ci_np["lo"]), 4),
    ci_hi = round(c(ci_all["hi"], ci_p["hi"], ci_np["hi"]), 4),
    med_cv_cr = round(c(med_cr, median(plasma_df$cv_cr), median(nonplasma_df$cv_cr)), 2),
    med_cv_h = round(c(med_h, median(plasma_df$cv_h), median(nonplasma_df$cv_h)), 2),
    wilcox_p = c(w_test$p.value, NA, NA),
    delta_cv_wilcox_p = c(NA, w_plasma$p.value, NA)
  ),
  file.path(DAT, "audit_baseline_cv_cor.csv")
)

message("F02/CRvH Panel F done")
