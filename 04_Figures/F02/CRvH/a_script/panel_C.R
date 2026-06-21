# Figure 2 (CRvH) — Panel C: PCA Biplot + PERMANOVA
# CRvH sample space: all 5 groups (CRE_T1, CRE_T2, PLA_T1, PLA_T2, H_T1).
# PERMANOVA with Group_Time factor. 80% confidence ellipses. Bootstrap PC variance CIs.
# Outputs: pC (ggplot object), panel_C_pca_MAIN.pdf/.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F02/a_script/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(vegan)
})

PC_W <- 145; PC_H <- 100

RPT_DIR <- "04_Figures/F02/CRvH/b_reports"
DAT_DIR <- "04_Figures/F02/CRvH/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Load data & metadata ──
# Imputed matrix from the proteoDA DAList (annotation + per-sample columns).
.dal <- readRDS("02_Normalization/imputation/c_data/DAList_imputed_imp4p.rds")
imp_df <- tibble::as_tibble(cbind(
  .dal$annotation[, c("uniprot_id", "protein", "gene", "description")],
  as.data.frame(.dal$data)
))
meta_full <- read_csv("00_input/CvH_meta.csv", show_col_types = FALSE)

ann_cols   <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(imp_df), ann_cols)

# CRvH sample space
meta <- meta_full |>
  filter(Group %in% c("CR_CRE", "CR_PLA", "PPS"), Col_ID %in% samp_names)

samp_crvh <- meta$Col_ID

meta$Group_Time  <- factor(meta$Group_Time,
                           levels = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1"))
meta$Timepoint   <- factor(meta$Timepoint, levels = c("T1", "T2"))
meta$Subject_ID  <- factor(meta$Subject_ID)

imp_mat <- as.matrix(imp_df[, samp_crvh])
rownames(imp_mat) <- imp_df$gene

pdf_device <- get_pdf_device()

# ── PCA ──
pca <- prcomp(t(imp_mat), center = TRUE, scale. = TRUE)
var_pct <- round(100 * summary(pca)$importance[2, 1:2], 1)

# Bootstrap PC variance CIs
set.seed(42)
n_boot <- 1000
n_prot <- nrow(imp_mat)
boot_var <- matrix(NA_real_, nrow = n_boot, ncol = 2)
for (b in seq_len(n_boot)) {
  idx <- sample(n_prot, replace = TRUE)
  bp  <- prcomp(t(imp_mat[idx, ]), center = TRUE, scale. = TRUE)
  boot_var[b, ] <- 100 * summary(bp)$importance[2, 1:2]
}
var_ci <- data.frame(
  PC      = c("PC1", "PC2"),
  var_pct = var_pct,
  ci_lo   = apply(boot_var, 2, quantile, 0.025),
  ci_hi   = apply(boot_var, 2, quantile, 0.975)
)

pca_df <- as.data.frame(pca$x[, 1:2]) |>
  mutate(sample_id = rownames(pca$x)) |>
  left_join(meta, by = c("sample_id" = "Col_ID"))

# ── PERMANOVA — Group_Time factor ──
dist_mat <- dist(scale(t(imp_mat)))

# Reorder meta to match dist_mat row ordering (samp_crvh)
meta_ordered <- meta[match(samp_crvh, meta$Col_ID), ]

set.seed(42)
perm_res <- adonis2(dist_mat ~ Group_Time, data = meta_ordered,
                    permutations = 999, by = "terms")

gt_r2 <- perm_res["Group_Time", "R2"]
gt_pv <- perm_res["Group_Time", "Pr(>F)"]
perm_label <- sprintf("PERMANOVA\nGroup_Time  R\u00b2 = %.3f,  %s",
                       gt_r2, fmt_p(gt_pv))

# Betadisper homogeneity check
bd_gt  <- betadisper(dist_mat, meta_ordered$Group_Time)
bd_gt_p <- permutest(bd_gt, pairwise = FALSE, permutations = 999)$tab$`Pr(>F)`[1]
if (bd_gt_p < 0.05)
  warning("Heterogeneous dispersions detected -- interpret PERMANOVA with caution")

# ── PCA labels ──
pca_labels <- c(CRE_T1 = "CRE T1", CRE_T2 = "CRE T2",
                PLA_T1 = "PLA T1", PLA_T2 = "PLA T2",
                H_T1   = "Healthy T1")

pC <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group_Time, shape = Group_Time)) +
  stat_ellipse(aes(fill = Group_Time), geom = "polygon",
               alpha = 0.10, level = 0.80, show.legend = FALSE) +
  stat_ellipse(aes(group = Group_Time), level = 0.80, linewidth = 0.4,
               linetype = "dashed", show.legend = FALSE) +
  geom_point(size = 2.5, alpha = 0.85) +
  annotate("text", x = Inf, y = Inf, label = perm_label,
           hjust = 1.05, vjust = 1.15,
           size = scale_text(BASE_STAT - 1.2, PC_W), color = "grey30",
           fontface = "bold") +
  scale_color_manual(values = PCA_COLORS, labels = pca_labels,
                     guide = guide_legend(override.aes = list(size = 2.5))) +
  scale_fill_manual(values = PCA_COLORS, guide = "none") +
  scale_shape_manual(values = PCA_SHAPES, labels = pca_labels) +
  labs(title = "Principal Component Analysis (PCA)",
       subtitle = sprintf("%s proteins (imputed), %d samples | PERMANOVA R\u00b2 = %.2f, %s",
                          format(nrow(imp_df), big.mark = ","), nrow(meta),
                          gt_r2, fmt_p(gt_pv)),
       x = sprintf("PC1 (%.1f%% [%.1f, %.1f])", var_pct[1], var_ci$ci_lo[1], var_ci$ci_hi[1]),
       y = sprintf("PC2 (%.1f%% [%.1f, %.1f])", var_pct[2], var_ci$ci_lo[2], var_ci$ci_hi[2]),
       tag = "C") +
  FIG_THEME + theme(legend.position = c(0.88, 0.12),
                    legend.background = element_rect(fill = alpha("white", 0.8),
                                                     color = NA),
                    legend.title = element_blank(),
                    legend.text = element_text(size = FIG_LEGEND_TEXT),
                    legend.key.size = unit(3, "mm"),
                    legend.spacing.y = unit(0.5, "mm"))

# ── Save ──
write.csv(var_ci, file.path(DAT_DIR, "audit_panel_C_pca_variance_ci.csv"),
          row.names = FALSE)
betadisper_results <- data.frame(
  factor      = "Group_Time",
  p_value     = bd_gt_p,
  significant = bd_gt_p < 0.05
)
write.csv(betadisper_results, file.path(DAT_DIR, "audit_panel_C_betadisper.csv"),
          row.names = FALSE)

ggsave(file.path(RPT_DIR, "panel_C_pca_MAIN.pdf"), pC,
       width = PC_W, height = PC_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_C_pca_MAIN.png"), pC,
       width = PC_W, height = PC_H, units = "mm", dpi = 300)
