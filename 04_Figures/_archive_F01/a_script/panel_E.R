# Figure 1 — Panel E: PCA Biplot + PERMANOVA
# Outputs: pE (ggplot object), panel_E_pca.pdf/.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(vegan)
})

PE_W <- 145; PE_H <- 100

RPT_DIR <- "04_Figures/F01/b_reports"
DAT_DIR <- "04_Figures/F01/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Load imputed data ---
imp_df   <- read_csv("02_Imputation/c_data/01_imputed.csv", show_col_types = FALSE)
ann_cols <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(imp_df), ann_cols)

imp_mat <- as.matrix(imp_df[, samp_names])
rownames(imp_mat) <- imp_df$gene

# --- Load metadata, filter to samples in imputed data ---
meta <- read_csv("00_input/CvH_meta.csv", show_col_types = FALSE) |>
  filter(Col_ID %in% samp_names) |>
  mutate(Group_Time = factor(Group_Time,
    levels = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1")))

# Align matrix columns to metadata order
imp_mat <- imp_mat[, meta$Col_ID]

pdf_device <- get_pdf_device()

# --- PCA ---
pca     <- prcomp(t(imp_mat), center = TRUE, scale. = TRUE)
var_pct <- round(100 * summary(pca)$importance[2, 1:2], 1)

# --- Bootstrap PC variance CIs ---
set.seed(42)
n_boot <- 1000; n_prot <- nrow(imp_mat)
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
  ci_hi   = apply(boot_var, 2, quantile, 0.975))

# --- Join PCA scores with metadata ---
pca_df <- as.data.frame(pca$x[, 1:2]) |>
  mutate(Col_ID = rownames(pca$x)) |>
  left_join(meta, by = "Col_ID")

# --- PERMANOVA: Group_Time as single factor (unbalanced design) ---
dist_mat <- dist(scale(t(imp_mat)))
set.seed(42)
perm_res <- adonis2(dist_mat ~ Group_Time, data = meta,
                    permutations = 999, by = "terms")

perm_r2 <- perm_res["Group_Time", "R2"]
perm_pv <- perm_res["Group_Time", "Pr(>F)"]

# --- CR-only PERMANOVA: Supplement x Timepoint with paired design ---
cr_idx <- meta$Group %in% c("CR_CRE", "CR_PLA")
cr_meta <- meta[cr_idx, ] |>
  mutate(Supplement = factor(Supplement, levels = c("PLA", "CRE")),
         Timepoint  = factor(Timepoint,  levels = c("T1", "T2")))
cr_dist <- dist(scale(t(imp_mat[, cr_meta$Col_ID])))
set.seed(42)
perm_cr <- adonis2(cr_dist ~ Supplement * Timepoint, data = cr_meta,
                   permutations = how(nperm = 999, blocks = cr_meta$Subject_ID),
                   by = "terms")

cr_terms <- c("Supplement", "Timepoint", "Supplement:Timepoint")
cr_r2 <- perm_cr[cr_terms, "R2"]
cr_pv <- perm_cr[cr_terms, "Pr(>F)"]

# --- Format PERMANOVA annotation ---
perm_label <- sprintf(
  "PERMANOVA\nGroup  R\u00b2 = %.3f,  p %s\nCR-only: Suppl R\u00b2 = %.3f  p %s\n         Time R\u00b2 = %.3f  p %s\n         S\u00d7T R\u00b2 = %.3f  p %s",
  perm_r2, fmt_p(perm_pv),
  cr_r2[1], fmt_p(cr_pv[1]),
  cr_r2[2], fmt_p(cr_pv[2]),
  cr_r2[3], fmt_p(cr_pv[3]))

# --- Betadisper ---
bd_grp   <- betadisper(dist_mat, meta$Group_Time)
bd_grp_p <- permutest(bd_grp, pairwise = FALSE, permutations = 999)$tab$`Pr(>F)`[1]
if (bd_grp_p < 0.05)
  warning("Heterogeneous dispersions detected — interpret PERMANOVA with caution")

# --- Legend labels ---
pca_labels <- c(CRE_T1 = "CRE T1", CRE_T2 = "CRE T2",
                PLA_T1 = "PLA T1", PLA_T2 = "PLA T2",
                H_T1   = "Healthy")

# --- Plot ---
pE <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group_Time, shape = Group_Time)) +
  stat_ellipse(aes(fill = Group_Time), geom = "polygon",
               alpha = 0.08, level = 0.80, show.legend = FALSE) +
  stat_ellipse(aes(group = Group_Time), level = 0.80, linewidth = 0.4,
               linetype = "dashed", show.legend = FALSE) +
  geom_point(size = 2.5, alpha = 0.85) +
  annotate("text", x = Inf, y = Inf, label = perm_label,
           hjust = 1.05, vjust = 1.10,
           size = scale_text(BASE_STAT - 1.4, PE_W), color = "grey30",
           fontface = "bold") +
  scale_color_manual(values = PCA_COLORS, labels = pca_labels,
                     guide = guide_legend(override.aes = list(size = 2.5))) +
  scale_fill_manual(values = PCA_COLORS, guide = "none") +
  scale_shape_manual(values = PCA_SHAPES, labels = pca_labels) +
  labs(title = "Principal Component Analysis (PCA)",
       subtitle = sprintf("%s proteins | missForest-imputed | n = %d",
                          format(nrow(imp_df), big.mark = ","), nrow(meta)),
       x = sprintf("PC1 (%.1f%% [%.1f, %.1f])", var_pct[1],
                    var_ci$ci_lo[1], var_ci$ci_hi[1]),
       y = sprintf("PC2 (%.1f%% [%.1f, %.1f])", var_pct[2],
                    var_ci$ci_lo[2], var_ci$ci_hi[2]),
       tag = "E") +
  FIG_THEME + theme(legend.position = c(0.12, 0.15),
                    legend.background = element_rect(fill = alpha("white", 0.8),
                                                     color = NA),
                    legend.title = element_blank(),
                    legend.text  = element_text(size = FIG_LEGEND_TEXT),
                    legend.key.size  = unit(3, "mm"),
                    legend.spacing.y = unit(0.5, "mm"))

# --- Exports ---
write.csv(var_ci, file.path(DAT_DIR, "panel_E_pca_variance_ci.csv"),
          row.names = FALSE)

betadisper_results <- data.frame(
  factor      = "Group_Time",
  p_value     = bd_grp_p,
  significant = bd_grp_p < 0.05)
write.csv(betadisper_results, file.path(DAT_DIR, "panel_E_betadisper.csv"),
          row.names = FALSE)

ggsave(file.path(RPT_DIR, "panel_E_pca.pdf"), pE,
       width = PE_W, height = PE_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_E_pca.png"), pE,
       width = PE_W, height = PE_H, units = "mm", dpi = 300)

cat(sprintf("Panel E done: %s\n", file.path(RPT_DIR, "panel_E_pca.png")))
