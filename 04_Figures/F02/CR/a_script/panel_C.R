# Figure 2 CR — Panel C: PCA Biplot + PERMANOVA
# CR-only samples (no Healthy). PERMANOVA: Supplement * Timepoint, blocks=Subject_ID.
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

RPT_DIR <- "04_Figures/F02/CR/b_reports"
DAT_DIR <- "04_Figures/F02/CR/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

imp_df <- read_csv("02_Imputation/c_data/01_imputed.csv",
                   show_col_types = FALSE)
meta_full <- read_csv("00_input/CvH_meta.csv", show_col_types = FALSE)

# Filter to CR-only
meta_cr <- meta_full |>
  filter(Group %in% c("CR_CRE", "CR_PLA")) |>
  mutate(
    Supplement = factor(Supplement, levels = c("CRE", "PLA")),
    Timepoint  = factor(Timepoint, levels = c("T1", "T2")),
    Group_Time = factor(Group_Time, levels = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2")),
    Subject_ID = factor(Subject_ID)
  )

ann_cols   <- c("uniprot_id", "protein", "gene", "description")
cr_samples <- intersect(meta_cr$Col_ID, setdiff(names(imp_df), ann_cols))
meta_cr    <- meta_cr |> filter(Col_ID %in% cr_samples)

imp_mat <- as.matrix(imp_df[, cr_samples])
rownames(imp_mat) <- imp_df$gene

pdf_device <- get_pdf_device()

# Drop zero-variance proteins before scaling (blood filter may leave constant columns in CR subset)
zero_var <- apply(imp_mat, 1, var, na.rm = TRUE) == 0
if (any(zero_var)) message(sprintf("  Dropped %d zero-variance proteins for PCA", sum(zero_var)))
pca_mat <- imp_mat[!zero_var, , drop = FALSE]
pca <- prcomp(t(pca_mat), center = TRUE, scale. = TRUE)
var_pct <- round(100 * summary(pca)$importance[2, 1:2], 1)

set.seed(42)
n_boot <- 1000
n_prot <- nrow(pca_mat)
boot_var <- matrix(NA_real_, nrow = n_boot, ncol = 2)
for (b in seq_len(n_boot)) {
  idx <- sample(n_prot, replace = TRUE)
  bp  <- prcomp(t(pca_mat[idx, ]), center = TRUE, scale. = TRUE)
  boot_var[b, ] <- 100 * summary(bp)$importance[2, 1:2]
}
var_ci <- data.frame(
  PC      = c("PC1", "PC2"),
  var_pct = var_pct,
  ci_lo   = apply(boot_var, 2, quantile, 0.025),
  ci_hi   = apply(boot_var, 2, quantile, 0.975)
)

pca_df <- as.data.frame(pca$x[, 1:2]) |>
  mutate(Col_ID = rownames(pca$x)) |>
  left_join(meta_cr, by = "Col_ID")

# PERMANOVA — paired permutations within subjects
dist_mat <- dist(scale(t(imp_mat)))
set.seed(42)
perm_res <- adonis2(dist_mat ~ Supplement * Timepoint, data = meta_cr,
                    permutations = how(nperm = 999, blocks = meta_cr$Subject_ID),
                    by = "terms")

perm_terms <- c("Supplement", "Timepoint", "Supplement:Timepoint")
perm_r2 <- perm_res[perm_terms, "R2"]
perm_pv <- perm_res[perm_terms, "Pr(>F)"]
perm_label <- sprintf(
  "PERMANOVA\nSupp  R\u00b2 = %.3f,  %s\nTime  R\u00b2 = %.3f,  %s\nSupp\u00d7Time  R\u00b2 = %.3f,  %s",
  perm_r2[1], fmt_p(perm_pv[1]),
  perm_r2[2], fmt_p(perm_pv[2]),
  perm_r2[3], fmt_p(perm_pv[3]))

bd_supp <- betadisper(dist_mat, meta_cr$Supplement)
bd_time <- betadisper(dist_mat, meta_cr$Timepoint)
bd_grp  <- betadisper(dist_mat, meta_cr$Group_Time)
bd_supp_p <- permutest(bd_supp, pairwise = FALSE, permutations = 999)$tab$`Pr(>F)`[1]
bd_time_p <- permutest(bd_time, pairwise = FALSE, permutations = 999)$tab$`Pr(>F)`[1]
bd_grp_p  <- permutest(bd_grp,  pairwise = FALSE, permutations = 999)$tab$`Pr(>F)`[1]
if (bd_supp_p < 0.05 || bd_time_p < 0.05)
  warning("Heterogeneous dispersions detected -- interpret PERMANOVA with caution")

# CR-only PCA palette (4 groups, no H_T1)
cr_pca_colors <- PCA_COLORS[c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2")]
cr_pca_shapes <- PCA_SHAPES[c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2")]

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
  scale_color_manual(values = cr_pca_colors,
                     labels = c("CRE T1", "CRE T2", "PLA T1", "PLA T2"),
                     guide = guide_legend(override.aes = list(size = 2.5))) +
  scale_fill_manual(values = cr_pca_colors, guide = "none") +
  scale_shape_manual(values = cr_pca_shapes,
                     labels = c("CRE T1", "CRE T2", "PLA T1", "PLA T2")) +
  labs(title = "Principal Component Analysis (PCA)",
       subtitle = sprintf("%s proteins (imputed), %d samples | PERMANOVA R\u00b2 = %.2f, %s",
                          format(nrow(imp_df), big.mark = ","), nrow(meta_cr),
                          sum(perm_r2), fmt_p(min(perm_pv))),
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

write.csv(var_ci, file.path(DAT_DIR, "audit_panel_C_pca_variance_ci.csv"),
          row.names = FALSE)
betadisper_results <- data.frame(
  factor      = c("Supplement", "Timepoint", "Group_Time"),
  p_value     = c(bd_supp_p, bd_time_p, bd_grp_p),
  significant = c(bd_supp_p < 0.05, bd_time_p < 0.05, bd_grp_p < 0.05)
)
write.csv(betadisper_results, file.path(DAT_DIR, "audit_panel_C_betadisper.csv"),
          row.names = FALSE)

ggsave(file.path(RPT_DIR, "panel_C_pca_MAIN.pdf"), pC,
       width = PC_W, height = PC_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_C_pca_MAIN.png"), pC,
       width = PC_W, height = PC_H, units = "mm", dpi = 300)
