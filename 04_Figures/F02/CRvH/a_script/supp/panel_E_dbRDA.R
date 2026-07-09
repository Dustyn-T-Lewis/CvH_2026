# Figure 2 (CRvH) — Supp Panel E (db-RDA): Constrained Ordination
# Distance-based RDA on CRvH sample space.
# Full model: Group_Time. Within-subject model: Timepoint with Condition(Subject_ID).
# Outputs: pE_dbrda (ggplot object), panel_E_dbRDA_SUPP.pdf/.png

setwd(here::here())
source("04_Figures/F02/a_script/style.R")

pacman::p_load(dplyr, tidyr, stringr, readr, ggplot2, vegan)

PE_RDA_W <- 145; PE_RDA_H <- 100

RPT_DIR <- "04_Figures/F02/CRvH/b_reports/supp"
DAT_DIR <- "04_Figures/F02/CRvH/c_data/supp"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Load data & metadata ──
.dal <- readRDS("02_Normalization/imputation/c_data/DAList_imputed_imp4p.rds")
imp_df <- tibble::as_tibble(cbind(
  .dal$annotation[, c("uniprot_id", "protein", "gene", "description")],
  as.data.frame(.dal$data)))
meta_full <- read_csv("00_input/CvH_meta.csv", show_col_types = FALSE)

ann_cols   <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(imp_df), ann_cols)

# CRvH sample space: CR_CRE, CR_PLA, PPS
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

# Reorder meta to match samp_crvh
meta_ordered <- meta[match(samp_crvh, meta$Col_ID), ]

# ── db-RDA: constrained ordination ──
dist_mat <- vegdist(scale(t(imp_mat)), method = "euclidean")

# Full model: Group_Time
rda_full <- dbrda(dist_mat ~ Group_Time, data = meta_ordered)

# Within-subject model: Timepoint + Condition(Subject_ID)
# Only for paired CR subjects (PPS has no T2)
meta_paired <- meta_ordered |> filter(Group %in% c("CR_CRE", "CR_PLA"))
subj_counts <- meta_paired |> count(Subject_ID)
paired_subj <- subj_counts$Subject_ID[subj_counts$n == 2]
meta_within <- meta_paired |> filter(Subject_ID %in% paired_subj)

samp_within <- meta_within$Col_ID
imp_mat_within <- imp_mat[, samp_within]
dist_within <- vegdist(scale(t(imp_mat_within)), method = "euclidean")

meta_within_ordered <- meta_within[match(samp_within, meta_within$Col_ID), ]
meta_within_ordered$Subject_ID <- droplevels(meta_within_ordered$Subject_ID)

rda_within <- dbrda(dist_within ~ Timepoint + Condition(Subject_ID),
                    data = meta_within_ordered)

# Use full model for visualization (shows all groups)
site_scores  <- as.data.frame(scores(rda_full, display = "sites", choices = 1:2))
site_scores$sample_id <- rownames(site_scores)
site_scores <- left_join(site_scores, meta_ordered, by = c("sample_id" = "Col_ID"))

# Variance from full model
eigenvals_full <- eigenvals(rda_full)
constrained_eig <- eigenvals_full[grepl("^dbRDA", names(eigenvals_full))]
total_inertia   <- rda_full$tot.chi
constrained_var <- rda_full$CCA$tot.chi / total_inertia * 100

ax_var <- 100 * constrained_eig[1:2] / total_inertia

# ANOVA: full model and within-subject model
set.seed(42)
anova_full   <- anova(rda_full, by = "terms", permutations = 999)
anova_within <- anova(rda_within, by = "terms", permutations = 999)

# Column name depends on vegan version
var_col <- intersect(c("Variance", "SumOfSqs"), colnames(anova_full))[1]

get_anova_val <- function(av, term, col) {
  rn <- rownames(av)
  idx <- grep(paste0("^", term, "$"), rn)
  if (length(idx) == 0) return(NA_real_)
  av[idx, col]
}

total_ss_full   <- sum(anova_full[[var_col]], na.rm = TRUE)
total_ss_within <- sum(anova_within[[var_col]], na.rm = TRUE)

gt_r2  <- get_anova_val(anova_full, "Group_Time", var_col) / total_ss_full
gt_pv  <- get_anova_val(anova_full, "Group_Time", "Pr(>F)")
tp_r2  <- get_anova_val(anova_within, "Timepoint", var_col) / total_ss_within
tp_pv  <- get_anova_val(anova_within, "Timepoint", "Pr(>F)")

adj_r2 <- RsquareAdj(rda_full)$adj.r.squared

# Safe p-value formatter
safe_fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  fmt_p(p)
}

anova_label <- sprintf(
  "db-RDA (adj. R\u00b2 = %.3f)\nGroup_Time  R\u00b2 = %.3f,  %s\nTimepoint  R\u00b2 = %.3f,  %s  (cond.)",
  adj_r2,
  gt_r2, safe_fmt_p(gt_pv),
  tp_r2, safe_fmt_p(tp_pv)
)

subtitle_txt <- sprintf(
  "Constrained R\u00b2 = %.1f%% | Timepoint tested with Condition(Subject_ID) | %s proteins",
  constrained_var, format(nrow(imp_df), big.mark = ",")
)

# ── PCA labels ──
pca_labels <- c(CRE_T1 = "CRE T1", CRE_T2 = "CRE T2",
                PLA_T1 = "PLA T1", PLA_T2 = "PLA T2",
                H_T1   = "Healthy T1")

pE_dbrda <- ggplot(site_scores, aes(x = dbRDA1, y = dbRDA2,
                                     color = Group_Time, shape = Group_Time)) +
  stat_ellipse(aes(fill = Group_Time), geom = "polygon",
               alpha = 0.10, level = 0.80, show.legend = FALSE) +
  stat_ellipse(aes(group = Group_Time), level = 0.80, linewidth = 0.4,
               linetype = "dashed", show.legend = FALSE) +
  geom_point(size = 2.5, alpha = 0.85) +
  annotate("text", x = Inf, y = Inf, label = anova_label,
           hjust = 1.05, vjust = 1.15,
           size = scale_text(BASE_STAT - 1.2, PE_RDA_W), color = "grey30",
           fontface = "bold") +
  scale_color_manual(values = PCA_COLORS, labels = pca_labels,
                     guide = guide_legend(override.aes = list(size = 2.5))) +
  scale_fill_manual(values = PCA_COLORS, guide = "none") +
  scale_shape_manual(values = PCA_SHAPES, labels = pca_labels) +
  labs(title = "Constrained Ordination (db-RDA)",
       subtitle = subtitle_txt,
       x = sprintf("dbRDA1 (%.1f%%)", ax_var[1]),
       y = sprintf("dbRDA2 (%.1f%%)", ax_var[2]),
       tag = "E'") +
  FIG_THEME + theme(legend.position = c(0.88, 0.12),
                    legend.background = element_rect(fill = alpha("white", 0.8),
                                                     color = NA),
                    legend.title = element_blank(),
                    legend.text = element_text(size = FIG_LEGEND_TEXT),
                    legend.key.size = unit(3, "mm"),
                    legend.spacing.y = unit(0.5, "mm"))

# ── Save ──
audit_df <- data.frame(
  term         = c("Group_Time", "Timepoint_conditioned"),
  R2           = as.numeric(c(gt_r2, tp_r2)),
  p_value      = as.numeric(c(gt_pv, tp_pv)),
  constrained_R2_pct = rep(constrained_var, 2),
  adj_R2       = rep(adj_r2, 2),
  dbRDA1_pct   = rep(as.numeric(ax_var[1]), 2),
  dbRDA2_pct   = rep(as.numeric(ax_var[2]), 2)
)
write.csv(audit_df, file.path(DAT_DIR, "panel_E_dbRDA.csv"),
          row.names = FALSE)

ggsave(file.path(RPT_DIR, "panel_E_dbRDA_SUPP.pdf"), pE_dbrda,
       width = PE_RDA_W, height = PE_RDA_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_E_dbRDA_SUPP.png"), pE_dbrda,
       width = PE_RDA_W, height = PE_RDA_H, units = "mm", dpi = 300)
