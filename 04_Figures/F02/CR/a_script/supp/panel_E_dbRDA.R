# Figure 2 CR — Supp Panel E (db-RDA): Constrained Ordination
# CR-only. Full: Supplement * Timepoint. Within: Timepoint + Condition(Subject_ID).
# Outputs: pE_dbrda (ggplot object), panel_E_dbRDA_SUPP.pdf/.png

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

PE_RDA_W <- 145; PE_RDA_H <- 100

RPT_DIR <- "04_Figures/F02/CR/b_reports/supp"
DAT_DIR <- "04_Figures/F02/CR/c_data/supp"
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

# db-RDA: constrained ordination
dist_mat <- vegdist(scale(t(imp_mat)), method = "euclidean")

# Full model: Supplement + Timepoint + Supplement:Timepoint (no conditioning)
rda_full <- dbrda(dist_mat ~ Supplement * Timepoint, data = meta_cr)

# Within-subject model: Timepoint + Condition(Subject_ID)
# Condition(Subject_ID) absorbs Supplement since subject is nested within supplement
rda_within <- dbrda(dist_mat ~ Timepoint + Condition(Subject_ID), data = meta_cr)

# Use full model for visualization
site_scores  <- as.data.frame(scores(rda_full, display = "sites", choices = 1:2))
site_scores$Col_ID <- rownames(site_scores)
site_scores <- left_join(site_scores, meta_cr, by = "Col_ID")

# Variance from full model
eigenvals_full  <- eigenvals(rda_full)
constrained_eig <- eigenvals_full[grepl("^dbRDA", names(eigenvals_full))]
total_inertia   <- rda_full$tot.chi
constrained_var <- rda_full$CCA$tot.chi / total_inertia * 100

ax_var <- 100 * constrained_eig[1:2] / total_inertia

# ANOVA: full model terms and within-subject model for time
set.seed(42)
anova_full   <- anova(rda_full, by = "terms", permutations = 999)
anova_within <- anova(rda_within, by = "terms", permutations = 999)

# Extract from data frames by row name matching
get_anova_val <- function(av, term, col) {
  rn <- rownames(av)
  idx <- grep(paste0("^", term, "$"), rn)
  if (length(idx) == 0) return(NA_real_)
  av[idx, col]
}

# Column name depends on vegan version
var_col <- intersect(c("Variance", "SumOfSqs"), colnames(anova_full))[1]

total_ss_full   <- sum(anova_full[[var_col]], na.rm = TRUE)
total_ss_within <- sum(anova_within[[var_col]], na.rm = TRUE)

supp_r2 <- get_anova_val(anova_full, "Supplement", var_col) / total_ss_full
supp_pv <- get_anova_val(anova_full, "Supplement", "Pr(>F)")
time_r2 <- get_anova_val(anova_within, "Timepoint", var_col) / total_ss_within
time_pv <- get_anova_val(anova_within, "Timepoint", "Pr(>F)")
int_r2  <- get_anova_val(anova_full, "Supplement:Timepoint", var_col) / total_ss_full
int_pv  <- get_anova_val(anova_full, "Supplement:Timepoint", "Pr(>F)")

adj_r2 <- RsquareAdj(rda_full)$adj.r.squared

# Safe p-value formatter
safe_fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  fmt_p(p)
}

# Build annotation
anova_label <- sprintf(
  "db-RDA (adj. R\u00b2 = %.3f)\nSupp  R\u00b2 = %.3f,  %s\nTime  R\u00b2 = %.3f,  %s  (cond.)\nSupp\u00d7Time  R\u00b2 = %.3f,  %s",
  adj_r2,
  supp_r2, safe_fmt_p(supp_pv),
  time_r2, safe_fmt_p(time_pv),
  int_r2,  safe_fmt_p(int_pv)
)

subtitle_txt <- sprintf(
  "Constrained R\u00b2 = %.1f%% | Time tested with Condition(Subject_ID) | %s proteins",
  constrained_var, format(nrow(imp_df), big.mark = ",")
)

# CR-only palette (4 groups)
cr_pca_colors <- PCA_COLORS[c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2")]
cr_pca_shapes <- PCA_SHAPES[c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2")]

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
  scale_color_manual(values = cr_pca_colors,
                     labels = c("CRE T1", "CRE T2", "PLA T1", "PLA T2"),
                     guide = guide_legend(override.aes = list(size = 2.5))) +
  scale_fill_manual(values = cr_pca_colors, guide = "none") +
  scale_shape_manual(values = cr_pca_shapes,
                     labels = c("CRE T1", "CRE T2", "PLA T1", "PLA T2")) +
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

audit_df <- data.frame(
  term         = c("Supplement", "Timepoint_conditioned", "Supplement:Timepoint"),
  R2           = as.numeric(c(supp_r2, time_r2, int_r2)),
  p_value      = as.numeric(c(supp_pv, time_pv, int_pv)),
  constrained_R2_pct = rep(constrained_var, 3),
  adj_R2       = rep(adj_r2, 3),
  dbRDA1_pct   = rep(as.numeric(ax_var[1]), 3),
  dbRDA2_pct   = rep(as.numeric(ax_var[2]), 3)
)
write.csv(audit_df, file.path(DAT_DIR, "panel_E_dbRDA.csv"),
          row.names = FALSE)

ggsave(file.path(RPT_DIR, "panel_E_dbRDA_SUPP.pdf"), pE_dbrda,
       width = PE_RDA_W, height = PE_RDA_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_E_dbRDA_SUPP.png"), pE_dbrda,
       width = PE_RDA_W, height = PE_RDA_H, units = "mm", dpi = 300)
