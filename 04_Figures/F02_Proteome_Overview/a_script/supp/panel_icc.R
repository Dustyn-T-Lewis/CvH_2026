# Figure 2 (CRvH) — Supp Panel A (ICC): Intraclass Correlation Coefficient
# ICC(3,1) per protein within each group (CRE, PLA, Healthy).
# Healthy has only T1 -- ICC not computable; reported for CRE and PLA only.
# High ICC = stable subject trait; low ICC = responsive to training.
# Outputs: pA_ICC (ggplot object), icc.pdf/.png

setwd(here::here())
source("04_Figures/shared/style.R")

pacman::p_load(dplyr, tidyr, stringr, readr, ggplot2, psych)

PA_ICC_W <- 110; PA_ICC_H <- 120

RPT_DIR <- "04_Figures/F02_Proteome_Overview/b_reports/supp/panels"
DAT_DIR <- "04_Figures/F02_Proteome_Overview/c_data/supp"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Load data & metadata ──
.dal <- readRDS("02_Normalization/c_data/DAList_normalized.rds")
norm_df <- tibble::as_tibble(cbind(
  .dal$annotation[, c("uniprot_id", "protein", "gene", "description")],
  as.data.frame(.dal$data)))
meta_full <- read_csv("00_input/CvH_meta.csv", show_col_types = FALSE)

ann_cols   <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(norm_df), ann_cols)

# Only CR subjects with paired T1/T2 (CRE, PLA)
meta <- meta_full |>
  filter(Group %in% c("CR_CRE", "CR_PLA"), Col_ID %in% samp_names)

subj_counts <- meta |> count(Subject_ID)
paired_subj <- subj_counts$Subject_ID[subj_counts$n == 2]
meta <- meta |> filter(Subject_ID %in% paired_subj)

log_mat <- as.matrix(norm_df[, samp_names])
rownames(log_mat) <- norm_df$gene

pdf_device <- get_pdf_device()

# ── ICC(3,1) per protein within each supplement group ──
compute_icc_per_group <- function(supp_group) {
  grp_meta <- meta |> filter(Supplement == supp_group)
  subjects <- unique(grp_meta$Subject_ID)

  t1_ids <- grp_meta$Col_ID[grp_meta$Timepoint == "T1"]
  t2_ids <- grp_meta$Col_ID[grp_meta$Timepoint == "T2"]

  # Match subjects for paired structure
  t1_subj <- grp_meta$Subject_ID[grp_meta$Timepoint == "T1"]
  t2_subj <- grp_meta$Subject_ID[grp_meta$Timepoint == "T2"]
  common  <- intersect(t1_subj, t2_subj)

  t1_idx <- match(common, t1_subj)
  t2_idx <- match(common, t2_subj)

  t1_mat <- log_mat[, t1_ids[t1_idx], drop = FALSE]
  t2_mat <- log_mat[, t2_ids[t2_idx], drop = FALSE]

  n_prot <- nrow(log_mat)
  icc_vals <- numeric(n_prot)

  for (i in seq_len(n_prot)) {
    t1_vals <- t1_mat[i, ]
    t2_vals <- t2_mat[i, ]

    ok <- !is.na(t1_vals) & !is.na(t2_vals)
    if (sum(ok) < 3) { icc_vals[i] <- NA_real_; next }

    rating_mat <- cbind(t1_vals[ok], t2_vals[ok])
    icc_res <- tryCatch(
      psych::ICC(rating_mat, missing = FALSE, lmer = FALSE),
      error = function(e) NULL
    )
    if (is.null(icc_res)) { icc_vals[i] <- NA_real_; next }
    icc_vals[i] <- icc_res$results$ICC[3]
  }

  tibble(gene = norm_df$gene, group = supp_group, icc = icc_vals)
}

icc_cre <- compute_icc_per_group("CRE")
icc_pla <- compute_icc_per_group("PLA")
icc_df  <- bind_rows(icc_cre, icc_pla) |> filter(!is.na(icc))
icc_df$group <- factor(icc_df$group, levels = c("CRE", "PLA"))

# ── Bootstrap 95% CI on median ICC per group ──
set.seed(42)
boot_median_ci <- function(x, R = 2000, conf = 0.95) {
  meds <- replicate(R, median(sample(x, replace = TRUE)))
  qs   <- quantile(meds, c((1 - conf) / 2, (1 + conf) / 2))
  c(lower = unname(qs[1]), upper = unname(qs[2]))
}

icc_summary <- icc_df |>
  group_by(group) |>
  summarise(
    n        = n(),
    med_icc  = median(icc),
    ci_lo    = boot_median_ci(icc)[["lower"]],
    ci_hi    = boot_median_ci(icc)[["upper"]],
    pct_high = 100 * mean(icc > 0.75),
    pct_low  = 100 * mean(icc < 0.40),
    .groups  = "drop"
  )

wt <- wilcox.test(icc ~ group, data = icc_df)

sub_txt <- sprintf(
  "ICC(3,1) | CRE: %.2f [%.2f, %.2f] | PLA: %.2f [%.2f, %.2f] | Wilcoxon %s",
  icc_summary$med_icc[1], icc_summary$ci_lo[1], icc_summary$ci_hi[1],
  icc_summary$med_icc[2], icc_summary$ci_lo[2], icc_summary$ci_hi[2],
  fmt_p(wt$p.value)
)

SUPP_FILL <- c(CRE = "#2166AC", PLA = "#D6604D")

pA_ICC <- ggplot(icc_df, aes(x = group, y = icc, fill = group)) +
  geom_violin(alpha = 0.5, linewidth = 0.3, color = "black", scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.3,
               color = "black", fill = "white", coef = 0) +
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "grey50",
             linewidth = 0.4) +
  annotate("text", x = 2.4, y = 0.77, label = "Good (>0.75)",
           hjust = 0, size = scale_text(BASE_STAT - 1, PA_ICC_W),
           color = "grey50", fontface = "italic") +
  geom_label(data = icc_summary,
             aes(x = group, y = 1.05,
                 label = sprintf("%.2f [%.2f, %.2f]", med_icc, ci_lo, ci_hi)),
             size = scale_text(BASE_COUNT + 0.5, PA_ICC_W),
             fontface = "bold", fill = scales::alpha("white", 0.8),
             linewidth = 0.2, label.padding = unit(1.5, "pt")) +
  scale_fill_manual(values = SUPP_FILL) +
  coord_cartesian(ylim = c(-0.2, 1.1)) +
  labs(title = "Test-Retest Reliability (ICC)",
       subtitle = sub_txt,
       x = NULL, y = "ICC(3,1)",
       tag = "A'") +
  FIG_THEME + theme(legend.position = "none")

# ── Save ──
write.csv(icc_summary, file.path(DAT_DIR, "icc.csv"),
          row.names = FALSE)

ggsave(file.path(RPT_DIR, "icc.pdf"), pA_ICC,
       width = PA_ICC_W, height = PA_ICC_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "icc.png"), pA_ICC,
       width = PA_ICC_W, height = PA_ICC_H, units = "mm", dpi = 300)
