# Figure 2 CR — Panel E: Intra-Individual Proteomic Variability (Imputed)
# One boxplot per CR subject, faceted by Supplement. Ordered by median log2FC.
# Annotated with per-subject imputation fractions.
# Outputs: pE (ggplot object), panel_E_imputed_SUPP.pdf/.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F02/a_script/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
})

PE_W <- 160; PE_H <- 90

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
    Timepoint  = factor(Timepoint, levels = c("T1", "T2"))
  )

ann_cols   <- c("uniprot_id", "protein", "gene", "description")
cr_samples <- intersect(meta_cr$Col_ID, setdiff(names(imp_df), ann_cols))
meta_cr    <- meta_cr |> filter(Col_ID %in% cr_samples)

imp_mat <- as.matrix(imp_df[, cr_samples])
n_proteins <- nrow(imp_mat)

# Load imputation mask (TRUE = was missing, i.e. imputed)
mask_df <- read_csv("02_Imputation/c_data/07_imputation_mask.csv",
                    show_col_types = FALSE)
mask_mat <- as.matrix(mask_df[, intersect(cr_samples, names(mask_df))])
rownames(mask_mat) <- mask_df$gene

# Load MAR/MNAR classification
mnar_df <- read_csv("02_Imputation/c_data/02_mar_mnar_classification.csv",
                    show_col_types = FALSE)
mnar_genes <- mnar_df$gene[mnar_df$classification == "MNAR"]

pdf_device <- get_pdf_device()
subjects <- unique(meta_cr$Subject_ID)

lfc_list <- lapply(subjects, function(s) {
  t1_id <- meta_cr$Col_ID[meta_cr$Subject_ID == s & meta_cr$Timepoint == "T1"]
  t2_id <- meta_cr$Col_ID[meta_cr$Subject_ID == s & meta_cr$Timepoint == "T2"]
  if (length(t1_id) != 1 || length(t2_id) != 1) return(NULL)

  lfc <- imp_mat[, t2_id] - imp_mat[, t1_id]  # log2(T2) - log2(T1)

  # Count logFC values depending on at least one imputed value
  t1_imputed <- if (t1_id %in% colnames(mask_mat)) mask_mat[, t1_id] else rep(FALSE, n_proteins)
  t2_imputed <- if (t2_id %in% colnames(mask_mat)) mask_mat[, t2_id] else rep(FALSE, n_proteins)
  either_imputed <- t1_imputed | t2_imputed
  n_imputed_lfc  <- sum(either_imputed)
  pct_imputed    <- 100 * n_imputed_lfc / n_proteins

  # Count MNAR-dependent logFC values
  genes <- imp_df$gene
  mnar_and_imputed <- either_imputed & (genes %in% mnar_genes)
  n_mnar_lfc  <- sum(mnar_and_imputed)
  pct_mnar    <- 100 * n_mnar_lfc / n_proteins

  supp <- as.character(meta_cr$Supplement[meta_cr$Subject_ID == s][1])

  tibble(
    subject       = s,
    supplement    = supp,
    lfc           = as.numeric(lfc),
    n_imputed_lfc = n_imputed_lfc,
    pct_imputed   = pct_imputed,
    n_mnar_lfc    = n_mnar_lfc,
    pct_mnar      = pct_mnar
  )
})

lfc_long <- bind_rows(lfc_list)
lfc_long$supplement <- factor(lfc_long$supplement, levels = c("CRE", "PLA"))

subj_summary <- lfc_long |>
  group_by(subject, supplement, n_imputed_lfc, pct_imputed,
           n_mnar_lfc, pct_mnar) |>
  summarise(
    median_lfc = median(lfc, na.rm = TRUE),
    mad_lfc    = mad(lfc, na.rm = TRUE),
    sd_lfc     = sd(lfc, na.rm = TRUE),
    iqr_lfc    = IQR(lfc, na.rm = TRUE),
    q25        = quantile(lfc, 0.25, na.rm = TRUE),
    q75        = quantile(lfc, 0.75, na.rm = TRUE),
    n_proteins = n(),
    .groups    = "drop"
  )

subj_summary <- subj_summary |>
  arrange(supplement, median_lfc) |>
  mutate(subj_order = factor(subject, levels = unique(subject)))

lfc_long <- lfc_long |>
  mutate(subj_order = factor(subject, levels = levels(subj_summary$subj_order)))

group_summary <- subj_summary |>
  group_by(supplement) |>
  summarise(
    mean_median = mean(median_lfc),
    sd_median   = sd(median_lfc),
    n           = n(),
    .groups     = "drop"
  )

wt <- wilcox.test(median_lfc ~ supplement, data = subj_summary)
n1 <- sum(subj_summary$supplement == "CRE")
n2 <- sum(subj_summary$supplement == "PLA")
r_rb <- 1 - 2 * wt$statistic / (n1 * n2)

mean_pct_imp  <- mean(subj_summary$pct_imputed)
mean_pct_mnar <- mean(subj_summary$pct_mnar)
subtitle_text <- sprintf(
  "%s proteins (imputed) | %.0f%% imputed, %.0f%% MNAR | Wilcoxon %s",
  format(n_proteins, big.mark = ","), mean_pct_imp, mean_pct_mnar, fmt_p(wt$p.value)
)

# Supplement-based fill colors
supp_fill <- c(CRE = scales::alpha("#2166AC", 0.5),
               PLA = scales::alpha("#D6604D", 0.5))

pE <- ggplot(lfc_long, aes(x = subj_order, y = lfc, fill = supplement)) +
  geom_boxplot(width = 0.5, linewidth = 0.3, color = "black",
               outlier.shape = NA, alpha = 0.5) +
  facet_grid(~ supplement, scales = "free_x", space = "free_x") +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  scale_fill_manual(values = supp_fill) +
  labs(x = "Subject",
       y = expression(bold(Delta~log[2]*"FC (T2/T1)")),
       title = "Intra-Individual Proteomic Variability",
       subtitle = subtitle_text,
       tag = "E") +
  FIG_THEME +
  theme(legend.position = "none",
        panel.spacing = unit(3, "mm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = FIG_AXIS_TEXT - 1.5))

write.csv(subj_summary |>
            select(subject, supplement, median_lfc, mad_lfc, sd_lfc,
                   iqr_lfc, q25, q75, n_proteins,
                   n_imputed_lfc, pct_imputed, n_mnar_lfc, pct_mnar),
          file.path(DAT_DIR, "audit_panel_E_imputed.csv"),
          row.names = FALSE)

write.csv(group_summary,
          file.path(DAT_DIR, "audit_panel_E_wilcoxon.csv"),
          row.names = FALSE)

ggsave(file.path(RPT_DIR, "panel_E_imputed_SUPP.pdf"), pE,
       width = PE_W, height = PE_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_E_imputed_SUPP.png"), pE,
       width = PE_W, height = PE_H, units = "mm", dpi = 300)
