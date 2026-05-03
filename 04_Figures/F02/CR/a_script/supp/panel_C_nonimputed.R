# Figure 2 CR — Supp Panel C (Non-Imputed): Intra-Individual Proteomic Variability
# Complete-case logFC only (both T1 and T2 observed). No imputation.
# Annotated with n_complete and %MNAR excluded per subject.
# Outputs: pC_nonimp (ggplot object), panel_C_nonimputed_SUPP.pdf/.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F02/a_script/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
})

PC_NI_W <- 160; PC_NI_H <- 100

RPT_DIR <- "04_Figures/F02/CR/b_reports/supp"
DAT_DIR <- "04_Figures/F02/CR/c_data/supp"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

# Non-imputed normalized data
norm_df <- read_csv("01_normalization/c_data/02_normalized.csv",
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
cr_samples <- intersect(meta_cr$Col_ID, setdiff(names(norm_df), ann_cols))
meta_cr    <- meta_cr |> filter(Col_ID %in% cr_samples)

norm_mat <- as.matrix(norm_df[, cr_samples])
n_total  <- nrow(norm_mat)

# Load MAR/MNAR classification for annotation
mnar_df <- read_csv("02_Imputation/c_data/02_mar_mnar_classification.csv",
                    show_col_types = FALSE)
mnar_genes <- mnar_df$gene[mnar_df$classification == "MNAR"]

pdf_device <- get_pdf_device()
subjects <- unique(meta_cr$Subject_ID)

lfc_list <- lapply(subjects, function(s) {
  t1_id <- meta_cr$Col_ID[meta_cr$Subject_ID == s & meta_cr$Timepoint == "T1"]
  t2_id <- meta_cr$Col_ID[meta_cr$Subject_ID == s & meta_cr$Timepoint == "T2"]
  if (length(t1_id) != 1 || length(t2_id) != 1) return(NULL)

  t1_vals <- norm_mat[, t1_id]
  t2_vals <- norm_mat[, t2_id]

  # Complete case: both observed (not NA)
  complete <- !is.na(t1_vals) & !is.na(t2_vals)
  lfc <- t2_vals[complete] - t1_vals[complete]

  # Count MNAR proteins excluded (missing in at least one timepoint)
  missing_genes <- norm_df$gene[!complete]
  n_mnar_excluded <- sum(missing_genes %in% mnar_genes)

  supp <- as.character(meta_cr$Supplement[meta_cr$Subject_ID == s][1])

  tibble(
    subject         = s,
    supplement      = supp,
    lfc             = as.numeric(lfc),
    n_complete      = sum(complete),
    n_missing       = sum(!complete),
    n_mnar_excluded = n_mnar_excluded,
    pct_mnar_excl   = 100 * n_mnar_excluded / n_total
  )
})

lfc_long <- bind_rows(lfc_list)
lfc_long$supplement <- factor(lfc_long$supplement, levels = c("CRE", "PLA"))

subj_summary <- lfc_long |>
  group_by(subject, supplement, n_complete, n_missing,
           n_mnar_excluded, pct_mnar_excl) |>
  summarise(
    median_lfc = median(lfc, na.rm = TRUE),
    mad_lfc    = mad(lfc, na.rm = TRUE),
    sd_lfc     = sd(lfc, na.rm = TRUE),
    iqr_lfc    = IQR(lfc, na.rm = TRUE),
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
    mean_median    = mean(median_lfc),
    sd_median      = sd(median_lfc),
    mean_complete  = mean(n_complete),
    mean_pct_mnar  = mean(pct_mnar_excl),
    n              = n(),
    .groups        = "drop"
  )

wt <- wilcox.test(median_lfc ~ supplement, data = subj_summary)
n1 <- sum(subj_summary$supplement == "CRE")
n2 <- sum(subj_summary$supplement == "PLA")
r_rb <- 1 - 2 * wt$statistic / (n1 * n2)

subtitle_text <- sprintf(
  "Complete-case logFC (non-imputed) | n varies by subject (mean %.0f) | Wilcoxon %s",
  mean(subj_summary$n_complete), fmt_p(wt$p.value)
)

# Per-subject n_complete label at bottom of each box
n_label_df <- subj_summary |>
  mutate(label = as.character(n_complete))

supp_fill <- c(CRE = scales::alpha("#2166AC", 0.5),
               PLA = scales::alpha("#D6604D", 0.5))

pC_nonimp <- ggplot(lfc_long, aes(x = subj_order, y = lfc, fill = supplement)) +
  geom_boxplot(width = 0.5, linewidth = 0.3, color = "black",
               outlier.shape = NA, alpha = 0.5) +
  geom_text(data = n_label_df,
            aes(x = subj_order, y = -1.45, label = label),
            inherit.aes = FALSE, size = scale_text(BASE_COUNT - 1.5, PC_NI_W),
            color = "grey50", angle = 0) +
  facet_grid(~ supplement, scales = "free_x", space = "free_x") +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  scale_fill_manual(values = supp_fill) +
  labs(x = "Subject",
       y = expression(bold(Delta~log[2]*"FC (T2/T1)")),
       title = "Intra-Individual Variability (Non-Imputed)",
       subtitle = subtitle_text,
       tag = "C'") +
  FIG_THEME +
  theme(legend.position = "none",
        panel.spacing = unit(3, "mm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = FIG_AXIS_TEXT - 1.5))

write.csv(subj_summary |>
            select(subject, supplement, n_complete, n_missing,
                   n_mnar_excluded, pct_mnar_excl, median_lfc,
                   mad_lfc, sd_lfc, iqr_lfc),
          file.path(DAT_DIR, "panel_C_nonimputed.csv"),
          row.names = FALSE)

ggsave(file.path(RPT_DIR, "panel_C_nonimputed_SUPP.pdf"), pC_nonimp,
       width = PC_NI_W, height = PC_NI_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_C_nonimputed_SUPP.png"), pC_nonimp,
       width = PC_NI_W, height = PC_NI_H, units = "mm", dpi = 300)
