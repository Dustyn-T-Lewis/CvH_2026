# Figure 1 — Panel C: Intra-Individual Proteomic Variability
# One boxplot per subject, faceted by Supplement (CRE/PLA), ordered by median log2FC.
# Outputs: pC (ggplot object), panel_C_intra_variability.pdf/.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
})

PC_W <- 160; PC_H <- 90

RPT_DIR <- "04_Figures/F01/b_reports"
DAT_DIR <- "04_Figures/F01/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Load data ---
imp_df <- read_csv("02_Imputation/c_data/01_imputed.csv", show_col_types = FALSE)
meta   <- read_csv("00_input/CvH_meta.csv", show_col_types = FALSE)

ann_cols   <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(imp_df), ann_cols)

# --- Identify paired CR subjects (both T1 and T2 in imputed data) ---
cr_meta <- meta |>
  filter(Group != "PPS", Col_ID %in% samp_names)

paired_subjects <- cr_meta |>
  group_by(Subject_ID) |>
  filter(all(c("T1", "T2") %in% Timepoint)) |>
  ungroup()

cat(sprintf("Paired CR subjects: %d\n", n_distinct(paired_subjects$Subject_ID)))

# --- Compute log2FC per protein per subject ---
imp_mat    <- as.matrix(imp_df[, samp_names])
n_proteins <- nrow(imp_mat)

subjects <- unique(paired_subjects$Subject_ID)

lfc_list <- lapply(subjects, function(s) {
  s_meta  <- paired_subjects |> filter(Subject_ID == s)
  t1_id   <- s_meta$Col_ID[s_meta$Timepoint == "T1"]
  t2_id   <- s_meta$Col_ID[s_meta$Timepoint == "T2"]
  supp    <- s_meta$Supplement[1]

  lfc <- imp_mat[, t2_id] - imp_mat[, t1_id]  # already log2
  tibble(
    subject    = s,
    supplement = supp,
    lfc        = as.numeric(lfc)
  )
})

lfc_long <- bind_rows(lfc_list)
lfc_long$supplement <- factor(lfc_long$supplement, levels = c("CRE", "PLA"))

# --- Subject-level summary ---
subj_summary <- lfc_long |>
  group_by(subject, supplement) |>
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

# Order subjects by median logFC within each supplement facet
subj_summary <- subj_summary |>
  arrange(supplement, median_lfc) |>
  mutate(subj_order = factor(subject, levels = unique(subject)))

lfc_long <- lfc_long |>
  mutate(subj_order = factor(subject, levels = levels(subj_summary$subj_order)))

# --- Wilcoxon test: CRE vs PLA median logFC ---
group_summary <- subj_summary |>
  group_by(supplement) |>
  summarise(
    mean_median = mean(median_lfc),
    sd_median   = sd(median_lfc),
    n           = n(),
    .groups     = "drop"
  )

wt <- wilcox.test(median_lfc ~ supplement, data = subj_summary)
wt_label <- fmt_p(wt$p.value)

subtitle_text <- sprintf(
  "%d proteins per subject | Wilcoxon %s",
  n_proteins, wt_label
)

# --- Fill colors by supplement ---
supp_fill <- c(CRE = "#2166AC", PLA = "#D6604D")

pdf_device <- get_pdf_device()

# --- Plot ---
pC <- ggplot(lfc_long, aes(x = subj_order, y = lfc, fill = supplement)) +
  geom_boxplot(width = 0.5, linewidth = 0.3, color = "black",
               outlier.shape = NA, alpha = 0.5) +
  facet_grid(~ supplement, scales = "free_x", space = "free_x",
             labeller = labeller(supplement = SUPP_LABELS)) +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  scale_fill_manual(values = supp_fill) +
  labs(x = "Subject",
       y = expression(bold(Delta~log[2]*"FC (T2/T1)")),
       title = "Intra-Individual Proteomic Variability",
       subtitle = subtitle_text,
       tag = "C") +
  FIG_THEME +
  theme(legend.position = "none",
        panel.spacing = unit(8, "mm"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = FIG_AXIS_TEXT - 1.5))

# --- Export audit CSVs ---
write.csv(subj_summary |>
            select(subject, supplement, median_lfc, mad_lfc, sd_lfc,
                   iqr_lfc, q25, q75, n_proteins),
          file.path(DAT_DIR, "panel_C_intra_variability.csv"),
          row.names = FALSE)

write.csv(group_summary,
          file.path(DAT_DIR, "panel_C_wilcoxon.csv"),
          row.names = FALSE)

# --- Save figures ---
ggsave(file.path(RPT_DIR, "panel_C_intra_variability.pdf"), pC,
       width = PC_W, height = PC_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_C_intra_variability.png"), pC,
       width = PC_W, height = PC_H, units = "mm", dpi = 300)

cat("Panel C done.\n")
