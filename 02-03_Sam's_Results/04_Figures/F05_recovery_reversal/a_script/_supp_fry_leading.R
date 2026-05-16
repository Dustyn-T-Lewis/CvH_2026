#!/usr/bin/env Rscript
# Sourced by 02_supp_panels.R — expects style.R already loaded.
# F05 Supplementary Panel F: fry Leading-Edge Proteins
# Supports main Panel C — dotplot of top 20-30 driving proteins ranked by |t-stat|
# in Training CR, colored by cancer direction (Cancer_vs_Healthy up/down).

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
})

BASE    <- "02-03_Sam's_Results/04_Figures/F05_recovery_reversal"
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
DAT     <- file.path(BASE, "c_data", "panel_supp")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT,     recursive = TRUE, showWarnings = FALSE)
pdf_device <- grDevices::pdf  # force base pdf (cairo DLL fails)

# Data
# Try xlsx sheet first, fall back to computing from DEP results
xlsx_path <- file.path(BASE, "c_data", "F05_supplementary.xlsx")
driving_df <- NULL
if (file.exists(xlsx_path)) {
  tryCatch({
    driving_df <- read_excel(xlsx_path, sheet = "panel_C_fry_driving")
    message("  Read fry driving proteins from xlsx")
  }, error = function(e) {
    message("  Could not read panel_C_fry_driving sheet: ", e$message)
  })
}

if (is.null(driving_df) || nrow(driving_df) == 0) {
  # Fall back: identify driving proteins from DEP results
  dep <- read_csv("02-03_Sam's_Results/03_DEP/c_data/F05_combined_CvHvTCR.csv",
                  show_col_types = FALSE)

  # Cancer-significant DEPs by pi-score
  cancer_up <- dep |>
    filter(pi_score_Cancer_vs_Healthy < 0.05 & logFC_Cancer_vs_Healthy > 0) |>
    pull(gene)
  cancer_dn <- dep |>
    filter(pi_score_Cancer_vs_Healthy < 0.05 & logFC_Cancer_vs_Healthy < 0) |>
    pull(gene)

  # For driving proteins: those in cancer sets whose Training_CR t-stat
  # opposes the cancer direction (reversal drivers).
  # Use base R to avoid AnnotationDbi/dplyr masking issues from GO Slim panel.
  driving_df <- dep[dep$gene %in% c(cancer_up, cancer_dn), ]
  driving_df$set       <- ifelse(driving_df$gene %in% cancer_up, "cancer_up", "cancer_dn")
  driving_df$logFC_TCR <- driving_df$logFC_Training_CR
  driving_df$t_TCR     <- driving_df$t_Training_CR
  driving_df$is_driving <- (driving_df$set == "cancer_up" & driving_df$t_TCR < 0) |
                            (driving_df$set == "cancer_dn" & driving_df$t_TCR > 0)
  driving_df <- driving_df[driving_df$is_driving,
                            c("gene", "set", "logFC_Cancer_vs_Healthy",
                              "logFC_TCR", "t_TCR", "is_driving")]
  message("  Computed driving proteins from DEP results")
}

# Ensure required columns exist (base R — avoid dplyr masking)
if (!"t_TCR" %in% names(driving_df) && "t_Training_CR" %in% names(driving_df)) {
  names(driving_df)[names(driving_df) == "t_Training_CR"] <- "t_TCR"
}
if (!"logFC_TCR" %in% names(driving_df) && "logFC_Training_CR" %in% names(driving_df)) {
  names(driving_df)[names(driving_df) == "logFC_Training_CR"] <- "logFC_TCR"
}
# Map any variant column names from xlsx
t_candidates <- c("t_Training_CR", "t_test", "t_training_cr", "t")
for (tc in t_candidates) {
  if (!"t_TCR" %in% names(driving_df) && tc %in% names(driving_df)) {
    names(driving_df)[names(driving_df) == tc] <- "t_TCR"; break
  }
}
lfc_candidates <- c("logFC_Training_CR", "logFC_training_cr", "logFC")
for (lc in lfc_candidates) {
  if (!"logFC_TCR" %in% names(driving_df) && lc %in% names(driving_df)) {
    names(driving_df)[names(driving_df) == lc] <- "logFC_TCR"; break
  }
}

# Top 25 by |t-stat| (base R to avoid dplyr masking)
driving_df$abs_t <- abs(driving_df$t_TCR)
driving_df <- driving_df[order(-driving_df$abs_t), ]
top_driving <- head(driving_df, 25)

write_csv(top_driving, file.path(DAT, "SUPP_fry_leading_edge.csv"))

# Ensure set column exists
if (!"set" %in% names(top_driving)) {
  top_driving <- top_driving |>
    mutate(set = ifelse(logFC_Cancer_vs_Healthy > 0, "cancer_up", "cancer_dn"))
}

# Plot
dir_colors <- c("cancer_up" = "#D6604D", "cancer_dn" = "#4393C3")
dir_labels <- c("cancer_up" = "Cancer Up (reversed down)",
                "cancer_dn" = "Cancer Down (reversed up)")

pS_fry_lead <- ggplot(top_driving,
                      aes(x = t_TCR,
                          y = reorder(gene, abs(t_TCR)),
                          color = set)) +
  geom_point(size = 2) +
  geom_segment(aes(xend = 0, yend = reorder(gene, abs(t_TCR))),
               linewidth = 0.4) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  scale_color_manual(values = dir_colors, labels = dir_labels,
                     name = "Cancer DEP set") +
  labs(title    = "fry Leading-Edge Proteins (Reversal Drivers)",
       subtitle = sprintf("Top %d by |t(Training CR)| | colored by cancer direction",
                          nrow(top_driving)),
       x = "t-statistic (Training CR)",
       y = NULL) +
  FIG_THEME +
  theme(legend.position  = "bottom",
        legend.direction = "horizontal",
        axis.text.y      = element_text(size = FIG_AXIS_TEXT, face = "italic"))

PW <- 89; PH <- 100
ggsave(file.path(RPT_PNG, "SUPP_fry_leading.png"), pS_fry_lead,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_fry_leading.pdf"), pS_fry_lead,
       width = PW, height = PH, units = "mm", device = pdf_device)

message("F05 SUPP Panel F (fry leading edge) saved")

# Expose for composite
pS_lead_title    <- "fry Leading-Edge Proteins (Reversal Drivers)"
pS_lead_subtitle <- sprintf("Top %d by |t(Training CR)| | colored by cancer direction",
                             nrow(top_driving))
pS_fry_lead      <- strip_for_composite(pS_fry_lead)
