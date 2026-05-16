#!/usr/bin/env Rscript
# SUPP Panel E: fry Leading-Edge Proteins
# Supports main Panel C — shows the top driving proteins from the fry
# rotation test, ranked by |t-stat| and colored by concordant direction.
# Sourced by 02_supp_panels.R — expects style.R already loaded.
# Exports: pS_fry_lead (ggplot)

library(readxl)

BASE <- "02-03_Sam's_Results/04_Figures/F04_training_concordance"
DAT  <- file.path(BASE, "c_data", "panel_supp")
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

# Data
xlsx_path <- file.path(BASE, "c_data", "F04_supplementary.xlsx")
driving_df <- tryCatch({
  df <- read_excel(xlsx_path, sheet = "panel_C_fry_driving")
  message("Read fry driving proteins from xlsx")
  as.data.frame(df)
}, error = function(e) {
  message("xlsx read failed: ", e$message)
  NULL
})

# Fallback: read from DEP results and reconstruct (base R to avoid dplyr masking)
if (is.null(driving_df)) {
  message("Reconstructing driving proteins from DEP results")
  dep_df <- as.data.frame(readr::read_csv(
    "02-03_Sam's_Results/03_DEP/c_data/F04_combined_CREvPLA.csv",
    show_col_types = FALSE))
  sig <- dep_df[!is.na(dep_df$pi_score_Training_CRE) & dep_df$pi_score_Training_CRE < 0.05, ]
  driving_df <- data.frame(
    gene      = sig$gene,
    t_test    = sig$t_Training_PLA,
    logFC_TY  = sig$logFC_Training_CRE,
    logFC_TO  = sig$logFC_Training_PLA,
    direction = ifelse(sig$logFC_Training_CRE > 0, "Up", "Down"),
    stringsAsFactors = FALSE
  )
  driving_df <- driving_df[!is.na(driving_df$t_test), ]
}

# Standardise column names (base R)
for (old_nm in c("t_stat", "t", "t_Training_PLA", "t_training_pla", "t_training_old")) {
  if (!"t_test" %in% names(driving_df) && old_nm %in% names(driving_df))
    names(driving_df)[names(driving_df) == old_nm] <- "t_test"
}
if (!"direction" %in% names(driving_df)) {
  lfc_col <- intersect(c("logFC_TY", "logFC_Training_CRE"), names(driving_df))[1]
  if (!is.na(lfc_col)) driving_df$direction <- ifelse(driving_df[[lfc_col]] > 0, "Up", "Down")
}
for (old_nm in c("logFC_Training_CRE", "logFC_training_cre", "logFC_Training_Young")) {
  if (!"logFC_TY" %in% names(driving_df) && old_nm %in% names(driving_df))
    names(driving_df)[names(driving_df) == old_nm] <- "logFC_TY"
}
for (old_nm in c("logFC_Training_PLA", "logFC_training_pla", "logFC_Training_Old")) {
  if (!"logFC_TO" %in% names(driving_df) && old_nm %in% names(driving_df))
    names(driving_df)[names(driving_df) == old_nm] <- "logFC_TO"
}

# Top 20 by |t-stat| (base R)
driving_df$abs_t <- abs(driving_df$t_test)
driving_df <- driving_df[order(-driving_df$abs_t), ]
top_df <- head(driving_df, 20)
top_df$gene <- factor(top_df$gene, levels = rev(top_df$gene))
top_df$dir_label <- ifelse(top_df$direction == "Up", "Concordant Up", "Concordant Down")

write_csv(top_df, file.path(DAT, "SUPP_fry_leading_edge.csv"))
message(sprintf("Top 20 fry driving proteins (of %d total)", nrow(driving_df)))

# Plot: lollipop of top 20 by |t-stat|
DIR_COLS <- c("Concordant Up" = "#D6604D", "Concordant Down" = "#4393C3")

pS_fry_lead <- ggplot(top_df, aes(x = t_test, y = gene, color = dir_label)) +
  geom_segment(aes(x = 0, xend = t_test, y = gene, yend = gene),
               linewidth = 0.6) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = DIR_COLS, name = "Direction") +
  labs(title = "fry Leading-Edge Proteins",
       subtitle = "Top 20 driving proteins ranked by |t-stat| in Training PLA",
       x = "t-statistic (Training PLA)", y = NULL) +
  FIG_THEME

# Standalone save
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)

PW <- 89; PH <- 70
ggsave(file.path(RPT_PNG, "SUPP_fry_leading.png"), pS_fry_lead,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_fry_leading.pdf"), pS_fry_lead,
       width = PW, height = PH, units = "mm", device = pdf_device)

pS_fry_lead <- strip_for_composite(pS_fry_lead)

message("SUPP Panel E (fry leading edge) done")
