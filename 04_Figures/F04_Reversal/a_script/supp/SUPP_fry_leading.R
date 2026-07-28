# SUPP Fry Leading-Edge Dotplot
# Top 25 fry driving proteins ranked by |t_CR_Training|
setwd(here::here())
source("04_Figures/shared/style.R")

pacman::p_load(tidyverse, openxlsx)

RPT_PNG <- "04_Figures/F04_Reversal/b_reports/supp/png/panels"
RPT_PDF <- "04_Figures/F04_Reversal/b_reports/supp/pdf/panels"
DAT <- "04_Figures/F04_Reversal/c_data/panel_supp"

dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# The workbook, not panel_D_fry/driving_proteins.csv: 90_stitch_F04.R sources this
# script after cleanup_after_workbook() has already deleted that c_data subdir.
xlsx_path <- "04_Figures/F04_Reversal/c_data/F04_supplementary.xlsx"
stopifnot("F04 supplementary workbook missing" = file.exists(xlsx_path))
drivers <- read.xlsx(xlsx_path, sheet = "panel_C_fry_driving")

# Load DEP results for t-statistics
source("04_Figures/F04_Reversal/a_script/f04_data.R")
dep <- dep_df %>%
  dplyr::select(gene,
    t_TR = t_CR_Training, t_CvH = t_CRvH_Baseline,
    logFC_CvH = logFC_CRvH_Baseline
  )

# Merge and rank
drivers_ann <- drivers %>%
  left_join(dep, by = "gene") %>%
  filter(!is.na(t_TR)) %>%
  arrange(desc(abs(t_TR))) %>%
  slice_head(n = 25) %>%
  mutate(
    cancer_dir = ifelse(logFC_CvH > 0, "Cancer Up", "Cancer Down"),
    gene       = fct_reorder(gene, abs(t_TR))
  )

# Cancer direction colors
dir_cols <- if (exists("CANCER_DIR_COLORS")) {
  CANCER_DIR_COLORS
} else {
  c("Cancer Up" = "#E57373", "Cancer Down" = "#64B5F6")
}

write.csv(drivers_ann %>% dplyr::select(gene, t_TR, t_CvH, logFC_CvH, cancer_dir),
  file.path(DAT, "SUPP_fry_leading_edge.csv"),
  row.names = FALSE
)

pS_fry_lead <- ggplot(
  drivers_ann,
  aes(
    x = t_TR, y = gene, colour = cancer_dir,
    size = abs(t_CvH)
  )
) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = dir_cols) +
  scale_size_continuous(range = c(1.5, 5), name = "|t| Cancer") +
  labs(
    x = expression(italic(t) ~ "(Training CR)"),
    y = NULL, colour = "Cancer direction",
    title = "Top 25 fry driving proteins",
    subtitle = "Ranked by |t| Training CR"
  ) +
  FIG_THEME +
  theme(legend.position = "right")

save_fig(pS_fry_lead, "SUPP_fry_leading", RPT_PDF, RPT_PNG,
  width = 140, height = 100,
  bg = "white"
)
message("Done: SUPP_fry_leading  [", nrow(drivers_ann), " proteins plotted]")
