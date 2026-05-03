# SUPP Fry Leading-Edge Dotplot
# Top 25 fry driving proteins ranked by |t_Training_CR|
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

library(tidyverse)
library(openxlsx)

RPT_PNG <- "04_Figures/Reversal/b_reports/supp/png/panels"
RPT_PDF <- "04_Figures/Reversal/b_reports/supp/pdf/panels"
DAT     <- "04_Figures/Reversal/c_data/panel_supp"

dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT,     recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# -- Load driving proteins (CSV with xlsx fallback) ----------------------------
csv_path  <- "04_Figures/Reversal/c_data/panel_D_fry/driving_proteins.csv"
xlsx_path <- "04_Figures/Reversal/c_data/Reversal_supplementary.xlsx"

if (file.exists(csv_path)) {
  drivers <- read_csv(csv_path, show_col_types = FALSE)
  message("Read driving proteins from CSV")
} else if (file.exists(xlsx_path)) {
  drivers <- read.xlsx(xlsx_path, sheet = "panel_C_fry_driving")
  message("Read driving proteins from xlsx fallback")
} else {
  stop("Cannot find driving proteins at CSV or xlsx path")
}

# -- Load DEP results for t-statistics ----------------------------------------
dep <- read_csv("03_DEP/c_data/03_combined_results_CRvH.csv",
                show_col_types = FALSE) %>%
  dplyr::select(gene, t_TR = t_Training_CR, t_CvH = t_Cancer_vs_Healthy,
         logFC_CvH = logFC_Cancer_vs_Healthy)

# -- Merge and rank ------------------------------------------------------------
drivers_ann <- drivers %>%
  left_join(dep, by = "gene") %>%
  filter(!is.na(t_TR)) %>%
  arrange(desc(abs(t_TR))) %>%
  slice_head(n = 25) %>%
  mutate(
    cancer_dir = ifelse(logFC_CvH > 0, "Cancer Up", "Cancer Down"),
    gene       = fct_reorder(gene, abs(t_TR))
  )

# -- Cancer direction colors ---------------------------------------------------
dir_cols <- if (exists("CANCER_DIR_COLORS")) {
  CANCER_DIR_COLORS
} else {
  c("Cancer Up" = "#E57373", "Cancer Down" = "#64B5F6")
}

# -- Export CSV ----------------------------------------------------------------
write.csv(drivers_ann %>% dplyr::select(gene, t_TR, t_CvH, logFC_CvH, cancer_dir),
          file.path(DAT, "SUPP_fry_leading_edge.csv"), row.names = FALSE)

# -- Plot ----------------------------------------------------------------------
pS_fry_lead <- ggplot(drivers_ann,
                      aes(x = t_TR, y = gene, colour = cancer_dir,
                          size = abs(t_CvH))) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = dir_cols) +
  scale_size_continuous(range = c(1.5, 5), name = "|t| Cancer") +
  labs(x = expression(italic(t) ~ "(Training CR)"),
       y = NULL, colour = "Cancer direction",
       title = "Top 25 fry driving proteins",
       subtitle = "Ranked by |t| Training CR") +
  FIG_THEME +
  theme(legend.position = "right")

ggsave(file.path(RPT_PNG, "SUPP_fry_leading.png"), pS_fry_lead,
       width = 140, height = 100, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(RPT_PDF, "SUPP_fry_leading.pdf"), pS_fry_lead,
       width = 140, height = 100, units = "mm", device = pdf_device)

message("Done: SUPP_fry_leading  [", nrow(drivers_ann), " proteins plotted]")
