# F03 Panel A: Training Response (Creatine) Volcano Ring
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/volcano_ring.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

PA_W <- 190
PH   <- 180

RPT <- "04_Figures/F03/b_reports"
DAT <- "04_Figures/F03/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_A"), recursive = TRUE, showWarnings = FALSE)

dep_df <- read_csv("03_DEP/c_data/03_combined_results_CR.csv", show_col_types = FALSE)

fgsea_all <- read_csv("04_Figures/F02/c_data/06_panel_H_CR_fgsea_results.csv",
                       show_col_types = FALSE)

# Select ring terms — use Hallmark + GO:BP (no GO Slim in this cache)
avail_dbs <- unique(fgsea_all$database)
ring_dbs <- intersect(c("Hallmark", "GO Slim", "GO:BP"), avail_dbs)
if (length(ring_dbs) == 0) ring_dbs <- avail_dbs[1:2]

top_terms_A <- select_ring_terms(fgsea_all, "Training_CRE",
                                  databases = ring_dbs)
ring_A      <- build_ring_with_gaps(top_terms_A, "Training_CRE", fgsea_all,
                                     databases = ring_dbs)

pA <- make_volcano_ring(
  de_df = dep_df, go_df = fgsea_all, contrast = "Training_CRE",
  title = NULL,
  contrast_title    = "Training Response (Creatine)",
  contrast_subtitle = "CRE_T2 \u2212 CRE_T1  |  2,582 proteins \u00b7 limma + dupCor \u00b7 missForest",
  ring_data_override = ring_A,
  databases     = ring_dbs,
  label_size    = scale_text(BASE_PATHWAY, PA_W),
  title_size    = scale_text(BASE_TAG, PA_W),
  point_size    = 1.2,
  point_alpha   = 0.55,
  count_label_size = scale_text(BASE_COUNT, PA_W)
)

ggsave(file.path(RPT, "panel_A_volcano_CRE.pdf"), pA,
       width = PA_W, height = PH, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT, "panel_A_volcano_CRE.png"), pA,
       width = PA_W, height = PH, units = "mm", dpi = 300)

ring_data_A <- attr(pA, "ring_data")
if (!is.null(ring_data_A) && nrow(ring_data_A) > 0) {
  write_csv(ring_data_A %>% dplyr::select(-gene_list),
            file.path(DAT, "panel_A", "ring_terms.csv"))
}

dep_df %>%
  transmute(
    gene,
    log2_fold_change = round(logFC_Training_CRE, 4),
    neg_log10_pvalue = round(-log10(P.Value_Training_CRE), 4),
    pi_score         = round(pi_score_Training_CRE, 6),
    adjusted_pvalue  = round(adj.P.Val_Training_CRE, 6),
    direction = case_when(
      pi_score_Training_CRE < 0.05 & logFC_Training_CRE > 0 ~ "Up",
      pi_score_Training_CRE < 0.05 & logFC_Training_CRE < 0 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  filter(!is.na(log2_fold_change), !is.na(neg_log10_pvalue)) %>%
  arrange(pi_score) %>%
  write_csv(file.path(DAT, "panel_A", "volcano_CRE.csv"))

message("F03 Panel A done")
