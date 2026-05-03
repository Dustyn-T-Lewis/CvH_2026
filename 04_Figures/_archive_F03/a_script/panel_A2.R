# F03 Panel A2: Training Response (Placebo) Volcano Ring
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
dir.create(file.path(DAT, "panel_A2"), recursive = TRUE, showWarnings = FALSE)

dep_df <- read_csv("03_DEP/c_data/03_combined_results_CR.csv", show_col_types = FALSE)

fgsea_all <- read_csv("04_Figures/F02/c_data/06_panel_H_CR_fgsea_results.csv",
                       show_col_types = FALSE)

# Select ring terms — use Hallmark + GO:BP (no GO Slim in this cache)
avail_dbs <- unique(fgsea_all$database)
ring_dbs <- intersect(c("Hallmark", "GO Slim", "GO:BP"), avail_dbs)
if (length(ring_dbs) == 0) ring_dbs <- avail_dbs[1:2]

top_terms_A2 <- select_ring_terms(fgsea_all, "Training_PLA",
                                  databases = ring_dbs)
ring_A2      <- build_ring_with_gaps(top_terms_A2, "Training_PLA", fgsea_all,
                                     databases = ring_dbs)

pA2 <- make_volcano_ring(
  de_df = dep_df, go_df = fgsea_all, contrast = "Training_PLA",
  title = NULL,
  contrast_title    = "Training Response (Placebo)",
  contrast_subtitle = "PLA_T2 \u2212 PLA_T1  |  2,582 proteins \u00b7 limma + dupCor \u00b7 missForest",
  ring_data_override = ring_A2,
  databases     = ring_dbs,
  label_size    = scale_text(BASE_PATHWAY, PA_W),
  title_size    = scale_text(BASE_TAG, PA_W),
  point_size    = 1.2,
  point_alpha   = 0.55,
  count_label_size = scale_text(BASE_COUNT, PA_W)
)

ggsave(file.path(RPT, "panel_A2_volcano_PLA.pdf"), pA2,
       width = PA_W, height = PH, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT, "panel_A2_volcano_PLA.png"), pA2,
       width = PA_W, height = PH, units = "mm", dpi = 300)

ring_data_A2 <- attr(pA2, "ring_data")
if (!is.null(ring_data_A2) && nrow(ring_data_A2) > 0) {
  write_csv(ring_data_A2 %>% dplyr::select(-gene_list),
            file.path(DAT, "panel_A2", "ring_terms.csv"))
}

dep_df %>%
  transmute(
    gene,
    log2_fold_change = round(logFC_Training_PLA, 4),
    neg_log10_pvalue = round(-log10(P.Value_Training_PLA), 4),
    pi_score         = round(pi_score_Training_PLA, 6),
    adjusted_pvalue  = round(adj.P.Val_Training_PLA, 6),
    direction = case_when(
      pi_score_Training_PLA < 0.05 & logFC_Training_PLA > 0 ~ "Up",
      pi_score_Training_PLA < 0.05 & logFC_Training_PLA < 0 ~ "Down",
      TRUE ~ "NS"
    )
  ) %>%
  filter(!is.na(log2_fold_change), !is.na(neg_log10_pvalue)) %>%
  arrange(pi_score) %>%
  write_csv(file.path(DAT, "panel_A2", "volcano_PLA.csv"))

message("F03 Panel A2 done")
