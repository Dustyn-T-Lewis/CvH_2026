# Panel A — per-sample stacked bar of blood-marker percentages.
# Bars sorted ascending by B_M_ratio (left = cleanest, right = most contaminated).
# Stack 5 markers only (HBB/HBA1/MB/ALB/CKM); muscle/other left implicit.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(RColorBrewer)
})

setwd(rprojroot::find_rstudio_root_file())
source("02-03_Sam's_Results/04_Figures/build_data_index.R")
source("02-03_Sam's_Results/04_Figures/shared/style.R")
source("02-03_Sam's_Results/04_Figures/F02_blood_contamination/a_script/_shared_keys.R")

sam <- readRDS(sam_idx$sam$dalist_rds)
meta <- as.data.frame(sam$metadata)
meta$sample_id <- normalize_sample_id(meta$sample_id)

# Build Group_Time from supp_time, falling back to cancer_time for control
# samples (Sam's CTL rows have supp = NA, so supp_time is NA_T1; map to
# the cancer_time value instead, e.g. CTL_T1).
meta$Group_Time <- ifelse(
  is.na(meta$supp) | meta$supp == "" | grepl("^NA_", meta$supp_time),
  as.character(meta$cancer_time),
  as.character(meta$supp_time)
)

panel_A_data <- meta |>
  select(sample_id, Group_Time, HBB_pct, HBA1_pct, MB_pct,
         ALB_pct, CKM_pct, B_M_ratio) |>
  pivot_longer(c(HBB_pct, HBA1_pct, MB_pct, ALB_pct, CKM_pct),
               names_to = "marker", values_to = "pct") |>
  mutate(marker = sub("_pct$", "", marker))

bm_order <- meta |> arrange(B_M_ratio) |> pull(sample_id)
panel_A_data$sample_id <- factor(panel_A_data$sample_id, levels = bm_order)

gt_levels <- intersect(c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2",
                         "CTL_T1", "CTL_T2"), unique(panel_A_data$Group_Time))
panel_A_data$Group_Time <- factor(panel_A_data$Group_Time, levels = gt_levels)

write_csv(panel_A_data, "02-03_Sam's_Results/04_Figures/F02_blood_contamination/c_data/panel_A_stacked_data.csv")

marker_palette <- setNames(brewer.pal(5, "Set2"),
                           c("HBB", "HBA1", "MB", "ALB", "CKM"))

p_A <- ggplot(panel_A_data, aes(sample_id, pct, fill = marker)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = marker_palette, name = "Marker") +
  labs(x = "Sample (sorted by B_M_ratio, ascending)",
       y = "% of total signal",
       title = "Per-sample blood-marker share") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        legend.position = "right")

ggsave("02-03_Sam's_Results/04_Figures/F02_blood_contamination/b_reports/main/png/panels/panel_A.png",
       p_A, width = 7, height = 4, dpi = 300, bg = "white")

panel_A <- p_A
