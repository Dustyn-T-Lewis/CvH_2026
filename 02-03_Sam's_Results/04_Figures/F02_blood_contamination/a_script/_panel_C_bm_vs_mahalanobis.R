# Panel C — cross-pipeline scatter: Sam's B_M_ratio vs our Mahalanobis distance.
# Color = consensus_outlier from our 4-method 01_normalization outlier_diag.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(ggrepel)
})

setwd(rprojroot::find_rstudio_root_file())
source("02-03_Sam's_Results/04_Figures/build_data_index.R")
source("02-03_Sam's_Results/04_Figures/shared/style.R")
source("02-03_Sam's_Results/04_Figures/F02_blood_contamination/a_script/_shared_keys.R")

sam <- readRDS(sam_idx$sam$dalist_rds)
sam_meta <- as.data.frame(sam$metadata) |>
  mutate(sample_id = normalize_sample_id(sample_id))

our_intermediates <- readRDS("01_normalization/c_data/00_report_intermediates.rds")
our_outlier <- our_intermediates$outlier_diag |>
  mutate(sample_id = normalize_sample_id(Col_ID))

joined <- inner_join(
  sam_meta |> select(sample_id, Group_Time = supp_time, B_M_ratio),
  our_outlier |> select(sample_id, mahal_dist, n_flags, consensus_outlier),
  by = "sample_id"
)

bm_cut <- quantile(joined$B_M_ratio, 0.90, na.rm = TRUE)
joined$label_this <- joined$consensus_outlier | joined$B_M_ratio >= bm_cut

write_csv(joined,
          "02-03_Sam's_Results/04_Figures/F02_blood_contamination/c_data/panel_C_join_data.csv")

sp <- suppressWarnings(cor.test(joined$B_M_ratio, joined$mahal_dist,
                                method = "spearman", exact = FALSE))
pe <- cor.test(joined$B_M_ratio, joined$mahal_dist, method = "pearson")
stats_label <- sprintf("Spearman ρ = %.2f (p = %.2g)\nPearson r = %.2f (p = %.2g)\nn = %d",
                      sp$estimate, sp$p.value, pe$estimate, pe$p.value, nrow(joined))

p_C <- ggplot(joined, aes(B_M_ratio, mahal_dist,
                          color = consensus_outlier,
                          shape = consensus_outlier)) +
  geom_point(size = 2.5) +
  geom_text_repel(data = filter(joined, label_this),
                  aes(label = sample_id), size = 2.6, max.overlaps = Inf,
                  box.padding = 0.5) +
  scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey40"),
                     name = "Our 4-method\nconsensus_outlier") +
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16),
                     name = "Our 4-method\nconsensus_outlier") +
  annotate("text", x = -Inf, y = Inf, label = stats_label,
           hjust = -0.1, vjust = 1.2, size = 3) +
  labs(x = "Sam's B_M_ratio",
       y = "Our Mahalanobis distance (PC1–2)",
       title = "Cross-pipeline contamination agreement") +
  theme_minimal(base_size = 10)

ggsave("02-03_Sam's_Results/04_Figures/F02_blood_contamination/b_reports/main/png/panels/panel_C.png",
       p_C, width = 5.5, height = 4, dpi = 300, bg = "white")

panel_C <- p_C
