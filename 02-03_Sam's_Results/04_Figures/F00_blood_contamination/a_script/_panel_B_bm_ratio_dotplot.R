# Panel B — B_M_ratio dotplot by Group_Time, cohort top-decile (top 4 of 35)
# colored red and labeled with sample_id.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(ggrepel)
})

setwd(rprojroot::find_rstudio_root_file())
source("02-03_Sam's_Results/04_Figures/build_data_index.R")
source("02-03_Sam's_Results/04_Figures/shared/style.R")
source("02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_shared_keys.R")

sam <- readRDS(sam_idx$sam$dalist_rds)
meta <- as.data.frame(sam$metadata)
meta$sample_id <- normalize_sample_id(meta$sample_id)

# Group_Time derivation mirrors panel A: supp_time for SURV samples,
# cancer_time for CTL fallback.
meta$Group_Time <- ifelse(
  is.na(meta$supp) | meta$supp == "" | grepl("^NA_", meta$supp_time),
  as.character(meta$cancer_time),
  as.character(meta$supp_time)
)

bm_cut <- quantile(meta$B_M_ratio, 0.90, na.rm = TRUE)
gt_levels <- intersect(c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2",
                         "CTL_T1", "CTL_T2"), unique(meta$Group_Time))

panel_B_data <- meta |>
  transmute(sample_id,
            Group_Time = factor(Group_Time, levels = gt_levels),
            B_M_ratio,
            top_decile_flag = B_M_ratio >= bm_cut)

write_csv(panel_B_data,
          "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/panel_B_dotplot_data.csv")

kw_p <- kruskal.test(B_M_ratio ~ Group_Time, data = panel_B_data)$p.value
kw_label <- sprintf("Kruskal–Wallis p = %.3g", kw_p)

# Extend GROUP_COLORS with CTL_T1/CTL_T2 = the same green as H_T1.
group_palette <- c(GROUP_COLORS,
                   CTL_T1 = unname(GROUP_COLORS["H_T1"]),
                   CTL_T2 = unname(GROUP_COLORS["H_T1"]))

p_B <- ggplot(panel_B_data, aes(Group_Time, B_M_ratio)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.2, color = "grey50") +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "grey30") +
  geom_jitter(aes(color = Group_Time, alpha = top_decile_flag,
                  size = top_decile_flag), width = 0.15, height = 0) +
  geom_text_repel(data = filter(panel_B_data, top_decile_flag),
                  aes(label = sample_id), size = 2.6, max.overlaps = Inf,
                  box.padding = 0.5, color = "firebrick") +
  scale_color_manual(values = group_palette, guide = "none") +
  scale_alpha_manual(values = c(`FALSE` = 0.7, `TRUE` = 1), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 1.8, `TRUE` = 2.6), guide = "none") +
  annotate("text", x = 0.7, y = max(panel_B_data$B_M_ratio, na.rm = TRUE),
           label = kw_label, hjust = 0, size = 3) +
  labs(x = "Group × Timepoint",
       y = "B_M_ratio (blood-to-muscle signal ratio)",
       title = "B_M_ratio by Group_Time (top-decile labeled)") +
  theme_minimal(base_size = 10)

ggsave("02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/png/panels/panel_B.png",
       p_B, width = 5.5, height = 4, dpi = 300, bg = "white")

panel_B <- p_B
