# Figure 2 (CRvH) — Panel A: CV% Violins (Inter-Individual Variability)
# Faceted by CR pooled, CRE, PLA, Healthy. T1/T2 x-axis for CR groups.
# Median labels with bootstrap CIs, Wilcoxon + Cliff's delta for paired T1/T2.
# Outputs: pA (ggplot object), cv_violin.pdf/.png

setwd(here::here())
source("04_Figures/shared/style.R")

pacman::p_load(dplyr, tidyr, stringr, readr, ggplot2, ggbeeswarm)

PA_W <- 210; PA_H <- 120

RPT_DIR <- "04_Figures/F02_Proteome_Overview/b_reports/supp/panels"
DAT_DIR <- "04_Figures/F02_Proteome_Overview/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

# Load data & metadata
# Normalized (non-imputed) matrix from the proteoDA DAList.
.dal <- readRDS("02_Normalization/c_data/DAList_normalized.rds")
norm_df <- tibble::as_tibble(cbind(
  .dal$annotation[, c("uniprot_id", "protein", "gene", "description")],
  as.data.frame(.dal$data)
))
meta_full <- read_csv("00_input/CvH_meta.csv", show_col_types = FALSE)

ann_cols   <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(norm_df), ann_cols)

# CRvH sample space: CR_CRE, CR_PLA, PPS
meta <- meta_full |>
  filter(Group %in% c("CR_CRE", "CR_PLA", "PPS"), Col_ID %in% samp_names)

samp_crvh <- meta$Col_ID

pdf_device <- get_pdf_device()

# CV on linear (not log) scale per Brenes 2024
lin_mat <- 2^as.matrix(norm_df[, samp_crvh])

# Group_Time levels: CRE_T1, CRE_T2, PLA_T1, PLA_T2, H_T1
group_levels <- c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1")

cv_list <- lapply(group_levels, function(g) {
  idx <- meta$Col_ID[meta$Group_Time == g]
  sub <- lin_mat[, idx, drop = FALSE]
  cv_pct <- apply(sub, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) return(NA_real_)
    sd(x) / mean(x) * 100
  })
  tibble(group_time = g, cv = cv_pct)
})
cv_df <- bind_rows(cv_list) |> filter(!is.na(cv))
cv_df$group_time <- factor(cv_df$group_time, levels = group_levels)

# Facet: CR pooled (CRE+PLA), CRE, PLA, Healthy
cv_df <- cv_df |>
  mutate(
    timepoint = ifelse(grepl("T2$", group_time), "T2", "T1"),
    facet = case_when(
      grepl("^CRE", group_time) ~ "CRE",
      grepl("^PLA", group_time) ~ "PLA",
      grepl("^H",   group_time) ~ "Healthy",
      TRUE ~ NA_character_
    )
  )
cv_df$timepoint <- factor(cv_df$timepoint, levels = c("T1", "T2"))
cv_df$facet     <- factor(cv_df$facet, levels = c("CRE", "PLA", "Healthy"))

# Also build CR pooled data
cv_pooled <- cv_df |>
  filter(facet %in% c("CRE", "PLA")) |>
  mutate(facet = "CR (pooled)")
cv_all <- bind_rows(cv_df, cv_pooled)
cv_all$facet <- factor(cv_all$facet,
                       levels = c("CR (pooled)", "CRE", "PLA", "Healthy"))

# Bootstrap 95% CI on median CV per group
set.seed(42)
boot_median_ci <- function(x, R = 2000, conf = 0.95) {
  meds <- replicate(R, median(sample(x, replace = TRUE)))
  qs   <- quantile(meds, c((1 - conf) / 2, (1 + conf) / 2))
  c(lower = unname(qs[1]), upper = unname(qs[2]))
}

# Wilcoxon tests: T1 vs T2 within each facet (paired where possible)
cliffs_delta <- function(x, y) {
  nx <- length(x); ny <- length(y)
  d <- outer(x, y, function(a, b) sign(a - b))
  sum(d) / (nx * ny)
}

bracket_facets <- c("CR (pooled)", "CRE", "PLA")
wilcox_results <- lapply(bracket_facets, function(f) {
  t1 <- cv_all$cv[cv_all$facet == f & cv_all$timepoint == "T1"]
  t2 <- cv_all$cv[cv_all$facet == f & cv_all$timepoint == "T2"]
  wt <- wilcox.test(t1, t2)
  data.frame(
    facet    = f,
    p_raw    = wt$p.value,
    cliffs_d = cliffs_delta(t1, t2),
    stringsAsFactors = FALSE
  )
})
wilcox_df <- bind_rows(wilcox_results)
wilcox_df$p_bh <- p.adjust(wilcox_df$p_raw, method = "BH")

# Summary stats per facet x timepoint
cv_ci <- cv_all |>
  group_by(facet, timepoint) |>
  summarise(
    med    = median(cv),
    ci_lo  = boot_median_ci(cv)[["lower"]],
    ci_hi  = boot_median_ci(cv)[["upper"]],
    cv_max = max(cv),
    .groups = "drop"
  )

n_prot     <- nrow(norm_df)
grand_med  <- median(cv_df$cv)
grand_ci   <- boot_median_ci(cv_df$cv)

# Delta median CV per facet (T2 - T1) for facets with paired timepoints
delta_cv <- cv_ci |>
  filter(facet != "Healthy") |>
  select(facet, timepoint, med) |>
  pivot_wider(names_from = timepoint, values_from = med) |>
  mutate(delta = T2 - T1,
         arrow_label = sprintf("%+.1f%%", delta))

arrow_df <- delta_cv |>
  mutate(x = 1, xend = 2,
         y_mid = (T1 + T2) / 2)

sub_txt <- sprintf(
  "%s proteins | CV %.0f%% [%.0f\u2013%.0f] | CR %+.1f%%, CRE %+.1f%%, PLA %+.1f%%",
  format(n_prot, big.mark = ","), grand_med, grand_ci[1], grand_ci[2],
  delta_cv$delta[delta_cv$facet == "CR (pooled)"],
  delta_cv$delta[delta_cv$facet == "CRE"],
  delta_cv$delta[delta_cv$facet == "PLA"]
)

# Plot
pA <- ggplot(cv_all, aes(x = timepoint, y = cv, fill = group_time)) +
  geom_violin(alpha = 0.5, linewidth = 0.3, color = "black", scale = "width") +
  geom_quasirandom(aes(color = group_time), alpha = 0.15, size = 0.5,
                   width = 0.25, groupOnX = TRUE, show.legend = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.3,
               color = "black", fill = "white", coef = 0) +
  geom_hline(yintercept = 25, linetype = "dashed", color = "grey50",
             linewidth = 0.4) +
  # Median labels with bootstrap CI
  geom_label(data = cv_ci,
             aes(x = timepoint, y = cv_max + 3,
                 label = sprintf("%.0f%% [%.0f\u2013%.0f]", med, ci_lo, ci_hi)),
             size = scale_text(BASE_COUNT + 0.5, PA_W / 4),
             fontface = "bold", fill = alpha("white", 0.8),
             linewidth = 0.2, label.padding = unit(1.5, "pt"),
             hjust = 0, nudge_x = -0.4) +
  # Directional arrows between T1/T2 medians
  geom_segment(data = arrow_df,
               aes(x = x, xend = xend, y = T1, yend = T2),
               inherit.aes = FALSE, color = "grey30",
               arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
               linewidth = 0.6) +
  geom_label(data = arrow_df,
             aes(x = 1.5, y = y_mid, label = arrow_label),
             inherit.aes = FALSE, size = scale_text(BASE_COUNT, PA_W / 4),
             fontface = "bold.italic", fill = alpha("white", 0.85),
             label.padding = unit(1.5, "pt"), linewidth = 0.2,
             color = "grey30") +
  facet_wrap(~ facet, nrow = 1) +
  scale_fill_manual(values = GROUP_FILL) +
  scale_color_manual(values = GROUP_FILL) +
  coord_cartesian(ylim = c(0, max(cv_ci$cv_max) + 12)) +
  labs(title = "Inter-Individual Variability (CV%)",
       subtitle = sub_txt,
       x = NULL, y = "CV (%)",
       tag = "A") +
  FIG_THEME + theme(legend.position = "none",
                    panel.spacing = unit(8, "mm"))

# Save
write.csv(as.data.frame(cv_ci),
          file.path(DAT_DIR, "audit_cv_violin_median_ci.csv"), row.names = FALSE)
write.csv(wilcox_df,
          file.path(DAT_DIR, "audit_cv_violin_wilcoxon.csv"), row.names = FALSE)

ggsave(file.path(RPT_DIR, "cv_violin.pdf"), pA,
       width = PA_W, height = PA_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "cv_violin.png"), pA,
       width = PA_W, height = PA_H, units = "mm", dpi = 300)
