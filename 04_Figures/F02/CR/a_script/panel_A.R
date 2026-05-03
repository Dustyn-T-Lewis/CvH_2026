# Figure 2 CR — Panel A: CV% Violins (Inter-Individual Variability)
# Faceted by Supplement (CRE | PLA), T1/T2 x-axis. Median labels with bootstrap CIs.
# Wilcoxon T1 vs T2 per supplement.
# Outputs: pA (ggplot object), panel_A_cv_SUPP.pdf/.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F02/a_script/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(ggbeeswarm)
})

PA_W <- 110; PA_H <- 120

RPT_DIR <- "04_Figures/F02/CR/b_reports"
DAT_DIR <- "04_Figures/F02/CR/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

norm_df <- read_csv("01_normalization/c_data/02_normalized.csv",
                    show_col_types = FALSE)
meta_full <- read_csv("00_input/CvH_meta.csv", show_col_types = FALSE)

# Filter to CR-only (CRE + PLA, no PPS/Healthy)
meta_cr <- meta_full |>
  filter(Group %in% c("CR_CRE", "CR_PLA")) |>
  mutate(
    Supplement = factor(Supplement, levels = c("CRE", "PLA")),
    Timepoint  = factor(Timepoint, levels = c("T1", "T2")),
    Group_Time = factor(Group_Time, levels = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2"))
  )

ann_cols   <- c("uniprot_id", "protein", "gene", "description")
cr_samples <- intersect(meta_cr$Col_ID, setdiff(names(norm_df), ann_cols))
meta_cr    <- meta_cr |> filter(Col_ID %in% cr_samples)

pdf_device <- get_pdf_device()

# CV on linear (not log) scale per Brenes 2024
lin_mat <- 2^as.matrix(norm_df[, cr_samples])

cv_list <- lapply(levels(meta_cr$Group_Time), function(g) {
  idx <- meta_cr$Col_ID[meta_cr$Group_Time == g]
  sub <- lin_mat[, idx, drop = FALSE]
  cv_pct <- apply(sub, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) return(NA_real_)
    sd(x) / mean(x) * 100
  })
  tibble(group = g, cv = cv_pct)
})
cv_df <- bind_rows(cv_list) |> filter(!is.na(cv))
cv_df$group <- factor(cv_df$group,
                      levels = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2"))

cv_df$supplement <- factor(ifelse(grepl("CRE", cv_df$group), "CRE", "PLA"),
                           levels = c("CRE", "PLA"))
cv_df$timepoint  <- factor(ifelse(grepl("T1", cv_df$group), "T1", "T2"),
                           levels = c("T1", "T2"))

# Bootstrap 95% CI on median CV per group
set.seed(42)
boot_median_ci <- function(x, R = 2000, conf = 0.95) {
  meds <- replicate(R, median(sample(x, replace = TRUE)))
  qs   <- quantile(meds, c((1 - conf) / 2, (1 + conf) / 2))
  c(lower = unname(qs[1]), upper = unname(qs[2]))
}

# Pairwise Wilcoxon tests: T1 vs T2 per supplement, BH corrected
bracket_comps <- list(c("CRE_T1", "CRE_T2"), c("PLA_T1", "PLA_T2"),
                      c("CRE_T1", "PLA_T1"), c("CRE_T2", "PLA_T2"))
bracket_pvals_raw <- sapply(bracket_comps, function(pair)
  wilcox.test(cv_df$cv[cv_df$group == pair[1]],
              cv_df$cv[cv_df$group == pair[2]])$p.value)
bracket_pvals <- p.adjust(bracket_pvals_raw, method = "BH")

cliffs_delta <- function(x, y) {
  nx <- length(x); ny <- length(y)
  d <- outer(x, y, function(a, b) sign(a - b))
  sum(d) / (nx * ny)
}
cliff_results <- data.frame(
  comparison = sapply(bracket_comps, paste, collapse = " vs "),
  p_raw      = bracket_pvals_raw,
  p_bh       = bracket_pvals,
  cliffs_d   = sapply(bracket_comps, function(pair)
    cliffs_delta(cv_df$cv[cv_df$group == pair[1]],
                 cv_df$cv[cv_df$group == pair[2]]))
)

cv_ci <- cv_df |>
  group_by(supplement, timepoint, group) |>
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

# Delta median CV per supplement (T2 - T1)
delta_cv <- cv_ci |>
  select(supplement, timepoint, med) |>
  pivot_wider(names_from = timepoint, values_from = med) |>
  mutate(delta = T2 - T1,
         arrow_label = sprintf("%+.1f%%", delta))

# Arrow annotation data
arrow_df <- delta_cv |>
  mutate(x = 1, xend = 2,
         y_mid = (T1 + T2) / 2)

sub_txt <- sprintf(
  "%s proteins | CV %.0f%% [%.0f\u2013%.0f] | CRE %+.1f%%, PLA %+.1f%%",
  format(n_prot, big.mark = ","), grand_med, grand_ci[1], grand_ci[2],
  delta_cv$delta[delta_cv$supplement == "CRE"],
  delta_cv$delta[delta_cv$supplement == "PLA"]
)

pA <- ggplot(cv_df, aes(x = timepoint, y = cv, fill = group)) +
  geom_violin(alpha = 0.5, linewidth = 0.3, color = "black", scale = "width") +
  geom_quasirandom(aes(color = group), alpha = 0.15, size = 0.5,
                   width = 0.25, groupOnX = TRUE, show.legend = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.3,
               color = "black", fill = "white", coef = 0) +
  geom_hline(yintercept = 25, linetype = "dashed", color = "grey50",
             linewidth = 0.4) +
  geom_label(data = cv_ci,
             aes(x = timepoint, y = cv_max + 3,
                 label = sprintf("%.0f%% [%.0f\u2013%.0f]", med, ci_lo, ci_hi)),
             size = scale_text(BASE_COUNT + 0.5, PA_W / 2),
             fontface = "bold", fill = alpha("white", 0.8),
             linewidth = 0.2, label.padding = unit(1.5, "pt"),
             hjust = 0, nudge_x = -0.4) +
  geom_segment(data = arrow_df,
               aes(x = x, xend = xend, y = T1, yend = T2),
               inherit.aes = FALSE, color = "grey30",
               arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
               linewidth = 0.6) +
  geom_label(data = arrow_df,
             aes(x = 1.5, y = y_mid, label = arrow_label),
             inherit.aes = FALSE, size = scale_text(BASE_COUNT, PA_W / 2),
             fontface = "bold.italic", fill = alpha("white", 0.85),
             label.padding = unit(1.5, "pt"), linewidth = 0.2,
             color = "grey30") +
  facet_wrap(~ supplement, nrow = 1) +
  scale_fill_manual(values = GROUP_FILL) +
  scale_color_manual(values = GROUP_FILL) +
  coord_cartesian(ylim = c(0, max(cv_ci$cv_max) + 12)) +
  labs(title = "Inter-Individual Variability (CV%)",
       subtitle = sub_txt,
       x = NULL, y = "CV (%)",
       tag = "A") +
  FIG_THEME + theme(legend.position = "none",
                    panel.spacing = unit(8, "mm"))

write.csv(as.data.frame(cv_ci),
          file.path(DAT_DIR, "audit_panel_A_median_cv_ci.csv"), row.names = FALSE)
write.csv(cliff_results,
          file.path(DAT_DIR, "audit_panel_A_wilcoxon_effects.csv"), row.names = FALSE)

ggsave(file.path(RPT_DIR, "panel_A_cv_SUPP.pdf"), pA,
       width = PA_W, height = PA_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_A_cv_SUPP.png"), pA,
       width = PA_W, height = PA_H, units = "mm", dpi = 300)
