# Figure 1 — Panel A: CV% Violins (Inter-Individual Variability)
# 4 facets: Cancer Recovery (pooled CRE+PLA, T1/T2), Creatine (T1/T2),
#           Placebo (T1/T2), Healthy (T1 only).
# Wilcoxon + Cliff's delta for CR, CRE, and PLA paired comparisons.
# Outputs: pA (ggplot object), panel_A_cv.pdf/.png, audit CSVs.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggbeeswarm)
})

PA_W <- 210; PA_H <- 120

RPT_DIR <- "04_Figures/F01/b_reports"
DAT_DIR <- "04_Figures/F01/c_data"
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Load data and metadata ---
norm_df <- read_csv("01_normalization/c_data/02_normalized.csv",
                    show_col_types = FALSE)
ann_cols   <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(norm_df), ann_cols)

meta <- read_csv("00_input/CvH_meta.csv", show_col_types = FALSE) |>
  filter(Col_ID %in% samp_names)

meta$group_time <- factor(meta$Group_Time,
                          levels = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1"))

# Build CR pooled metadata rows (CRE+PLA pooled by timepoint)
meta_cr <- meta |>
  filter(Supplement %in% c("CRE", "PLA")) |>
  mutate(group_time = factor(
    ifelse(Timepoint == "T1", "CR_T1", "CR_T2"),
    levels = c("CR_T1", "CR_T2")))

# Facet variable: supplement group label
meta$facet <- factor(
  case_when(
    meta$Supplement == "CRE" ~ "Creatine",
    meta$Supplement == "PLA" ~ "Placebo",
    TRUE                     ~ "Healthy"
  ),
  levels = c("Cancer Recovery", "Creatine", "Placebo", "Healthy")
)
meta_cr$facet <- factor("Cancer Recovery",
  levels = c("Cancer Recovery", "Creatine", "Placebo", "Healthy"))

meta$timepoint <- factor(meta$Timepoint, levels = c("T1", "T2"))

pdf_device <- get_pdf_device()

# --- CV% on linear scale per Brenes 2024 ---
lin_mat <- 2^as.matrix(norm_df[, samp_names])

# CR pooled fill colors
CR_FILL <- c(CR_T1 = scales::alpha("#8B5E3C", 0.7),
             CR_T2 = scales::alpha("#C49A6C", 0.7))
ALL_FILL <- c(CR_FILL, GROUP_FILL)

compute_cv_group <- function(mat, sample_ids) {
  sub <- mat[, sample_ids, drop = FALSE]
  apply(sub, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) return(NA_real_)
    sd(x) / mean(x) * 100
  })
}

# Compute CV for the 5 original groups + 2 CR pooled groups
all_groups <- c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1")
cv_list <- lapply(all_groups, function(g) {
  idx <- meta$Col_ID[meta$group_time == g]
  tibble(group_time = g, cv = compute_cv_group(lin_mat, idx))
})

# CR pooled (CRE+PLA at each timepoint)
cv_list_cr <- lapply(c("CR_T1", "CR_T2"), function(g) {
  idx <- meta_cr$Col_ID[meta_cr$group_time == g]
  tibble(group_time = g, cv = compute_cv_group(lin_mat, idx))
})

cv_df <- bind_rows(c(cv_list_cr, cv_list)) |> filter(!is.na(cv))
cv_df$group_time <- factor(cv_df$group_time,
  levels = c("CR_T1", "CR_T2", "CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1"))

# Add facet and timepoint columns
facet_levels <- c("Cancer Recovery", "Creatine", "Placebo", "Healthy")
cv_df <- cv_df |>
  mutate(
    facet = factor(
      case_when(
        grepl("^CR_", group_time)  ~ "Cancer Recovery",
        grepl("CRE", group_time)   ~ "Creatine",
        grepl("PLA", group_time)   ~ "Placebo",
        TRUE                       ~ "Healthy"
      ),
      levels = facet_levels
    ),
    timepoint = factor(
      ifelse(grepl("T2", group_time), "T2", "T1"),
      levels = c("T1", "T2")
    )
  )

# --- Bootstrap 95% CI on median CV ---
set.seed(42)
boot_median_ci <- function(x, R = 2000, conf = 0.95) {
  meds <- replicate(R, median(sample(x, replace = TRUE)))
  qs   <- quantile(meds, c((1 - conf) / 2, (1 + conf) / 2))
  c(lower = unname(qs[1]), upper = unname(qs[2]))
}

cv_ci <- cv_df |>
  group_by(facet, timepoint, group_time) |>
  summarise(
    med    = median(cv),
    ci_lo  = boot_median_ci(cv)[["lower"]],
    ci_hi  = boot_median_ci(cv)[["upper"]],
    cv_max = max(cv),
    .groups = "drop"
  )

# --- Wilcoxon tests: CR, CRE, PLA T1 vs T2 (BH-corrected) ---
bracket_comps <- list(c("CR_T1", "CR_T2"), c("CRE_T1", "CRE_T2"), c("PLA_T1", "PLA_T2"))
bracket_pvals_raw <- sapply(bracket_comps, function(pair)
  wilcox.test(cv_df$cv[cv_df$group_time == pair[1]],
              cv_df$cv[cv_df$group_time == pair[2]])$p.value)
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
    cliffs_delta(cv_df$cv[cv_df$group_time == pair[1]],
                 cv_df$cv[cv_df$group_time == pair[2]]))
)

# --- Plot ---
# Drop unused factor levels so Healthy facet only shows T1
cv_df$timepoint <- droplevels(cv_df$timepoint)

pA <- ggplot(cv_df, aes(x = timepoint, y = cv, fill = group_time)) +
  geom_violin(alpha = 0.5, linewidth = 0.3, color = "black", scale = "width") +
  geom_quasirandom(aes(color = group_time), alpha = 0.15, size = 0.4,
                   width = 0.25, groupOnX = TRUE, show.legend = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.3,
               color = "black", fill = "white", coef = 0) +
  geom_hline(yintercept = 25, linetype = "dashed", color = "grey50",
             linewidth = 0.4) +
  geom_label(data = cv_ci,
             aes(x = timepoint, y = cv_max + 5,
                 label = sprintf("%.0f%% [%.0f-%.0f]", med, ci_lo, ci_hi)),
             size = scale_text(BASE_STAT - 0.5, PA_W),
             fontface = "bold", fill = alpha("white", 0.8),
             linewidth = 0.2, label.padding = unit(1.5, "pt"),
             hjust = 0.5) +
  facet_wrap(~ facet, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = ALL_FILL) +
  scale_color_manual(values = ALL_FILL) +
  coord_cartesian(ylim = c(0, max(cv_ci$cv_max) + max(cv_ci$cv_max) * 0.15)) +
  labs(title = "Inter-Individual Variability (CV%)",
       subtitle = "Protein-level CV% (cycloess-normalized)",
       x = NULL, y = "CV (%)",
       tag = "A") +
  FIG_THEME + theme(legend.position = "none",
                    panel.spacing = unit(8, "mm"))

# --- Save outputs ---
write.csv(as.data.frame(cv_ci),
          file.path(DAT_DIR, "panel_A_median_cv_ci.csv"), row.names = FALSE)
write.csv(cliff_results,
          file.path(DAT_DIR, "panel_A_wilcoxon_effects.csv"), row.names = FALSE)

ggsave(file.path(RPT_DIR, "panel_A_cv.pdf"), pA,
       width = PA_W, height = PA_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_DIR, "panel_A_cv.png"), pA,
       width = PA_W, height = PA_H, units = "mm", dpi = 300)

cat("Panel A done.\n")
