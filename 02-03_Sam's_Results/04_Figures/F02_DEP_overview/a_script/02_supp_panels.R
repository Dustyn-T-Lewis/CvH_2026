#!/usr/bin/env Rscript
# F02_DEP_overview Supp — SA: CV scatter  SB: CV violin  SC: per-subject variability
#
# Sam-parallel of YvO F02 supp. Adapted for CvH 5-group design:
#   CRE_T1 / CRE_T2 / PLA_T1 / PLA_T2 / H_T1
#
# SA: Per-protein CV% scatter: group-pair scatter (T1 vs T2 per arm)
# SB: Per-group CV% violin: inter-individual variability across all 5 groups
# SC: Per-subject intra-individual variability (CRE/PLA arms: T1→T2 logFC)
#
# Run from A_CvH_2026/ (rprojroot resolves to A_CvH_2026/).

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(ggbeeswarm)
  library(patchwork)
  library(cowplot)
  library(purrr)
})

source("02-03_Sam's_Results/04_Figures/shared/style.R")

# cairo_pdf requires X11 (unavailable on this system); fall back to base pdf().
get_pdf_device <- function() grDevices::pdf

# ── Local helpers ──────────────────────────────────────────────────────────
strip_for_composite <- function(p) {
  p + labs(title = NULL, subtitle = NULL, tag = NULL) +
    theme(legend.position = "none")
}

composite_text_sizes <- function(comp_h_mm) {
  list(
    title    = pmax(6, pmin(8, round(5 + comp_h_mm / 80))),
    subtitle = pmax(4, pmin(6, round(3 + comp_h_mm / 100))),
    tag      = 8
  )
}

boot_median_ci <- function(x, R = 2000, conf = 0.95) {
  meds <- replicate(R, median(sample(x, replace = TRUE)))
  qs   <- quantile(meds, c((1 - conf) / 2, (1 + conf) / 2))
  c(lower = unname(qs[1]), upper = unname(qs[2]))
}

# ── Paths ────────────────────────────────────────────────────────────────────
SAM_ROOT <- "02-03_Sam's_Results"
BASE     <- file.path(SAM_ROOT, "04_Figures", "F02_DEP_overview")
RPT_PNG  <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF  <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
DAT      <- file.path(BASE, "c_data")
for (d in c(RPT_PNG, RPT_PDF, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

pdf_dev <- get_pdf_device()

HEATMAP_LO <- "#2166AC"; HEATMAP_HI <- "#B2182B"

# ── Load normalized DAList ────────────────────────────────────────────────────
dal_norm <- readRDS(file.path(SAM_ROOT, "00_input",
                               "01_normalized_DAList_SURV_stringent_muscle.RDS"))
norm_mat  <- as.matrix(dal_norm$data)
# Metadata: sample_id, supp, cancer (SURV/CTL), supp_time, timepoint
norm_meta <- as_tibble(dal_norm$metadata) |>
  mutate(
    group_label = case_when(
      cancer == "CTL" ~ "CTL",
      supp_time != "" ~ supp_time,
      TRUE            ~ paste0(supp, "_", timepoint)
    ),
    group     = factor(group_label,
                       levels = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "CTL")),
    arm       = factor(if_else(cancer == "CTL", "H",
                               if_else(supp != "", supp, "Unknown")),
                       levels = c("CRE", "PLA", "H")),
    subject   = sub("_T[12]$", "", sample_id)
  )

ann_cols  <- c("uniprot_id", "protein", "gene", "description")
ann_df    <- as_tibble(dal_norm$annotation) |>
  select(any_of(ann_cols))
norm_df   <- bind_cols(ann_df, as_tibble(norm_mat))
samp_names <- norm_meta$sample_id

# ── Panel SA: CV% Scatter (T1 vs T2 per arm) ─────────────────────────────────
PA_W <- 178; PA_H <- 70
PA_SUB <- 55

compute_cv <- function(mat, idx) {
  sub <- mat[, idx, drop = FALSE]
  apply(sub, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) return(NA_real_)
    sd(x) / mean(x) * 100
  })
}

lin_mat <- 2^norm_mat

# Pairs: CRE T1 vs T2, PLA T1 vs T2, H vs CRE_T1 (cross-group comparison)
pair_defs <- list(
  list(label = "CRE: T1 vs T2",  g1 = "CRE_T1", g2 = "CRE_T2"),
  list(label = "PLA: T1 vs T2",  g1 = "PLA_T1", g2 = "PLA_T2")
)

scatter_df <- bind_rows(lapply(pair_defs, function(pd) {
  idx1 <- norm_meta$sample_id[norm_meta$group == pd$g1]
  idx2 <- norm_meta$sample_id[norm_meta$group == pd$g2]
  cv1  <- compute_cv(lin_mat, idx1)
  cv2  <- compute_cv(lin_mat, idx2)
  tibble(gene    = ann_df$gene,
         cv_g1   = cv1,
         cv_g2   = cv2,
         pair    = pd$label,
         g1_name = pd$g1,
         g2_name = pd$g2)
})) |>
  filter(!is.na(cv_g1), !is.na(cv_g2)) |>
  mutate(
    max_cv        = pmax(cv_g1, cv_g2),
    pair          = factor(pair, levels = sapply(pair_defs, `[[`, "label"))
  )

max_cv_cap         <- quantile(scatter_df$max_cv, 0.98, na.rm = TRUE)
scatter_df$max_cv_capped <- pmin(scatter_df$max_cv, max_cv_cap)

top_cv_labels <- scatter_df |>
  slice_max(max_cv, n = 12, with_ties = FALSE, by = pair)

r_annots <- scatter_df |>
  group_by(pair) |>
  summarise(
    r_val = cor(cv_g1, cv_g2, use = "complete.obs"),
    n_val = sum(!is.na(cv_g1) & !is.na(cv_g2)),
    .groups = "drop"
  ) |>
  mutate(label = sprintf("r = %.2f (n = %d)", r_val, n_val))

axis_max_cv <- 250

pSA <- ggplot(scatter_df, aes(x = cv_g1, y = cv_g2)) +
  facet_wrap(~pair, nrow = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey50", linewidth = 0.4) +
  geom_point(aes(color = max_cv_capped), alpha = 0.35, size = 0.9) +
  geom_label_repel(data = top_cv_labels,
                   aes(label = gene, fill = max_cv_capped),
                   color = "white", fontface = "bold",
                   size = scale_text(BASE_GENE, PA_SUB),
                   label.padding = unit(1, "pt"), label.size = 0.3,
                   max.overlaps = 20, segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, seed = 42, show.legend = FALSE) +
  geom_label(data = r_annots, aes(label = label),
             x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4,
             size = scale_text(BASE_STAT, PA_SUB),
             color = "grey30", fontface = "bold",
             fill = scales::alpha("white", 0.85), linewidth = 0,
             label.padding = unit(2, "pt"), inherit.aes = FALSE) +
  scale_color_viridis_c(option = "inferno", direction = -1, begin = 0.1, end = 0.85,
                        name = "CV%",
                        guide = guide_colorbar(barwidth = unit(2, "mm"),
                                               barheight = unit(12, "mm"),
                                               title.position = "top", title.hjust = 0.5)) +
  scale_fill_viridis_c(option = "inferno", direction = -1, begin = 0.1, end = 0.85,
                       name = "CV%", guide = "none") +
  coord_fixed(xlim = c(0, axis_max_cv), ylim = c(0, axis_max_cv)) +
  labs(title    = "Per-Protein Variability (CV%)",
       subtitle = sprintf("%s proteins | CRE arm r = %.2f | PLA arm r = %.2f",
                          format(nrow(ann_df), big.mark = ","),
                          r_annots$r_val[r_annots$pair == "CRE: T1 vs T2"],
                          r_annots$r_val[r_annots$pair == "PLA: T1 vs T2"]),
       x = expression(bold(CV * "%" [T1])),
       y = expression(bold(CV * "%" [T2])),
       tag = "a") +
  FIG_THEME +
  theme(legend.position      = c(0.97, 0.02), legend.justification = c(1, 0),
        legend.background    = element_rect(fill = scales::alpha("white", 0.8), color = NA),
        strip.text           = element_text(face = "bold", size = FIG_STRIP_SIZE),
        plot.margin          = margin(t = 0, r = 5.5, b = 0, l = 5.5))

write.csv(scatter_df |> select(gene, cv_g1, cv_g2, max_cv, pair),
          file.path(DAT, "SUPP_panel_SA_cv_scatter.csv"), row.names = FALSE)
ggsave(file.path(RPT_PNG, "SUPP_panel_SA_cv_scatter.png"), pSA,
       width = PA_W, height = PA_H, units = "mm", dpi = 300)

pSA_title    <- "Per-Protein Variability (CV%)"
pSA_subtitle <- sprintf("%s proteins (cycloess-normalized) | T1–T2 per arm",
                         format(nrow(ann_df), big.mark = ","))
pSA <- strip_for_composite(pSA)

# ── Panel SB: CV% Violin (Inter-Individual Variability) ──────────────────────
PB_W <- 110; PB_H <- 80

cv_list <- lapply(levels(norm_meta$group), function(g) {
  idx    <- norm_meta$sample_id[norm_meta$group == g]
  if (length(idx) == 0) return(NULL)
  sub    <- lin_mat[, idx, drop = FALSE]
  cv_pct <- apply(sub, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 2) return(NA_real_)
    sd(x) / mean(x) * 100
  })
  tibble(group = g, cv = cv_pct)
})
cv_df       <- bind_rows(cv_list) |> filter(!is.na(cv))
cv_df$group <- factor(cv_df$group,
                      levels = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "CTL"))
cv_df$arm   <- factor(
  dplyr::case_when(
    grepl("^CRE", cv_df$group) ~ "CRE",
    grepl("^PLA", cv_df$group) ~ "PLA",
    TRUE                       ~ "H"
  ),
  levels = c("CRE", "PLA", "H")
)

set.seed(42)
cv_ci <- cv_df |>
  group_by(group) |>
  summarise(
    med    = median(cv),
    ci_lo  = boot_median_ci(cv)[["lower"]],
    ci_hi  = boot_median_ci(cv)[["upper"]],
    cv_max = max(cv),
    .groups = "drop"
  )

# Delta median for CRE and PLA arms (T1 -> T2)
delta_cv <- cv_ci |>
  filter(as.character(group) %in% c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2")) |>
  mutate(arm = sub("_T[12]$", "", as.character(group)),
         tp  = sub("^.*_", "", as.character(group))) |>
  select(arm, tp, med) |>
  pivot_wider(names_from = tp, values_from = med) |>
  mutate(delta = T2 - T1, arrow_label = sprintf("%+.1f%%", delta))

arm_xpos   <- c(CRE = 1, PLA = 3)
arm_xpos_e <- c(CRE = 2, PLA = 4)
arrow_df <- delta_cv |>
  filter(arm %in% c("CRE", "PLA")) |>
  mutate(
    x     = arm_xpos[arm],
    xend  = arm_xpos_e[arm],
    Pre   = T1, Post = T2,
    y_mid = (T1 + T2) / 2
  )

n_prot    <- nrow(norm_df)
grand_med <- median(cv_df$cv)
set.seed(42)
grand_ci  <- boot_median_ci(cv_df$cv)

GROUP_LABELS <- c(CRE_T1 = "Pre", CRE_T2 = "Post",
                  PLA_T1 = "Pre", PLA_T2 = "Post", CTL = "Healthy")
ARM_LABELS   <- c("CRE", "PLA", "Healthy")
ARM_X_POS    <- c(1.5, 3.5, 5)
# Extend GROUP_FILL to include CTL
GROUP_FILL_F2 <- c(GROUP_FILL, CTL = scales::alpha("#4DAF4A", 0.7))

sub_txt <- sprintf(
  "Inter-individual CV%% by group | %s proteins | grand median %.0f%% [%.0f–%.0f]",
  format(n_prot, big.mark = ","), grand_med, grand_ci[1], grand_ci[2]
)

pSB <- ggplot(cv_df, aes(x = group, y = cv, fill = group)) +
  annotate("text", x = ARM_X_POS[1], y = Inf, label = "CRE", vjust = 1.3,
           size = 2.0, fontface = "bold", color = "grey25") +
  annotate("text", x = ARM_X_POS[2], y = Inf, label = "PLA", vjust = 1.3,
           size = 2.0, fontface = "bold", color = "grey25") +
  annotate("text", x = ARM_X_POS[3], y = Inf, label = "Healthy", vjust = 1.3,
           size = 2.0, fontface = "bold", color = "grey25") +
  geom_violin(alpha = 0.5, linewidth = 0.3, color = "black", scale = "width") +
  geom_quasirandom(aes(color = group), alpha = 0.15, size = 0.5,
                   width = 0.25, groupOnX = TRUE, show.legend = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.3,
               color = "black", fill = "white", coef = 0) +
  geom_hline(yintercept = 25, linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  geom_label(data = cv_ci,
             aes(x = group, y = cv_max + 3,
                 label = sprintf("%.0f%% [%.0f–%.0f]", med, ci_lo, ci_hi)),
             size = scale_text(BASE_COUNT - 1.5, PB_W),
             fontface = "bold", fill = scales::alpha("white", 0.8),
             linewidth = 0.2, label.padding = unit(1.0, "pt"), hjust = 0.5) +
  {if (nrow(arrow_df) > 0)
    geom_segment(data = arrow_df,
                 aes(x = x, xend = xend, y = Pre, yend = Post),
                 inherit.aes = FALSE, color = "grey30",
                 arrow = arrow(length = unit(1.5, "mm"), type = "closed"),
                 linewidth = 0.6)} +
  {if (nrow(arrow_df) > 0)
    geom_label(data = arrow_df,
               aes(x = (x + xend) / 2, y = y_mid, label = arrow_label),
               inherit.aes = FALSE, size = scale_text(BASE_COUNT - 1.5, PB_W),
               fontface = "bold.italic", fill = scales::alpha("white", 0.85),
               label.padding = unit(1.5, "pt"), linewidth = 0.2, color = "grey30")} +
  scale_fill_manual(values  = GROUP_FILL_F2) +
  scale_color_manual(values = GROUP_FILL_F2) +
  scale_x_discrete(labels = GROUP_LABELS) +
  coord_cartesian(ylim = c(0, max(cv_ci$cv_max) + 20)) +
  labs(title    = "Inter-Individual Variability (CV%)",
       subtitle = sub_txt,
       x = NULL, y = "CV (%)", tag = "b") +
  FIG_THEME +
  theme(legend.position = "none",
        plot.title      = element_text(margin = margin(b = 0)),
        plot.subtitle   = element_text(size = FIG_SUBTITLE_SIZE - 1.0,
                                       face = "bold.italic", color = "grey40",
                                       margin = margin(t = 0, b = 1)),
        plot.margin     = margin(t = 0, r = 5.5, b = 5.5, l = 5.5))

write.csv(as.data.frame(cv_ci),
          file.path(DAT, "SUPP_panel_SB_median_cv_ci.csv"), row.names = FALSE)
ggsave(file.path(RPT_PNG, "SUPP_panel_SB_cv_violin.png"), pSB,
       width = PB_W, height = PB_H, units = "mm", dpi = 300)

pSB_title    <- "Inter-Individual Variability (CV%)"
pSB_subtitle <- sub_txt
pSB <- strip_for_composite(pSB)

# ── Panel SC: Per-Subject Intra-Individual Variability (CR arm only) ──────────
PE_W <- 110; PE_H <- 55

# Use normalized data as proxy for intra-individual variability (T1 -> T2 logFC)
subjects_cr <- norm_meta |>
  filter(arm %in% c("CRE", "PLA")) |>
  select(sample_id, subject, arm, timepoint, group) |>
  distinct()

lfc_list <- lapply(unique(subjects_cr$subject), function(s) {
  t1_id <- subjects_cr |> filter(subject == s, timepoint == "T1") |> pull(sample_id)
  t2_id <- subjects_cr |> filter(subject == s, timepoint == "T2") |> pull(sample_id)
  if (length(t1_id) != 1 || length(t2_id) != 1) return(NULL)
  lfc <- norm_mat[, t2_id] - norm_mat[, t1_id]
  arm <- subjects_cr |> filter(subject == s) |> pull(arm) |> unique()
  tibble(subject = s, arm = arm[1], lfc = as.numeric(lfc))
})
lfc_long     <- bind_rows(compact(lfc_list))
lfc_long$arm <- factor(lfc_long$arm, levels = c("CRE", "PLA"))

subj_summary <- lfc_long |>
  group_by(subject, arm) |>
  summarise(
    median_lfc = median(lfc, na.rm = TRUE),
    sd_lfc     = sd(lfc, na.rm = TRUE),
    iqr_lfc    = IQR(lfc, na.rm = TRUE),
    q25        = quantile(lfc, 0.25, na.rm = TRUE),
    q75        = quantile(lfc, 0.75, na.rm = TRUE),
    n_proteins = n(),
    .groups    = "drop"
  ) |>
  arrange(arm, median_lfc) |>
  mutate(subj_order = factor(subject, levels = unique(subject)))

lfc_long <- lfc_long |>
  mutate(subj_order = factor(subject, levels = levels(subj_summary$subj_order)))

wt   <- wilcox.test(median_lfc ~ arm, data = subj_summary)
n1   <- sum(subj_summary$arm == "CRE")
n2   <- sum(subj_summary$arm == "PLA")

n_proteins_sc <- nrow(ann_df)
subtitle_sc   <- sprintf(
  "Per-subject Delta log2FC (T2 - T1) | %s proteins | Wilcoxon CRE vs PLA %s",
  format(n_proteins_sc, big.mark = ","), fmt_p(wt$p.value)
)

arm_label_df <- data.frame(
  arm   = factor(c("CRE", "PLA"), levels = c("CRE", "PLA")),
  x_mid = c((n1 + 1) / 2, (n2 + 1) / 2),
  label = c("CRE", "PLA")
)

ARM_COLORS <- c(CRE = "#2166AC", PLA = "#D6604D")

pSC <- ggplot(lfc_long, aes(x = subj_order, y = lfc, fill = arm)) +
  geom_boxplot(width = 0.5, linewidth = 0.3, color = "black",
               outlier.shape = NA, alpha = 0.5) +
  geom_text(data = arm_label_df, aes(x = x_mid, y = Inf, label = label),
            inherit.aes = FALSE, vjust = 1.3, hjust = 0.5,
            size = FIG_STRIP_SIZE / .pt, fontface = "bold", color = "grey20") +
  facet_grid(~ arm, scales = "free_x", space = "free_x") +
  coord_cartesian(ylim = c(-1.5, 1.5), clip = "off") +
  scale_fill_manual(values = ARM_COLORS) +
  labs(x = "Subject",
       y = expression(bold(Delta~log[2]*"FC (T2/T1)")),
       title    = "Intra-Individual Proteomic Variability (T1→T2)",
       subtitle = subtitle_sc,
       tag = "c") +
  FIG_THEME +
  theme(legend.position = "none",
        panel.spacing   = unit(3, "mm"),
        axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                       size = FIG_AXIS_TEXT - 1.5),
        strip.text      = element_blank(),
        plot.margin     = margin(t = 0, r = 5.5, b = 5.5, l = 5.5))

write.csv(subj_summary |>
            select(subject, arm, median_lfc, sd_lfc, iqr_lfc, q25, q75, n_proteins),
          file.path(DAT, "SUPP_panel_SC_intra_variability.csv"),
          row.names = FALSE)
ggsave(file.path(RPT_PNG, "SUPP_panel_SC_intra_variability.png"), pSC,
       width = PE_W, height = PE_H, units = "mm", dpi = 300)

pSC_title    <- "Intra-Individual Proteomic Variability"
pSC_subtitle <- subtitle_sc
pSC <- strip_for_composite(pSC)

# ── Supp composite (SA full-width top; SB + SC side by side) ─────────────────
SUPP_PNG <- file.path(BASE, "b_reports", "supp", "png")
SUPP_PDF <- file.path(BASE, "b_reports", "supp", "pdf")
dir.create(SUPP_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_PDF, recursive = TRUE, showWarnings = FALSE)

COMP_W <- 178; COMP_H <- 115
txt    <- composite_text_sizes(COMP_H)

composite <- (wrap_elements(pSA) / (pSB | pSC)) +
  plot_layout(heights = c(1, 0.7))

Y_TOP <- 0.985; Y_BOT <- 0.500
composite <- ggdraw(composite) +
  draw_label("A",          x = 0.01, y = Y_TOP, size = txt$tag, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSA_title,    x = 0.04, y = Y_TOP, size = txt$title, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("B",          x = 0.01, y = Y_BOT, size = txt$tag, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSB_title,    x = 0.04, y = Y_BOT, size = txt$title, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("C",          x = 0.52, y = Y_BOT, size = txt$tag, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(pSC_title,    x = 0.55, y = Y_BOT, size = txt$title, fontface = "bold", hjust = 0, vjust = 1)

ggsave(file.path(SUPP_PDF, "SUPP_F02_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_dev)
ggsave(file.path(SUPP_PNG, "SUPP_F02_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message("F02 supp composite done")
