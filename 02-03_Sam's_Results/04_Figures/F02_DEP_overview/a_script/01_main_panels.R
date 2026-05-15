#!/usr/bin/env Rscript
# F02_DEP_overview Main — Sam-parallel DEP overview (6-panel 3×2 grid)
# A: PCA biplot  B: logFC density  C: DEPs per contrast
# D: UpSet overlap  E: fGSEA pathways  F: Barcode rank
#
# Run from A_CvH_2026/ (rprojroot resolves to A_CvH_2026/).

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(cowplot)
  library(vegan)
  library(ComplexHeatmap)
  library(purrr)
})

source("02-03_Sam's_Results/04_Figures/shared/style.R")

# cairo_pdf requires X11 (unavailable on this system); fall back to base pdf().
get_pdf_device <- function() grDevices::pdf

# ── Local helpers (not in CvH style.R) ──────────────────────────────────────
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

# ── Paths ────────────────────────────────────────────────────────────────────
SAM_ROOT <- "02-03_Sam's_Results"
BASE     <- file.path(SAM_ROOT, "04_Figures", "F02_DEP_overview")
RPT_PNG  <- file.path(BASE, "b_reports", "main", "png")
RPT_PDF  <- file.path(BASE, "b_reports", "main", "pdf")
PNL_PNG  <- file.path(RPT_PNG, "panels")
PNL_PDF  <- file.path(RPT_PDF, "panels")
DAT      <- file.path(BASE, "c_data")
for (d in c(PNL_PNG, PNL_PDF, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

DEP_ROOT     <- file.path(SAM_ROOT, "03_DEP", "c_data")
PER_CONTRAST <- file.path(DEP_ROOT, "04_per_contrast_results")
FGSEA_CACHE  <- file.path(SAM_ROOT, "04_Figures", "shared", "fgsea_cache")

pdf_dev <- get_pdf_device()

# ── Load normalized DAList (for PCA) ─────────────────────────────────────────
dal_norm <- readRDS(file.path(SAM_ROOT, "00_input",
                               "01_normalized_DAList_SURV_stringent_muscle.RDS"))
norm_mat  <- as.matrix(dal_norm$data)
# Metadata: sample_id, supp (CRE/PLA/""), cancer (SURV/CTL), supp_time, timepoint
norm_meta <- as_tibble(dal_norm$metadata) |>
  mutate(
    # Derive display group label: use supp_time for CR participants, "CTL" for controls
    group_label = case_when(
      cancer == "CTL"  ~ "CTL",
      supp_time != ""  ~ supp_time,
      TRUE             ~ paste0(supp, "_", timepoint)
    ),
    group     = factor(group_label,
                       levels = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "CTL")),
    arm       = factor(cancer, levels = c("SURV", "CTL")),
    subject   = sub("_T[12]$", "", sample_id)
  )

# ── Load combined results + merge pi_score from per-contrast files ────────────
dep_crvh <- read_csv(file.path(DEP_ROOT, "03_combined_results_CRvH.csv"),
                     show_col_types = FALSE)
dep_cr   <- read_csv(file.path(DEP_ROOT, "03_combined_results_CR.csv"),
                     show_col_types = FALSE)

# Merge all contrasts into one wide table (annotation cols + stat cols)
ann_cols <- c("uniprot_id", "protein", "gene", "description")

# Per-contrast pi_score and sig_pi are in the per-contrast CSVs
load_pi <- function(contrast_name) {
  f <- file.path(PER_CONTRAST, paste0(contrast_name, ".csv"))
  read_csv(f, show_col_types = FALSE) |>
    select(uniprot_id, gene,
           !!paste0("pi_score_", contrast_name) := pi_score,
           !!paste0("sig_pi_",   contrast_name) := sig_pi)
}

pi_tabs <- lapply(c("Cancer_vs_Healthy", "Training_CR",
                    "Baseline_Supplement", "Training_CRE",
                    "Training_PLA", "Supplement_Interaction"),
                  load_pi)

# Build a single wide dep_df
dep_df <- dep_crvh |>
  select(all_of(ann_cols), starts_with("logFC_"), starts_with("t_"),
         starts_with("P.Value_"), starts_with("adj.P.Val_"),
         starts_with("sig.PVal_"), starts_with("sig.FDR_")) |>
  left_join(
    dep_cr |>
      select(-all_of(setdiff(ann_cols, "uniprot_id"))),
    by = "uniprot_id"
  ) |>
  left_join(reduce(pi_tabs, full_join, by = c("uniprot_id", "gene")),
            by = c("uniprot_id", "gene"))

all_genes <- unique(dep_df$gene[!is.na(dep_df$gene)])

# ── Contrast configuration ───────────────────────────────────────────────────
# 6 contrasts for panels B, E, F; 4 for UpSet (D)
CONTRASTS_ALL <- c("Cancer_vs_Healthy", "Training_CR",
                   "Baseline_Supplement", "Training_CRE",
                   "Training_PLA", "Supplement_Interaction")

CONTRASTS_UPSET <- c("Cancer_vs_Healthy", "Training_CR",
                     "Training_CRE", "Training_PLA")

SET_LABELS_UPSET <- c(
  Cancer_vs_Healthy = "CR vs H",
  Training_CR       = "Tr.(CR)",
  Training_CRE      = "Tr.(CRE)",
  Training_PLA      = "Tr.(PLA)"
)

# ── Panel A: PCA biplot + PERMANOVA ─────────────────────────────────────────
PC_W <- 67; PC_H <- 55

# Filter to complete cases (no NA) for PCA; log2 norm data may have missing values
complete_rows <- apply(norm_mat, 1, function(x) all(!is.na(x)))
pca_mat <- norm_mat[complete_rows, , drop = FALSE]
message(sprintf("PCA: using %d / %d proteins (complete cases)", nrow(pca_mat), nrow(norm_mat)))

pca     <- prcomp(t(pca_mat), center = TRUE, scale. = TRUE)
var_pct <- round(100 * summary(pca)$importance[2, 1:2], 1)

set.seed(42)
boot_var <- replicate(1000, {
  idx <- sample(nrow(pca_mat), replace = TRUE)
  100 * summary(prcomp(t(pca_mat[idx, ]), center = TRUE, scale. = TRUE))$importance[2, 1:2]
})
var_ci <- data.frame(PC = c("PC1", "PC2"), var_pct = var_pct,
                     ci_lo = apply(boot_var, 1, quantile, 0.025),
                     ci_hi = apply(boot_var, 1, quantile, 0.975))

pca_df <- as.data.frame(pca$x[, 1:2]) |>
  mutate(sample_id = rownames(pca$x)) |>
  left_join(norm_meta, by = "sample_id")

dist_mat <- dist(scale(t(pca_mat)))
set.seed(42)
# PERMANOVA: cancer (SURV/CTL) + timepoint (T1/T2) for SURV participants only
# Cannot use repeated-measures blocking across CTL (single timepoint) — use strata instead
perm_res <- adonis2(dist_mat ~ cancer + timepoint, data = norm_meta,
                    permutations = how(nperm = 999),
                    by = "terms")
perm_terms <- c("cancer", "timepoint")
perm_r2 <- perm_res[perm_terms, "R2"]
perm_pv <- perm_res[perm_terms, "Pr(>F)"]
perm_label <- sprintf(
  " PERMANOVA\nGroup     R² = %.3f,  %s\nTime      R² = %.3f,  %s",
  perm_r2[1], fmt_p(perm_pv[1]),
  perm_r2[2], fmt_p(perm_pv[2]))

PCA_LABEL_MAP <- c(
  CRE_T1 = "CRE Pre", CRE_T2 = "CRE Post",
  PLA_T1 = "PLA Pre", PLA_T2 = "PLA Post",
  CTL    = "Healthy"
)
# Extend PCA_COLORS for CTL
PCA_COLORS_F2 <- c(PCA_COLORS, CTL = "#4DAF4A")
PCA_SHAPES_F2 <- c(PCA_SHAPES, CTL = 15L)

pA <- ggplot(pca_df, aes(PC1, PC2, color = group, shape = group)) +
  stat_ellipse(aes(fill = group), geom = "polygon",
               alpha = 0.10, level = 0.80, show.legend = FALSE) +
  stat_ellipse(aes(group = group), level = 0.80, linewidth = 0.4,
               linetype = "dashed", show.legend = FALSE) +
  geom_point(size = 1.8, alpha = 0.85) +
  annotate("label", x = -Inf, y = Inf, label = perm_label,
           hjust = -0.02, vjust = 1.05, lineheight = 0.9,
           size = scale_text(BASE_COUNT, PC_W) - 0.5, color = "grey20", fontface = "bold",
           fill = scales::alpha("white", 0.85), linewidth = 0.2,
           label.padding = unit(0.12, "lines")) +
  scale_color_manual(values = PCA_COLORS_F2, labels = PCA_LABEL_MAP,
                     guide = guide_legend(override.aes = list(size = 1.6))) +
  scale_fill_manual(values  = PCA_COLORS_F2, guide = "none") +
  scale_shape_manual(values = PCA_SHAPES_F2, labels = PCA_LABEL_MAP) +
  labs(title = "Sample PCA",
       subtitle = sprintf("n = %d, %s proteins (normalized)",
                          nrow(norm_meta),
                          format(nrow(norm_mat), big.mark = ",")),
       x = sprintf("PC1 (%.1f%% [%.1f, %.1f])",
                   var_pct[1], var_ci$ci_lo[1], var_ci$ci_hi[1]),
       y = sprintf("PC2 (%.1f%% [%.1f, %.1f])",
                   var_pct[2], var_ci$ci_lo[2], var_ci$ci_hi[2]),
       tag = "a") +
  FIG_THEME +
  theme(plot.subtitle    = element_text(size = FIG_SUBTITLE_SIZE, face = "bold.italic",
                                         color = "grey30"),
        legend.position  = c(0.88, 0.15), legend.background = element_blank(),
        legend.key       = element_blank(), legend.title = element_blank(),
        legend.text      = element_text(size = FIG_LEGEND_TEXT + 0.5),
        legend.key.size  = unit(3, "mm"),
        plot.margin      = margin(6, 6, 2, 8))

write.csv(var_ci, file.path(DAT, "panel_A_pca_variance_ci.csv"), row.names = FALSE)
ggsave(file.path(PNL_PNG, "MAIN_panel_A_pca.png"), pA,
       width = PC_W, height = PC_H, units = "mm", dpi = 300)

pA_title    <- "Sample PCA"
pA_subtitle <- sprintf("n = %d, %s proteins (normalized)",
                        nrow(norm_meta), format(nrow(norm_mat), big.mark = ","))
pA <- strip_for_composite(pA)

# ── Panel B: logFC Density Histograms (4 contrasts) ─────────────────────────
PD_W <- 48; PD_H <- 55

# Use all 4 "main story" contrasts for density
CONTRASTS_DENS <- c("Cancer_vs_Healthy", "Training_CR",
                    "Training_CRE", "Training_PLA")

lfc_long_all <- dep_df |>
  select(any_of(c("uniprot_id", "gene")),
         starts_with("logFC_")) |>
  pivot_longer(starts_with("logFC_"), names_to = "contrast", values_to = "logFC") |>
  mutate(contrast = str_remove(contrast, "logFC_")) |>
  filter(!is.na(logFC))
write_csv(lfc_long_all, file.path(DAT, "panel_B_logfc_long.csv"))

lfc_long <- lfc_long_all |>
  filter(contrast %in% CONTRASTS_DENS) |>
  mutate(contrast = factor(contrast, levels = CONTRASTS_DENS))

lfc_stats <- lfc_long |>
  summarise(
    med_abs_lfc = median(abs(logFC)),
    n_above_05  = sum(abs(logFC) > 0.5),
    .by = contrast
  )
write_csv(lfc_stats, file.path(DAT, "panel_B_stats.csv"))

lfc_stats$annotation <- sprintf("Med.|logFC| = %.2f\nn(>0.5) = %d",
                                 lfc_stats$med_abs_lfc, lfc_stats$n_above_05)

lfc_binwidth <- 4 / 50
pB <- ggplot(lfc_long, aes(logFC, fill = contrast)) +
  geom_histogram(bins = 50, color = "black", linewidth = 0.2, alpha = 0.85) +
  geom_density(aes(y = after_stat(count) * lfc_binwidth),
               alpha = 0.15, linewidth = 0.5, color = "grey20") +
  geom_text(data = lfc_stats,
            aes(x = 0, y = 695, label = CTR_SHORT[as.character(contrast)]),
            inherit.aes = FALSE, hjust = 0.5, vjust = 1,
            size = scale_text(BASE_COUNT, PD_W) + 0.2,
            color = "grey30", fontface = "bold") +
  geom_label(data = lfc_stats,
             aes(x = -1.08, y = 695, label = annotation),
             inherit.aes = FALSE, hjust = 0, vjust = 1,
             size = scale_text(BASE_COUNT, PD_W) - 0.5,
             color = "grey20", fontface = "bold",
             lineheight = 0.9, fill = scales::alpha("white", 0.85),
             linewidth = 0.2, label.padding = unit(0.12, "lines")) +
  facet_wrap(~contrast, ncol = 1,
             labeller = labeller(contrast = CTR_SHORT)) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(0, 700)) +
  scale_fill_manual(values = CONTRAST_COLORS[CONTRASTS_DENS]) +
  labs(title = "Effect Size Distribution", subtitle = NULL,
       x = expression(bold(log[2]~FC)), y = " ", tag = "b") +
  FIG_THEME +
  theme(legend.position = "none", strip.text = element_blank(),
        strip.background = element_blank(),
        panel.spacing.y = unit(0, "pt"),
        axis.text.y  = element_text(size = FIG_AXIS_TEXT - 1.5, color = "grey40"),
        axis.ticks.y = element_blank(),
        plot.margin  = margin(6, 4, 0, 4))
ggsave(file.path(PNL_PNG, "MAIN_panel_B_logfc_density.png"), pB,
       width = PD_W, height = PD_H, units = "mm", dpi = 300)

pB_title    <- "Effect Size Distribution"
pB_subtitle <- sprintf("%s proteins | 4 contrasts",
                        format(length(unique(lfc_long$gene)), big.mark = ","))
pB <- strip_for_composite(pB)

# ── Panel C: DEPs per Contrast (stacked bar, 3 thresholds) ─────────────────
PA_W <- 67; PA_H <- 55
n_total <- length(all_genes)

SET_LABELS_C <- c(
  Cancer_vs_Healthy      = "CR vs H",
  Training_CR            = "Tr.(CR)",
  Baseline_Supplement    = "BL(Supp)",
  Training_CRE           = "Tr.(CRE)",
  Training_PLA           = "Tr.(PLA)",
  Supplement_Interaction = "CRE×PLA"
)

frac_df <- bind_rows(lapply(CONTRASTS_ALL, function(ctr) {
  pi_col  <- paste0("pi_score_", ctr)
  fdr_col <- paste0("adj.P.Val_", ctr)
  p_col   <- paste0("P.Value_",   ctr)
  tibble(
    contrast  = SET_LABELS_C[ctr],
    threshold = c("p < 0.05", "q < 0.05", "Pi < 0.05"),
    n = c(
      if (p_col   %in% names(dep_df)) sum(!is.na(dep_df[[p_col]])   & dep_df[[p_col]]   < 0.05) else 0L,
      if (fdr_col %in% names(dep_df)) sum(!is.na(dep_df[[fdr_col]]) & dep_df[[fdr_col]] < 0.10) else 0L,
      if (pi_col  %in% names(dep_df)) sum(!is.na(dep_df[[pi_col]])  & dep_df[[pi_col]]  < 0.05) else 0L
    )
  )
})) |>
  mutate(
    contrast  = factor(contrast,
                       levels = rev(unname(SET_LABELS_C[CONTRASTS_ALL]))),
    threshold = factor(threshold,
                       levels = c("p < 0.05", "q < 0.05", "Pi < 0.05")),
    pct       = 100 * n / n_total,
    fill_key  = paste(contrast, threshold, sep = "___")
  ) |>
  filter(n > 0)

FRAC_FILL <- c()
for (cname in unique(frac_df$contrast)) {
  col_idx <- match(cname, unname(SET_LABELS_C[CONTRASTS_ALL]))
  col     <- CONTRAST_COLORS[CONTRASTS_ALL[col_idx]]
  FRAC_FILL[paste(cname, "p < 0.05",  sep = "___")] <- adjustcolor(col, alpha.f = 0.15)
  FRAC_FILL[paste(cname, "q < 0.05",  sep = "___")] <- adjustcolor(col, alpha.f = 0.40)
  FRAC_FILL[paste(cname, "Pi < 0.05", sep = "___")] <- col
}

THRESH_LABEL <- c("p < 0.05" = "p", "q < 0.05" = "FDR", "Pi < 0.05" = "Pi")
label_df <- frac_df |>
  arrange(contrast, threshold) |>
  mutate(
    next_pct  = lead(pct, default = 0),
    seg_width = pct - next_pct,
    label_y   = (next_pct + pct) / 2,
    label     = THRESH_LABEL[as.character(threshold)],
    text_col  = if_else(threshold == "p < 0.05", "grey20", "white"),
    .by = contrast
  ) |>
  filter(seg_width > 0.3)

pi_total  <- sum(sapply(CONTRASTS_ALL, \(ctr) {
  pc <- paste0("pi_score_", ctr)
  if (pc %in% names(dep_df)) sum(!is.na(dep_df[[pc]]) & dep_df[[pc]] < 0.05) else 0
}))
fdr_total <- sum(sapply(CONTRASTS_ALL, \(ctr) {
  fc <- paste0("adj.P.Val_", ctr)
  if (fc %in% names(dep_df)) sum(!is.na(dep_df[[fc]]) & dep_df[[fc]] < 0.10) else 0
}))

pC <- ggplot(frac_df, aes(contrast, pct, fill = fill_key)) +
  geom_col(position = "identity", width = 0.75,
           color = "black", linewidth = 0.3) +
  geom_text(data = label_df,
            aes(x = contrast, y = label_y, label = label, color = I(text_col)),
            inherit.aes = FALSE, hjust = 0.5, size = 2.2, fontface = "bold") +
  scale_fill_manual(values = FRAC_FILL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     breaks = seq(0, 30, by = 5), limits = c(0, 30)) +
  coord_flip() +
  labs(title = "DEPs per Contrast",
       subtitle = sprintf("%s proteins | Pi %d | FDR(0.10) %d",
                          format(n_total, big.mark = ","), pi_total, fdr_total),
       x = NULL, y = "% of proteome", tag = "c") +
  FIG_THEME +
  theme(plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE,
                                      face = "bold.italic", color = "grey40"),
        legend.position = "none",
        axis.text.y = element_text(face = "bold", size = FIG_AXIS_TEXT - 0.5))
ggsave(file.path(PNL_PNG, "MAIN_panel_C_dep_counts.png"), pC,
       width = PA_W, height = PA_H, units = "mm", dpi = 300)

pC_title    <- "DEPs per Contrast"
pC_subtitle <- sprintf("%s proteins | Pi %d | FDR(0.10) %d",
                        format(n_total, big.mark = ","), pi_total, fdr_total)
pC <- strip_for_composite(pC)

# ── Panel D: UpSet overlap plot (4 contrasts) ────────────────────────────────
PC_W_D <- 67; PC_H_D <- 55

UPSET_COLORS <- CONTRAST_COLORS[CONTRASTS_UPSET]
names(UPSET_COLORS) <- SET_LABELS_UPSET[CONTRASTS_UPSET]

# Build sig sets and direction maps
sig_sets <- list(); dir_map <- list()
for (ctr in CONTRASTS_UPSET) {
  sig_pi_col  <- paste0("sig_pi_",   ctr)
  lfc_col     <- paste0("logFC_",    ctr)
  if (!sig_pi_col %in% names(dep_df)) next
  is_sig <- !is.na(dep_df[[sig_pi_col]]) & dep_df[[sig_pi_col]] != 0
  sig_sets[[ctr]] <- dep_df$gene[is_sig]
  dir_map[[ctr]]  <- setNames(
    ifelse(dep_df[[lfc_col]][is_sig] > 0, "Up", "Down"),
    dep_df$gene[is_sig]
  )
}

bin_mat           <- sapply(sig_sets, function(s) as.integer(all_genes %in% s))
rownames(bin_mat) <- all_genes
colnames(bin_mat) <- SET_LABELS_UPSET[colnames(bin_mat)]
label_to_contrast <- setNames(names(SET_LABELS_UPSET), SET_LABELS_UPSET)

cm      <- make_comb_mat(bin_mat, mode = "distinct")
cs      <- comb_size(cm)
keep_cm <- cs > 0 & comb_degree(cm) > 0
cm_sub  <- cm[keep_cm]

comb_names_vec    <- comb_name(cm_sub)
n_comb            <- length(comb_names_vec)
set_names_ordered <- set_name(cm_sub)
up_counts    <- numeric(n_comb)
down_counts  <- numeric(n_comb)

for (i in seq_len(n_comb)) {
  members <- extract_comb(cm_sub, comb_names_vec[i])
  if (length(members) == 0) next
  bits             <- as.logical(as.integer(strsplit(comb_names_vec[i], "")[[1]]))
  active_contrasts <- label_to_contrast[set_names_ordered[bits]]
  gene_dirs <- sapply(members, function(g) {
    dirs <- sapply(active_contrasts, function(ctr) {
      if (g %in% names(dir_map[[ctr]])) dir_map[[ctr]][g] else NA
    })
    dirs <- dirs[!is.na(dirs)]
    if (all(dirs == "Up")) "Up" else if (all(dirs == "Down")) "Down" else "Mixed"
  })
  up_counts[i]   <- sum(gene_dirs == "Up")
  down_counts[i] <- sum(gene_dirs == "Down")
}

# Pairwise overlap Fisher tests
overlap_tests  <- list()
contrast_pairs <- combn(CONTRASTS_UPSET, 2, simplify = FALSE)
n_bg           <- length(all_genes)
for (pair in contrast_pairs) {
  a <- sig_sets[[pair[1]]]; b <- sig_sets[[pair[2]]]
  n_both <- length(intersect(a, b))
  n_a    <- length(a); n_b <- length(b)
  mat <- matrix(c(n_both, n_a - n_both, n_b - n_both,
                  n_bg - n_a - n_b + n_both), nrow = 2)
  ft <- fisher.test(mat, alternative = "greater")
  overlap_tests[[paste(pair, collapse = " & ")]] <- data.frame(
    set_A = pair[1], set_B = pair[2],
    n_A = n_a, n_B = n_b, overlap = n_both,
    expected    = round(n_a * n_b / n_bg, 1),
    odds_ratio  = round(ft$estimate, 2),
    p_value     = ft$p.value
  )
}
overlap_df       <- bind_rows(overlap_tests)
overlap_df$p_bh  <- p.adjust(overlap_df$p_value, method = "BH")
n_sig_overlaps   <- sum(overlap_df$p_bh < 0.05, na.rm = TRUE)
write.csv(overlap_df, file.path(DAT, "panel_D_upset_overlap_enrich.csv"),
          row.names = FALSE)

# Sort by total size
display_total <- up_counts + down_counts
keep_display  <- display_total > 0
comb_ord      <- which(keep_display)[order(-display_total[keep_display])]
up_ord        <- up_counts[comb_ord]
down_ord      <- down_counts[comb_ord]
n_unique_deps <- sum(comb_size(cm_sub))

set_display_order <- c("CR vs H", "Tr.(CR)", "Tr.(CRE)", "Tr.(PLA)")
n_int             <- length(comb_ord)
comb_names_ord    <- comb_names_vec[comb_ord]
set_order_ch      <- set_name(cm_sub)
set_y_levels      <- rev(set_display_order)

comb_deg_ord <- vapply(comb_names_ord, function(cn) {
  sum(as.integer(strsplit(cn, "")[[1]]))
}, integer(1))

bar_long <- tibble(
  x         = rep(seq_len(n_int), 2),
  direction = factor(rep(c("Up", "Down"), each = n_int), levels = c("Down", "Up")),
  count     = c(up_ord, down_ord),
  is_single = rep(comb_deg_ord == 1, 2)
)

dot_df        <- expand_grid(x = seq_len(n_int), set = set_display_order)
dot_df$active <- vapply(seq_len(nrow(dot_df)), function(r) {
  bits <- as.integer(strsplit(comb_names_ord[dot_df$x[r]], "")[[1]])
  as.logical(bits[match(dot_df$set[r], set_order_ch)])
}, logical(1))
dot_df$set  <- factor(dot_df$set, levels = set_y_levels)
dot_df$ynum <- as.numeric(dot_df$set)

seg_list <- vector("list", n_int)
for (i in seq_len(n_int)) {
  bits <- as.logical(as.integer(strsplit(comb_names_ord[i], "")[[1]]))
  ypos <- match(set_order_ch[bits], set_y_levels)
  if (length(ypos) > 1)
    seg_list[[i]] <- tibble(x = i, ymin = min(ypos), ymax = max(ypos))
}
seg_df <- bind_rows(seg_list)

stripe_fills <- adjustcolor(unname(UPSET_COLORS[set_y_levels]), alpha.f = 0.20)

bar_bg_list <- lapply(seq_len(n_int), function(i) {
  bits <- as.logical(as.integer(strsplit(comb_names_ord[i], "")[[1]]))
  if (sum(bits) != 1) return(NULL)
  tibble(xmin = i - 0.5, xmax = i + 0.5,
         fill = unname(UPSET_COLORS[set_order_ch[bits]]))
})
bar_bg  <- bind_rows(compact(bar_bg_list))
lbl_sz  <- scale_text(BASE_COUNT, PC_W_D)
y_max   <- ceiling(max(c(up_ord, down_ord)) * 1.15 / 5) * 5

pD_bars <- ggplot(bar_long, aes(x, count, fill = direction)) +
  {if (nrow(bar_bg) > 0)
    geom_rect(data = bar_bg,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = bar_bg$fill, alpha = 0.20,
              color = "grey70", linewidth = 0.2, inherit.aes = FALSE)} +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.3) +
  geom_text(data = \(d) d |> filter(count > 4, is_single),
            aes(label = as.integer(count), y = count / 2),
            position = position_dodge(width = 0.7), vjust = 0.5,
            size = lbl_sz - 0.9, color = "white", fontface = "bold") +
  geom_text(data = \(d) d |> filter(count > 0, count <= 4, is_single),
            aes(label = as.integer(count), y = count + 1),
            position = position_dodge(width = 0.7), vjust = 0,
            size = lbl_sz - 0.9, color = "black", fontface = "bold") +
  geom_text(data = \(d) d |> filter(count == 0, is_single),
            aes(label = "0", y = 1),
            position = position_dodge(width = 0.7), vjust = 0,
            size = lbl_sz - 0.9, color = "black", fontface = "bold") +
  geom_text(data = \(d) d |> filter(count > 0, !is_single),
            aes(label = as.integer(count), y = count + 1),
            position = position_dodge(width = 0.7), vjust = 0,
            size = lbl_sz - 0.9, color = "black", fontface = "bold",
            check_overlap = TRUE) +
  scale_fill_manual(values = c(Up = unname(DIR_COLORS["Up"]),
                                Down = unname(DIR_COLORS["Down"]))) +
  scale_x_continuous(expand = expansion(add = 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)),
                     breaks = scales::breaks_pretty(n = 4)) +
  coord_cartesian(xlim = c(0.5, n_int + 0.5), ylim = c(0, y_max)) +
  labs(title = "Contrast Overlap (UpSet)",
       subtitle = sprintf("%d unique Pi DEPs | %d/%d sig pairwise overlaps",
                          n_unique_deps, n_sig_overlaps, nrow(overlap_df)),
       y = NULL) +
  FIG_THEME +
  theme(plot.subtitle      = element_text(size = FIG_SUBTITLE_SIZE,
                                          face = "bold.italic", color = "grey40"),
        axis.text.x        = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x       = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey80", linewidth = 0.3),
        panel.grid.minor   = element_blank(),
        legend.position    = "none",
        plot.margin        = margin(2, 0, 0, 0))

pD_dots <- ggplot() +
  annotate("rect",
           xmin = 0.5, xmax = n_int + 0.5,
           ymin = seq_along(set_y_levels) - 0.45,
           ymax = seq_along(set_y_levels) + 0.45,
           fill = stripe_fills) +
  {if (nrow(seg_df) > 0)
    geom_segment(data = seg_df,
                 aes(x = x, xend = x, y = ymin, yend = ymax),
                 linewidth = 0.6, color = "grey25")} +
  geom_point(data = dot_df |> filter(!active),
             aes(x = x, y = ynum), color = "grey78", size = 1.6) +
  geom_point(data = dot_df |> filter(active),
             aes(x = x, y = ynum), color = "grey15", size = 1.6) +
  scale_x_continuous(expand = expansion(add = 0)) +
  scale_y_continuous(breaks = seq_along(set_y_levels), labels = set_y_levels,
                     expand = expansion(add = 0)) +
  coord_cartesian(xlim = c(0.5, n_int + 0.5), ylim = c(0.50, 4.50)) +
  labs(x = NULL, y = NULL) +
  FIG_THEME +
  theme(axis.text.x  = element_blank(), axis.ticks = element_blank(),
        panel.grid   = element_blank(),
        panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.3),
        axis.text.y  = element_text(size = FIG_AXIS_TEXT - 0.5, face = "bold",
                                     margin = margin(r = 1)),
        plot.margin  = margin(0, 0, 0, 0))

dir_key_df_D <- tibble(
  label = c("Up", "Down"), y = c(0, -0.002),
  fill  = c(unname(DIR_COLORS["Up"]), unname(DIR_COLORS["Down"]))
)
p_key_dir_D <- ggplot(dir_key_df_D) +
  geom_point(aes(x = 0, y = y), shape = 22, size = 1.8,
             fill = dir_key_df_D$fill, color = "grey30", stroke = 0.3) +
  geom_text(aes(x = 0.20, y = y, label = label),
            size = 1.5, color = "grey20", hjust = 0) +
  scale_x_continuous(limits = c(-0.2, 1.2)) +
  coord_cartesian(ylim = c(-0.006, 0.003), clip = "off") +
  theme_void() + theme(plot.margin = margin(0, 0, 0, 0))

pD_pw_standalone <- (pD_bars / pD_dots) +
  plot_layout(heights = c(0.78, 0.22)) +
  plot_annotation(theme = theme(plot.margin = margin(t = 2, r = 2, b = 4, l = 0)))

pD_bars_clean <- pD_bars + labs(title = NULL, subtitle = NULL)
pD_pw <- (pD_bars_clean / pD_dots) +
  plot_layout(heights = c(0.78, 0.22)) +
  plot_annotation(theme = theme(plot.margin = margin(t = 4, r = 2, b = 4, l = 0)))

pD_standalone <- ggdraw(pD_pw_standalone) +
  draw_label("Intersection size", x = 0.02, y = 0.58, angle = 90,
             size = 5, fontface = "bold") +
  draw_plot(p_key_dir_D, x = 0.83, y = 0.70, width = 0.14, height = 0.22)

ggsave(file.path(PNL_PNG, "MAIN_panel_D_upset.png"), pD_standalone,
       width = PC_W_D, height = PC_H_D, units = "mm", dpi = 300)

pD <- ggdraw(pD_pw) +
  draw_label("Intersection size", x = 0.02, y = 0.58, angle = 90,
             size = 5, fontface = "bold") +
  draw_plot(p_key_dir_D, x = 0.85, y = 0.92, width = 0.13, height = 0.14)

pD_title    <- "Contrast Overlap (UpSet)"
pD_subtitle <- sprintf("%d unique Pi DEPs | %d/%d sig overlaps",
                        n_unique_deps, n_sig_overlaps, nrow(overlap_df))

# ── Panel E: fGSEA Stacked Bar (all 6 contrasts) ────────────────────────────
PE_W <- 52; PE_H <- 55

DB_ORDER          <- c("GO:BP", "Reactome", "Hallmark", "KEGG", "GO Slim")
DISPLAY_CONTRASTS <- CONTRASTS_ALL   # all 6

# Load and stack fgsea cache
fgsea_list <- lapply(DISPLAY_CONTRASTS, function(ctr) {
  f <- file.path(FGSEA_CACHE, paste0(ctr, "_fgsea.rds"))
  if (!file.exists(f)) return(NULL)
  df <- readRDS(f)
  df$contrast <- ctr
  df
})
fgsea_raw <- bind_rows(compact(fgsea_list))

sig_pathways <- fgsea_raw |>
  filter(!is.na(padj), padj < 0.05,
         database %in% DB_ORDER,
         contrast %in% DISPLAY_CONTRASTS) |>
  select(contrast, database, pathway, pval, padj, ES, NES, size) |>
  arrange(contrast, database, padj)

write_csv(sig_pathways, file.path(DAT, "panel_E_fgsea_sig.csv"))

count_df <- sig_pathways |>
  mutate(direction = ifelse(NES > 0, "Up", "Down")) |>
  group_by(contrast, direction, database) |>
  summarise(count = n(), .groups = "drop")

sig_counts_wide <- count_df |>
  pivot_wider(names_from = direction, values_from = count, values_fill = 0L) |>
  arrange(contrast, database)
write_csv(sig_counts_wide, file.path(DAT, "panel_E_fgsea_counts.csv"))

full_grid <- expand_grid(
  contrast  = DISPLAY_CONTRASTS,
  direction = c("Up", "Down"),
  database  = DB_ORDER
)
count_df <- full_grid |>
  left_join(count_df, by = c("contrast", "direction", "database")) |>
  mutate(count = replace_na(count, 0))

count_df$database <- factor(count_df$database, levels = DB_ORDER)

red_shades  <- colorRampPalette(c("#B2182B", "#D6604D", "#F4A582"))(length(DB_ORDER))
blue_shades <- colorRampPalette(c("#2166AC", "#4393C3", "#92C5DE"))(length(DB_ORDER))
names(red_shades)  <- DB_ORDER
names(blue_shades) <- DB_ORDER

BAR_W     <- 0.28
DODGE_GAP <- 0.07
ctr_centers <- setNames(seq_along(DISPLAY_CONTRASTS), DISPLAY_CONTRASTS)
count_df <- count_df |>
  mutate(
    x_center = ctr_centers[contrast] +
      ifelse(direction == "Up",
             -(BAR_W / 2 + DODGE_GAP / 2),
               BAR_W / 2 + DODGE_GAP / 2)
  )

count_df <- count_df |>
  arrange(contrast, direction, factor(database, levels = DB_ORDER)) |>
  group_by(contrast, direction) |>
  mutate(ymax = cumsum(count), ymin = ymax - count) |>
  ungroup() |>
  mutate(fill = ifelse(direction == "Up",
                       red_shades[as.character(database)],
                       blue_shades[as.character(database)]))

bar_tops <- count_df |>
  group_by(contrast, direction, x_center) |>
  summarise(total = sum(count), .groups = "drop")

bg_rects <- tibble(
  xmin = seq_along(DISPLAY_CONTRASTS) - 0.5,
  xmax = seq_along(DISPLAY_CONTRASTS) + 0.5,
  fill = CONTRAST_COLORS[DISPLAY_CONTRASTS]
)

lbl_sz_e <- scale_text(BASE_COUNT, PE_W)
y_max_e  <- max(bar_tops$total, na.rm = TRUE) * 1.15

pE_base <- ggplot() +
  geom_rect(data = bg_rects,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = bg_rects$fill, alpha = 0.20,
            color = "grey70", linewidth = 0.2) +
  geom_rect(data = count_df,
            aes(xmin = x_center - BAR_W / 2, xmax = x_center + BAR_W / 2,
                ymin = ymin, ymax = ymax),
            fill = count_df$fill, color = "white", linewidth = 0.25) +
  geom_text(data = bar_tops |> filter(total > 0),
            aes(x = x_center, y = total, label = total),
            vjust = -0.3, size = lbl_sz_e, fontface = "bold", color = "black") +
  scale_x_continuous(breaks = seq_along(DISPLAY_CONTRASTS),
                     labels = CTR_SHORT[DISPLAY_CONTRASTS],
                     expand = expansion(mult = 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)),
                     limits = c(0, max(y_max_e, 50))) +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "Significant pathways",
       title = "Pathway Enrichment (Up/Down)",
       subtitle = sprintf("fGSEA | 5 databases | per-db BH | %d sig / %d tested",
                          sum(count_df$count),
                          nrow(fgsea_raw |>
                               filter(database %in% DB_ORDER,
                                      contrast %in% DISPLAY_CONTRASTS)))) +
  FIG_THEME +
  theme(axis.title.y       = element_text(hjust = 0.54),
        axis.text.x        = element_text(angle = 35, hjust = 1,
                                          size = FIG_AXIS_TEXT - 0.5),
        legend.position    = "none",
        panel.grid.major.x = element_blank(),
        plot.margin        = margin(0, 0, 0, 0))

# Legends
grey_shades <- colorRampPalette(c("grey30", "grey75"))(length(DB_ORDER))
names(grey_shades) <- DB_ORDER
DB_LEGEND_ORDER <- rev(DB_ORDER)
shared_ylim <- c(-0.025, 0.005)

db_df <- tibble(
  label     = DB_LEGEND_ORDER,
  y         = -cumsum(c(0, rep(0.005, length(DB_ORDER) - 1))),
  fill      = grey_shades[DB_LEGEND_ORDER],
  is_header = rep(FALSE, length(DB_ORDER))
)
make_key_plot <- function(kdf) {
  ggplot(kdf) +
    geom_point(aes(x = 0, y = y), shape = 22, size = 1.8,
               fill = kdf$fill, color = "grey30", stroke = 0.3) +
    geom_text(aes(x = 0.35, y = y, label = label),
              size = 1.5, color = "grey20", hjust = 0) +
    scale_x_continuous(limits = c(-0.2, 1.5)) +
    coord_cartesian(ylim = shared_ylim, clip = "off") +
    theme_void() + theme(plot.margin = margin(0, 0, 0, 0))
}
dir_df <- tibble(
  label = c("Up", "Down"), y = c(0, -0.004),
  fill  = c(unname(DIR_COLORS["Up"]), unname(DIR_COLORS["Down"]))
)
p_key_db  <- make_key_plot(db_df)
p_key_dir <- ggplot(dir_df) +
  geom_point(aes(x = 0, y = y), shape = 22, size = 1.8,
             fill = dir_df$fill, color = "grey30", stroke = 0.3) +
  geom_text(aes(x = 0.35, y = y, label = label),
            size = 1.5, color = "grey20", hjust = 0) +
  scale_x_continuous(limits = c(-0.2, 1.5)) +
  coord_cartesian(ylim = c(-0.010, 0.003), clip = "off") +
  theme_void() + theme(plot.margin = margin(0, 0, 0, 0))

pE <- (pE_base +
  inset_element(p_key_db,  left = 0.60, right = 0.78,
                top = 1.00, bottom = 0.72) +
  inset_element(p_key_dir, left = 0.82, right = 1.00,
                top = 0.98, bottom = 0.85)) +
  plot_annotation(theme = theme(plot.margin = margin(t = 6, r = 3, b = 4, l = 3)))

ggsave(file.path(PNL_PNG, "MAIN_panel_E_fgsea.png"), pE,
       width = PE_W, height = PE_H, units = "mm", dpi = 300)

pE_title    <- "Pathway Enrichment (Up/Down)"
pE_subtitle <- sprintf("fGSEA | 5 databases | per-db BH | %d sig / %d tested",
                        sum(count_df$count),
                        nrow(fgsea_raw |>
                             filter(database %in% DB_ORDER,
                                    contrast %in% DISPLAY_CONTRASTS)))
pE <- pE & labs(title = NULL, subtitle = NULL, tag = NULL)

# ── Panel F: DEP Rank Barcode (4 contrasts) ──────────────────────────────────
PD_W_F <- 67; PD_H_F <- 55

rank_list <- lapply(CONTRASTS_UPSET, function(ctr) {
  t_col   <- paste0("t_",        ctr)
  pi_col  <- paste0("pi_score_", ctr)
  lfc_col <- paste0("logFC_",    ctr)
  if (!t_col %in% names(dep_df)) return(NULL)
  dep_df |>
    filter(!is.na(.data[[t_col]])) |>
    arrange(.data[[t_col]]) |>
    mutate(
      rank_frac = seq_len(n()) / n(),
      is_dep    = !is.na(.data[[pi_col]]) & .data[[pi_col]] < 0.05,
      direction = case_when(
        !is_dep              ~ NA_character_,
        .data[[lfc_col]] > 0 ~ "Up",
        TRUE                 ~ "Down"
      ),
      contrast = ctr
    ) |>
    select(gene, contrast, rank_frac, is_dep, direction)
})
rank_df          <- bind_rows(compact(rank_list))
rank_df$contrast <- factor(rank_df$contrast, levels = CONTRASTS_UPSET)

dep_only          <- rank_df |> filter(is_dep)
dep_only$direction <- factor(dep_only$direction, levels = c("Up", "Down"))

dep_counts <- dep_only |>
  group_by(contrast) |>
  summarise(n_up = sum(direction == "Up"), n_down = sum(direction == "Down"),
            n_total = n(), .groups = "drop") |>
  mutate(label = sprintf("n = %d  (%d Up  %d Dn)",
                          n_total, n_up, n_down))
write.csv(dep_counts, file.path(DAT, "panel_F_barcode_counts.csv"), row.names = FALSE)

DENS_PAD  <- 0.06
dens_list <- lapply(split(dep_only, dep_only$contrast), function(ctr_df) {
  lapply(split(ctr_df, ctr_df$direction, drop = TRUE), function(dir_df) {
    if (nrow(dir_df) < 2) return(NULL)
    d <- density(dir_df$rank_frac, adjust = 1.8,
                 from = -DENS_PAD, to = 1 + DENS_PAD, n = 512)
    tibble(x = d$x, y = d$y, direction = dir_df$direction[1],
           contrast = dir_df$contrast[1])
  }) |> bind_rows()
}) |> bind_rows()

dens_list <- dens_list |>
  group_by(contrast) |>
  mutate(y_norm = y / max(y)) |>
  ungroup()
dens_list$direction <- factor(dens_list$direction, levels = c("Up", "Down"))
dens_list$contrast  <- factor(dens_list$contrast, levels = CONTRASTS_UPSET)

TICK_DEPTH <- -0.25
ANNOT_SZ   <- scale_text(BASE_STAT - 0.5, PD_W_F)

peak_pos <- dens_list |>
  group_by(contrast, direction) |>
  slice_max(y_norm, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(contrast, direction, peak_x = x, peak_y = y_norm)

n_down <- dep_only |> filter(direction == "Down") |> count(contrast) |> tibble::deframe()
n_up   <- dep_only |> filter(direction == "Up")   |> count(contrast) |> tibble::deframe()

DESC_DOWN <- c(
  Cancer_vs_Healthy = "proteins lower in cancer recovery",
  Training_CR       = "proteins dec. with training",
  Training_CRE      = "proteins dec. (CRE arm)",
  Training_PLA      = "proteins dec. (PLA arm)"
)
DESC_UP <- c(
  Cancer_vs_Healthy = "proteins higher in cancer recovery",
  Training_CR       = "proteins inc. with training",
  Training_CRE      = "proteins inc. (CRE arm)",
  Training_PLA      = "proteins inc. (PLA arm)"
)

LABEL_NUDGE <- 0.06
ad_all <- peak_pos |>
  filter(direction == "Down") |>
  mutate(
    ctr     = as.character(contrast),
    label_x = peak_x + LABEL_NUDGE,
    label_y = peak_y * 0.78,
    label   = paste(n_down[ctr], DESC_DOWN[ctr])
  ) |> filter(!is.na(label))

au_all <- peak_pos |>
  filter(direction == "Up") |>
  mutate(
    ctr     = as.character(contrast),
    label_x = peak_x - LABEL_NUDGE,
    label_y = peak_y * 0.85,
    label   = paste(n_up[ctr], DESC_UP[ctr])
  ) |> filter(!is.na(label))

cd_all <- ad_all |>
  mutate(x_start = peak_x, y_start = peak_y, x_end = label_x, y_end = label_y)
cu_all <- au_all |>
  mutate(x_start = peak_x, y_start = peak_y, x_end = label_x, y_end = label_y)

bg_wash <- tibble(
  contrast = factor(CONTRASTS_UPSET, levels = CONTRASTS_UPSET),
  fill     = unname(CONTRAST_COLORS[CONTRASTS_UPSET]),
  xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
)

CTR_FACET_F <- setNames(CTR_SHORT[CONTRASTS_UPSET], CONTRASTS_UPSET)

pF <- ggplot() +
  geom_rect(data = bg_wash,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = bg_wash$fill, alpha = 0.18, inherit.aes = FALSE) +
  geom_ribbon(data = dens_list,
              aes(x = x, ymin = 0, ymax = y_norm, fill = direction),
              alpha = 0.30, outline.type = "upper") +
  geom_line(data = dens_list,
            aes(x = x, y = y_norm, color = direction), linewidth = 0.5) +
  geom_segment(data = dep_only,
               aes(x = rank_frac, xend = rank_frac,
                   y = 0, yend = TICK_DEPTH, color = direction),
               linewidth = 0.35, alpha = 0.8) +
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey50") +
  {if (nrow(cd_all) > 0)
    geom_segment(data = cd_all,
                 aes(x = x_start, xend = x_end, y = y_start, yend = y_end),
                 linewidth = 0.3, color = unname(DIR_COLORS["Down"]),
                 alpha = 0.4, inherit.aes = FALSE)} +
  {if (nrow(ad_all) > 0)
    geom_label(data = ad_all, aes(x = label_x, y = label_y, label = label),
               hjust = 0, vjust = 0.5, size = ANNOT_SZ,
               fill = unname(DIR_COLORS["Down"]), color = "white",
               fontface = "bold", linewidth = 0,
               label.padding = unit(0.08, "lines"), inherit.aes = FALSE)} +
  {if (nrow(cu_all) > 0)
    geom_segment(data = cu_all,
                 aes(x = x_start, xend = x_end, y = y_start, yend = y_end),
                 linewidth = 0.3, color = unname(DIR_COLORS["Up"]),
                 alpha = 0.4, inherit.aes = FALSE)} +
  {if (nrow(au_all) > 0)
    geom_label(data = au_all, aes(x = label_x, y = label_y, label = label),
               hjust = 1, vjust = 0.5, size = ANNOT_SZ,
               fill = unname(DIR_COLORS["Up"]), color = "white",
               fontface = "bold", linewidth = 0,
               label.padding = unit(0.08, "lines"), inherit.aes = FALSE)} +
  scale_fill_manual(values  = c(Up = unname(DIR_COLORS["Up"]),
                                 Down = unname(DIR_COLORS["Down"]))) +
  scale_color_manual(values = c(Up = unname(DIR_COLORS["Up"]),
                                 Down = unname(DIR_COLORS["Down"]))) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(-DENS_PAD, 1 + DENS_PAD), ylim = c(TICK_DEPTH, 1.02)) +
  facet_grid(contrast ~ ., switch = "y",
             labeller = labeller(contrast = CTR_FACET_F)) +
  labs(title = "DEP Rank Location",
       subtitle = sprintf("%s proteins | %d Pi DEPs (t-ranked)",
                          format(length(unique(rank_df$gene)), big.mark = ","),
                          sum(dep_counts$n_total)),
       x = "Rank position (by t-statistic)", y = NULL, tag = "f") +
  FIG_THEME +
  theme(legend.position    = "none",
        axis.text.y        = element_blank(), axis.ticks.y = element_blank(),
        strip.text.y.left  = element_text(face = "bold", size = FIG_AXIS_TEXT - 0.5,
                                           angle = 0, hjust = 1),
        strip.background   = element_blank(), strip.placement = "outside",
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        panel.spacing.y    = unit(2, "pt"))

ggsave(file.path(PNL_PNG, "MAIN_panel_F_barcode.png"), pF,
       width = PD_W_F, height = PD_H_F, units = "mm", dpi = 300)

pF_title    <- "DEP Rank Location"
pF_subtitle <- sprintf("%s proteins | %d Pi DEPs (t-ranked)",
                        format(length(unique(rank_df$gene)), big.mark = ","),
                        sum(dep_counts$n_total))
pF <- strip_for_composite(pF)

# ── Composite (3×2) ──────────────────────────────────────────────────────────
layout  <- "ABC\n###\nDEF"
ROW_TOP <- 0.46
SPACER  <- 0.04   # breathing room between rows

pA <- pA + theme(plot.margin = margin(14,  6, 18,  6))
pB <- pB + theme(plot.margin = margin(14,  6, 18,  6))
pC <- pC + theme(plot.margin = margin(14, 18, 18,  4))   # extra right margin so "p" label clears edge
pF <- pF + theme(plot.margin = margin( 6,  6,  6,  6))

composite <- wrap_elements(full = pA) + pB + pC +
             wrap_elements(full = pD) +
             wrap_elements(full = pE) +
             pF +
  plot_layout(
    design  = layout,
    widths  = c(155, 120, 150),   # widen col C, trim col B slightly
    heights = c(ROW_TOP, SPACER, 1 - ROW_TOP - SPACER)
  )

COMP_W <- 220; COMP_H <- 150
txt    <- composite_text_sizes(COMP_H)
TAG_SZ <- txt$tag
TTL_SZ <- txt$title
SUB_SZ <- txt$subtitle

TOP_Y      <- 0.985
BOT_Y      <- 1 - ROW_TOP - SPACER + 0.025
X_LEFT     <- 0.002
X_MID      <- 0.365   # col B: 155/425
X_RIGHT    <- 0.647   # col C: 275/425
X_TTL      <- 0.04
SUB_OFFSET <- 0.022

composite <- ggdraw(composite) +
  # Panel A
  draw_label("A",         x = X_LEFT,          y = TOP_Y,              size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pA_title,    x = X_LEFT + X_TTL,  y = TOP_Y,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pA_subtitle, x = X_LEFT + X_TTL,  y = TOP_Y - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30") +
  # Panel B
  draw_label("B",         x = X_MID,           y = TOP_Y,              size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pB_title,    x = X_MID + X_TTL,   y = TOP_Y,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pB_subtitle, x = X_MID + X_TTL,   y = TOP_Y - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30") +
  # Panel C
  draw_label("C",         x = X_RIGHT,         y = TOP_Y,              size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pC_title,    x = X_RIGHT + X_TTL, y = TOP_Y,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pC_subtitle, x = X_RIGHT + X_TTL, y = TOP_Y - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30") +
  # Panel D
  draw_label("D",         x = X_LEFT,          y = BOT_Y,              size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pD_title,    x = X_LEFT + X_TTL,  y = BOT_Y,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pD_subtitle, x = X_LEFT + X_TTL,  y = BOT_Y - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30") +
  # Panel E
  draw_label("E",         x = X_MID,           y = BOT_Y,              size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pE_title,    x = X_MID + X_TTL,   y = BOT_Y,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pE_subtitle, x = X_MID + X_TTL,   y = BOT_Y - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30") +
  # Panel F
  draw_label("F",         x = X_RIGHT,         y = BOT_Y,              size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pF_title,    x = X_RIGHT + X_TTL, y = BOT_Y,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(pF_subtitle, x = X_RIGHT + X_TTL, y = BOT_Y - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey30")

# Reset path vars (panels may have overwritten them)
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")

ggsave(file.path(RPT_PDF, "MAIN_F02_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_dev)
ggsave(file.path(RPT_PNG, "MAIN_F02_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message("F02 main composite done")
