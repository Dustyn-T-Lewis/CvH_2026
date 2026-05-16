#!/usr/bin/env Rscript
# F04 Training Concordance — Supplementary (6-panel composite)
# Sam's CvH parallel of YvO F04 SUPP diagnostics
#
# Panels:
#   A: ORA Dedup Sensitivity (concordant quadrant genes)
#   B: Pearson r Bootstrap (logFC CRE vs PLA)
#   C: Circularity Diagnostic (permuted null for r)
#   D: Threshold Sensitivity (concordant/discordant % by |logFC| cutoff)
#   E: GO Slim Category Distribution by concordance quadrant
#   F: Leading-Edge Proteins (top by |t_PLA|, from CRE leading edge)
#
# NOTE: 0 DEPs at Pi<0.05 in both arms. Panels A and F use CRE fGSEA
# leading-edge genes as proxy for "directed" proteins. Annotated in subtitles.
#
# Output: b_reports/supp/pdf/SUPP_F04_training_concordance_diagnostics.{pdf,png}

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(readr)
  library(stringr)
  library(boot)
})

# ── Style + helpers ──────────────────────────────────────────────────────────
source("02-03_Sam's_Results/04_Figures/shared/style.R")
source("02-03_Sam's_Results/04_Figures/shared/pathway_utils.R")

pdf_device <- grDevices::pdf

if (!exists("strip_for_composite")) {
  strip_for_composite <- function(p) {
    p + labs(title = NULL, subtitle = NULL, tag = NULL) +
      theme(legend.position = "none")
  }
}
if (!exists("composite_text_sizes")) {
  composite_text_sizes <- function(comp_h_mm) {
    list(title    = pmax(6, pmin(8, round(5 + comp_h_mm / 80))),
         subtitle = pmax(4, pmin(6, round(3 + comp_h_mm / 100))),
         tag      = 8)
  }
}

# ── Paths ────────────────────────────────────────────────────────────────────
SAM_DEP_DIR <- "02-03_Sam's_Results/03_DEP/c_data/04_per_contrast_results"
FGSEA_CACHE <- "02-03_Sam's_Results/04_Figures/shared/fgsea_cache"
BASE        <- "02-03_Sam's_Results/04_Figures/F04_training_concordance"
RPT_PDF     <- file.path(BASE, "b_reports", "supp", "pdf")
RPT_PNG     <- file.path(BASE, "b_reports", "supp", "png")
PNL_PDF     <- file.path(RPT_PDF, "panels")
PNL_PNG     <- file.path(RPT_PNG, "panels")
DAT         <- file.path(BASE, "c_data", "panel_supp")
for (d in c(RPT_PDF, RPT_PNG, PNL_PDF, PNL_PNG, DAT))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ── Load data ────────────────────────────────────────────────────────────────
cre_dep <- read_csv(file.path(SAM_DEP_DIR, "Training_CRE.csv"), show_col_types = FALSE)
pla_dep <- read_csv(file.path(SAM_DEP_DIR, "Training_PLA.csv"), show_col_types = FALSE)

dep_df <- cre_dep |>
  select(gene, logFC_CRE = logFC, t_CRE = t, pi_CRE = pi_score) |>
  inner_join(pla_dep |> select(gene, logFC_PLA = logFC, t_PLA = t, pi_PLA = pi_score),
             by = "gene") |>
  filter(!is.na(logFC_CRE), !is.na(logFC_PLA))

scatter_df <- dep_df |>
  mutate(quadrant = case_when(
    logFC_CRE > 0 & logFC_PLA > 0 ~ "Concordant Up",
    logFC_CRE < 0 & logFC_PLA < 0 ~ "Concordant Down",
    logFC_CRE > 0 & logFC_PLA < 0 ~ "Discordant (CRE up / PLA down)",
    TRUE                           ~ "Discordant (CRE down / PLA up)"
  ))

universe <- dep_df$gene

# fGSEA caches
cre_fgsea <- readRDS(file.path(FGSEA_CACHE, "Training_CRE_fgsea.rds")) |>
  filter(!is.na(NES))

# CRE leading-edge gene pools (proxy for directed gene sets, 0 DEPs at Pi<0.05)
cre_up_genes <- cre_fgsea |>
  filter(NES > 0, !is.na(padj), padj < 0.05) |>
  arrange(desc(NES)) |>
  head(10) |>
  pull(leadingEdge) |> strsplit(";") |> unlist() |> unique()
cre_dn_genes <- cre_fgsea |>
  filter(NES < 0, !is.na(padj), padj < 0.05) |>
  arrange(NES) |>
  head(10) |>
  pull(leadingEdge) |> strsplit(";") |> unlist() |> unique()

# ── Panel A: ORA Dedup Sensitivity ──────────────────────────────────────────
message("=== SUPP Panel A: ORA dedup sensitivity ===")

pw_collection <- build_pathway_collection(min_size = 15, max_size = 500,
                                           include_goslim = FALSE,
                                           exclude_variants = TRUE)
cutoffs  <- c(0.3, 0.5, 0.7, 1.0)
# Use concordant quadrant genes (all proteins, no pi filter) + CRE leading-edge
quad_sets <- list(
  "Concordant Up"   = scatter_df$gene[scatter_df$quadrant == "Concordant Up"],
  "Concordant Down" = scatter_df$gene[scatter_df$quadrant == "Concordant Down"],
  "CRE Leading Up"  = cre_up_genes,
  "CRE Leading Down" = cre_dn_genes
)

sens_results <- list()
for (qs_name in names(quad_sets)) {
  genes_q <- quad_sets[[qs_name]]
  if (length(genes_q) < 5) next
  for (jc in cutoffs) {
    ora_res <- tryCatch(
      run_ora_deduplicated(genes = genes_q, universe = universe,
                           pathways = pw_collection, jaccard_cutoff = jc,
                           min_size = 15, max_size = 500, padj_cutoff = 1),
      error = function(e) tibble()
    )
    n_sig <- if (nrow(ora_res) > 0) sum(ora_res$padj < 0.05) else 0L
    sens_results[[length(sens_results) + 1]] <- tibble(
      gene_set       = qs_name,
      jaccard_cutoff = jc,
      n_enriched     = n_sig,
      n_total_tested = nrow(ora_res)
    )
  }
}
sens_df <- bind_rows(sens_results) |>
  mutate(cutoff_label = factor(sprintf("J = %.1f", jaccard_cutoff)))

write_csv(sens_df, file.path(DAT, "SUPP_ora_dedup_sensitivity.csv"))

pS_ora_dedup <- ggplot(sens_df, aes(x = gene_set, y = n_enriched, fill = cutoff_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "grey30", linewidth = 0.3) +
  scale_fill_brewer(palette = "Blues", name = "Jaccard cutoff") +
  labs(title    = "ORA Dedup Sensitivity",
       subtitle = sprintf("Hypergeometric ORA | %d universe | concordant/leading-edge sets",
                          length(universe)),
       x = NULL, y = "Enriched pathways (FDR < 0.05)") +
  FIG_THEME +
  theme(legend.position = "right",
        axis.text.x = element_text(size = FIG_AXIS_TEXT - 1, face = "bold",
                                    angle = 20, hjust = 1))

PW <- 89; PH <- 70
ggsave(file.path(PNL_PNG, "SUPP_ora_dedup.png"), pS_ora_dedup,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(PNL_PDF, "SUPP_ora_dedup.pdf"), pS_ora_dedup,
       width = PW, height = PH, units = "mm", device = pdf_device)

pS_ora_title    <- "ORA Dedup Sensitivity"
pS_ora_subtitle <- sprintf("Hypergeometric ORA | concordant/leading-edge sets | %d universe",
                            length(universe))
pS_ora_dedup    <- strip_for_composite(pS_ora_dedup)
message("SUPP Panel A done")

# ── Panel B: Pearson r Bootstrap ─────────────────────────────────────────────
message("=== SUPP Panel B: r bootstrap ===")

set.seed(42)
boot_r_fn <- function(data, indices) {
  d <- data[indices, ]
  cor(d$logFC_CRE, d$logFC_PLA, use = "complete.obs")
}
b      <- boot(dep_df, statistic = boot_r_fn, R = 1000)
obs_r  <- b$t0
ci     <- boot.ci(b, type = "perc", conf = 0.95)
ci_lo  <- ci$percent[4]
ci_hi  <- ci$percent[5]

boot_df <- tibble(replicate = seq_len(1000), r = as.numeric(b$t))
write_csv(boot_df, file.path(DAT, "SUPP_r_bootstrap.csv"))

pS_r_boot <- ggplot(boot_df, aes(x = r)) +
  geom_histogram(bins = 40, fill = "#5DA5DA", color = "white", linewidth = 0.3) +
  geom_vline(xintercept = obs_r, color = "#D6604D", linewidth = 0.7) +
  geom_vline(xintercept = c(ci_lo, ci_hi), color = "#D6604D",
             linewidth = 0.5, linetype = "dashed") +
  annotate("label", x = obs_r, y = Inf, vjust = 1.5, hjust = 1.1,
           label = sprintf("r = %.3f\n95%% CI [%.3f, %.3f]", obs_r, ci_lo, ci_hi),
           size = BASE_STAT, fontface = "bold", fill = alpha("white", 0.9),
           label.padding = unit(1.5, "pt")) +
  labs(title    = "Concordance Pearson r Bootstrap",
       subtitle = sprintf("boot::boot() R = 1000 | n = %d proteins | percentile CI",
                          nrow(dep_df)),
       x = "Pearson r (logFC Training CRE vs logFC Training PLA)",
       y = "Count") +
  FIG_THEME

ggsave(file.path(PNL_PNG, "SUPP_r_bootstrap.png"), pS_r_boot,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(PNL_PDF, "SUPP_r_bootstrap.pdf"), pS_r_boot,
       width = PW, height = PH, units = "mm", device = pdf_device)

pS_rboot_title    <- "Concordance Pearson r Bootstrap"
pS_rboot_subtitle <- sprintf("boot::boot() R = 1000 | n = %d | r = %.3f [%.3f, %.3f]",
                               nrow(dep_df), obs_r, ci_lo, ci_hi)
pS_r_boot         <- strip_for_composite(pS_r_boot)
message("SUPP Panel B done")

# ── Panel C: Circularity Diagnostic ─────────────────────────────────────────
message("=== SUPP Panel C: circularity ===")

obs_r_perm <- cor(dep_df$logFC_CRE, dep_df$logFC_PLA, use = "complete.obs")
set.seed(42)
n_perm  <- 1000
perm_r  <- numeric(n_perm)
to_vec  <- dep_df$logFC_PLA
for (i in seq_len(n_perm)) {
  perm_r[i] <- cor(sample(dep_df$logFC_CRE), to_vec, use = "complete.obs")
}

null_mean <- mean(perm_r)
null_sd   <- sd(perm_r)
p_perm    <- mean(abs(perm_r) >= abs(obs_r_perm))

circ_df <- tibble(replicate = seq_len(n_perm), perm_r = perm_r)
write_csv(circ_df, file.path(DAT, "SUPP_concordance_circularity.csv"))

sub_circ <- sprintf("Observed r = %.3f | Null mean = %.4f | p_perm = %.3f",
                    obs_r_perm, null_mean, p_perm)

pS_circ <- ggplot(circ_df, aes(x = perm_r)) +
  geom_histogram(bins = 50, fill = "grey70", color = "white", linewidth = 0.3) +
  geom_vline(xintercept = obs_r_perm, color = "#D6604D", linewidth = 0.8) +
  geom_vline(xintercept = null_mean, color = "#4393C3", linewidth = 0.5,
             linetype = "dashed") +
  coord_cartesian(clip = "off") +
  annotate("label", x = obs_r_perm, y = Inf, vjust = 1.3, hjust = 1.1,
           label = sprintf("Observed\nr = %.3f", obs_r_perm),
           size = 2.5, fontface = "bold", color = "#D6604D",
           fill = alpha("white", 0.9), label.padding = unit(2, "pt")) +
  annotate("label", x = null_mean, y = Inf, vjust = 1.3, hjust = -0.1,
           label = sprintf("Null mean\n= %.4f", null_mean),
           size = 2.2, fontface = "bold", color = "#4393C3",
           fill = alpha("white", 0.9), label.padding = unit(2, "pt")) +
  labs(title    = "Circularity Diagnostic: Protein-Permuted Null",
       subtitle = sub_circ,
       x = "Permuted Pearson r",
       y = "Count") +
  FIG_THEME

ggsave(file.path(PNL_PNG, "SUPP_circularity.png"), pS_circ,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(PNL_PDF, "SUPP_circularity.pdf"), pS_circ,
       width = PW, height = PH, units = "mm", device = pdf_device)

pS_circ         <- pS_circ + theme(plot.margin = margin(18, 8, 6, 10))
pS_circ_title    <- "Circularity Diagnostic: Protein-Permuted Null"
pS_circ_subtitle <- sub_circ
pS_circ          <- strip_for_composite(pS_circ)
message("SUPP Panel C done")

# ── Panel D: Concordance Threshold Sensitivity ────────────────────────────────
message("=== SUPP Panel D: concordance threshold ===")

thresholds <- c(0.05, 0.1, 0.2, 0.3)
thresh_results <- lapply(thresholds, function(thr) {
  classified <- dep_df |>
    mutate(category = case_when(
      (logFC_CRE > thr  & logFC_PLA > thr)  |
        (logFC_CRE < -thr & logFC_PLA < -thr) ~ "Concordant",
      (logFC_CRE > thr  & logFC_PLA < -thr) |
        (logFC_CRE < -thr & logFC_PLA > thr)  ~ "Discordant",
      TRUE ~ "Negligible"
    ))
  n_total <- nrow(classified)
  classified |>
    count(category) |>
    mutate(threshold = thr, pct = 100 * n / n_total)
})
thresh_df <- bind_rows(thresh_results)
write_csv(thresh_df, file.path(DAT, "SUPP_concordance_threshold.csv"))

cat_colors <- c("Concordant" = "#D6604D", "Discordant" = "#4393C3", "Negligible" = "grey60")
thresh_df <- thresh_df |>
  mutate(category = factor(category, levels = c("Concordant", "Discordant", "Negligible")))

pS_thresh <- ggplot(thresh_df, aes(x = threshold, y = pct, color = category)) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2.5, shape = 16) +
  scale_color_manual(values = cat_colors, name = "Category") +
  scale_x_continuous(breaks = thresholds,
                     labels = sprintf("%.2f", thresholds)) +
  labs(title    = "Concordance Classification Threshold Sensitivity",
       subtitle = sprintf("logFC thresholds: %s | stable concordance %%",
                          paste(thresholds, collapse = ", ")),
       x = "|logFC| threshold",
       y = "Percentage of proteins") +
  FIG_THEME +
  theme(legend.position = "right")

ggsave(file.path(PNL_PNG, "SUPP_concordance_threshold.png"), pS_thresh,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(PNL_PDF, "SUPP_concordance_threshold.pdf"), pS_thresh,
       width = PW, height = PH, units = "mm", device = pdf_device)

pS_thresh_title    <- "Concordance Classification Threshold Sensitivity"
pS_thresh_subtitle <- sprintf("logFC thresholds: %s | stable concordance %%",
                               paste(thresholds, collapse = ", "))
pS_thresh          <- strip_for_composite(pS_thresh)
message("SUPP Panel D done")

# ── Panel E: GO Slim Category Distribution ────────────────────────────────────
message("=== SUPP Panel E: GO Slim bars ===")

# Use all proteins split by concordance quadrant (no pi filter — 0 DEPs)
dep_sig <- scatter_df |>
  mutate(quadrant = case_when(
    quadrant %in% c("Concordant Up", "Concordant Down") ~ quadrant,
    TRUE                                                 ~ "Discordant"
  ))

yvo_slim_path <- file.path(dirname(rprojroot::find_rstudio_root_file()),
                            "A_YvO_2026", "04_Figures", "shared", "go_slim_categories.R")

slim_merged <- tryCatch({
  if (file.exists(yvo_slim_path)) {
    source(yvo_slim_path)
    fg_genes  <- dep_sig$gene
    all_genes <- dep_df$gene
    slim_result <- assign_go_slim_consolidated(fg_genes, all_genes)
    dep_sig |> left_join(slim_result, by = "gene") |> filter(!is.na(consolidated))
  } else {
    stop("go_slim not available")
  }
}, error = function(e) {
  message("  GO Slim fallback: using pathway category classifier")
  dep_sig |>
    mutate(consolidated = classify_pathway_func(gene)) |>
    filter(!is.na(consolidated), consolidated != "Other")
})

slim_counts_raw <- slim_merged |>
  count(consolidated, quadrant, name = "count")
present_cats <- if (exists("CONSOLIDATED_PATHWAY_ORDER")) {
  CONSOLIDATED_PATHWAY_ORDER[CONSOLIDATED_PATHWAY_ORDER %in% unique(slim_counts_raw$consolidated)]
} else {
  unique(slim_counts_raw$consolidated)
}
slim_counts <- slim_counts_raw |>
  mutate(consolidated = factor(consolidated, levels = rev(present_cats)))

write_csv(slim_counts, file.path(DAT, "SUPP_goslim_distribution.csv"))

quad_colors_slim <- c("Concordant Up" = "#E57373", "Concordant Down" = "#64B5F6",
                      "Discordant"    = "grey60")
slim_counts <- slim_counts |>
  mutate(quadrant = factor(quadrant,
                           levels = c("Concordant Up", "Concordant Down", "Discordant")))

pS_goslim <- ggplot(slim_counts,
                    aes(x = count, y = consolidated, fill = quadrant)) +
  geom_col(position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = quad_colors_slim, name = "Quadrant") +
  labs(title    = "GO Slim Category Distribution (Training Concordance)",
       subtitle = sprintf("All proteins | %d proteins | %d categories (no pi filter; 0 DEPs)",
                          nrow(slim_merged),
                          n_distinct(slim_merged$consolidated)),
       x = "Protein count", y = NULL) +
  FIG_THEME +
  theme(legend.position    = "right",
        axis.text.y        = element_text(size = FIG_AXIS_TEXT),
        panel.grid.major.y = element_blank())

ggsave(file.path(PNL_PNG, "SUPP_goslim_bars.png"), pS_goslim,
       width = PW, height = 85, units = "mm", dpi = 300)
ggsave(file.path(PNL_PDF, "SUPP_goslim_bars.pdf"), pS_goslim,
       width = PW, height = 85, units = "mm", device = pdf_device)

pS_goslim_title    <- "GO Slim Category Distribution (Training Concordance)"
pS_goslim_subtitle <- sprintf("All proteins | %d proteins | %d categories (no pi filter; 0 DEPs)",
                               nrow(slim_merged), n_distinct(slim_merged$consolidated))
pS_goslim          <- strip_for_composite(pS_goslim)
message("SUPP Panel E done")

# ── Panel F: Leading-Edge Proteins ──────────────────────────────────────────
message("=== SUPP Panel F: leading edge proteins ===")

# Use CRE leading-edge genes and show their PLA t-stats (concordance signal)
t_PLA_named <- setNames(dep_df$t_PLA, dep_df$gene)
t_PLA_named <- t_PLA_named[!is.na(t_PLA_named)]

driving_df <- dep_df |>
  filter(gene %in% c(cre_up_genes, cre_dn_genes), !is.na(t_PLA)) |>
  mutate(
    set        = if_else(gene %in% cre_up_genes, "cre_up", "cre_dn"),
    set_label  = if_else(set == "cre_up", "CRE Leading Up", "CRE Leading Down"),
    is_conc    = (set == "cre_up" & t_PLA > 0) | (set == "cre_dn" & t_PLA < 0)
  ) |>
  arrange(desc(abs(t_PLA))) |>
  slice_head(n = 25)

write_csv(driving_df, file.path(DAT, "SUPP_fry_leading_edge.csv"))

dir_colors_f <- c("cre_up" = "#D6604D", "cre_dn" = "#4393C3")
dir_labels_f <- c("cre_up" = "CRE Leading Up (concordant = PLA up)",
                  "cre_dn" = "CRE Leading Down (concordant = PLA down)")

pS_fry_lead <- ggplot(driving_df,
                       aes(x = t_PLA, y = reorder(gene, abs(t_PLA)),
                           color = set)) +
  geom_point(size = 2) +
  geom_segment(aes(xend = 0, yend = reorder(gene, abs(t_PLA))),
               linewidth = 0.4) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  scale_color_manual(values = dir_colors_f, labels = dir_labels_f,
                     name = "CRE direction") +
  labs(title    = "Leading-Edge Proteins (Concordance Drivers)",
       subtitle = sprintf("Top %d CRE leading-edge genes by |t(Training PLA)|",
                          nrow(driving_df)),
       x = "t-statistic (Training PLA)",
       y = NULL) +
  FIG_THEME +
  theme(legend.position  = "bottom",
        legend.direction = "horizontal",
        axis.text.y      = element_text(size = FIG_AXIS_TEXT, face = "italic"))

ggsave(file.path(PNL_PNG, "SUPP_fry_leading.png"), pS_fry_lead,
       width = PW, height = 100, units = "mm", dpi = 300)
ggsave(file.path(PNL_PDF, "SUPP_fry_leading.pdf"), pS_fry_lead,
       width = PW, height = 100, units = "mm", device = pdf_device)

pS_lead_title    <- "Leading-Edge Proteins (Concordance Drivers)"
pS_lead_subtitle <- sprintf("Top %d CRE leading-edge genes by |t(Training PLA)|",
                              nrow(driving_df))
pS_fry_lead      <- strip_for_composite(pS_fry_lead)
message("SUPP Panel F done")

# ── Composite ────────────────────────────────────────────────────────────────
COMP_W <- 260
COMP_H <- 310
CTS    <- composite_text_sizes(COMP_H)
CTS$title    <- CTS$title + 1
CTS$subtitle <- CTS$subtitle + 0.5

axis_fix <- theme(axis.title.y = element_text(margin = margin(0, 2, 0, 0)),
                  axis.title.x = element_text(margin = margin(2, 0, 0, 0)))

grid <- (pS_ora_dedup | pS_r_boot) /
        (pS_circ      | pS_thresh) /
        (pS_goslim    | pS_fry_lead) &
  theme(plot.margin = margin(18, 8, 6, 5),
        axis.title  = element_text(size = 7, face = "bold")) &
  axis_fix

X_LEFT     <- 0.015
X_RIGHT    <- 0.525
X_TTL      <- 0.028
SUB_OFFSET <- 0.013
Y_R1 <- 0.993
Y_R2 <- 0.663
Y_R3 <- 0.333

composite_supp <- ggdraw(grid) +
  draw_label("A", x = X_LEFT, y = Y_R1,
             fontface = "bold", size = CTS$tag, hjust = 0, vjust = 1) +
  draw_label(pS_ora_title, x = X_LEFT + X_TTL, y = Y_R1,
             fontface = "bold", size = CTS$title, hjust = 0, vjust = 1) +
  draw_label(pS_ora_subtitle, x = X_LEFT + X_TTL, y = Y_R1 - SUB_OFFSET,
             fontface = "bold.italic", size = CTS$subtitle, colour = "grey30",
             hjust = 0, vjust = 1) +
  draw_label("B", x = X_RIGHT, y = Y_R1,
             fontface = "bold", size = CTS$tag, hjust = 0, vjust = 1) +
  draw_label(pS_rboot_title, x = X_RIGHT + X_TTL, y = Y_R1,
             fontface = "bold", size = CTS$title, hjust = 0, vjust = 1) +
  draw_label(pS_rboot_subtitle, x = X_RIGHT + X_TTL, y = Y_R1 - SUB_OFFSET,
             fontface = "bold.italic", size = CTS$subtitle, colour = "grey30",
             hjust = 0, vjust = 1) +
  draw_label("C", x = X_LEFT, y = Y_R2,
             fontface = "bold", size = CTS$tag, hjust = 0, vjust = 1) +
  draw_label(pS_circ_title, x = X_LEFT + X_TTL, y = Y_R2,
             fontface = "bold", size = CTS$title, hjust = 0, vjust = 1) +
  draw_label(pS_circ_subtitle, x = X_LEFT + X_TTL, y = Y_R2 - SUB_OFFSET,
             fontface = "bold.italic", size = CTS$subtitle, colour = "grey30",
             hjust = 0, vjust = 1) +
  draw_label("D", x = X_RIGHT, y = Y_R2,
             fontface = "bold", size = CTS$tag, hjust = 0, vjust = 1) +
  draw_label(pS_thresh_title, x = X_RIGHT + X_TTL, y = Y_R2,
             fontface = "bold", size = CTS$title, hjust = 0, vjust = 1) +
  draw_label(pS_thresh_subtitle, x = X_RIGHT + X_TTL, y = Y_R2 - SUB_OFFSET,
             fontface = "bold.italic", size = CTS$subtitle, colour = "grey30",
             hjust = 0, vjust = 1) +
  draw_label("E", x = X_LEFT, y = Y_R3,
             fontface = "bold", size = CTS$tag, hjust = 0, vjust = 1) +
  draw_label(pS_goslim_title, x = X_LEFT + X_TTL, y = Y_R3,
             fontface = "bold", size = CTS$title, hjust = 0, vjust = 1) +
  draw_label(pS_goslim_subtitle, x = X_LEFT + X_TTL, y = Y_R3 - SUB_OFFSET,
             fontface = "bold.italic", size = CTS$subtitle, colour = "grey30",
             hjust = 0, vjust = 1) +
  draw_label("F", x = X_RIGHT, y = Y_R3,
             fontface = "bold", size = CTS$tag, hjust = 0, vjust = 1) +
  draw_label(pS_lead_title, x = X_RIGHT + X_TTL, y = Y_R3,
             fontface = "bold", size = CTS$title, hjust = 0, vjust = 1) +
  draw_label(pS_lead_subtitle, x = X_RIGHT + X_TTL, y = Y_R3 - SUB_OFFSET,
             fontface = "bold.italic", size = CTS$subtitle, colour = "grey30",
             hjust = 0, vjust = 1)

ggsave(file.path(RPT_PDF, "SUPP_F04_training_concordance_diagnostics.pdf"),
       composite_supp, width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_PNG, "SUPP_F04_training_concordance_diagnostics.png"),
       composite_supp, width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message(sprintf("F04 SUPP composite saved -> %s / %s", RPT_PDF, RPT_PNG))
