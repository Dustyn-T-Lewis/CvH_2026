#!/usr/bin/env Rscript
# F05 Recovery Reversal — Main (5-panel composite)
# Sam's CvH parallel of YvO F05 "Aging Reversal"
#
# Question: does resistance training REVERSE the Cancer-vs-Healthy proteome?
# Contrast mapping:
#   YvO Aging        -> Cancer_vs_Healthy (SURV - CTL)
#   YvO Training_Old -> Training_CR       (combined pre->post)
#
# Panels (inline — no per-panel helper files):
#   A: Quadrant ORA scatter + flanking bars
#   B: Pattern heatmap (per-protein reversal classification)
#   C: fry rotation test (Cancer DEP sets vs Training_CR t-stats)
#   D: NES scatter (pathway-level reversal)
#   E: RRHO2 (threshold-free rank-rank overlap)
#
# Layout (3-col):
#   Top row:    A (Quadrant ORA) | B (Pattern heatmap, full height)
#   Bottom row: C (fry barcode)  | D (NES scatter) | E (RRHO2)

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(ggrepel)
  library(fgsea)
  library(RRHO2)
  library(boot)
})

# ── Style + helpers ──────────────────────────────────────────────────────────
source("02-03_Sam's_Results/04_Figures/shared/style.R")
source("02-03_Sam's_Results/04_Figures/shared/pathway_utils.R")

# cairo_pdf requires X11 (unavailable on this system); fall back to base pdf()
pdf_device <- grDevices::pdf

# Print-scale for 380mm-wide composite (matches YvO panel_A_ORA sizing)
PRINT_SCALE <- 380 / 178

# Helpers not present in CvH style.R
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
SAM_DEP_DIR  <- "02-03_Sam's_Results/03_DEP/c_data/04_per_contrast_results"
FGSEA_CACHE  <- "02-03_Sam's_Results/04_Figures/shared/fgsea_cache"
BASE         <- "02-03_Sam's_Results/04_Figures/F05_recovery_reversal"
RPT_PDF      <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG      <- file.path(BASE, "b_reports", "main", "png")
PNL_PNG      <- file.path(RPT_PNG, "panels")
PNL_PDF      <- file.path(RPT_PDF, "panels")
DAT          <- file.path(BASE, "c_data")
SUPP_PNL_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
SUPP_PNL_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
for (d in c(RPT_PDF, RPT_PNG, PNL_PNG, PNL_PDF, DAT, SUPP_PNL_PNG, SUPP_PNL_PDF,
            file.path(DAT, "panel_A"), file.path(DAT, "panel_B_heatmap"),
            file.path(DAT, "panel_C_fry"), file.path(DAT, "panel_D"),
            file.path(DAT, "panel_E")))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ── Data loading ─────────────────────────────────────────────────────────────
cvh_dep <- read_csv(file.path(SAM_DEP_DIR, "Cancer_vs_Healthy.csv"),
                    show_col_types = FALSE)
tcr_dep <- read_csv(file.path(SAM_DEP_DIR, "Training_CR.csv"),
                    show_col_types = FALSE)

# Merge on gene: bring Cancer and Training columns together
dep_df <- cvh_dep |>
  select(gene, logFC_CvH = logFC, t_CvH = t, pi_CvH = pi_score) |>
  inner_join(
    tcr_dep |>
      select(gene, logFC_TCR = logFC, t_TCR = t, pi_TCR = pi_score),
    by = "gene"
  ) |>
  filter(!is.na(logFC_CvH), !is.na(logFC_TCR))

message(sprintf("Merged: %d proteins", nrow(dep_df)))

# ── Panel A — Quadrant ORA scatter + flanking bars ──────────────────────────
message("=== Panel A: Quadrant ORA ===")

N_SHOW     <- 5
COMP_RED   <- unname(DIR_COLORS["Up"])
COMP_BLUE  <- unname(DIR_COLORS["Down"])

# Significance classification (2-axis pi-score)
dep_df <- dep_df |>
  mutate(
    sig_class = case_when(
      pi_CvH < 0.05 & pi_TCR < 0.05 ~ "Sig Both",
      pi_CvH < 0.05                  ~ "Sig Cancer only",
      pi_TCR < 0.05                  ~ "Sig Training only",
      TRUE                           ~ "NS"
    ),
    is_sig   = sig_class != "NS",
    quadrant = case_when(
      logFC_CvH > 0 & logFC_TCR < 0 ~ "Reversed (Cancer Up / Training Down)",
      logFC_CvH < 0 & logFC_TCR > 0 ~ "Reversed (Cancer Down / Training Up)",
      logFC_CvH > 0 & logFC_TCR > 0 ~ "Exacerbated Up",
      TRUE                           ~ "Exacerbated Down"
    )
  )

# sig colors — adapt SIG_COLORS_F3 labels to CvH biology
SIG_COLS_F5 <- c(
  "Sig Both"          = unname(SIG_COLORS_F3["Sig Both"]),
  "Sig Cancer only"   = unname(SIG_COLORS_F3["Sig Aging only"]),
  "Sig Training only" = unname(SIG_COLORS_F3["Sig Training only"]),
  "NS"                = "grey75"
)

universe     <- dep_df$gene
pw_collection <- build_pathway_collection(min_size = 15, max_size = 500,
                                           include_goslim = FALSE,
                                           exclude_variants = TRUE)

run_set_ora <- function(genes, set_name) {
  if (length(genes) < 5) return(tibble())
  res <- tryCatch(
    run_ora_deduplicated(genes = genes, universe = universe,
                          pathways = pw_collection, jaccard_cutoff = 0.5,
                          min_size = 15, max_size = 500, padj_cutoff = 1),
    error = function(e) { message("  ORA error: ", e$message); tibble() })
  if (nrow(res) == 0) return(tibble())
  res |>
    mutate(set = set_name,
           pathway_label    = clean_pathway_name(pathway),
           neg_log10_padj   = -log10(padj + 1e-16),
           significant      = padj < 0.05) |>
    arrange(desc(neg_log10_padj)) |>
    slice_head(n = N_SHOW)
}

ora_tl <- run_set_ora(dep_df$gene[dep_df$quadrant == "Reversed (Cancer Down / Training Up)"],
                       "Reversed (Cancer Down / Training Up)")
ora_tr <- run_set_ora(dep_df$gene[dep_df$quadrant == "Exacerbated Up"],
                       "Exacerbated Up")
ora_bl <- run_set_ora(dep_df$gene[dep_df$quadrant == "Exacerbated Down"],
                       "Exacerbated Down")
ora_br <- run_set_ora(dep_df$gene[dep_df$quadrant == "Reversed (Cancer Up / Training Down)"],
                       "Reversed (Cancer Up / Training Down)")

all_quad_ora <- bind_rows(ora_tl, ora_tr, ora_bl, ora_br)
if (nrow(all_quad_ora) > 0)
  write_csv(all_quad_ora, file.path(DAT, "panel_A", "ora_quadrant.csv"))

xlim_range <- c(-3.1, 3.1)
ylim_range <- c(-2.8, 2.8)

ns_df  <- filter(dep_df, sig_class == "NS")
sig_df <- filter(dep_df, sig_class != "NS")

q_df     <- dep_df |>
  mutate(q = case_when(
    logFC_CvH > 0 & logFC_TCR < 0 ~ "BR",
    logFC_CvH < 0 & logFC_TCR > 0 ~ "TL",
    logFC_CvH > 0 & logFC_TCR > 0 ~ "TR",
    TRUE ~ "BL"))
q_counts <- q_df |> count(q) |> deframe()
q_sig    <- q_df |> filter(sig_class != "NS") |> count(q) |> deframe()
for (qq in c("BR","TL","TR","BL")) if (is.na(q_sig[qq])) q_sig[qq] <- 0

label_df <- sig_df |>
  group_by(sig_class) |>
  arrange(desc(abs(logFC_CvH) + abs(logFC_TCR))) |>
  slice_head(n = 5) |>
  ungroup()

txt_gene <- scale_text(BASE_GENE, 190) * 0.82 + 1
txt_quad <- scale_text(BASE_QUADRANT, 190) * 0.88

x_breaks <- seq(-3, 3, 1)
y_breaks <- seq(-2, 2, 1)
x_tick_df <- tibble(x = x_breaks[x_breaks != 0], y = 0,
                    label = as.character(x_breaks[x_breaks != 0]))
y_tick_df <- tibble(x = 0, y = y_breaks[y_breaks != 0],
                    label = as.character(y_breaks[y_breaks != 0]))

p_scatter <- ggplot(mapping = aes(x = logFC_CvH, y = logFC_TCR)) +
  annotate("rect", xmin = 0, xmax = Inf,  ymin = -Inf, ymax = 0,
           fill = "#DCEEFF", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#DCEEFF", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 0, xmax = Inf,  ymin = 0, ymax = Inf,
           fill = "#FFE0E0", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#FFE0E0", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  geom_text(data = x_tick_df, aes(x = x, y = y, label = label),
            vjust = 1.5, size = 1.3 * PRINT_SCALE, color = "grey40", fontface = "bold") +
  geom_text(data = y_tick_df, aes(x = x, y = y, label = label),
            hjust = -0.5, size = 1.3 * PRINT_SCALE, color = "grey40", fontface = "bold") +
  geom_point(data = ns_df, aes(x = logFC_CvH, y = logFC_TCR),
             color = "grey80", fill = "grey85", shape = 21,
             size = 0.35, alpha = 0.3, stroke = 0.10) +
  geom_point(data = sig_df, aes(x = logFC_CvH, y = logFC_TCR, fill = sig_class),
             shape = 21,
             size = ifelse(sig_df$sig_class == "NS", 0.6, 0.9),
             color = "grey75",
             alpha = ifelse(sig_df$sig_class == "Sig Both", 0.75, 0.85),
             stroke = 0.6) +
  scale_fill_manual(values = SIG_COLS_F5, name = "Significance") +
  geom_label_repel(data = label_df, aes(x = logFC_CvH, y = logFC_TCR, label = gene),
                   fill = "white", color = "grey15",
                   size = txt_gene, fontface = "italic", max.overlaps = 40,
                   segment.size = 0.3, segment.color = "grey25",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.3, point.padding = 0.3,
                   force = 6, force_pull = 0.3,
                   label.padding = unit(1, "pt"), label.r = unit(0.5, "pt"),
                   linewidth = 0.15, seed = 42,
                   xlim = c(-3, 3) * 0.85, ylim = c(-2.7, 2.7) * 0.85) +
  annotate("label", x = xlim_range[1], y = ylim_range[2],
           label = sprintf("Reversed (Cav Tr^)\n%s/%s", q_sig["TL"], q_counts["TL"]),
           hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
           color = COMP_BLUE, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt"), lineheight = 0.9) +
  annotate("label", x = xlim_range[2], y = ylim_range[2],
           label = sprintf("Exacerbated Up\n%s/%s", q_sig["TR"], q_counts["TR"]),
           hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
           color = COMP_RED, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt"), lineheight = 0.9) +
  annotate("label", x = xlim_range[1], y = ylim_range[1],
           label = sprintf("%s/%s\nExacerbated Down", q_sig["BL"], q_counts["BL"]),
           hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
           color = COMP_RED, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt"), lineheight = 0.9) +
  annotate("label", x = xlim_range[2], y = ylim_range[1],
           label = sprintf("%s/%s\nReversed (Ca^ Trv)", q_sig["BR"], q_counts["BR"]),
           hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
           color = COMP_BLUE, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt"), lineheight = 0.9) +
  annotate("text", x = 2.5, y = 0,
           label = expression(log[2]*FC ~ "(Cancer vs Healthy)"),
           hjust = 0.5, vjust = -0.4, size = 1.3 * PRINT_SCALE, color = "grey30",
           fontface = "bold") +
  annotate("text", x = 0, y = 2.0,
           label = expression(log[2]*FC ~ "(Training CR)"),
           hjust = 0.5, vjust = -0.4, size = 1.3 * PRINT_SCALE, color = "grey30",
           fontface = "bold", angle = 90) +
  coord_cartesian(xlim = xlim_range, ylim = ylim_range, expand = FALSE) +
  labs(x = NULL, y = NULL) +
  FIG_THEME +
  theme(plot.title       = element_blank(),
        plot.subtitle    = element_blank(),
        axis.text        = element_blank(),
        axis.ticks       = element_blank(),
        axis.title       = element_blank(),
        plot.margin      = margin(2, 0, 0, 0, "mm"),
        legend.position  = "none")

# Significance key
key_lvls    <- c("Sig Both", "Sig Cancer only", "Sig Training only")
key_display <- c("Sig Both", "Sig Cancer", "Sig Training")
key_df <- tibble(
  category = factor(key_lvls, levels = key_lvls),
  display  = key_display,
  fill_col = unname(SIG_COLS_F5[key_lvls]),
  x        = c(1.25, 1.90, 2.55),
  y        = 0
)
p_key <- ggplot(key_df, aes(x = x, y = y)) +
  geom_point(aes(fill = category), shape = 21, size = 2.5 * PRINT_SCALE,
             color = "grey50", stroke = 0.6, alpha = 0.85, show.legend = FALSE) +
  geom_text(aes(label = display), nudge_x = 0.06, hjust = 0,
            size = 2.0 * PRINT_SCALE, fontface = "bold", color = "grey25") +
  scale_fill_manual(values = setNames(key_df$fill_col, key_df$category)) +
  scale_x_continuous(limits = c(0.2, 4.0), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.15, 0.15), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(-24, 0, 0, 0, "mm"))

# Half-bar builder
make_half_bars <- function(df, fill_color, side, ylim) {
  bar_h  <- 0.42
  n_bars <- if (is.null(df) || nrow(df) == 0) 0L else min(nrow(df), 5L)
  if (n_bars == 0)
    return(ggplot() + theme_void() +
           scale_y_continuous(limits = ylim, expand = c(0, 0)))

  y_pos <- if (ylim[1] >= 0) {
    rev(seq(0.3, 2.3, length.out = 5))[seq_len(n_bars)]
  } else {
    seq(-0.3, -2.5, length.out = 5)[seq_len(n_bars)]
  }

  bars <- df |>
    arrange(desc(neg_log10_padj)) |>
    slice_head(n = 5) |>
    mutate(
      y        = y_pos,
      bar_fill = ifelse(significant, scales::alpha(fill_color, 0.85),
                        scales::alpha(fill_color, 0.30)),
      star     = sig_stars(padj)
    )
  x_max         <- max(bars$neg_log10_padj)
  x_display_max <- x_max * 1.18
  is_upper      <- ylim[1] >= 0
  brk_fn        <- function(limits) { b <- scales::pretty_breaks(n = 3)(limits); b[b != 0] }

  bars <- bars |>
    mutate(
      label_inside = neg_log10_padj >= x_max * 0.10,
      label_x      = ifelse(label_inside, neg_log10_padj * 0.5,
                            neg_log10_padj + x_max * 0.03),
      label_hjust  = ifelse(label_inside, 0.5, 0),
      label_color  = ifelse(label_inside,
                            ifelse(significant, "white", "grey15"), "grey20"),
      text_size    = scale_text(BASE_PATHWAY, 190) * 0.80
    )

  star_x_mult <- if (side == "left") 0.12 else 0.035

  p <- ggplot(bars, aes(y = y)) +
    geom_rect(aes(xmin = 0, xmax = neg_log10_padj,
                  ymin = y - bar_h / 2, ymax = y + bar_h / 2),
              fill = bars$bar_fill, color = "black", linewidth = 0.3) +
    geom_text(aes(x = label_x, y = y, label = pathway_label),
              hjust = bars$label_hjust, size = bars$text_size,
              fontface = "bold", color = bars$label_color, lineheight = 0.85) +
    geom_text(aes(x = neg_log10_padj + x_max * star_x_mult, label = star),
              hjust = 0, vjust = 0.5,
              size = 2.2 * PRINT_SCALE, fontface = "bold", color = "black") +
    labs(x = if (!is_upper) expression(-log[10](p[adj])) else NULL, y = NULL) +
    theme_minimal(base_size = 9) +
    theme(panel.grid   = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x  = element_text(size = FIG_AXIS_TEXT, face = "bold",
                                       margin = margin(t = 0, unit = "mm")),
          axis.title.x = if (!is_upper)
            element_text(size = 5 * PRINT_SCALE, face = "bold",
                         margin = margin(t = 0, unit = "mm"))
            else element_blank(),
          axis.line.x  = element_line(color = "grey50", linewidth = 0.3),
          axis.ticks.x = element_line(color = "grey50", linewidth = 0.3),
          plot.margin  = if (is_upper && side == "left") margin(4, 0, 0, 3, "mm")
                         else if (is_upper) margin(4, 3, 0, 0, "mm")
                         else if (side == "left") margin(2, 0, 0, 3, "mm")
                         else margin(2, 3, 0, 0, "mm"))

  if (side == "left") {
    p + scale_x_reverse(limits = c(x_display_max, 0),
                         breaks = brk_fn, expand = expansion(mult = c(0, 0))) +
        scale_y_continuous(limits = ylim, expand = c(0, 0)) +
        coord_cartesian(clip = "off")
  } else {
    p + scale_x_continuous(limits = c(0, x_display_max),
                            breaks = brk_fn, expand = expansion(mult = c(0, 0))) +
        scale_y_continuous(limits = ylim, expand = c(0, 0)) +
        coord_cartesian(clip = "off")
  }
}

p_ul <- make_half_bars(ora_tl, scales::alpha(COMP_BLUE, 0.30), "left",  c(0, 2.8))
p_ll <- make_half_bars(ora_bl, scales::alpha(COMP_RED, 0.30),  "left",  c(-2.8, 0))
p_ur <- make_half_bars(ora_tr, scales::alpha(COMP_RED, 0.30),  "right", c(0, 2.8))
p_lr <- make_half_bars(ora_br, scales::alpha(COMP_BLUE, 0.30), "right", c(-2.8, 0))

design_A <- c(
  area(1, 1), area(1, 2, 2, 2), area(1, 3),
  area(2, 1), area(2, 3),       area(3, 1, 3, 3)
)
n_total_A  <- nrow(dep_df)
n_sig_A    <- sum(dep_df$is_sig)
n_enrich_A <- if (nrow(all_quad_ora) > 0) sum(all_quad_ora$significant) else 0L
r_spear_A  <- cor(dep_df$logFC_CvH, dep_df$logFC_TCR, use = "complete.obs",
                  method = "spearman")

composite_A <- p_ul + p_scatter + p_ur + p_ll + p_lr + p_key +
  plot_layout(design = design_A,
              widths  = c(70, 100, 70) / 240,
              heights = c(85, 85, 8) / 178) +
  plot_annotation(
    title    = "Recovery Reversal: Quadrant ORA",
    subtitle = sprintf(
      "Threshold-free ORA (hypergeometric) | N = %d | %d DEPs (Pi < 0.05) | %d enriched (FDR < 0.05) | rho(all) = %.2f",
      n_total_A, n_sig_A, n_enrich_A, r_spear_A),
    theme = theme(plot.title    = element_text(size = FIG_TITLE_SIZE, face = "bold"),
                  plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, hjust = 0, color = "grey30")))

png(file.path(PNL_PNG, "MAIN_panel_A_ORA_composite.png"),
    width = 200, height = 120, units = "mm", res = 300)
print(composite_A)
dev.off()
pdf(file.path(PNL_PDF, "MAIN_panel_A_ORA_composite.pdf"),
    width = 200 / 25.4, height = 120 / 25.4)
print(composite_A)
dev.off()
message("Panel A done")

# Strip title/subtitle/legend for final composite
composite_A <- composite_A &
  labs(title = NULL, subtitle = NULL, tag = NULL) &
  theme(legend.position = "none")
composite_A <- composite_A +
  plot_annotation(title = NULL, subtitle = NULL,
                  theme = theme(plot.title = element_blank(), plot.subtitle = element_blank()))

# ── Panel D — fGSEA NES scatter (pathway-level reversal) ────────────────────
message("=== Panel D: NES scatter ===")

cvh_fgsea <- readRDS(file.path(FGSEA_CACHE, "Cancer_vs_Healthy_fgsea.rds")) |>
  filter(!is.na(NES))
tcr_fgsea <- readRDS(file.path(FGSEA_CACHE, "Training_CR_fgsea.rds")) |>
  filter(!is.na(NES))

fgsea_wide <- cvh_fgsea |>
  select(pathway, NES_CvH = NES, padj_CvH = padj) |>
  inner_join(
    tcr_fgsea |> select(pathway, NES_TCR = NES, padj_TCR = padj),
    by = "pathway"
  ) |>
  mutate(
    sig_CvH   = !is.na(padj_CvH) & padj_CvH < 0.05,
    sig_TCR   = !is.na(padj_TCR) & padj_TCR < 0.05,
    sig_class = case_when(
      sig_CvH & sig_TCR ~ "Sig Both",
      sig_CvH           ~ "Sig Cancer only",
      sig_TCR           ~ "Sig Training only",
      TRUE              ~ "NS"
    ),
    direction = case_when(
      sig_CvH & sig_TCR & (NES_CvH * NES_TCR < 0) ~ "Reversed",
      sig_CvH & sig_TCR & (NES_CvH * NES_TCR > 0) ~ "Exacerbated",
      TRUE ~ "NS"
    )
  )

write_csv(fgsea_wide, file.path(DAT, "panel_D", "nes_scatter.csv"))

# Spearman + Fisher Z CI for all / sig-both pathways
nes_cor_all <- cor.test(fgsea_wide$NES_CvH, fgsea_wide$NES_TCR,
                         method = "spearman", exact = FALSE)
rho_D     <- as.numeric(nes_cor_all$estimate)
n_pw_D    <- nrow(fgsea_wide)
n_sig_pw_D <- sum(fgsea_wide$sig_CvH | fgsea_wide$sig_TCR)

sig_both  <- filter(fgsea_wide, sig_CvH & sig_TCR)
pw_rev_D  <- if (nrow(sig_both) > 0) mean(sig_both$direction == "Reversed") else 0

# Labels: top sig-both pathways by distance from diagonal
label_pw <- sig_both |>
  mutate(score = abs(NES_CvH) + abs(NES_TCR)) |>
  arrange(desc(score)) |>
  slice_head(n = 10) |>
  mutate(label = clean_pathway_name(pathway))

pD <- ggplot(fgsea_wide, aes(x = NES_CvH, y = NES_TCR)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "#DCEEFF", alpha = 0.45) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#DCEEFF", alpha = 0.45) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "#FFE0E0", alpha = 0.45) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#FFE0E0", alpha = 0.45) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_abline(slope = -1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  geom_point(data = filter(fgsea_wide, sig_class == "NS"),
             color = "grey80", fill = "grey85", shape = 21,
             size = 0.5, alpha = 0.3, stroke = 0.1) +
  geom_point(data = filter(fgsea_wide, sig_class != "NS"),
             aes(fill = sig_class), shape = 21,
             size = 1.2, alpha = 0.8, stroke = 0.4) +
  scale_fill_manual(values = c(
    "Sig Both"          = unname(SIG_COLORS_F3["Sig Both"]),
    "Sig Cancer only"   = unname(SIG_COLORS_F3["Sig Aging only"]),
    "Sig Training only" = unname(SIG_COLORS_F3["Sig Training only"])
  ), name = "Significance") +
  geom_label_repel(data = label_pw, aes(x = NES_CvH, y = NES_TCR, label = label),
                   size = scale_text(BASE_PATHWAY, 146) * 0.80,
                   fill = alpha("white", 0.9), color = "grey15",
                   fontface = "bold", max.overlaps = 30,
                   segment.size = 0.25, box.padding = 0.2, seed = 42,
                   label.padding = unit(1, "pt"), label.r = unit(0.5, "pt"),
                   linewidth = 0.15) +
  labs(title    = "Pathway Reversal (fGSEA)",
       subtitle = sprintf("rho = %.2f | %.0f%% reversed (sig-both)", rho_D, pw_rev_D * 100),
       x        = "NES (Cancer vs Healthy)",
       y        = "NES (Training CR)") +
  FIG_THEME +
  theme(legend.position = "none")

ggsave(file.path(PNL_PNG, "MAIN_panel_D_nes_scatter.png"), pD,
       width = 80, height = 80, units = "mm", dpi = 300)
ggsave(file.path(PNL_PDF, "MAIN_panel_D_nes_scatter.pdf"), pD,
       width = 80, height = 80, units = "mm", device = pdf_device)
message("Panel D done")

# ── Panel C — fry rotation test ──────────────────────────────────────────────
message("=== Panel C: fry ===")

# Cancer-significant DEP sets
cancer_up <- dep_df |> filter(pi_CvH < 0.05 & logFC_CvH > 0) |> pull(gene)
cancer_dn <- dep_df |> filter(pi_CvH < 0.05 & logFC_CvH < 0) |> pull(gene)

# t-stats for Training_CR (for fry input)
t_TCR <- setNames(dep_df$t_TCR, dep_df$gene)
t_TCR <- t_TCR[!is.na(t_TCR)]

message(sprintf("  Cancer-Up DEPs: %d | Cancer-Down DEPs: %d", length(cancer_up), length(cancer_dn)))

# fry rotation test
fry_sets <- list(
  cancer_up = cancer_up[cancer_up %in% names(t_TCR)],
  cancer_dn = cancer_dn[cancer_dn %in% names(t_TCR)]
)
fry_sets <- fry_sets[lengths(fry_sets) >= 5]

fry_res <- tryCatch({
  requireNamespace("limma", quietly = TRUE)
  t_vec   <- t_TCR[!is.na(t_TCR)]
  t_mat   <- matrix(t_vec, ncol = 1, dimnames = list(names(t_vec), "t"))
  # Build index list for limma::fry
  idx_list <- lapply(fry_sets, function(g) which(rownames(t_mat) %in% g))
  idx_list <- idx_list[lengths(idx_list) >= 5]
  if (length(idx_list) < 1) stop("No valid sets for fry")
  limma::fry(t_mat, index = idx_list, geneid = rownames(t_mat))
}, error = function(e) {
  message("  fry error: ", e$message, " — building stub from t-stats")
  NULL
})

# Build barcode plot — ranked t-stat plot for cancer-up and cancer-down sets
t_df <- tibble(gene = names(t_TCR), t = as.numeric(t_TCR)) |>
  arrange(t) |>
  mutate(rank = seq_len(n()),
         set  = case_when(gene %in% cancer_up ~ "cancer_up",
                          gene %in% cancer_dn ~ "cancer_dn",
                          TRUE ~ "NS"))

n_all_C   <- nrow(t_df)
n_up_C    <- length(cancer_up)
n_dn_C    <- length(cancer_dn)

# Enrichment score for annotation
ks_up <- if (n_up_C >= 5) {
  ks.test(t_df$rank[t_df$set == "cancer_up"], t_df$rank[t_df$set == "NS"])
} else { NULL }
ks_dn <- if (n_dn_C >= 5) {
  ks.test(t_df$rank[t_df$set == "cancer_dn"], t_df$rank[t_df$set == "NS"])
} else { NULL }

# fry result text
if (!is.null(fry_res)) {
  write_csv(as.data.frame(fry_res) |> tibble::rownames_to_column("set"),
            file.path(DAT, "panel_C_fry", "fry_results_all.csv"))
  fry_up_p <- if ("cancer_up" %in% rownames(fry_res)) fry_res["cancer_up", "PValue"] else NA
  fry_dn_p <- if ("cancer_dn" %in% rownames(fry_res)) fry_res["cancer_dn", "PValue"] else NA
  fry_up_dir <- if ("cancer_up" %in% rownames(fry_res)) fry_res["cancer_up", "Direction"] else NA
  fry_dn_dir <- if ("cancer_dn" %in% rownames(fry_res)) fry_res["cancer_dn", "Direction"] else NA
} else {
  fry_up_p <- fry_dn_p <- fry_up_dir <- fry_dn_dir <- NA
}

make_barcode <- function(set_label, set_name, color, subtitle_extra = "") {
  set_genes <- if (set_name == "cancer_up") cancer_up else cancer_dn
  n_set     <- length(set_genes)
  df_set    <- t_df |> filter(gene %in% set_genes)

  # Running enrichment score (simple GSEA-style)
  all_ranks <- t_df$rank
  hit_ranks <- sort(df_set$rank)
  n_total   <- length(all_ranks)
  n_hit     <- length(hit_ranks)
  if (n_hit < 1) return(ggplot() + theme_void())

  step_up   <- 1 / n_hit
  step_dn   <- 1 / (n_total - n_hit)
  es_x      <- c(0, rep(hit_ranks, each = 2), n_total)
  es_y      <- numeric(length(es_x))
  cur_es    <- 0
  prev_rank <- 0
  k         <- 1
  for (i in seq_along(hit_ranks)) {
    gap_steps <- hit_ranks[i] - prev_rank - 1
    if (gap_steps > 0) {
      cur_es   <- cur_es - step_dn * gap_steps
      k        <- k + 1
      es_x[k]  <- hit_ranks[i] - 0.5
      es_y[k]  <- cur_es
      k        <- k + 1
    }
    cur_es    <- cur_es + step_up
    k         <- k + 1
    es_x[k]   <- hit_ranks[i]
    es_y[k]   <- cur_es
    prev_rank <- hit_ranks[i]
  }
  if (prev_rank < n_total) {
    cur_es <- cur_es - step_dn * (n_total - prev_rank)
    k <- k + 1
    es_x[k] <- n_total
    es_y[k] <- cur_es
  }
  es_df <- tibble(x = es_x[1:k], y = es_y[1:k])
  max_es <- max(abs(es_df$y))

  tick_df <- tibble(x = hit_ranks, y = 0)

  stat_text <- paste0(
    set_label, " (n=", n_hit, ")",
    if (!is.na(fry_up_p) && set_name == "cancer_up")
      sprintf("\nfry p=%.3f (%s)", fry_up_p, fry_up_dir)
    else if (!is.na(fry_dn_p) && set_name == "cancer_dn")
      sprintf("\nfry p=%.3f (%s)", fry_dn_p, fry_dn_dir)
    else ""
  )

  ggplot() +
    geom_line(data = es_df, aes(x = x, y = y), color = color, linewidth = 0.5) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
    geom_segment(data = tick_df, aes(x = x, xend = x, y = -0.02, yend = 0.02),
                 color = color, linewidth = 0.2, alpha = 0.6) +
    annotate("text", x = n_total * 0.02, y = max_es * 0.85,
             label = stat_text,
             hjust = 0, vjust = 1, size = scale_text(BASE_STAT, 178) * 0.85,
             fontface = "bold", color = color) +
    scale_x_continuous(limits = c(1, n_total), expand = c(0, 0),
                       labels = NULL, breaks = NULL) +
    labs(x = "Proteins ranked by Training CR t-stat (Low -> High)",
         y = "Enrichment score") +
    FIG_THEME
}

pC_up <- make_barcode("Cancer-Up DEPs", "cancer_up", COMP_RED)
pC_dn <- make_barcode("Cancer-Down DEPs", "cancer_dn", COMP_BLUE)
pC_fry <- pC_up / pC_dn +
  plot_annotation(
    title    = "fry: Cancer Reversal",
    subtitle = sprintf("Cancer DEP sets vs Training_CR t-stats | n = %d proteins", n_all_C),
    theme    = theme(plot.title    = element_text(size = FIG_TITLE_SIZE, face = "bold"),
                     plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, color = "grey30")))

ggsave(file.path(PNL_PNG, "MAIN_panel_C_fry.png"), pC_fry,
       width = 130, height = 80, units = "mm", dpi = 300)
ggsave(file.path(PNL_PDF, "MAIN_panel_C_fry.pdf"), pC_fry,
       width = 130, height = 80, units = "mm", device = pdf_device)

# Save driving proteins
if (nrow(filter(dep_df, gene %in% cancer_up | gene %in% cancer_dn)) > 0) {
  driving_df <- dep_df |>
    filter(gene %in% c(cancer_up, cancer_dn)) |>
    mutate(set = if_else(gene %in% cancer_up, "cancer_up", "cancer_dn"),
           is_driving = (set == "cancer_up" & t_TCR < 0) |
                        (set == "cancer_dn" & t_TCR > 0)) |>
    filter(is_driving) |>
    select(gene, set, logFC_CvH, logFC_TCR, t_TCR)
  write_csv(driving_df, file.path(DAT, "panel_C_fry", "driving_proteins.csv"))
}

n_all_C  <- n_all_C
cor_imp_C <- cor(dep_df$logFC_CvH, dep_df$logFC_TCR, use = "complete.obs")
message("Panel C done")

# ── Panel E — RRHO2 ──────────────────────────────────────────────────────────
message("=== Panel E: RRHO2 ===")

# Build input data frame for RRHO2: gene, t_CvH, t_TCR
rrho_df <- dep_df |>
  filter(!is.na(t_CvH), !is.na(t_TCR)) |>
  transmute(gene, t_CvH, t_TCR)

n_shared <- nrow(rrho_df)
message(sprintf("  RRHO2 input: %d proteins", n_shared))

rrho_dir <- file.path(DAT, "panel_E")
dir.create(rrho_dir, recursive = TRUE, showWarnings = FALSE)

pE_heat <- tryCatch({
  # RRHO2 requires two data frames: gene + numeric rank score
  df1 <- data.frame(gene = rrho_df$gene, rank = rrho_df$t_CvH)
  df2 <- data.frame(gene = rrho_df$gene, rank = rrho_df$t_TCR)

  rrho_obj <- RRHO2_initialize(df1, df2,
                                labels = c("Cancer vs Healthy", "Training CR"),
                                boundary = 0.05,
                                log10.ind = TRUE,
                                method = "hyper")

  # Extract the overlap matrix
  mat <- rrho_obj$hypermat

  # Summarize quadrant maxima
  nr <- nrow(mat); nc <- ncol(mat)
  h1  <- floor(nr / 2); h2 <- floor(nc / 2)
  max_UU <- max(mat[1:h1, 1:h2], na.rm = TRUE)
  max_DD <- max(mat[(h1+1):nr, (h2+1):nc], na.rm = TRUE)
  max_UD <- max(mat[1:h1, (h2+1):nc], na.rm = TRUE)
  max_DU <- max(mat[(h1+1):nr, 1:h2], na.rm = TRUE)

  # Hotspot genes
  get_hotspot <- function(i_range, j_range) {
    sub <- mat[i_range, j_range]
    idx <- which(sub == max(sub, na.rm = TRUE), arr.ind = TRUE)
    if (nrow(idx) == 0) return(character(0))
    # Reconstruct gene sets corresponding to row/col of hotspot
    # Upper-left of RRHO2 = Up in both (sorted descending)
    character(0)  # genes extracted separately below
  }

  rrho_summary <- tibble(
    quadrant    = c("UU", "UD", "DU", "DD"),
    label       = c("Exacerbated Up", "Reversed (CvH Up / Tr Down)",
                    "Reversed (CvH Down / Tr Up)", "Exacerbated Down"),
    max_overlap = c(max_UU, max_UD, max_DU, max_DD)
  )
  write_csv(rrho_summary, file.path(rrho_dir, "rrho2_summary.csv"))

  # Plot the heatmap
  mat_df <- as.data.frame(mat) |>
    tibble::rownames_to_column("row_idx") |>
    pivot_longer(-row_idx, names_to = "col_idx", values_to = "log10p") |>
    mutate(row_i = as.integer(row_idx),
           col_j = as.integer(str_extract(col_idx, "[0-9]+")))

  ggplot(mat_df, aes(x = col_j, y = row_i, fill = log10p)) +
    geom_raster() +
    scale_fill_gradientn(
      colors = c("#313695", "#4575B4", "#74ADD1", "#ABD9E9",
                 "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61",
                 "#F46D43", "#D73027", "#A50026"),
      name = "-log10(p)",
      na.value = "grey95"
    ) +
    geom_hline(yintercept = nr / 2, color = "grey30", linewidth = 0.4,
               linetype = "dashed") +
    geom_vline(xintercept = nc / 2, color = "grey30", linewidth = 0.4,
               linetype = "dashed") +
    annotate("text", x = nc * 0.25, y = nr * 0.12,
             label = sprintf("Exacerbated Up\n%.1f", max_UU),
             color = "white", size = scale_text(BASE_STAT, 80) * 0.85, fontface = "bold") +
    annotate("text", x = nc * 0.75, y = nr * 0.12,
             label = sprintf("Reversed\n(Ca^ Trv)\n%.1f", max_UD),
             color = "white", size = scale_text(BASE_STAT, 80) * 0.85, fontface = "bold") +
    annotate("text", x = nc * 0.25, y = nr * 0.88,
             label = sprintf("Reversed\n(Cav Tr^)\n%.1f", max_DU),
             color = "white", size = scale_text(BASE_STAT, 80) * 0.85, fontface = "bold") +
    annotate("text", x = nc * 0.75, y = nr * 0.88,
             label = sprintf("Exacerbated Down\n%.1f", max_DD),
             color = "white", size = scale_text(BASE_STAT, 80) * 0.85, fontface = "bold") +
    scale_x_continuous(name = expression("Cancer vs Healthy rank"~(Up %->% Down)),
                       expand = c(0, 0)) +
    scale_y_continuous(name = expression("Training CR rank"~(Up %->% Down)),
                       expand = c(0, 0)) +
    labs(title    = "RRHO2: Threshold-Free Reversal",
         subtitle = sprintf("%d proteins | warm off-diagonal = training reverses cancer signature",
                            n_shared)) +
    FIG_THEME +
    theme(axis.text     = element_blank(),
          axis.ticks    = element_blank(),
          legend.key.height = unit(3, "mm"),
          legend.title  = element_text(size = FIG_LEGEND_TITLE, face = "bold"),
          legend.text   = element_text(size = FIG_LEGEND_TEXT))

}, error = function(e) {
  message("  RRHO2 error: ", e$message, " — building stub panel")
  rrho_summary <- tibble(quadrant = character(), max_overlap = numeric())
  write_csv(rrho_summary, file.path(rrho_dir, "rrho2_summary.csv"))

  # Count quadrant membership using t-stat rank
  dep_df2 <- dep_df |> filter(!is.na(t_CvH), !is.na(t_TCR)) |>
    mutate(
      q = case_when(
        t_CvH > 0 & t_TCR < 0 ~ "Reversed (Ca^ Trv)",
        t_CvH < 0 & t_TCR > 0 ~ "Reversed (Cav Tr^)",
        t_CvH > 0 & t_TCR > 0 ~ "Exacerbated Up",
        TRUE                   ~ "Exacerbated Down"
      ))
  q_cnt <- dep_df2 |> count(q)

  ggplot(q_cnt, aes(x = reorder(q, n), y = n, fill = q)) +
    geom_col(color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c(
      "Exacerbated Up"              = "#D6604D",
      "Exacerbated Down"            = "#F4A261",
      "Reversed (Ca^ Trv)" = "#4393C3",
      "Reversed (Cav Tr^)" = "#74ADD1"
    )) +
    geom_text(aes(label = n), hjust = -0.2, size = scale_text(BASE_STAT, 80)) +
    coord_flip() +
    labs(title    = "RRHO2 Stub — Quadrant Protein Counts",
         subtitle = sprintf("RRHO2 pending | %d proteins by t-stat quadrant", n_shared),
         x = NULL, y = "Protein count") +
    FIG_THEME +
    theme(legend.position = "none")
})

max_UD <- max_DU <- n_UD <- n_DU <- NA_real_

# Attempt to read summary if RRHO2 succeeded
rrho_sum_path <- file.path(rrho_dir, "rrho2_summary.csv")
if (file.exists(rrho_sum_path)) {
  rrho_sum <- read_csv(rrho_sum_path, show_col_types = FALSE)
  if (nrow(rrho_sum) > 0 && "max_overlap" %in% names(rrho_sum)) {
    max_UD <- rrho_sum$max_overlap[rrho_sum$quadrant == "UD"]
    max_DU <- rrho_sum$max_overlap[rrho_sum$quadrant == "DU"]
    if (length(max_UD) == 0) max_UD <- NA_real_
    if (length(max_DU) == 0) max_DU <- NA_real_
  }
}
n_shared_E <- n_shared
n_rev_E    <- if (!is.na(max_UD) && !is.na(max_DU)) max(max_UD, max_DU) else 0

ggsave(file.path(PNL_PNG, "MAIN_panel_E_rrho2.png"), pE_heat,
       width = 90, height = 90, units = "mm", dpi = 300)
ggsave(file.path(PNL_PDF, "MAIN_panel_E_rrho2.pdf"), pE_heat,
       width = 90, height = 90, units = "mm", device = pdf_device)
message("Panel E done")

# ── Panel B — Pattern heatmap ────────────────────────────────────────────────
message("=== Panel B: Pattern heatmap ===")

ROW_H <- 0.078

# Classify proteins into reversal patterns
heatmap_df <- dep_df |>
  filter(!is.na(logFC_CvH), !is.na(logFC_TCR)) |>
  filter(pi_CvH < 0.05 | pi_TCR < 0.05) |>
  mutate(
    quadrant = case_when(
      logFC_CvH > 0 & logFC_TCR < 0 ~ "Reversed Up",
      logFC_CvH < 0 & logFC_TCR > 0 ~ "Reversed Down",
      TRUE                            ~ "Non-reversed"
    ),
    sig_cat = case_when(
      pi_CvH < 0.05 & pi_TCR < 0.05 ~ "Both",
      pi_CvH < 0.05                  ~ "Cancer",
      pi_TCR < 0.05                  ~ "Tr.(CR)",
      TRUE                           ~ "NS"
    )
  ) |>
  arrange(match(quadrant, c("Reversed Up", "Reversed Down", "Non-reversed")),
          desc(logFC_CvH))

# GO Slim annotation for heatmap rows (via pathway_utils classify_pathway_func)
# Use consolidated pathway categories as a proxy for functional annotation
write_csv(heatmap_df |> select(gene, quadrant, sig_cat, logFC_CvH, logFC_TCR, pi_CvH, pi_TCR),
          file.path(DAT, "panel_B_heatmap", "pattern_classification.csv"))

n_total <- nrow(heatmap_df)
n_rev_up  <- sum(heatmap_df$quadrant == "Reversed Up")
n_rev_dn  <- sum(heatmap_df$quadrant == "Reversed Down")
n_non_rev <- sum(heatmap_df$quadrant == "Non-reversed")

QUAD_COLORS <- c("Reversed Up" = "#B2182B", "Reversed Down" = "#2166AC",
                 "Non-reversed" = "#1B7837")
QUAD_BG     <- c("Reversed Up" = "#F4D9D2", "Reversed Down" = "#D5DEEF",
                 "Non-reversed" = "#C8E0CD")
SIG_COLORS_B <- c("Both" = "#2E7D32", "Cancer" = "#E05A4E",
                   "Tr.(CR)" = "#5DA5DA", "NS" = "grey70")

# Simple dot-heatmap: sorted proteins x 2 contrasts, colored by direction
n_pw <- 3  # Reversed Up, Reversed Down, Non-reversed (used in subtitle)

Y_MAX <- n_total * ROW_H

heatmap_df2 <- heatmap_df |>
  mutate(
    row_idx  = seq_len(n()),
    y        = n_total - row_idx + 1,
    lfc_fill_CvH = pmax(-3, pmin(3, logFC_CvH)),
    lfc_fill_TCR = pmax(-3, pmin(3, logFC_TCR)),
    bg_color = QUAD_BG[quadrant]
  )

# Color scale: diverging
make_lfc_color <- function(lfc) {
  cols <- colorRampPalette(c("#08306B", "#4393C3", "white", "#D6604D", "#67000D"))(101)
  idx  <- round((pmax(-3, pmin(3, lfc)) + 3) / 6 * 100) + 1
  cols[idx]
}

heatmap_df2 <- heatmap_df2 |>
  mutate(fill_CvH = make_lfc_color(lfc_fill_CvH),
         fill_TCR = make_lfc_color(lfc_fill_TCR))

x_CvH <- 0.5; x_TCR <- 1.5
x_bar_start <- 2.3
X_BAR_MAX   <- x_bar_start + max(abs(dep_df$logFC_CvH), abs(dep_df$logFC_TCR), na.rm = TRUE) * 0.4
BAR_YMAX    <- max(heatmap_df2$y) + ROW_H

# Build heatmap as geom_tile (two columns side by side)
tile_df <- bind_rows(
  heatmap_df2 |> transmute(gene, y, x = x_CvH, lfc = lfc_fill_CvH, contrast = "Cancer"),
  heatmap_df2 |> transmute(gene, y, x = x_TCR, lfc = lfc_fill_TCR, contrast = "Training CR")
)

sig_dot_df <- heatmap_df2 |>
  mutate(x_sig = 2.0, sig_fill = SIG_COLORS_B[sig_cat])

pB <- ggplot() +
  # Quadrant background (use pre-computed fill column to avoid scale conflict)
  geom_rect(data = heatmap_df2,
            aes(xmin = -0.2, xmax = X_BAR_MAX + 1.5,
                ymin = y - ROW_H / 2, ymax = y + ROW_H / 2),
            fill = heatmap_df2$bg_color, color = NA) +
  # Heatmap tiles with continuous fill
  geom_tile(data = tile_df,
            aes(x = x, y = y, fill = lfc),
            width = 0.6, height = ROW_H * 0.85, color = NA) +
  scale_fill_gradientn(
    colors = c("#08306B", "#4393C3", "white", "#D6604D", "#67000D"),
    values = scales::rescale(c(-3, -1, 0, 1, 3)),
    limits = c(-3, 3), oob = scales::squish,
    name = expression(log[2]*FC)) +
  # Sig dots (use color not fill to avoid scale conflict)
  geom_point(data = sig_dot_df,
             aes(x = x_sig, y = y),
             color = sig_dot_df$sig_fill, size = 0.3, shape = 16, alpha = 0.8) +
  # Quadrant separation lines
  geom_hline(yintercept = n_rev_up + 0.5, color = "grey30", linewidth = 0.4) +
  geom_hline(yintercept = n_rev_dn + 0.5, color = "grey30", linewidth = 0.4) +
  # Column headers
  annotate("text", x = x_CvH, y = BAR_YMAX + ROW_H * 1,
           label = "Cancer", hjust = 0.5, size = scale_text(BASE_STAT, 90),
           fontface = "bold", color = unname(CONTRAST_COLORS["Cancer_vs_Healthy"])) +
  annotate("text", x = x_TCR, y = BAR_YMAX + ROW_H * 1,
           label = "Tr.(CR)", hjust = 0.5, size = scale_text(BASE_STAT, 90),
           fontface = "bold", color = unname(CONTRAST_COLORS["Training_CR"])) +
  # Quadrant count labels
  annotate("text", x = 0.5, y = n_total - n_rev_up / 2 + 0.5,
           label = sprintf("Rev. Up\n(n=%d)", n_rev_up),
           hjust = 0, size = scale_text(BASE_STAT, 90) * 0.85,
           color = unname(QUAD_COLORS["Reversed Up"]), fontface = "bold") +
  scale_x_continuous(limits = c(-0.25, X_BAR_MAX + 1.75), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-ROW_H, BAR_YMAX + ROW_H * 3), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(title    = "Protein-to-Pathway",
       subtitle = sprintf("%d proteins | %d patterns", n_total, n_pw)) +
  FIG_THEME +
  theme(axis.text    = element_blank(),
        axis.ticks   = element_blank(),
        axis.title   = element_blank(),
        panel.border = element_blank(),
        panel.grid   = element_blank(),
        plot.margin  = margin(1, -28, 8, -14, "mm"))

ggsave(file.path(PNL_PNG, "MAIN_panel_B_heatmap.png"), pB,
       width = 80, height = 120, units = "mm", dpi = 300)
ggsave(file.path(PNL_PDF, "MAIN_panel_B_heatmap.pdf"), pB,
       width = 80, height = 120, units = "mm", device = pdf_device)
message("Panel B done")

# Restore paths
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf")
RPT_PNG <- file.path(BASE, "b_reports", "main", "png")

# ── Quadrant legend (for composite) ─────────────────────────────────────────
inset_quad_df <- tibble(
  quadrant  = factor(c("Reversed Up", "Reversed Down", "Non-reversed"),
                     levels = c("Reversed Up", "Reversed Down", "Non-reversed")),
  bg_color  = unname(QUAD_BG[c("Reversed Up", "Reversed Down", "Non-reversed")]),
  bar_color = unname(QUAD_COLORS[c("Reversed Up", "Reversed Down", "Non-reversed")])
)
nudge_idx3 <- 0.15
quad_legend <- ggplot(inset_quad_df) +
  geom_rect(aes(xmin = (as.integer(quadrant) - 1) * 3.5 +
                       (as.integer(quadrant) == 3) * nudge_idx3,
                xmax = (as.integer(quadrant) - 1) * 3.5 + 0.7 +
                       (as.integer(quadrant) == 3) * nudge_idx3,
                ymin = -0.35, ymax = 0.35),
            fill = inset_quad_df$bg_color, color = "black", linewidth = 0.5) +
  geom_rect(aes(xmin = (as.integer(quadrant) - 1) * 3.5 + 0.10 +
                       (as.integer(quadrant) == 3) * nudge_idx3,
                xmax = (as.integer(quadrant) - 1) * 3.5 + 0.60 +
                       (as.integer(quadrant) == 3) * nudge_idx3,
                ymin = -0.15, ymax = 0.15),
            fill = inset_quad_df$bar_color, color = "black", linewidth = 0.3) +
  geom_text(aes(x = (as.integer(quadrant) - 1) * 3.5 + 0.85 +
                    (as.integer(quadrant) == 3) * nudge_idx3,
                y = 0, label = as.character(quadrant)),
            hjust = 0, size = 3.5, fontface = "bold", color = "grey15") +
  coord_cartesian(xlim = c(0, 10.5), ylim = c(-0.7, 0.7), clip = "off") +
  theme_void() +
  theme(plot.background = element_blank(), panel.background = element_blank(),
        plot.margin = margin(0, 0, 0, 0, "mm"))

# ── Composite ────────────────────────────────────────────────────────────────
COMP_W      <- 420
COMP_H      <- 310
PRINT_SCALE2 <- 380 / 178
TAG_SZ      <- round(10 * PRINT_SCALE2 * 0.85)
TTL_SZ      <- round(10 * PRINT_SCALE2 * 0.85)
SUB_SZ      <- round(7  * PRINT_SCALE2 * 0.85)

ttl_A <- "Quadrant ORA (Reversal)"
sub_A <- sprintf("N = %d | %d DEPs (Pi) | %d enriched (FDR) | rho = %.2f",
                 n_total_A, n_sig_A, n_enrich_A, r_spear_A)
ttl_B <- "Protein-to-Pathway"
sub_B <- sprintf("%d proteins | %d patterns", n_total, n_pw)
ttl_C <- "fry: Reversal"
sub_C <- sprintf("n = %d | r = %.3f", n_all_C, cor_imp_C)
ttl_D <- "Pathway Reversal"
sub_D <- sprintf("rho = %.2f | %.0f%% reversed (sig-both)", rho_D, pw_rev_D * 100)
ttl_E <- "RRHO2 Reversal"
sub_E <- sprintf("%d proteins | max %.1f", n_shared_E, n_rev_E)

layout <- paste(
  "##############",
  "AAAAAAAABBBBBB",
  "AAAAAAAABBBBBB",
  "AAAAAAAABBBBBB",
  "AAAAAAAABBBBBB",
  "AAAAAAAABBBBBB",
  "AAAAAAAABBBBBB",
  "##############",
  "##############",
  "CCCCCCDDDDEEEE",
  "CCCCCCDDDDEEEE",
  "CCCCCCDDDDEEEE",
  "CCCCCCDDDDEEEE",
  "CCCCCCDDDDEEEE",
  "CCCCCCDDDDEEEE",
  sep = "\n"
)

composite_A_final <- composite_A + plot_annotation(
  theme = theme(plot.margin = margin(-2.5, -1, -2.5, -1, "mm")))
pC_fry_final <- pC_fry + plot_annotation(
  theme = theme(plot.margin = margin(3, 5, 0, 0, "mm")))
pD_final <- pD + theme(plot.margin = margin(-2.8, 5, 2.8, -5, "mm"))
pE_final <- pE_heat + theme(plot.margin = margin(-2.1, -0.2, 3.4, -3.5, "mm"),
                             axis.title = element_text(face = "bold", size = 8))
pB_final <- pB + coord_cartesian(xlim = c(-0.25, X_BAR_MAX + 1.75),
                                  ylim = c(BAR_YMAX + ROW_H * 6.5, -ROW_H * 0.05),
                                  expand = FALSE) +
                 theme(plot.margin = margin(1, -28, 8, -14, "mm"))

fig <- wrap_elements(full = composite_A_final) +
       wrap_elements(full = pB_final) +
       wrap_elements(full = pC_fry_final) +
       wrap_elements(full = pD_final) +
       wrap_elements(full = pE_final) +
       plot_layout(design = layout,
                   widths  = rep(1, 14),
                   heights = c(6.5, rep(10, 6), 4, 4.5, rep(12, 6)))

X_A <- 0.005; X_B <- 0.549; X_C <- 0.012; X_D <- 0.406; X_E <- 0.693
X_TTL      <- 0.030
TAG_DY     <- -0.002
SUB_OFFSET <- 0.020
Y_A <- 0.984; Y_B <- 0.984
Y_C <- 0.512; Y_D <- 0.511; Y_E <- 0.511

composite_final <- ggdraw(fig) +
  draw_label("A",   x = X_A,          y = Y_A - TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_A, x = X_A + X_TTL,  y = Y_A,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_A, x = X_A + X_TTL,  y = Y_A - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("B",   x = X_B,          y = Y_B - TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_B, x = X_B + X_TTL,  y = Y_B,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_B, x = X_B + X_TTL,  y = Y_B - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("C",   x = X_C,          y = Y_C - TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_C, x = X_C + X_TTL,  y = Y_C,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_C, x = X_C + X_TTL,  y = Y_C - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("D",   x = X_D,          y = Y_D - TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_D, x = X_D + X_TTL,  y = Y_D,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_D, x = X_D + X_TTL,  y = Y_D - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("E",   x = X_E,          y = Y_E - TAG_DY,     size = TAG_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(ttl_E, x = X_E + X_TTL,  y = Y_E,              size = TTL_SZ, fontface = "bold",        hjust = 0, vjust = 1) +
  draw_label(sub_E, x = X_E + X_TTL,  y = Y_E - SUB_OFFSET, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_plot(quad_legend, x = 0.64, y = 0.524, width = 0.30, height = 0.045)

ggsave(file.path(RPT_PDF, "MAIN_F05_composite.pdf"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_PNG, "MAIN_F05_composite.png"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message("F05 composite (5-panel) saved")

# ── Biological summary ───────────────────────────────────────────────────────
reversal_total <- q_counts["BR"] + q_counts["TL"]
exacer_total   <- q_counts["TR"] + q_counts["BL"]
message(sprintf("\n=== Biological Summary ==="))
message(sprintf("  Total proteins: %d", nrow(dep_df)))
message(sprintf("  Reversal quadrants (BR+TL): %d (%.1f%%)",
                reversal_total, 100 * reversal_total / nrow(dep_df)))
message(sprintf("  Exacerbation quadrants (TR+BL): %d (%.1f%%)",
                exacer_total, 100 * exacer_total / nrow(dep_df)))
message(sprintf("  Cancer-Up / Training-Down (BR): %d | Cancer-Down / Training-Up (TL): %d",
                q_counts["BR"], q_counts["TL"]))
message(sprintf("  DEPs: Cancer %d | Training CR %d | Sig Both %d",
                sum(dep_df$pi_CvH < 0.05, na.rm = TRUE),
                sum(dep_df$pi_TCR < 0.05, na.rm = TRUE),
                sum(dep_df$sig_class == "Sig Both")))
message(sprintf("  Reversal DEPs: %d/%d (%.1f%%)",
                sum(dep_df$is_sig & dep_df$quadrant %in%
                    c("Reversed (Cancer Up / Training Down)",
                      "Reversed (Cancer Down / Training Up)")),
                n_sig_A,
                100 * sum(dep_df$is_sig & dep_df$quadrant %in%
                    c("Reversed (Cancer Up / Training Down)",
                      "Reversed (Cancer Down / Training Up)")) / max(n_sig_A, 1)))

write_csv(dep_df |>
  select(gene, logFC_CvH, logFC_TCR, t_CvH, t_TCR, pi_CvH, pi_TCR,
         sig_class, is_sig, quadrant),
  file.path(DAT, "scatter_protein_quadrants.csv"))
