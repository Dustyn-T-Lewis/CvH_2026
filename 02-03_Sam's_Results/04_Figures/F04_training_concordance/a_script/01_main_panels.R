#!/usr/bin/env Rscript
# F04 Training Concordance — Main (5-panel composite)
# Sam's CvH parallel of YvO F04 "Training Concordance"
#
# Question: does resistance training produce the same proteomic signature
# in creatine-supplemented (CRE) vs placebo (PLA) participants?
# Contrast mapping:
#   YvO Training_Young -> Training_CRE  (creatine arm training response)
#   YvO Training_Old   -> Training_PLA  (placebo arm training response)
#
# PRE-WARNING: Both contrasts have 0 DEPs at Pi < 0.05.
# Panel D (fGSEA NES scatter) carries the primary signal:
#   Training_CRE: 44 sig pathways; Training_PLA: 1 sig pathway.
# Protein-level panels (A, B, C) use relaxed thresholds and are annotated
# as sparse in captions.
#
# Panels (inline):
#   A: Quadrant ORA scatter (logFC CRE vs logFC PLA) + flanking ORA bars
#   B: Pattern heatmap (per-protein concordance; relaxed pi < 0.10)
#   C: fry rotation test (CRE-direction genes -> PLA t-stats)
#   D: fGSEA NES scatter (pathway concordance) -- primary signal
#   E: RRHO2 threshold-free rank-rank overlap
#
# Layout (3-col, identical to YvO F04 and Sam F05):
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

pdf_device <- grDevices::pdf
PRINT_SCALE <- 380 / 178

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
BASE         <- "02-03_Sam's_Results/04_Figures/F04_training_concordance"
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
cre_dep <- read_csv(file.path(SAM_DEP_DIR, "Training_CRE.csv"),
                    show_col_types = FALSE)
pla_dep <- read_csv(file.path(SAM_DEP_DIR, "Training_PLA.csv"),
                    show_col_types = FALSE)

dep_df <- cre_dep |>
  select(gene, logFC_CRE = logFC, t_CRE = t, pi_CRE = pi_score) |>
  inner_join(
    pla_dep |> select(gene, logFC_PLA = logFC, t_PLA = t, pi_PLA = pi_score),
    by = "gene"
  ) |>
  filter(!is.na(logFC_CRE), !is.na(logFC_PLA))

message(sprintf("Merged: %d proteins", nrow(dep_df)))

# ── Panel A — Quadrant ORA scatter + flanking bars ──────────────────────────
message("=== Panel A: Quadrant ORA ===")

# NOTE: 0 DEPs at Pi < 0.05 in both contrasts.
# Scatter shows all proteins; quadrant ORA runs on all proteins in each quadrant
# (unsorted, no pi filter) to reveal pathway-level preferences.
# Caption annotates: "No proteins reach Pi < 0.05; ORA uses all quadrant proteins."

N_SHOW    <- 5
COMP_RED  <- unname(DIR_COLORS["Up"])
COMP_BLUE <- unname(DIR_COLORS["Down"])

dep_df <- dep_df |>
  mutate(
    sig_class = case_when(
      pi_CRE < 0.05 & pi_PLA < 0.05 ~ "Sig Both",
      pi_CRE < 0.05                  ~ "Sig CRE only",
      pi_PLA < 0.05                  ~ "Sig PLA only",
      TRUE                           ~ "NS"
    ),
    is_sig   = sig_class != "NS",
    quadrant = case_when(
      logFC_CRE > 0 & logFC_PLA > 0 ~ "Concordant Up",
      logFC_CRE < 0 & logFC_PLA < 0 ~ "Concordant Down",
      logFC_CRE > 0 & logFC_PLA < 0 ~ "Discordant (CRE up / PLA down)",
      TRUE                           ~ "Discordant (CRE down / PLA up)"
    )
  )

SIG_COLS_F4 <- c(
  "Sig Both"     = unname(SIG_COLORS_F3["Sig Both"]),
  "Sig CRE only" = unname(SIG_COLORS_F3["Sig CRE only"]),
  "Sig PLA only" = unname(SIG_COLORS_F3["Sig PLA only"]),
  "NS"           = "grey75"
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
           pathway_label  = clean_pathway_name(pathway),
           neg_log10_padj = -log10(padj + 1e-16),
           significant    = padj < 0.05) |>
    arrange(desc(neg_log10_padj)) |>
    slice_head(n = N_SHOW)
}

ora_conc_up <- run_set_ora(
  dep_df$gene[dep_df$quadrant == "Concordant Up"],   "Concordant Up")
ora_conc_dn <- run_set_ora(
  dep_df$gene[dep_df$quadrant == "Concordant Down"],  "Concordant Down")
ora_disc_tr <- run_set_ora(
  dep_df$gene[dep_df$quadrant == "Discordant (CRE up / PLA down)"],
  "Discordant (CRE up/PLA down)")
ora_disc_bl <- run_set_ora(
  dep_df$gene[dep_df$quadrant == "Discordant (CRE down / PLA up)"],
  "Discordant (CRE down/PLA up)")

all_quad_ora <- bind_rows(ora_conc_up, ora_conc_dn, ora_disc_tr, ora_disc_bl)
if (nrow(all_quad_ora) > 0)
  write_csv(all_quad_ora, file.path(DAT, "panel_A", "ora_quadrant.csv"))

xlim_range <- c(-3.1, 3.1)
ylim_range <- c(-2.8, 2.8)

ns_df  <- filter(dep_df, sig_class == "NS")
sig_df <- filter(dep_df, sig_class != "NS")

q_df     <- dep_df |>
  mutate(q = case_when(
    logFC_CRE > 0 & logFC_PLA > 0 ~ "TR",
    logFC_CRE < 0 & logFC_PLA < 0 ~ "BL",
    logFC_CRE > 0 & logFC_PLA < 0 ~ "BR",
    TRUE                           ~ "TL"))
q_counts <- q_df |> count(q) |> deframe()
q_sig    <- q_df |> filter(sig_class != "NS") |> count(q) |> deframe()
for (qq in c("TR","BL","BR","TL")) if (is.na(q_sig[qq])) q_sig[qq] <- 0

label_df <- dep_df |>
  arrange(desc(abs(logFC_CRE) + abs(logFC_PLA))) |>
  slice_head(n = 10)

txt_gene <- scale_text(BASE_GENE, 190) * 0.82 + 1
txt_quad <- scale_text(BASE_QUADRANT, 190) * 0.88

x_breaks <- seq(-3, 3, 1)
y_breaks <- seq(-2, 2, 1)
x_tick_df <- tibble(x = x_breaks[x_breaks != 0], y = 0,
                    label = as.character(x_breaks[x_breaks != 0]))
y_tick_df <- tibble(x = 0, y = y_breaks[y_breaks != 0],
                    label = as.character(y_breaks[y_breaks != 0]))

p_scatter <- ggplot(mapping = aes(x = logFC_CRE, y = logFC_PLA)) +
  # Concordant quadrants (upper-right, lower-left) red; discordant blue
  annotate("rect", xmin = 0, xmax = Inf,  ymin = 0, ymax = Inf,
           fill = "#FFE0E0", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#FFE0E0", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 0, xmax = Inf,  ymin = -Inf, ymax = 0,
           fill = "#DCEEFF", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#DCEEFF", alpha = 0.55, color = "grey70", linewidth = 0.2) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  geom_text(data = x_tick_df, aes(x = x, y = y, label = label),
            vjust = 1.5, size = 1.3 * PRINT_SCALE, color = "grey40", fontface = "bold") +
  geom_text(data = y_tick_df, aes(x = x, y = y, label = label),
            hjust = -0.5, size = 1.3 * PRINT_SCALE, color = "grey40", fontface = "bold") +
  geom_point(data = ns_df, aes(x = logFC_CRE, y = logFC_PLA),
             color = "grey80", fill = "grey85", shape = 21,
             size = 0.35, alpha = 0.3, stroke = 0.10) +
  geom_point(data = sig_df, aes(x = logFC_CRE, y = logFC_PLA, fill = sig_class),
             shape = 21,
             size = 0.9,
             color = "grey75",
             alpha = 0.85,
             stroke = 0.6) +
  scale_fill_manual(values = SIG_COLS_F4, name = "Significance") +
  geom_label_repel(data = label_df, aes(x = logFC_CRE, y = logFC_PLA, label = gene),
                   fill = "white", color = "grey15",
                   size = txt_gene, fontface = "italic", max.overlaps = 40,
                   segment.size = 0.3, segment.color = "grey25",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.3, point.padding = 0.3,
                   force = 6, force_pull = 0.3,
                   label.padding = unit(1, "pt"), label.r = unit(0.5, "pt"),
                   linewidth = 0.15, seed = 42,
                   xlim = c(-3, 3) * 0.85, ylim = c(-2.7, 2.7) * 0.85) +
  annotate("label", x = xlim_range[2], y = ylim_range[2],
           label = sprintf("Concordant Up\n%s/%s", q_sig["TR"], q_counts["TR"]),
           hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
           color = COMP_RED, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt"), lineheight = 0.9) +
  annotate("label", x = xlim_range[1], y = ylim_range[1],
           label = sprintf("%s/%s\nConcordant Down", q_sig["BL"], q_counts["BL"]),
           hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
           color = COMP_RED, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt"), lineheight = 0.9) +
  annotate("label", x = xlim_range[1], y = ylim_range[2],
           label = sprintf("Discordant\n(CRE↓ PLA↑)\n%s/%s", q_sig["TL"], q_counts["TL"]),
           hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
           color = COMP_BLUE, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt"), lineheight = 0.9) +
  annotate("label", x = xlim_range[2], y = ylim_range[1],
           label = sprintf("%s/%s\nDiscordant\n(CRE↑ PLA↓)", q_sig["BR"], q_counts["BR"]),
           hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
           color = COMP_BLUE, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt"), lineheight = 0.9) +
  annotate("text", x = 2.5, y = 0,
           label = expression(log[2]*FC ~ "(Training CRE)"),
           hjust = 0.5, vjust = -0.4, size = 1.3 * PRINT_SCALE, color = "grey30",
           fontface = "bold") +
  annotate("text", x = 0, y = 2.0,
           label = expression(log[2]*FC ~ "(Training PLA)"),
           hjust = 0.5, vjust = -0.4, size = 1.3 * PRINT_SCALE, color = "grey30",
           fontface = "bold", angle = 90) +
  annotate("text", x = 0, y = ylim_range[1] * 0.92,
           label = "NOTE: 0 DEPs at Π<0.05 in both arms; labels = top by |logFC|",
           hjust = 0.5, vjust = 0, size = 1.1 * PRINT_SCALE, color = "grey50",
           fontface = "italic") +
  coord_cartesian(xlim = xlim_range, ylim = ylim_range, expand = FALSE) +
  labs(x = NULL, y = NULL) +
  FIG_THEME +
  theme(plot.title      = element_blank(),
        plot.subtitle   = element_blank(),
        axis.text       = element_blank(),
        axis.ticks      = element_blank(),
        axis.title      = element_blank(),
        plot.margin     = margin(2, 0, 0, 0, "mm"),
        legend.position = "none")

# Significance key
key_lvls    <- c("Sig Both", "Sig CRE only", "Sig PLA only")
key_display <- c("Sig Both", "Sig CRE", "Sig PLA")
key_df <- tibble(
  category = factor(key_lvls, levels = key_lvls),
  display  = key_display,
  fill_col = unname(SIG_COLS_F4[key_lvls]),
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

# Half-bar builder (shared with F05 pattern)
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

p_ul <- make_half_bars(ora_conc_dn, scales::alpha(COMP_RED,  0.30), "left",  c(0, 2.8))
p_ll <- make_half_bars(ora_disc_bl, scales::alpha(COMP_BLUE, 0.30), "left",  c(-2.8, 0))
p_ur <- make_half_bars(ora_conc_up, scales::alpha(COMP_RED,  0.30), "right", c(0, 2.8))
p_lr <- make_half_bars(ora_disc_tr, scales::alpha(COMP_BLUE, 0.30), "right", c(-2.8, 0))

design_A <- c(
  area(1, 1), area(1, 2, 2, 2), area(1, 3),
  area(2, 1), area(2, 3),       area(3, 1, 3, 3)
)
n_total_A  <- nrow(dep_df)
n_sig_A    <- sum(dep_df$is_sig)
n_enrich_A <- if (nrow(all_quad_ora) > 0) sum(all_quad_ora$significant) else 0L
r_spear_A  <- cor(dep_df$logFC_CRE, dep_df$logFC_PLA, use = "complete.obs",
                  method = "spearman")

composite_A <- p_ul + p_scatter + p_ur + p_ll + p_lr + p_key +
  plot_layout(design = design_A,
              widths  = c(70, 100, 70) / 240,
              heights = c(85, 85, 8) / 178) +
  plot_annotation(
    title    = "Training Concordance: Quadrant ORA",
    subtitle = sprintf(
      "All-protein ORA (hypergeometric) | N = %d | %d DEPs (Π<0.05; SPARSE — see caption) | %d enriched (FDR) | ρ = %.2f",
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

# Strip for composite
composite_A <- composite_A &
  labs(title = NULL, subtitle = NULL, tag = NULL) &
  theme(legend.position = "none")
composite_A <- composite_A +
  plot_annotation(title = NULL, subtitle = NULL,
                  theme = theme(plot.title = element_blank(), plot.subtitle = element_blank()))

# ── Panel D — fGSEA NES scatter (pathway-level concordance) ─────────────────
message("=== Panel D: NES scatter ===")

cre_fgsea <- readRDS(file.path(FGSEA_CACHE, "Training_CRE_fgsea.rds")) |>
  filter(!is.na(NES))
pla_fgsea <- readRDS(file.path(FGSEA_CACHE, "Training_PLA_fgsea.rds")) |>
  filter(!is.na(NES))

fgsea_wide <- cre_fgsea |>
  select(pathway, NES_CRE = NES, padj_CRE = padj) |>
  inner_join(
    pla_fgsea |> select(pathway, NES_PLA = NES, padj_PLA = padj),
    by = "pathway"
  ) |>
  mutate(
    sig_CRE   = !is.na(padj_CRE) & padj_CRE < 0.05,
    sig_PLA   = !is.na(padj_PLA) & padj_PLA < 0.05,
    sig_class = case_when(
      sig_CRE & sig_PLA ~ "Sig Both",
      sig_CRE           ~ "Sig CRE only",
      sig_PLA           ~ "Sig PLA only",
      TRUE              ~ "NS"
    ),
    concordant = case_when(
      sig_CRE & sig_PLA & (NES_CRE * NES_PLA > 0) ~ "Concordant",
      sig_CRE & sig_PLA & (NES_CRE * NES_PLA < 0) ~ "Discordant",
      TRUE ~ "NS"
    )
  )

write_csv(fgsea_wide, file.path(DAT, "panel_D", "nes_scatter.csv"))

nes_cor_all <- cor.test(fgsea_wide$NES_CRE, fgsea_wide$NES_PLA,
                         method = "spearman", exact = FALSE)
rho_D      <- as.numeric(nes_cor_all$estimate)
n_pw_D     <- nrow(fgsea_wide)
n_sig_pw_D <- sum(fgsea_wide$sig_CRE | fgsea_wide$sig_PLA)

sig_both     <- filter(fgsea_wide, sig_CRE & sig_PLA)
pw_conc_D    <- if (nrow(sig_both) > 0) mean(sig_both$concordant == "Concordant") else 0

# Bootstrap CI for rho
set.seed(42)
boot_rho_fn <- function(data, indices) {
  d <- data[indices, ]
  cor(d$NES_CRE, d$NES_PLA, use = "complete.obs", method = "spearman")
}
b_rho  <- tryCatch(boot::boot(fgsea_wide, statistic = boot_rho_fn, R = 500),
                   error = function(e) NULL)
rho_ci <- if (!is.null(b_rho)) {
  ci <- tryCatch(boot::boot.ci(b_rho, type = "perc", conf = 0.95),
                 error = function(e) NULL)
  if (!is.null(ci)) ci$percent[4:5] else c(NA_real_, NA_real_)
} else c(NA_real_, NA_real_)

# Labels: sig-CRE pathways sorted by |NES_CRE|
label_pw <- filter(fgsea_wide, sig_CRE) |>
  mutate(score = abs(NES_CRE) + abs(NES_PLA)) |>
  arrange(desc(score)) |>
  slice_head(n = 12) |>
  mutate(label = clean_pathway_name(pathway))

pD <- ggplot(fgsea_wide, aes(x = NES_CRE, y = NES_PLA)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "#FFE0E0", alpha = 0.45) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#FFE0E0", alpha = 0.45) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "#DCEEFF", alpha = 0.45) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#DCEEFF", alpha = 0.45) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  geom_point(data = filter(fgsea_wide, sig_class == "NS"),
             color = "grey80", fill = "grey85", shape = 21,
             size = 0.5, alpha = 0.3, stroke = 0.1) +
  geom_point(data = filter(fgsea_wide, sig_class != "NS"),
             aes(fill = sig_class), shape = 21,
             size = 1.2, alpha = 0.8, stroke = 0.4) +
  scale_fill_manual(values = c(
    "Sig Both"     = unname(SIG_COLORS_F3["Sig Both"]),
    "Sig CRE only" = unname(SIG_COLORS_F3["Sig CRE only"]),
    "Sig PLA only" = unname(SIG_COLORS_F3["Sig PLA only"])
  ), name = "Significance") +
  geom_label_repel(data = label_pw, aes(x = NES_CRE, y = NES_PLA, label = label),
                   size = scale_text(BASE_PATHWAY, 146) * 0.80,
                   fill = alpha("white", 0.9), color = "grey15",
                   fontface = "bold", max.overlaps = 30,
                   segment.size = 0.25, box.padding = 0.2, seed = 42,
                   label.padding = unit(1, "pt"), label.r = unit(0.5, "pt"),
                   linewidth = 0.15) +
  labs(title    = "Pathway Concordance (fGSEA NES)",
       subtitle = sprintf(
         "ρ = %.2f [95%% CI %.2f–%.2f] | %d sig pathways (CRE+PLA combined) | primary signal carrier",
         rho_D, rho_ci[1], rho_ci[2], n_sig_pw_D),
       x = "NES (Training CRE)",
       y = "NES (Training PLA)") +
  FIG_THEME +
  theme(legend.position = "none")

ggsave(file.path(PNL_PNG, "MAIN_panel_D_nes_scatter.png"), pD,
       width = 80, height = 80, units = "mm", dpi = 300)
ggsave(file.path(PNL_PDF, "MAIN_panel_D_nes_scatter.pdf"), pD,
       width = 80, height = 80, units = "mm", device = pdf_device)
message("Panel D done")

# ── Panel C — fry rotation test ──────────────────────────────────────────────
message("=== Panel C: fry ===")

# NOTE: 0 DEPs at Pi < 0.05 for both CRE and PLA. Use sign(NES) of fGSEA
# leading edge genes from Training_CRE as the "source" gene set.
# Specifically: top 20 CRE up-genes and top 20 CRE down-genes by NES.
# Then test whether those sets are enriched in PLA t-stats.

cre_up_genes <- cre_fgsea |>
  filter(!is.na(NES), NES > 0, !is.na(padj), padj < 0.05) |>
  arrange(desc(NES)) |>
  head(10) |>
  pull(leadingEdge) |>
  unlist() |>
  unique()

cre_dn_genes <- cre_fgsea |>
  filter(!is.na(NES), NES < 0, !is.na(padj), padj < 0.05) |>
  arrange(NES) |>
  head(10) |>
  pull(leadingEdge) |>
  unlist() |>
  unique()

message(sprintf("  CRE leading-edge Up: %d genes | Down: %d genes",
                length(cre_up_genes), length(cre_dn_genes)))

# t-stats for Training_PLA (for fry input)
t_PLA_named <- setNames(dep_df$t_PLA, dep_df$gene)
t_PLA_named <- t_PLA_named[!is.na(t_PLA_named)]

t_df <- tibble(gene = names(t_PLA_named), t = as.numeric(t_PLA_named)) |>
  arrange(t) |>
  mutate(rank = seq_len(n()),
         set  = case_when(gene %in% cre_up_genes ~ "cre_up",
                          gene %in% cre_dn_genes ~ "cre_dn",
                          TRUE ~ "NS"))

n_all_C <- nrow(t_df)
n_up_C  <- length(intersect(cre_up_genes, names(t_PLA_named)))
n_dn_C  <- length(intersect(cre_dn_genes, names(t_PLA_named)))

# fry rotation test
fry_sets_c <- list(
  cre_up = cre_up_genes[cre_up_genes %in% names(t_PLA_named)],
  cre_dn = cre_dn_genes[cre_dn_genes %in% names(t_PLA_named)]
)
fry_sets_c <- fry_sets_c[lengths(fry_sets_c) >= 5]

fry_res <- tryCatch({
  requireNamespace("limma", quietly = TRUE)
  t_vec   <- t_PLA_named[!is.na(t_PLA_named)]
  t_mat   <- matrix(t_vec, ncol = 1, dimnames = list(names(t_vec), "t"))
  idx_list <- lapply(fry_sets_c, function(g) which(rownames(t_mat) %in% g))
  idx_list <- idx_list[lengths(idx_list) >= 5]
  if (length(idx_list) < 1) stop("No valid sets for fry")
  limma::fry(t_mat, index = idx_list, geneid = rownames(t_mat))
}, error = function(e) {
  message("  fry error: ", e$message)
  NULL
})

if (!is.null(fry_res)) {
  write_csv(as.data.frame(fry_res) |> tibble::rownames_to_column("set"),
            file.path(DAT, "panel_C_fry", "fry_results_all.csv"))
  fry_up_p   <- if ("cre_up" %in% rownames(fry_res)) fry_res["cre_up", "PValue"] else NA
  fry_dn_p   <- if ("cre_dn" %in% rownames(fry_res)) fry_res["cre_dn", "PValue"] else NA
  fry_up_dir <- if ("cre_up" %in% rownames(fry_res)) fry_res["cre_up", "Direction"] else NA
  fry_dn_dir <- if ("cre_dn" %in% rownames(fry_res)) fry_res["cre_dn", "Direction"] else NA
} else {
  fry_up_p <- fry_dn_p <- fry_up_dir <- fry_dn_dir <- NA
}

make_barcode_conc <- function(set_label, set_name, color) {
  set_genes <- if (set_name == "cre_up") cre_up_genes else cre_dn_genes
  df_set    <- t_df |> filter(gene %in% set_genes)
  n_set     <- nrow(df_set)

  all_ranks <- t_df$rank
  hit_ranks <- sort(df_set$rank)
  n_total   <- length(all_ranks)
  n_hit     <- length(hit_ranks)
  if (n_hit < 1) return(ggplot() + theme_void())

  step_up <- 1 / n_hit
  step_dn <- 1 / (n_total - n_hit)
  es_x    <- numeric(n_hit * 2 + 2)
  es_y    <- numeric(n_hit * 2 + 2)
  cur_es  <- 0
  prev_rank <- 0
  k <- 1
  for (i in seq_along(hit_ranks)) {
    gap_steps <- hit_ranks[i] - prev_rank - 1
    if (gap_steps > 0) {
      cur_es  <- cur_es - step_dn * gap_steps
      k       <- k + 1
      es_x[k] <- hit_ranks[i] - 0.5
      es_y[k] <- cur_es
      k       <- k + 1
    }
    cur_es    <- cur_es + step_up
    k         <- k + 1
    es_x[k]   <- hit_ranks[i]
    es_y[k]   <- cur_es
    prev_rank <- hit_ranks[i]
  }
  if (prev_rank < n_total) {
    cur_es  <- cur_es - step_dn * (n_total - prev_rank)
    k       <- k + 1
    es_x[k] <- n_total
    es_y[k] <- cur_es
  }
  es_df  <- tibble(x = es_x[1:k], y = es_y[1:k])
  max_es <- max(abs(es_df$y))

  tick_df   <- tibble(x = hit_ranks, y = 0)
  fry_p_val <- if (set_name == "cre_up") fry_up_p else fry_dn_p
  fry_dir   <- if (set_name == "cre_up") fry_up_dir else fry_dn_dir

  stat_text <- paste0(
    set_label, " (n=", n_hit, ")",
    if (!is.na(fry_p_val))
      sprintf("\nfry p=%.3f (%s)", fry_p_val, fry_dir)
    else "\n(using CRE leading-edge genes; 0 DEPs at Pi<0.05)"
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
    labs(x = "Proteins ranked by Training PLA t-stat (Low -> High)",
         y = "Enrichment score") +
    FIG_THEME
}

pC_up  <- make_barcode_conc("CRE Leading-Edge Up (top pathways)", "cre_up",  COMP_RED)
pC_dn  <- make_barcode_conc("CRE Leading-Edge Down (top pathways)", "cre_dn", COMP_BLUE)
pC_fry <- pC_up / pC_dn +
  plot_annotation(
    title    = "fry: Training Concordance",
    subtitle = sprintf(
      "CRE leading-edge genes vs Training_PLA t-stats | n = %d proteins | 0 DEPs at Π<0.05 (uses pathway leading-edges)",
      n_all_C),
    theme = theme(plot.title    = element_text(size = FIG_TITLE_SIZE, face = "bold"),
                   plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, color = "grey30")))

ggsave(file.path(PNL_PNG, "MAIN_panel_C_fry.png"), pC_fry,
       width = 130, height = 80, units = "mm", dpi = 300)
ggsave(file.path(PNL_PDF, "MAIN_panel_C_fry.pdf"), pC_fry,
       width = 130, height = 80, units = "mm", device = pdf_device)

# Save driving proteins
driving_df <- dep_df |>
  filter(gene %in% c(cre_up_genes, cre_dn_genes)) |>
  mutate(set = if_else(gene %in% cre_up_genes, "cre_up", "cre_dn"))
if (nrow(driving_df) > 0)
  write_csv(driving_df, file.path(DAT, "panel_C_fry", "driving_proteins.csv"))

n_all_C   <- n_all_C
cor_imp_C <- cor(dep_df$logFC_CRE, dep_df$logFC_PLA, use = "complete.obs")
message("Panel C done")

# ── Panel E — RRHO2 ──────────────────────────────────────────────────────────
message("=== Panel E: RRHO2 ===")

rrho_df <- dep_df |>
  filter(!is.na(t_CRE), !is.na(t_PLA)) |>
  transmute(gene, t_CRE, t_PLA)

n_shared <- nrow(rrho_df)
message(sprintf("  RRHO2 input: %d proteins", n_shared))

rrho_dir <- file.path(DAT, "panel_E")
dir.create(rrho_dir, recursive = TRUE, showWarnings = FALSE)

pE_heat <- tryCatch({
  df1 <- data.frame(gene = rrho_df$gene, rank = rrho_df$t_CRE)
  df2 <- data.frame(gene = rrho_df$gene, rank = rrho_df$t_PLA)

  rrho_obj <- RRHO2_initialize(df1, df2,
                                labels = c("Training CRE", "Training PLA"),
                                boundary = 0.05,
                                log10.ind = TRUE,
                                method = "hyper")
  mat <- rrho_obj$hypermat
  nr <- nrow(mat); nc <- ncol(mat)
  h1 <- floor(nr / 2); h2 <- floor(nc / 2)

  max_UU <- max(mat[1:h1, 1:h2], na.rm = TRUE)
  max_DD <- max(mat[(h1+1):nr, (h2+1):nc], na.rm = TRUE)
  max_UD <- max(mat[1:h1, (h2+1):nc], na.rm = TRUE)
  max_DU <- max(mat[(h1+1):nr, 1:h2], na.rm = TRUE)

  rrho_summary <- tibble(
    quadrant    = c("UU", "UD", "DU", "DD"),
    label       = c("Concordant Up", "Discordant (CRE Up / PLA Down)",
                    "Discordant (CRE Down / PLA Up)", "Concordant Down"),
    max_overlap = c(max_UU, max_UD, max_DU, max_DD)
  )
  write_csv(rrho_summary, file.path(rrho_dir, "rrho2_summary.csv"))

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
      name = "-log10(p)", na.value = "grey95"
    ) +
    geom_hline(yintercept = nr / 2, color = "grey30", linewidth = 0.4, linetype = "dashed") +
    geom_vline(xintercept = nc / 2, color = "grey30", linewidth = 0.4, linetype = "dashed") +
    annotate("text", x = nc * 0.25, y = nr * 0.12,
             label = sprintf("Concordant Up\n%.1f", max_UU),
             color = "white", size = scale_text(BASE_STAT, 80) * 0.85, fontface = "bold") +
    annotate("text", x = nc * 0.75, y = nr * 0.12,
             label = sprintf("Discordant\n(CRE↑ PLA↓)\n%.1f", max_UD),
             color = "white", size = scale_text(BASE_STAT, 80) * 0.85, fontface = "bold") +
    annotate("text", x = nc * 0.25, y = nr * 0.88,
             label = sprintf("Discordant\n(CRE↓ PLA↑)\n%.1f", max_DU),
             color = "white", size = scale_text(BASE_STAT, 80) * 0.85, fontface = "bold") +
    annotate("text", x = nc * 0.75, y = nr * 0.88,
             label = sprintf("Concordant Down\n%.1f", max_DD),
             color = "white", size = scale_text(BASE_STAT, 80) * 0.85, fontface = "bold") +
    scale_x_continuous(name = expression("Training CRE rank"~(Up %->% Down)),
                       expand = c(0, 0)) +
    scale_y_continuous(name = expression("Training PLA rank"~(Up %->% Down)),
                       expand = c(0, 0)) +
    labs(title    = "RRHO2: Threshold-Free Concordance",
         subtitle = sprintf("%d proteins | warm diagonal = shared training signature (CRE = PLA)",
                            n_shared)) +
    FIG_THEME +
    theme(axis.text         = element_blank(),
          axis.ticks        = element_blank(),
          legend.key.height = unit(3, "mm"),
          legend.title      = element_text(size = FIG_LEGEND_TITLE, face = "bold"),
          legend.text       = element_text(size = FIG_LEGEND_TEXT))

}, error = function(e) {
  message("  RRHO2 error: ", e$message, " — building stub panel")
  rrho_summary <- tibble(quadrant = character(), max_overlap = numeric())
  write_csv(rrho_summary, file.path(rrho_dir, "rrho2_summary.csv"))

  dep_df2 <- dep_df |> filter(!is.na(t_CRE), !is.na(t_PLA)) |>
    mutate(q = case_when(
      t_CRE > 0 & t_PLA > 0 ~ "Concordant Up",
      t_CRE < 0 & t_PLA < 0 ~ "Concordant Down",
      t_CRE > 0 & t_PLA < 0 ~ "Discordant (CRE↑ PLA↓)",
      TRUE                   ~ "Discordant (CRE↓ PLA↑)"
    ))
  q_cnt <- dep_df2 |> count(q)

  ggplot(q_cnt, aes(x = reorder(q, n), y = n, fill = q)) +
    geom_col(color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c(
      "Concordant Up"              = "#D6604D",
      "Concordant Down"            = "#F4A261",
      "Discordant (CRE↑ PLA↓)" = "#4393C3",
      "Discordant (CRE↓ PLA↑)" = "#74ADD1"
    )) +
    geom_text(aes(label = n), hjust = -0.2, size = scale_text(BASE_STAT, 80)) +
    coord_flip() +
    labs(title    = "RRHO2 Stub — Quadrant Protein Counts",
         subtitle = sprintf("RRHO2 pending | %d proteins by t-stat quadrant", n_shared),
         x = NULL, y = "Protein count") +
    FIG_THEME +
    theme(legend.position = "none")
})

max_UU_E <- max_DD_E <- NA_real_
rrho_sum_path <- file.path(rrho_dir, "rrho2_summary.csv")
if (file.exists(rrho_sum_path)) {
  rrho_sum <- read_csv(rrho_sum_path, show_col_types = FALSE)
  if (nrow(rrho_sum) > 0 && "max_overlap" %in% names(rrho_sum)) {
    max_UU_E <- rrho_sum$max_overlap[rrho_sum$quadrant == "UU"]
    max_DD_E <- rrho_sum$max_overlap[rrho_sum$quadrant == "DD"]
    if (length(max_UU_E) == 0) max_UU_E <- NA_real_
    if (length(max_DD_E) == 0) max_DD_E <- NA_real_
  }
}
n_shared_E <- n_shared
n_conc_E   <- if (!is.na(max_UU_E)) max(max_UU_E, max_DD_E, na.rm = TRUE) else 0

ggsave(file.path(PNL_PNG, "MAIN_panel_E_rrho2.png"), pE_heat,
       width = 90, height = 90, units = "mm", dpi = 300)
ggsave(file.path(PNL_PDF, "MAIN_panel_E_rrho2.pdf"), pE_heat,
       width = 90, height = 90, units = "mm", device = pdf_device)
message("Panel E done")

# ── Panel B — Pattern heatmap ────────────────────────────────────────────────
message("=== Panel B: Pattern heatmap ===")

ROW_H <- 0.078

# 0 DEPs at Pi < 0.05; relax to padj < 0.10 for protein inclusion
heatmap_df <- dep_df |>
  filter(!is.na(logFC_CRE), !is.na(logFC_PLA)) |>
  mutate(
    pi_CRE_adj = pi_CRE,
    pi_PLA_adj = pi_PLA,
    padj_flag  = pi_CRE < 0.10 | pi_PLA < 0.10
  ) |>
  filter(padj_flag) |>
  mutate(
    quadrant = case_when(
      logFC_CRE > 0 & logFC_PLA > 0 ~ "Concordant Up",
      logFC_CRE < 0 & logFC_PLA < 0 ~ "Concordant Down",
      TRUE                            ~ "Discordant"
    ),
    sig_cat = case_when(
      pi_CRE < 0.05 & pi_PLA < 0.05 ~ "Both",
      pi_CRE < 0.05                  ~ "CRE",
      pi_PLA < 0.05                  ~ "PLA",
      pi_CRE < 0.10                  ~ "CRE (q<0.10)",
      TRUE                           ~ "PLA (q<0.10)"
    )
  ) |>
  arrange(match(quadrant, c("Concordant Up", "Concordant Down", "Discordant")),
          desc(logFC_CRE))

# Warn if completely empty
if (nrow(heatmap_df) == 0) {
  message("  No proteins pass pi < 0.10 — heatmap will be stub")
  heatmap_df <- dep_df |>
    arrange(desc(abs(logFC_CRE))) |>
    slice_head(n = 40) |>
    mutate(
      quadrant = case_when(
        logFC_CRE > 0 & logFC_PLA > 0 ~ "Concordant Up",
        logFC_CRE < 0 & logFC_PLA < 0 ~ "Concordant Down",
        TRUE                            ~ "Discordant"
      ),
      sig_cat = "Top |logFC| (no pi<0.10)"
    ) |>
    arrange(match(quadrant, c("Concordant Up", "Concordant Down", "Discordant")),
            desc(logFC_CRE))
}

write_csv(heatmap_df |> select(gene, quadrant, sig_cat, logFC_CRE, logFC_PLA,
                                pi_CRE, pi_PLA),
          file.path(DAT, "panel_B_heatmap", "pattern_classification.csv"))

n_total   <- nrow(heatmap_df)
n_conc_up <- sum(heatmap_df$quadrant == "Concordant Up")
n_conc_dn <- sum(heatmap_df$quadrant == "Concordant Down")
n_discord <- sum(heatmap_df$quadrant == "Discordant")
n_pw      <- 3

QUAD_COLORS <- c("Concordant Up" = "#B2182B", "Concordant Down" = "#2166AC",
                 "Discordant"    = "#1B7837")
QUAD_BG     <- c("Concordant Up" = "#F4D9D2", "Concordant Down" = "#D5DEEF",
                 "Discordant"    = "#C8E0CD")
SIG_COLORS_B <- c("Both" = "#2E7D32", "CRE" = "#2166AC",
                   "PLA" = "#D6604D", "CRE (q<0.10)" = "#A8C8E8",
                   "PLA (q<0.10)" = "#F4B8A0", "Top |logFC| (no pi<0.10)" = "grey50",
                   "NS" = "grey70")

Y_MAX <- n_total * ROW_H
heatmap_df2 <- heatmap_df |>
  mutate(
    row_idx      = seq_len(n()),
    y            = n_total - row_idx + 1,
    lfc_fill_CRE = pmax(-3, pmin(3, logFC_CRE)),
    lfc_fill_PLA = pmax(-3, pmin(3, logFC_PLA)),
    bg_color     = QUAD_BG[quadrant]
  )

make_lfc_color <- function(lfc) {
  cols <- colorRampPalette(c("#08306B", "#4393C3", "white", "#D6604D", "#67000D"))(101)
  idx  <- round((pmax(-3, pmin(3, lfc)) + 3) / 6 * 100) + 1
  cols[idx]
}

heatmap_df2 <- heatmap_df2 |>
  mutate(fill_CRE = make_lfc_color(lfc_fill_CRE),
         fill_PLA = make_lfc_color(lfc_fill_PLA))

x_CRE       <- 0.5; x_PLA <- 1.5
x_bar_start <- 2.3
X_BAR_MAX   <- x_bar_start + max(abs(dep_df$logFC_CRE), abs(dep_df$logFC_PLA), na.rm = TRUE) * 0.4
BAR_YMAX    <- max(heatmap_df2$y) + ROW_H

tile_df <- bind_rows(
  heatmap_df2 |> transmute(gene, y, x = x_CRE, lfc = lfc_fill_CRE, contrast = "CRE"),
  heatmap_df2 |> transmute(gene, y, x = x_PLA, lfc = lfc_fill_PLA, contrast = "PLA")
)

sig_dot_df <- heatmap_df2 |>
  mutate(x_sig = 2.0, sig_fill = SIG_COLORS_B[sig_cat])

# Quadrant dividers — safe guard for 0-row quadrant
sep_y1 <- if (n_conc_up > 0 && n_conc_dn > 0) n_conc_dn + 0.5 else NULL
sep_y2 <- if (n_conc_dn > 0 && n_discord > 0) n_conc_dn + n_discord + 0.5 else NULL

pB <- ggplot() +
  geom_rect(data = heatmap_df2,
            aes(xmin = -0.2, xmax = X_BAR_MAX + 1.5,
                ymin = y - ROW_H / 2, ymax = y + ROW_H / 2),
            fill = heatmap_df2$bg_color, color = NA) +
  geom_tile(data = tile_df,
            aes(x = x, y = y, fill = lfc),
            width = 0.6, height = ROW_H * 0.85, color = NA) +
  scale_fill_gradientn(
    colors = c("#08306B", "#4393C3", "white", "#D6604D", "#67000D"),
    values = scales::rescale(c(-3, -1, 0, 1, 3)),
    limits = c(-3, 3), oob = scales::squish,
    name = expression(log[2]*FC)) +
  geom_point(data = sig_dot_df,
             aes(x = x_sig, y = y),
             color = sig_dot_df$sig_fill, size = 0.3, shape = 16, alpha = 0.8) +
  { if (!is.null(sep_y1))
      geom_hline(yintercept = sep_y1, color = "grey30", linewidth = 0.4)
    else list() } +
  { if (!is.null(sep_y2))
      geom_hline(yintercept = sep_y2, color = "grey30", linewidth = 0.4)
    else list() } +
  annotate("text", x = x_CRE, y = BAR_YMAX + ROW_H * 1,
           label = "Tr.(CRE)", hjust = 0.5, size = scale_text(BASE_STAT, 90),
           fontface = "bold",
           color = unname(CONTRAST_COLORS["Training_CRE"])) +
  annotate("text", x = x_PLA, y = BAR_YMAX + ROW_H * 1,
           label = "Tr.(PLA)", hjust = 0.5, size = scale_text(BASE_STAT, 90),
           fontface = "bold",
           color = unname(CONTRAST_COLORS["Training_PLA"])) +
  annotate("text", x = 0.5, y = n_total - n_conc_up / 2 + 0.5,
           label = sprintf("Conc. Up\n(n=%d)", n_conc_up),
           hjust = 0, size = scale_text(BASE_STAT, 90) * 0.85,
           color = unname(QUAD_COLORS["Concordant Up"]), fontface = "bold") +
  annotate("text", x = 0.5, y = n_conc_dn / 2 + 0.5,
           label = sprintf("Conc. Down\n(n=%d)", n_conc_dn),
           hjust = 0, size = scale_text(BASE_STAT, 90) * 0.85,
           color = unname(QUAD_COLORS["Concordant Down"]), fontface = "bold") +
  scale_x_continuous(limits = c(-0.25, X_BAR_MAX + 1.75), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-ROW_H, BAR_YMAX + ROW_H * 3), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(title    = "Protein-to-Pathway",
       subtitle = sprintf(
         "%d proteins (Π<0.10 relaxed threshold — see caption) | %d patterns",
         n_total, n_pw)) +
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

# ── Quadrant legend ──────────────────────────────────────────────────────────
inset_quad_df <- tibble(
  quadrant  = factor(c("Concordant Up", "Concordant Down", "Discordant"),
                     levels = c("Concordant Up", "Concordant Down", "Discordant")),
  bg_color  = unname(QUAD_BG[c("Concordant Up", "Concordant Down", "Discordant")]),
  bar_color = unname(QUAD_COLORS[c("Concordant Up", "Concordant Down", "Discordant")])
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

ttl_A <- "Quadrant ORA (Concordance)"
sub_A <- sprintf("N = %d | %d DEPs (Π SPARSE) | %d enriched (FDR) | ρ = %.2f",
                 n_total_A, n_sig_A, n_enrich_A, r_spear_A)
ttl_B <- "Protein-to-Pathway"
sub_B <- sprintf("%d proteins (Π<0.10) | %d patterns", n_total, n_pw)
ttl_C <- "fry: Concordance"
sub_C <- sprintf("n = %d | r = %.3f", n_all_C, cor_imp_C)
ttl_D <- "Pathway Concordance"
sub_D <- sprintf("ρ = %.2f | %.0f%% concordant (sig-both)", rho_D, pw_conc_D * 100)
ttl_E <- "RRHO2 Concordance"
sub_E <- sprintf("%d proteins | max %.1f", n_shared_E, n_conc_E)

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

ggsave(file.path(RPT_PDF, "MAIN_F04_composite.pdf"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_PNG, "MAIN_F04_composite.png"), composite_final,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message("F04 composite (5-panel) saved")

# ── Biological summary ───────────────────────────────────────────────────────
q_conc <- q_counts["TR"] + q_counts["BL"]
q_disc <- q_counts["BR"] + q_counts["TL"]
message("\n=== Biological Summary ===")
message(sprintf("  Total proteins: %d", nrow(dep_df)))
message(sprintf("  Concordant quadrants (TR+BL): %d (%.1f%%)",
                q_conc, 100 * q_conc / nrow(dep_df)))
message(sprintf("  Discordant quadrants (BR+TL): %d (%.1f%%)",
                q_disc, 100 * q_disc / nrow(dep_df)))
message(sprintf("  Spearman rho (all pathways): %.3f", rho_D))
sig_conc <- sum(fgsea_wide$sig_CRE & fgsea_wide$sig_PLA &
                  fgsea_wide$NES_CRE * fgsea_wide$NES_PLA > 0, na.rm = TRUE)
sig_disc <- sum(fgsea_wide$sig_CRE & fgsea_wide$sig_PLA &
                  fgsea_wide$NES_CRE * fgsea_wide$NES_PLA < 0, na.rm = TRUE)
message(sprintf("  Sig-both pathways: concordant = %d, discordant = %d", sig_conc, sig_disc))
message(sprintf("  CRE sig pathways: %d | PLA sig pathways: %d",
                sum(fgsea_wide$sig_CRE), sum(fgsea_wide$sig_PLA)))

write_csv(dep_df |>
  select(gene, logFC_CRE, logFC_PLA, t_CRE, t_PLA, pi_CRE, pi_PLA,
         sig_class, is_sig, quadrant),
  file.path(DAT, "scatter_protein_quadrants.csv"))
