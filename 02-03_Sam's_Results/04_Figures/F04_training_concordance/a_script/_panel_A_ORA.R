#!/usr/bin/env Rscript
# F04 Panel A ORA: Scatter + Flanking ORA Bars Composite
# Threshold-free ORA on all proteins per quadrant, displayed as bars flanking
# the concordance scatter. Red = same-direction (concordant), Blue = opposing.
#
# Sourced by 01_main_panels.R — expects style.R + pathway_utils.R already loaded.
# Exports: composite (patchwork), n_total, n_sig, n_enrich, r_spear

source("02-03_Sam's_Results/04_Figures/shared/print_scale_apply_380mm.R")
source("02-03_Sam's_Results/04_Figures/shared/pathway_utils.R")
library(fgsea)
library(ggrepel)

BASE <- "02-03_Sam's_Results/04_Figures/F04_training_concordance"
RPT_PNG <- file.path(BASE, "b_reports", "main", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "main", "pdf", "panels")
DAT <- file.path(BASE, "c_data")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_A"), recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

COMP_RED   <- unname(DIR_COLORS["Up"])
COMP_BLUE  <- unname(DIR_COLORS["Down"])
N_SHOW     <- 5

# F04-specific pathway label shortenings (used in the half-bar panels below)
DISPLAY_LABELS_F04 <- c(
  "Cargo Recognition For Clathrin Mediated Endocytosis" = "Clathrin Endocytosis",
  "The Role Of Gtse1 In G2 M Progression After G2 Checkpoint" = "GTSE1 G2/M Progression",
  "Reference Mitochondrial Complex Ucp1 In Thermogenesis" = "Mito Complex UCP1",
  "Reference Electron Transfer In Complex I" = "e- Transfer\nin Complex I",
  "Formation Of The Dystrophin Glycoprotein Complex Dgc" = "Dystrophin Complex",
  "Striated Muscle Contraction"        = "Striated Muscle\nContraction",
  "Metabolism Of Vitamins And Cofactors" = "Vitamins/Cofactors Metab.",
  "Asparagine N Linked Glycosylation"   = "Asparagine N-Linked\nGlycosylation",
  "Regulation Of Expression Of Slits And Robos" = "Slit/Robo Expression\nRegulation",
  "Ribosome Associated Quality Control" = "Ribosome QC",
  "Cellular Response To Starvation"     = "Starvation Response",
  "ATP Synthesis Coupled Electron Transport" = "ATP Synthesis e- Transport",
  "Respiratory Electron Transport"     = "Respiratory ETC",
  "Complex I Biogenesis"               = "Complex I Biogenesis",
  "Reference Translation Initiation"   = "Translation Initiation",
  "Non Integrin Membrane Ecm Interactions" = "Non-Integrin ECM",
  "Polyol Metabolic Process"            = "Polyol Metabolism",
  "Rac1 Gtpase Cycle"                   = "RAC1 GTPase Cycle",
  "Rac3 Gtpase Cycle"                   = "RAC3 GTPase Cycle"
)

# Data
dep_df <- read_csv("02-03_Sam's_Results/03_DEP/c_data/F04_combined_CREvPLA.csv",
                   show_col_types = FALSE)
imp_path <- "02-03_Sam's_Results/02_Imputation/c_data/02_mar_mnar_classification.csv"
imputation_df <- if (file.exists(imp_path)) {
  read_csv(imp_path, show_col_types = FALSE) |>
    transmute(gene, imputed = classification != "Complete")
} else {
  message("  Note: imputation classification not found, all marked non-imputed")
  tibble(gene = character(), imputed = logical())
}

scatter_df <- dep_df |>
  transmute(gene,
            logFC_TY = logFC_Training_CRE, logFC_TO = logFC_Training_PLA,
            pi_TY = pi_score_Training_CRE, pi_TO = pi_score_Training_PLA) |>
  filter(!is.na(logFC_TY), !is.na(logFC_TO)) |>
  left_join(imputation_df, by = "gene") |>
  mutate(
    imputed   = replace_na(imputed, FALSE),
    sig_class = classify_proteins_f2(pi_TY, pi_TO, NULL),
    is_sig    = sig_class != "NS",
    quadrant  = case_when(
      logFC_TY > 0 & logFC_TO > 0 ~ "Concordant Up",
      logFC_TY < 0 & logFC_TO < 0 ~ "Concordant Down",
      logFC_TY > 0 & logFC_TO < 0 ~ "Discordant (CRE Up / PLA Down)",
      TRUE                         ~ "Discordant (CRE Down / PLA Up)"))

universe <- scatter_df$gene
message(sprintf("  Total proteins: %d | Significant: %d",
                nrow(scatter_df), sum(scatter_df$is_sig)))

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
    mutate(set = set_name, pathway_label = clean_pathway_name(pathway),
           neg_log10_padj = -log10(padj), significant = padj < 0.05) |>
    arrange(desc(neg_log10_padj)) |>
    slice_head(n = N_SHOW)
}

message("\n--- Quadrant ORA (threshold-free) ---")
ora_q1 <- run_set_ora(scatter_df$gene[scatter_df$quadrant == "Concordant Up"],
                       "Concordant Up")
ora_q2 <- run_set_ora(scatter_df$gene[scatter_df$quadrant == "Discordant (CRE Down / PLA Up)"],
                       "Discordant (CRE Down / PLA Up)")
ora_q3 <- run_set_ora(scatter_df$gene[scatter_df$quadrant == "Concordant Down"],
                       "Concordant Down")
ora_q4 <- run_set_ora(scatter_df$gene[scatter_df$quadrant == "Discordant (CRE Up / PLA Down)"],
                       "Discordant (CRE Up / PLA Down)")

all_quad_ora <- bind_rows(ora_q1, ora_q2, ora_q3, ora_q4)
if (nrow(all_quad_ora) > 0)
  write_csv(all_quad_ora, file.path(DAT, "panel_A", "ora_quadrant.csv"))

# Scatter panel
xlim_range <- c(-3.1, 3.1)
ylim_range <- c(-2.8, 2.8)

ns_df  <- filter(scatter_df, sig_class == "NS")
sig_df <- filter(scatter_df, sig_class != "NS")

q_df <- scatter_df |>
  mutate(q = case_when(
    logFC_TY > 0 & logFC_TO > 0 ~ "Q1",
    logFC_TY < 0 & logFC_TO < 0 ~ "Q3",
    logFC_TY > 0 & logFC_TO < 0 ~ "Q4",
    TRUE ~ "Q2"))
q_counts <- q_df |> count(q) |> deframe()
q_sig    <- q_df |> filter(sig_class != "NS") |> count(q) |> deframe()
for (qq in c("Q1","Q2","Q3","Q4")) if (is.na(q_sig[qq])) q_sig[qq] <- 0

label_df <- sig_df |>
  group_by(sig_class) |>
  arrange(desc(abs(logFC_TY) + abs(logFC_TO))) |>
  slice_head(n = 5) |>
  ungroup() |>
  mutate(label_fill     = SIG_LABEL_FILL_F2[as.character(sig_class)],
         label_text_col = SIG_LABEL_TEXT_F2[as.character(sig_class)])

txt_gene <- scale_text(BASE_GENE, 190) * 0.82 + 1  # +1pt for print legibility
txt_quad <- scale_text(BASE_QUADRANT, 190) * 0.88

# Center-axis tick labels (behind points, at x=0 / y=0)
x_breaks <- seq(-3, 3, 1)
y_breaks <- seq(-2, 2, 1)
x_tick_df <- tibble(x = x_breaks[x_breaks != 0], y = 0,
                    label = as.character(x_breaks[x_breaks != 0]))
y_tick_df <- tibble(x = 0, y = y_breaks[y_breaks != 0],
                    label = as.character(y_breaks[y_breaks != 0]))

p_scatter <- ggplot(mapping = aes(x = logFC_TY, y = logFC_TO)) +
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
  geom_point(data = ns_df, color = "grey80", fill = "grey85", shape = 21,
             size = 0.35, alpha = 0.3, stroke = 0.10) +
  geom_point(data = sig_df, aes(fill = sig_class), shape = 21,
             size = ifelse(sig_df$sig_class == "NS", 0.6, 0.9),
             color = ifelse(sig_df$imputed, "black", "grey75"),
             alpha = case_when(
               sig_df$sig_class == "NS"          ~ 0.30,
               sig_df$sig_class == "Interaction" ~ 0.55,
               sig_df$sig_class == "Sig Both"    ~ 0.75,
               TRUE                               ~ 0.85),
             stroke = ifelse(sig_df$sig_class == "NS", 0.4, 0.6)) +
  scale_fill_manual(values = SIG_COLORS_F2, name = "Significance") +
  geom_label_repel(data = label_df, aes(label = gene),
                   fill = label_df$label_fill, color = label_df$label_text_col,
                   size = txt_gene, fontface = "italic", max.overlaps = 50,
                   segment.size = 0.3, segment.color = "grey25",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.3, point.padding = 0.2,
                   force = 6, force_pull = 0.3,
                   label.padding = unit(1, "pt"), label.r = unit(0.5, "pt"),
                   linewidth = 0.15, seed = 42,
                   xlim = c(-3, 3) * 0.85, ylim = c(-2.7, 2.7) * 0.85) +
  # Quadrant labels — stacked: title + counts, corner-aligned
  annotate("label", x = xlim_range[2], y = ylim_range[2],
           label = sprintf("Concordant Up\n%s/%s", q_sig["Q1"], q_counts["Q1"]),
           hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
           color = COMP_RED, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt"), lineheight = 0.9) +
  annotate("label", x = xlim_range[1], y = ylim_range[1],
           label = sprintf("%s/%s\nConcordant Down", q_sig["Q3"], q_counts["Q3"]),
           hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
           color = COMP_RED, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt"), lineheight = 0.9) +
  annotate("label", x = xlim_range[1], y = ylim_range[2],
           label = sprintf("Discordant (CREDn PLAUp)\n%s/%s", q_sig["Q2"], q_counts["Q2"]),
           hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
           color = COMP_BLUE, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt"), lineheight = 0.9) +
  annotate("label", x = xlim_range[2], y = ylim_range[1],
           label = sprintf("%s/%s\nDiscordant (CREUp PLADn)", q_sig["Q4"], q_counts["Q4"]),
           hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
           color = COMP_BLUE, fill = alpha("white", 0.92),
           label.padding = unit(2.5, "pt"), lineheight = 0.9) +
  # Axis titles
  annotate("text", x = 2.5, y = 0,
           label = expression(log[2]*FC ~ "(Training CRE)"),
           hjust = 0.5, vjust = -0.4, size = 1.3 * PRINT_SCALE, color = "grey30", fontface = "bold") +
  annotate("text", x = 0, y = 2.0,
           label = expression(log[2]*FC ~ "(Training PLA)"),
           hjust = 0.5, vjust = -0.4, size = 1.3 * PRINT_SCALE, color = "grey30", fontface = "bold",
           angle = 90) +
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

# Custom Significance key
key_lvls <- c("Sig Both", "Interaction", "Sig CRE only", "Sig PLA only")
key_display <- c("Sig Both", "Interaction", "Sig CRE", "Sig PLA")
key_df   <- tibble(
  category = factor(key_lvls, levels = key_lvls),
  display  = key_display,
  fill_col = unname(SIG_COLORS_F2[key_lvls]),
  x        = c(1.25, 1.80, 2.45, 3.05),
  y        = 0
)
p_key <- ggplot(key_df, aes(x = x, y = y)) +
  geom_point(aes(fill = category), shape = 21, size = 2.5 * PRINT_SCALE,
             color = "grey50", stroke = 0.6, alpha = 0.85,
             show.legend = FALSE) +
  geom_text(aes(label = display), nudge_x = 0.06, hjust = 0,
            size = 2.0 * PRINT_SCALE, fontface = "bold", color = "grey25") +
  scale_fill_manual(values = setNames(key_df$fill_col, key_df$category)) +
  scale_x_continuous(limits = c(0.2, 4.5),
                     expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.15, 0.15), expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(plot.margin = margin(-24, 0, 0, 0, "mm"))

# Half-bar builder
make_half_bars <- function(df, fill_color, side, ylim,
                            display_labels = character(0)) {
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
      display_name = ifelse(pathway_label %in% names(display_labels),
                            display_labels[pathway_label], pathway_label),
      star_raw = sig_stars(padj),
      star = star_raw)

  x_max <- max(bars$neg_log10_padj)
  x_display_max <- x_max * 1.18
  is_upper <- ylim[1] >= 0
  brk_fn <- function(limits) {
    b <- scales::pretty_breaks(n = 3)(limits)
    b[b != 0]
  }

  bars <- bars |>
    mutate(
      label_inside = neg_log10_padj >= x_max * 0.50,
      label_x      = ifelse(label_inside,
                            neg_log10_padj * 0.5,
                            neg_log10_padj + x_max * 0.03),
      label_hjust  = ifelse(label_inside, 0.5, 0),
      label_color  = ifelse(label_inside,
                            ifelse(significant, "white", "grey15"),
                            "grey20"),
      text_size    = scale_text(BASE_PATHWAY, 190) * 0.80)

  star_x_mult <- if (side == "left") 0.12 else 0.035

  p <- ggplot(bars, aes(y = y)) +
    geom_rect(aes(xmin = 0, xmax = neg_log10_padj,
                  ymin = y - bar_h / 2, ymax = y + bar_h / 2),
              fill = bars$bar_fill, color = "black", linewidth = 0.3) +
    geom_text(aes(x = label_x, y = y, label = display_name),
              hjust = bars$label_hjust, size = bars$text_size,
              fontface = "bold", color = bars$label_color, lineheight = 0.85) +
    geom_text(aes(x = neg_log10_padj + x_max * star_x_mult, label = star),
              hjust = 0, vjust = 0.5,
              size = 2.2 * PRINT_SCALE, fontface = "bold", color = "black",
              lineheight = 1.0) +
    labs(x = if (!is_upper) expression(-log[10](p[adj])) else NULL,
         y = NULL) +
    theme_minimal(base_size = 9) +
    theme(panel.grid   = element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x  = element_text(size = FIG_AXIS_TEXT, face = "bold",
                                       margin = margin(t = 0, unit = "mm")),
          axis.title.x = if (!is_upper) element_text(size = FIG_AXIS_TEXT, face = "bold",
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
                         breaks = brk_fn,
                         expand = expansion(mult = c(0, 0))) +
        scale_y_continuous(limits = ylim, expand = c(0, 0)) +
        coord_cartesian(clip = "off")
  } else {
    p + scale_x_continuous(limits = c(0, x_display_max),
                            breaks = brk_fn,
                            expand = expansion(mult = c(0, 0))) +
        scale_y_continuous(limits = ylim, expand = c(0, 0)) +
        coord_cartesian(clip = "off")
  }
}

# 4 half-bar panels (ylim matches scatter ±2.8)
p_ul <- make_half_bars(ora_q2, scales::alpha(COMP_BLUE, 0.30), "left",  c(0, 2.8),
                        display_labels = DISPLAY_LABELS_F04)
p_ll <- make_half_bars(ora_q3, scales::alpha(COMP_RED, 0.30),  "left",  c(-2.8, 0),
                        display_labels = DISPLAY_LABELS_F04)
p_ur <- make_half_bars(ora_q1, scales::alpha(COMP_RED, 0.30),  "right", c(0, 2.8),
                        display_labels = DISPLAY_LABELS_F04)
p_lr <- make_half_bars(ora_q4, scales::alpha(COMP_BLUE, 0.30), "right", c(-2.8, 0),
                        display_labels = DISPLAY_LABELS_F04)

# Composite
design <- c(
  area(1, 1),          # p_ul (top-left ORA bars)
  area(1, 2, 2, 2),   # p_scatter (rows 1-2, center)
  area(1, 3),          # p_ur (top-right ORA bars)
  area(2, 1),          # p_ll (bottom-left ORA bars)
  area(2, 3),          # p_lr (bottom-right ORA bars)
  area(3, 1, 3, 3)    # key spans full width below scatter, centered
)
n_total  <- nrow(scatter_df)
n_sig    <- sum(scatter_df$is_sig)
n_enrich <- if (nrow(all_quad_ora) > 0) sum(all_quad_ora$significant) else 0L
r_spear  <- cor(scatter_df$logFC_TY, scatter_df$logFC_TO, use = "complete.obs",
                method = "spearman")
sig_sub  <- scatter_df |> filter(is_sig)
r_pear   <- cor(sig_sub$logFC_TY, sig_sub$logFC_TO, use = "complete.obs")

composite <- p_ul + p_scatter + p_ur + p_ll + p_lr + p_key +
  plot_layout(design = design,
              widths  = c(70, 100, 70) / 240,
              heights = c(85, 85, 8) / 178) +
  plot_annotation(
    title    = "Training Concordance: Quadrant ORA",
    subtitle = sprintf("Threshold-free ORA (hypergeometric) | N = %d | %d DEPs (Pi < 0.05) | %d enriched (FDR < 0.05) | rho(all) = %.2f, r(sig) = %.2f",
                        n_total, n_sig, n_enrich, r_spear, r_pear),
    theme    = theme(plot.title    = element_text(size = FIG_TITLE_SIZE, face = "bold", hjust = 0),
                     plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, hjust = 0, color = "grey30"),
                     plot.title.position = "panel"))

COMP_W <- 200; COMP_H <- 120
ggsave(file.path(RPT_PNG, "MAIN_panel_A_ORA_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "MAIN_panel_A_ORA_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)

composite <- composite &
  labs(title = NULL, subtitle = NULL, tag = NULL) &
  theme(legend.position = "none")
composite <- composite +
  plot_annotation(title = NULL, subtitle = NULL,
                  theme = theme(plot.title = element_blank(),
                                plot.subtitle = element_blank()))

message("\nF04 Panel A ORA composite done")
