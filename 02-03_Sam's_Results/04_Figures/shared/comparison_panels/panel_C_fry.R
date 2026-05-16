# Shared Panel C — fry barcode + flanking ORA; used by F04 (concordance) and F05 (reversal)
# fry is an exact rotation-based gene-set test (Wu & Smyth 2010, Bioinformatics)

source("02-03_Sam's_Results/04_Figures/shared/style.R")
source("02-03_Sam's_Results/04_Figures/shared/print_scale_apply_380mm.R")
source("02-03_Sam's_Results/04_Figures/shared/pathway_utils.R")

library(tidyverse)
library(limma)
library(fgsea)
library(patchwork)

set.seed(42)

RPT_PNG     <- cfg$rpt_png
RPT_PDF     <- cfg$rpt_pdf
RPT_SUP_PNG <- cfg$rpt_sup_png
RPT_SUP_PDF <- cfg$rpt_sup_pdf
DAT         <- cfg$dat
PANEL_W     <- cfg$panel_w %||% 178

dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUP_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUP_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_C_fry"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- grDevices::pdf  # force base pdf (cairo DLL fails)

# Sam CvH: load normalized DAList + combined DEP results
# cfg$dal_rds: path to Sam's normalized DAList (has $data, $metadata)
# cfg$dep_combined_csv: path to combined results CSV for this figure
# cfg$test_contrast_col: name of the t-statistic column for the test contrast (e.g. "t_Training_PLA")
# cfg$test_contrast_design: the contrast level expression for makeContrasts (e.g. "supp_CRE_T2 - supp_CRE_T1")

sam_dal <- readRDS("02-03_Sam's_Results/00_input/01_normalized_DAList_SURV_stringent_muscle.RDS")
dep_df  <- read_csv(cfg$dep_combined_csv, show_col_types = FALSE)

meta <- sam_dal$metadata
mat_imp <- sam_dal$data   # already normalized matrix; rows = uniprot_id, cols = samples

# Sam design: supp_time (CRE_T1, CRE_T2, PLA_T1, PLA_T2) for SURV group
# Only SURV samples have T1+T2; CTL (healthy) have no pre/post
surv_meta <- meta[meta$cancer == "SURV", , drop = FALSE]
# Match matrix columns to SURV metadata
surv_samples <- intersect(colnames(mat_imp), surv_meta$sample_id)
surv_meta <- surv_meta[match(surv_samples, surv_meta$sample_id), , drop = FALSE]
# Keep only complete rows (fry requires finite, non-NA values)
mat_surv <- mat_imp[, surv_samples, drop = FALSE]
complete_rows <- rowSums(is.na(mat_surv) | !is.finite(mat_surv)) == 0
mat_surv <- mat_surv[complete_rows, , drop = FALSE]
message(sprintf("mat_surv: %d proteins x %d SURV samples (%d complete rows)",
                nrow(mat_surv), ncol(mat_surv), sum(complete_rows)))

surv_meta$supp_time <- factor(surv_meta$supp_time,
  levels = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2"))
design <- model.matrix(~ 0 + supp_time, data = surv_meta)
colnames(design) <- gsub("^supp_time", "supp_", colnames(design))

block_id <- surv_meta$pid

corfit_imp <- duplicateCorrelation(mat_surv, design, block = block_id)
cor_imp <- corfit_imp$consensus.correlation
message(sprintf("Within-subject cor: %.4f", cor_imp))

# cfg$test_contrast_expr: the contrast expression string, e.g. "supp_PLA_T2 - supp_PLA_T1"
# defaults to PLA training (F04) or CRE training (F05)
test_contrast_expr <- cfg$test_contrast_expr %||% "supp_PLA_T2 - supp_PLA_T1"
cm <- eval(parse(text = sprintf('makeContrasts("%s", levels = design)', test_contrast_expr)))
colnames(cm) <- cfg$contrast_test

circ_r <- NULL
if (isTRUE(cfg$has_circularity)) {
  circ_src_col  <- paste0("t_", cfg$contrast_source)
  circ_test_col <- paste0("t_", cfg$contrast_test)
  circ_r <- cor(dep_df[[circ_src_col]], dep_df[[circ_test_col]], use = "complete.obs")
  message(sprintf("Circularity: r(%s, %s) = %.3f", cfg$contrast_source, cfg$contrast_test, circ_r))
}

# Use SURV matrix rows for indexing (fry uses mat_surv which is SURV-only)
imp_ids <- rownames(mat_surv)
pi_col  <- paste0("pi_score_", cfg$contrast_source)
lfc_col <- paste0("logFC_", cfg$contrast_source)

sig_pi <- dep_df |>
  filter(.data[[pi_col]] < 0.05, uniprot_id %in% imp_ids)

sets_pi <- list(
  up       = match(sig_pi$uniprot_id[sig_pi[[lfc_col]] > 0], imp_ids),
  down     = match(sig_pi$uniprot_id[sig_pi[[lfc_col]] < 0], imp_ids),
  up_ids   = sig_pi$uniprot_id[sig_pi[[lfc_col]] > 0],
  down_ids = sig_pi$uniprot_id[sig_pi[[lfc_col]] < 0]
)

message(sprintf("Gene sets (Pi < 0.05): up = %d, down = %d",
                length(sets_pi$up), length(sets_pi$down)))

run_fry_set <- function(idx, set_name) {
  if (length(idx) < 3) return(tibble(set = set_name, n = length(idx),
                                      direction = NA_character_,
                                      PValue = NA_real_, PValue.Mixed = NA_real_))
  res <- fry(mat_surv, index = idx, design = design,
             contrast = cm[, cfg$contrast_test], block = block_id,
             correlation = cor_imp)
  tibble(set = set_name, n = length(idx), direction = res$Direction[1],
         PValue = res$PValue[1], PValue.Mixed = res$PValue.Mixed[1])
}

set_up_name <- paste0(cfg$set_prefix, "_up")
set_dn_name <- paste0(cfg$set_prefix, "_down")

fry_up <- run_fry_set(sets_pi$up, set_up_name) |>
  mutate(expected = cfg$expected_up, consistent = direction == expected)
fry_dn <- run_fry_set(sets_pi$down, set_dn_name) |>
  mutate(expected = cfg$expected_down, consistent = direction == expected)
fry_all <- bind_rows(fry_up, fry_dn) |>
  mutate(cor_within = cor_imp)

# Add circularity column if applicable
if (!is.null(circ_r)) fry_all$circularity_r <- circ_r

write_csv(fry_all, file.path(DAT, "panel_C_fry", "fry_results_all.csv"))

t_src_col  <- paste0("t_", cfg$contrast_source)
t_test_col <- paste0("t_", cfg$contrast_test)
lfc_src_col  <- paste0("logFC_", cfg$contrast_source)
lfc_test_col <- paste0("logFC_", cfg$contrast_test)
pi_src_col   <- paste0("pi_score_", cfg$contrast_source)

driving_up_filter <- if (cfg$driving_up_sign == "pos") {
  quote(t_test > 0)
} else {
  quote(t_test < 0)
}
driving_dn_filter <- if (cfg$driving_dn_sign == "pos") {
  quote(t_test > 0)
} else {
  quote(t_test < 0)
}

driving_up <- dep_df |>
  filter(uniprot_id %in% sets_pi$up_ids, uniprot_id %in% imp_ids) |>
  mutate(t_test = .data[[t_test_col]]) |>
  filter(eval(driving_up_filter)) |>
  transmute(gene, uniprot_id, set = set_up_name,
            t_source = .data[[t_src_col]], t_test = .data[[t_test_col]],
            !!lfc_src_col := .data[[lfc_src_col]],
            !!lfc_test_col := .data[[lfc_test_col]],
            !!pi_src_col := .data[[pi_src_col]])

driving_dn <- dep_df |>
  filter(uniprot_id %in% sets_pi$down_ids, uniprot_id %in% imp_ids) |>
  mutate(t_test = .data[[t_test_col]]) |>
  filter(eval(driving_dn_filter)) |>
  transmute(gene, uniprot_id, set = set_dn_name,
            t_source = .data[[t_src_col]], t_test = .data[[t_test_col]],
            !!lfc_src_col := .data[[lfc_src_col]],
            !!lfc_test_col := .data[[lfc_test_col]],
            !!pi_src_col := .data[[pi_src_col]])

driving_df <- bind_rows(driving_up, driving_dn)
write_csv(driving_df, file.path(DAT, "panel_C_fry", "driving_proteins.csv"))

pw_collection <- build_pathway_collection(min_size = 10, max_size = 500,
                                           include_goslim = TRUE,
                                           exclude_variants = TRUE)
all_genes <- dep_df$gene[dep_df$uniprot_id %in% imp_ids]

ora_leading_up <- if (nrow(driving_up) >= 5) {
  tryCatch(
    run_ora_deduplicated(genes = driving_up$gene, universe = all_genes,
                          pathways = pw_collection, jaccard_cutoff = 0.5,
                          min_size = 10, max_size = 500, padj_cutoff = 0.1) |>
      mutate(pathway_label = clean_pathway_name(pathway)) |>
      slice_head(n = 5),
    error = function(e) tibble())
} else tibble()

ora_leading_dn <- if (nrow(driving_dn) >= 5) {
  tryCatch(
    run_ora_deduplicated(genes = driving_dn$gene, universe = all_genes,
                          pathways = pw_collection, jaccard_cutoff = 0.5,
                          min_size = 10, max_size = 500, padj_cutoff = 0.1) |>
      mutate(pathway_label = clean_pathway_name(pathway)) |>
      slice_head(n = 5),
    error = function(e) tibble())
} else tibble()

t_rank <- dep_df |>
  filter(uniprot_id %in% imp_ids, !is.na(.data[[t_test_col]])) |>
  arrange(desc(.data[[t_test_col]])) |>
  mutate(rank = row_number(),
         in_up   = uniprot_id %in% sets_pi$up_ids,
         in_down = uniprot_id %in% sets_pi$down_ids)

running_es <- function(t_vals, in_set) {
  n <- length(t_vals); n_h <- sum(in_set)
  if (n_h == 0) return(rep(0, n))
  hit_w <- ifelse(in_set, abs(t_vals), 0)
  miss_w <- 1 / (n - n_h)
  cumsum(ifelse(in_set, hit_w / sum(hit_w), -miss_w))
}

t_rank$es_up   <- running_es(t_rank[[t_test_col]], t_rank$in_up)
t_rank$es_down <- running_es(t_rank[[t_test_col]], t_rank$in_down)

n_all <- nrow(t_rank)
txt_s <- scale_text(BASE_STAT, PANEL_W)

make_barcode <- function(t_df, in_col, es_col, fry_row, title, color,
                          stat_corner = "topright") {
  marks <- t_df |> filter(.data[[in_col]])
  is_sig <- !is.na(fry_row$PValue) && fry_row$PValue < 0.05
  line_color <- if (is_sig) color else scales::alpha(color, 0.4)

  consistent <- !is.na(fry_row$consistent) && fry_row$consistent
  p_label <- sprintf("fry %s, %s (n = %d)%s",
                      ifelse(is.na(fry_row$direction), "n/a", fry_row$direction),
                      fmt_p(fry_row$PValue),
                      fry_row$n,
                      if (consistent) "" else " \u2717")
  p_color <- if (consistent) "grey20" else "#DC2626"

  # Corner position for stat annotation
  stat_x <- if (stat_corner == "topright") Inf else -Inf
  stat_y <- if (stat_corner == "topright") Inf else -Inf
  stat_hjust <- if (stat_corner == "topright") 1.05 else -0.05
  stat_vjust <- if (stat_corner == "topright") 1.5 else -0.5

  p_es <- ggplot(t_df, aes(x = rank, y = .data[[es_col]])) +
    geom_area(fill = scales::alpha(line_color, 0.15), color = NA) +
    geom_line(color = line_color, linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60",
               linewidth = 0.3) +
    annotate("text", x = stat_x, y = stat_y, label = p_label,
             hjust = stat_hjust, vjust = stat_vjust,
             size = 1.6 * PRINT_SCALE, fontface = "bold", color = p_color) +
    labs(title = title, x = NULL, y = "ES") +
    scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    FIG_THEME +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = FIG_AXIS_TEXT, face = "bold",
                                       margin = margin(r = 1)),
          plot.title    = element_text(size = FIG_TITLE_SIZE, face = "bold", color = "grey15",
                                       margin = margin(b = 0.5, unit = "mm")),
          plot.margin   = margin(0, 1, 0, 0, "mm"))

  p_bc <- ggplot(marks, aes(x = rank, xend = rank, y = 0, yend = 1)) +
    geom_segment(color = line_color, linewidth = 0.3, alpha = 0.7) +
    scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    FIG_THEME +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey97"),
          plot.margin = margin(0, 1, 0, 2, "mm"))

  list(es = p_es, bc = p_bc)
}

up_title <- sprintf(cfg$up_title_fmt, length(sets_pi$up))
dn_ns_suffix <- if (is.na(fry_dn$PValue) || fry_dn$PValue > 0.05) "  (n.s.)" else ""
dn_title <- sprintf(cfg$dn_title_fmt, length(sets_pi$down), dn_ns_suffix)

FIG_COLOR <- cfg$fig_color

p1 <- make_barcode(t_rank, "in_up", "es_up", fry_up, up_title, FIG_COLOR,
                    stat_corner = cfg$stat_corner_up)
p2 <- make_barcode(t_rank, "in_down", "es_down", fry_dn, dn_title, FIG_COLOR,
                    stat_corner = cfg$stat_corner_dn)

T_STAT_COLOR <- unname(CONTRAST_COLORS[cfg$contrast_test] %||% CONTRAST_COLORS["Training_Old"])
p_t <- ggplot(t_rank, aes(x = rank, y = .data[[t_test_col]])) +
  geom_area(fill = scales::alpha(T_STAT_COLOR, 0.20),
            color = T_STAT_COLOR, linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  labs(x = sprintf("Rank (%s t-stat, n = %%d)", cfg$x_axis_label %||% cfg$contrast_test) |>
         sprintf(n_all),
       y = NULL) +
  scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
  FIG_THEME +
  theme(axis.title.x = element_text(size = FIG_AXIS_TEXT, face = "bold",
                                     margin = margin(t = 0)),
        axis.title.y = element_blank(),
        plot.margin = margin(0, 1, 1, 0, "mm"))

if (!is.null(circ_r)) {
  p_t <- p_t +
    annotate("text", x = n_all * 0.98, y = Inf,
             label = sprintf("Circularity: r(%s, %s) = %.3f",
                             cfg$contrast_source, cfg$contrast_test, circ_r),
             hjust = 1, vjust = 1.5, size = txt_s * 0.75, color = "grey50")
}

shorten_ora_label <- function(x, max_chars = 28, explicit_map = NULL,
                               mode = "wrap") {
  x <- gsub("Reference ", "", x)
  if (!is.null(explicit_map)) {
    idx <- x %in% names(explicit_map)
    x[idx] <- explicit_map[x[idx]]
  }
  x <- gsub("Regulation Of ", "Reg. ", x)
  x <- gsub("Negative Regulation Of ", "Neg. Reg. ", x)
  x <- gsub("Positive Regulation Of ", "Pos. Reg. ", x)
  x <- gsub("Epithelial Mesenchymal Transition", "EMT", x)
  x <- gsub("Catabolic Process", "Catabolism", x)
  x <- gsub("Metabolic Process", "Metabolism", x)
  x <- gsub(" Involved In ", " in ", x)
  x <- gsub("Regulated Microtubule Minus End Directed Transport", "MT Transport", x)
  if (mode == "wrap") {
    ifelse(nchar(x) > max_chars,
           stringr::str_wrap(x, width = max_chars),
           x)
  } else {
    ifelse(nchar(x) > max_chars,
           paste0(substr(x, 1, max_chars - 1), "\u2026"),
           x)
  }
}

make_flanking_ora <- function(ora_df, set_label, bar_color,
                               label_map = NULL, force_inside = NULL,
                               label_mode = "wrap") {
  if (is.null(ora_df) || nrow(ora_df) == 0) {
    return(ggplot() + theme_void() +
             annotate("text", x = 0.5, y = 0.5, label = "No sig. pathways",
                      size = 3, color = "grey60"))
  }
  bars <- ora_df |>
    slice_head(n = 5) |>
    mutate(neg_log_padj = -log10(pmax(padj, 1e-20)),
           significant  = padj < 0.05,
           bar_fill     = ifelse(significant,
                                 scales::alpha(bar_color, 0.85),
                                 scales::alpha(bar_color, 0.30)),
           short_label  = shorten_ora_label(pathway_label,
                                             explicit_map = label_map,
                                             mode = label_mode),
                    star_raw = sig_stars(padj),
           star = star_raw,
           y = rev(row_number()),
           bar_h = 0.85)

  x_max <- max(bars$neg_log_padj, na.rm = TRUE)
  x_display_max <- x_max * 1.15

  inside_check <- bars$neg_log_padj >= x_max * 0.55
  if (!is.null(force_inside)) {
    inside_check <- inside_check | bars$short_label %in% force_inside
  }
  bars <- bars |>
    mutate(label_inside = inside_check,
           label_x      = ifelse(label_inside,
                                 neg_log_padj * 0.5,
                                 neg_log_padj + x_max * 0.03),
           label_hjust  = ifelse(label_inside, 0.5, 0),
           label_color  = ifelse(label_inside,
                                 ifelse(significant, "white", "grey15"),
                                 "grey20"),
           text_size    = 1.6 * PRINT_SCALE)

  ggplot(bars, aes(y = y)) +
    geom_rect(aes(xmin = 0, xmax = neg_log_padj,
                  ymin = y - bar_h / 2, ymax = y + bar_h / 2),
              fill = bars$bar_fill, color = "black", linewidth = 0.3) +
    geom_text(aes(x = label_x, y = y, label = short_label),
              hjust = bars$label_hjust, size = bars$text_size, fontface = "bold",
              color = bars$label_color, lineheight = 0.85) +
    geom_text(aes(x = neg_log_padj + x_max * 0.01, label = star),
              hjust = 0, vjust = 0.5, size = 2.5 * PRINT_SCALE, fontface = "bold",
              color = "black", lineheight = 1.0) +
    labs(title = set_label, x = expression(-log[10](p[adj])), y = NULL) +
    scale_x_continuous(limits = c(0, x_display_max),
                       breaks = scales::pretty_breaks(n = 3),
                       expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = c(0.3, nrow(bars) + 0.7), expand = c(0, 0)) +
    theme_minimal(base_size = 9) +
    theme(panel.grid    = element_blank(),
          axis.text.y   = element_blank(),
          axis.ticks.y  = element_blank(),
          axis.title.y  = element_blank(),
          axis.text.x   = element_text(size = FIG_AXIS_TEXT),
          axis.title.x  = element_text(size = FIG_AXIS_TEXT),
          axis.line.x   = element_line(color = "grey40", linewidth = 0.3),
          plot.title    = element_text(face = "bold", size = FIG_TITLE_SIZE, hjust = 0.5),
          plot.margin   = margin(2, 1, 0, 0, "mm"))
}

p_ora_flank_up <- make_flanking_ora(ora_leading_up, cfg$ora_flank_up_label,
                                     unname(DIR_COLORS["Up"]),
                                     label_map = cfg$label_map,
                                     force_inside = cfg$force_inside_labels,
                                     label_mode = cfg$long_label_mode %||% "wrap")
p_ora_flank_dn <- make_flanking_ora(ora_leading_dn, cfg$ora_flank_dn_label,
                                     unname(DIR_COLORS["Down"]),
                                     label_map = cfg$label_map,
                                     force_inside = cfg$force_inside_labels,
                                     label_mode = cfg$long_label_mode %||% "wrap")

fry_design <- c(
  area(1, 1, 1, 1),  # p1$es  (row 1, col 1)
  area(2, 1, 2, 1),  # p1$bc  (row 2, col 1)
  area(3, 1, 3, 1),  # p2$es  (row 3, col 1)
  area(4, 1, 4, 1),  # p2$bc  (row 4, col 1)
  area(5, 1, 5, 1),  # p_t    (row 5, col 1)
  area(1, 2, 2, 2),  # ora_up (rows 1-2, col 2)
  area(3, 2, 4, 2)   # ora_dn (rows 3-4, col 2)
)

pD_subtitle_full <- if (!is.null(circ_r)) {
  sprintf(cfg$subtitle_fmt, circ_r, cor_imp, n_all)
} else {
  sprintf(cfg$subtitle_fmt, cor_imp, n_all)
}

pC_fry <- p1$es + p1$bc + p2$es + p2$bc + p_t +
  p_ora_flank_up + p_ora_flank_dn +
  plot_layout(design = fry_design,
              heights = c(2.0, 0.25, 2.0, 0.25, 0.6),
              widths  = c(1.73, 1.58)) +
  plot_annotation(
    title = cfg$title,
    subtitle = pD_subtitle_full,
    theme = theme(plot.title = element_text(size = FIG_TITLE_SIZE, face = "bold", hjust = 0,
                                            margin = margin(l = 12, unit = "mm")),
                  plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, color = "grey30", hjust = 0,
                                               margin = margin(l = 12, unit = "mm")),
                  plot.title.position = "panel")
  )

ggsave(file.path(RPT_PNG, "MAIN_panel_C_fry.png"), pC_fry,
       width = PANEL_W + 80, height = 175, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "MAIN_panel_C_fry.pdf"), pC_fry,
       width = PANEL_W + 80, height = 175, units = "mm", device = pdf_device)

make_ora_bars <- function(ora_df, set_label, bar_color, show_xaxis = FALSE) {
  if (nrow(ora_df) == 0) return(NULL)

  bars <- ora_df |>
    slice_head(n = 3) |>
    mutate(neg_log_padj = -log10(pmax(padj, 1e-20)),
           star = sig_stars(padj),
           y = rev(row_number()),
           bar_h = 0.85,
           text_size = case_when(
             nchar(pathway_label) > 35 ~ 2.8 * PRINT_SCALE,
             nchar(pathway_label) > 25 ~ 3.4 * PRINT_SCALE,
             TRUE ~ 4.0 * PRINT_SCALE
           ))

  x_max <- max(bars$neg_log_padj, na.rm = TRUE)
  x_display_max <- x_max * 1.15

  ggplot(bars, aes(y = y)) +
    geom_rect(aes(xmin = 0, xmax = neg_log_padj,
                  ymin = y - bar_h / 2, ymax = y + bar_h / 2),
              fill = bar_color, color = NA) +
    geom_text(aes(x = x_max * 0.03, y = y, label = pathway_label),
              hjust = 0, size = bars$text_size, fontface = "bold",
              color = "white", lineheight = 0.85) +
    geom_text(aes(x = neg_log_padj + x_max * 0.03, label = star),
              hjust = 0, vjust = 0.5, size = 3.5 * PRINT_SCALE, fontface = "bold",
              color = "black") +
    annotate("segment", x = 0, xend = x_display_max, y = -Inf, yend = -Inf,
             color = "grey40", linewidth = 0.3) +
    labs(title = set_label,
         x = if (show_xaxis) expression(-log[10](p[adj])) else NULL,
         y = NULL) +
    scale_x_continuous(limits = c(0, x_display_max),
                       breaks = scales::pretty_breaks(n = 3),
                       expand = expansion(mult = c(0, 0))) +
    scale_y_continuous(limits = c(0.5, nrow(bars) + 0.5), expand = c(0, 0)) +
    theme_minimal(base_size = 9) +
    theme(panel.grid    = element_blank(),
          axis.text.y   = element_blank(),
          axis.ticks.y  = element_blank(),
          axis.title.y  = element_blank(),
          axis.text.x   = if (show_xaxis) element_text(size = FIG_AXIS_TEXT, face = "bold")
                          else element_blank(),
          axis.title.x  = if (show_xaxis) element_text(size = FIG_AXIS_TEXT, face = "bold")
                          else element_blank(),
          axis.line.x   = element_blank(),
          axis.ticks.x  = if (show_xaxis) element_line(color = "grey40", linewidth = 0.3)
                          else element_blank(),
          plot.title    = element_text(face = "bold", size = FIG_TITLE_SIZE, hjust = 0.5),
          plot.margin   = margin(2, 6, 2, 2, "mm"))
}

p_ora_up <- make_ora_bars(ora_leading_up, cfg$ora_supp_up_label,
                           unname(DIR_COLORS["Up"]), show_xaxis = FALSE)
p_ora_dn <- make_ora_bars(ora_leading_dn, cfg$ora_supp_dn_label,
                           unname(DIR_COLORS["Down"]), show_xaxis = TRUE)

if (!is.null(p_ora_up) || !is.null(p_ora_dn)) {
  ora_panels <- Filter(Negate(is.null), list(p_ora_up, p_ora_dn))
  p_ora <- wrap_plots(ora_panels, ncol = 1) +
    plot_annotation(
      title = cfg$ora_supp_title,
      subtitle = cfg$ora_supp_subtitle,
      theme = theme(plot.title    = element_text(size = FIG_TITLE_SIZE, face = "bold", hjust = 0.5),
                    plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, hjust = 0.5, color = "grey30")))

  ggsave(file.path(RPT_SUP_PNG, "SUPP_panel_C_fry_ora.png"), p_ora,
         width = 160, height = 100, units = "mm", dpi = 300)
  ggsave(file.path(RPT_SUP_PDF, "SUPP_panel_C_fry_ora.pdf"), p_ora,
         width = 160, height = 100, units = "mm", device = pdf_device)
  message(sprintf("%s Panel C ORA (supp) done", cfg$fig_id))
}

pD_title    <- cfg$title
pD_subtitle <- pD_subtitle_full
pD_legend   <- NULL
pC_fry <- pC_fry &
  labs(title = NULL, subtitle = NULL, tag = NULL) &
  theme(legend.position = "none")
pC_fry <- pC_fry +
  plot_annotation(title = NULL, subtitle = NULL,
                  theme = theme(plot.title = element_blank(),
                                plot.subtitle = element_blank()))

message(sprintf("%s Panel C (fry) done", cfg$fig_id))
