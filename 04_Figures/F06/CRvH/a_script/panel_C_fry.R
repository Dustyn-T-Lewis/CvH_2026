# F06 CRvH Panel C -- fry Rotation Test: CvH Concordance Barcode
# Tests whether Cancer_vs_Healthy-significant protein sets respond concordantly
# in the Training_CR contrast, using limma's fry rotation framework.
#
# No circularity: Cancer_vs_Healthy = CR_T1 vs H_T1 (baseline comparison);
# Training_CR = CR_T2 - CR_T1 (within-CR training effect).
# No shared contrast terms.
#
# Reference: Wu & Smyth 2010, Bioinformatics -- ROAST/fry
# ---------------------------------------------------------------------------
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F06/a_script/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(patchwork)
})

set.seed(42)

RPT <- "04_Figures/F06/CRvH/b_reports"
DAT <- "04_Figures/F06/CRvH/c_data"
dir.create(file.path(DAT, "panel_C_fry"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()
PF_W <- 220

# -- Step 1: Load data --------------------------------------------------------

dal      <- readRDS("02_Imputation/c_data/01_DAList_imputed.rds")
dep_df   <- read_csv("03_DEP/c_data/03_combined_results_CRvH.csv",
                      show_col_types = FALSE)
imp_csv  <- read_csv("02_Imputation/c_data/01_imputed.csv",
                      show_col_types = FALSE)

meta <- dal$metadata
sample_cols <- meta$Col_ID

# -- Step 2: Build imputed matrix ----------------------------------------------

mat_imp <- imp_csv %>%
  select(uniprot_id, all_of(sample_cols)) %>%
  column_to_rownames("uniprot_id") %>%
  as.matrix()

n_imp <- nrow(mat_imp)
message(sprintf("Imputed matrix: %d proteins x %d samples", n_imp, ncol(mat_imp)))

# -- Step 3: Design matrix + duplicateCorrelation -----------------------------
# CRvH model: Group x Timepoint with Subject blocking

meta$subject <- sub("_T[12]$", "", meta$Col_ID)
meta$Group_Time <- factor(paste0(meta$Group, "_", meta$Timepoint))
design <- model.matrix(~ 0 + Group_Time, data = meta)
colnames(design) <- gsub("^Group_Time", "", colnames(design))

block_id <- meta$subject

corfit_imp <- duplicateCorrelation(mat_imp, design, block = block_id)
cor_imp <- corfit_imp$consensus.correlation
message(sprintf("Within-subject cor: %.4f", cor_imp))

# Training_CR contrast: average(CR_T2) - average(CR_T1)
available_cols <- colnames(design)
message(sprintf("Design columns: %s", paste(available_cols, collapse = ", ")))

t2_cols <- available_cols[grepl("T2", available_cols)]
t1_cols <- available_cols[grepl("T1", available_cols) & !grepl("PPS|H_", available_cols)]

if (length(t2_cols) == 0 || length(t1_cols) == 0) {
  stop("Could not identify T1/T2 columns for CR subjects")
}

cm <- matrix(0, nrow = ncol(design), ncol = 1)
rownames(cm) <- colnames(design)
colnames(cm) <- "Training_CR"
for (col in t2_cols) cm[col, 1] <-  1 / length(t2_cols)
for (col in t1_cols) cm[col, 1] <- -1 / length(t1_cols)

# -- Step 4: Define Cancer_vs_Healthy gene sets --------------------------------

imp_ids <- rownames(mat_imp)

sig_pi <- dep_df %>%
  filter(pi_score_Cancer_vs_Healthy < 0.05, uniprot_id %in% imp_ids)

sets_pi <- list(
  up      = match(sig_pi$uniprot_id[sig_pi$logFC_Cancer_vs_Healthy > 0], imp_ids),
  down    = match(sig_pi$uniprot_id[sig_pi$logFC_Cancer_vs_Healthy < 0], imp_ids),
  up_ids  = sig_pi$uniprot_id[sig_pi$logFC_Cancer_vs_Healthy > 0],
  down_ids = sig_pi$uniprot_id[sig_pi$logFC_Cancer_vs_Healthy < 0]
)

message(sprintf("Gene sets (Pi < 0.05): up = %d, down = %d",
                length(sets_pi$up), length(sets_pi$down)))

# -- Step 5: Run fry -----------------------------------------------------------

run_fry_set <- function(idx, set_name) {
  if (length(idx) < 3) return(tibble(set = set_name, n = length(idx),
                                      direction = NA_character_,
                                      PValue = NA_real_, PValue.Mixed = NA_real_))
  res <- fry(mat_imp, index = idx, design = design,
             contrast = cm[, "Training_CR"], block = block_id,
             correlation = cor_imp)
  tibble(set = set_name, n = length(idx), direction = res$Direction[1],
         PValue = res$PValue[1], PValue.Mixed = res$PValue.Mixed[1])
}

fry_up <- run_fry_set(sets_pi$up, "cvh_up") %>%
  mutate(expected = "Up", consistent = direction == expected)
fry_dn <- run_fry_set(sets_pi$down, "cvh_down") %>%
  mutate(expected = "Down", consistent = direction == expected)
fry_all <- bind_rows(fry_up, fry_dn) %>%
  mutate(cor_within = cor_imp)

write_csv(fry_all, file.path(DAT, "panel_C_fry", "fry_results_all.csv"))

# -- Step 6: Driving proteins + leading-edge ORA -------------------------------

driving_up <- dep_df %>%
  filter(uniprot_id %in% sets_pi$up_ids, uniprot_id %in% imp_ids,
         t_Training_CR > 0) %>%
  transmute(gene, uniprot_id, set = "cvh_up",
            t_cancer_vs_healthy = t_Cancer_vs_Healthy, t_training_cr = t_Training_CR,
            logFC_Cancer_vs_Healthy, logFC_Training_CR, pi_score_Cancer_vs_Healthy)

driving_dn <- dep_df %>%
  filter(uniprot_id %in% sets_pi$down_ids, uniprot_id %in% imp_ids,
         t_Training_CR < 0) %>%
  transmute(gene, uniprot_id, set = "cvh_down",
            t_cancer_vs_healthy = t_Cancer_vs_Healthy, t_training_cr = t_Training_CR,
            logFC_Cancer_vs_Healthy, logFC_Training_CR, pi_score_Cancer_vs_Healthy)

driving_df <- bind_rows(driving_up, driving_dn)
write_csv(driving_df, file.path(DAT, "panel_C_fry", "driving_proteins.csv"))

# Leading-edge ORA on driving proteins
pw_collection <- build_pathway_collection(min_size = 10, max_size = 500,
                                           include_goslim = TRUE,
                                           exclude_variants = TRUE)
all_genes <- dep_df$gene[dep_df$uniprot_id %in% imp_ids]

ora_leading_up <- if (nrow(driving_up) >= 5) {
  tryCatch(
    run_ora_deduplicated(genes = driving_up$gene, universe = all_genes,
                          pathways = pw_collection, jaccard_cutoff = 0.5,
                          min_size = 10, max_size = 500, padj_cutoff = 0.1) %>%
      mutate(pathway_label = clean_pathway_name(pathway)) %>%
      slice_head(n = 3),
    error = function(e) tibble())
} else tibble()

ora_leading_dn <- if (nrow(driving_dn) >= 5) {
  tryCatch(
    run_ora_deduplicated(genes = driving_dn$gene, universe = all_genes,
                          pathways = pw_collection, jaccard_cutoff = 0.5,
                          min_size = 10, max_size = 500, padj_cutoff = 0.1) %>%
      mutate(pathway_label = clean_pathway_name(pathway)) %>%
      slice_head(n = 3),
    error = function(e) tibble())
} else tibble()

# -- Step 7: Barcode data -----------------------------------------------------

t_rank <- dep_df %>%
  filter(uniprot_id %in% imp_ids, !is.na(t_Training_CR)) %>%
  arrange(desc(t_Training_CR)) %>%
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

t_rank$es_up   <- running_es(t_rank$t_Training_CR, t_rank$in_up)
t_rank$es_down <- running_es(t_rank$t_Training_CR, t_rank$in_down)

n_all <- nrow(t_rank)
txt_s <- scale_text(BASE_STAT, PF_W)

# -- Step 8: Barcode visualization -- both Up and Down -------------------------

FIG_COLOR <- unname(DIR_COLORS["Up"])  # #D6604D

make_barcode <- function(t_df, in_col, es_col, fry_row, title, color,
                          p_position = "top_right") {
  marks <- t_df %>% filter(.data[[in_col]])
  is_sig <- !is.na(fry_row$PValue) && fry_row$PValue < 0.05
  line_color <- if (is_sig) color else scales::alpha(color, 0.4)

  p_label <- sprintf("fry %s, %s (n = %d)%s",
                      fry_row$direction, fmt_p(fry_row$PValue),
                      fry_row$n,
                      if (fry_row$consistent) "" else " \u2717")
  p_color <- if (fry_row$consistent) "grey20" else "#DC2626"

  if (p_position == "bottom_left") {
    p_x <- n_all * 0.02; p_y <- -Inf; p_hjust <- 0; p_vjust <- -0.5
  } else {
    p_x <- n_all * 0.98; p_y <- Inf; p_hjust <- 1; p_vjust <- 1.3
  }

  p_es <- ggplot(t_df, aes(x = rank, y = .data[[es_col]])) +
    geom_area(fill = scales::alpha(line_color, 0.15), color = NA) +
    geom_line(color = line_color, linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60",
               linewidth = 0.3) +
    annotate("text", x = p_x, y = p_y,
             label = p_label,
             hjust = p_hjust, vjust = p_vjust, size = txt_s * 1.1,
             fontface = "bold", color = p_color) +
    labs(y = "ES", title = title) +
    scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    FIG_THEME +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(3, 4, 0, 4, "mm"),
          plot.title = element_text(size = 9.5, face = "bold"))

  p_bc <- ggplot(marks, aes(x = rank, xend = rank, y = 0, yend = 1)) +
    geom_segment(color = line_color, linewidth = 0.3, alpha = 0.7) +
    scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    FIG_THEME +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey97"),
          plot.margin = margin(0, 4, 0, 4, "mm"))

  list(es = p_es, bc = p_bc)
}

up_title <- sprintf("CvH-Up DEPs (Pi < 0.05, n = %d) \u2192 Tr.(CR) ranked t",
                     length(sets_pi$up))
dn_title <- sprintf("CvH-Down DEPs (Pi < 0.05, n = %d) \u2192 Tr.(CR) ranked t%s",
                     length(sets_pi$down),
                     if (fry_dn$PValue > 0.05) "  (n.s.)" else "")

p1 <- make_barcode(t_rank, "in_up", "es_up", fry_up, up_title, FIG_COLOR,
                    p_position = "top_right")
p2 <- make_barcode(t_rank, "in_down", "es_down", fry_dn, dn_title, FIG_COLOR,
                    p_position = "bottom_left")

p_t <- ggplot(t_rank, aes(x = rank, y = t_Training_CR)) +
  geom_area(fill = scales::alpha(CONTRAST_COLORS["Training_CR"], 0.20),
            color = CONTRAST_COLORS["Training_CR"], linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  labs(x = sprintf("Protein rank by t(Training CR)  [n = %d]", n_all),
       y = "t-stat") +
  scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
  FIG_THEME +
  theme(plot.margin = margin(2, 4, 4, 4, "mm"))

# -- Flanking ORA bar builder (compact, for right column) ---------------------
make_flanking_ora <- function(ora_df, set_label, bar_color) {
  if (is.null(ora_df) || nrow(ora_df) == 0) {
    return(ggplot() + theme_void() +
             annotate("text", x = 0.5, y = 0.5, label = "No sig. pathways",
                      size = 3, color = "grey60"))
  }
  bars <- ora_df %>%
    slice_head(n = 5) %>%
    mutate(neg_log_padj = -log10(pmax(padj, 1e-20)),
           star = sig_stars(padj),
           y = rev(row_number()),
           bar_h = 0.7,
           text_size = case_when(nchar(pathway_label) > 30 ~ 2.4,
                                 nchar(pathway_label) > 20 ~ 2.8,
                                 TRUE ~ 3.2))
  x_max <- max(bars$neg_log_padj, na.rm = TRUE)
  x_display_max <- x_max * 1.18

  ggplot(bars, aes(y = y)) +
    geom_rect(aes(xmin = 0, xmax = neg_log_padj,
                  ymin = y - bar_h / 2, ymax = y + bar_h / 2),
              fill = bar_color, color = NA) +
    geom_text(aes(x = neg_log_padj / 2, y = y, label = pathway_label),
              hjust = 0.5, size = bars$text_size, fontface = "bold",
              color = "white", lineheight = 0.85) +
    geom_text(aes(x = neg_log_padj + x_max * 0.04, label = star),
              hjust = 0, vjust = 0.5, size = 3, fontface = "bold",
              color = "black") +
    labs(title = set_label, x = expression(-log[10](p[adj])), y = NULL) +
    scale_x_continuous(limits = c(0, x_display_max),
                       breaks = scales::pretty_breaks(n = 3),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_y_continuous(limits = c(0.3, nrow(bars) + 0.7), expand = c(0, 0)) +
    theme_minimal(base_size = 8) +
    theme(panel.grid    = element_blank(),
          axis.text.y   = element_blank(),
          axis.ticks.y  = element_blank(),
          axis.title.y  = element_blank(),
          axis.text.x   = element_text(size = 7),
          axis.title.x  = element_text(size = 8),
          plot.title    = element_text(face = "bold", size = 8.5, hjust = 0.5),
          plot.margin   = margin(2, 4, 2, 2, "mm"))
}

p_ora_flank_up <- make_flanking_ora(ora_leading_up, "Concordant (Up)",
                                     unname(DIR_COLORS["Up"]))
p_ora_flank_dn <- make_flanking_ora(ora_leading_dn, "Concordant (Down)",
                                     unname(DIR_COLORS["Down"]))

# -- Composite: fry (left) + flanking ORA bars (right) ------------------------
fry_design <- c(
  area(1, 1, 1, 1),  # p1$es
  area(2, 1, 2, 1),  # p1$bc
  area(3, 1, 3, 1),  # p2$es
  area(4, 1, 4, 1),  # p2$bc
  area(5, 1, 5, 1),  # p_t
  area(1, 2, 2, 2),  # ora_up
  area(3, 2, 4, 2)   # ora_dn
)

pC_fry <- p1$es + p1$bc + p2$es + p2$bc + p_t +
  p_ora_flank_up + p_ora_flank_dn +
  plot_layout(design = fry_design,
              heights = c(2.5, 0.4, 2.5, 0.4, 1.2),
              widths  = c(3, 1.2)) +
  plot_annotation(
    title = "fry Gene-Set Rotation Test: CvH Concordance",
    subtitle = sprintf("Rotation-based set test (accounts for inter-gene correlation) | dupCor = %.3f | n = %d proteins",
                        cor_imp, n_all),
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 8.5, color = "grey30"))
  )

ggsave(file.path(RPT, "panel_C_fry_MAIN.pdf"), pC_fry,
       width = PF_W + 80, height = 175, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_C_fry_MAIN.png"), pC_fry,
       width = PF_W + 80, height = 175, units = "mm", dpi = 300)

# -- SUPP: Leading-edge ORA bar chart -----------------------------------------

make_ora_bars <- function(ora_df, set_label, bar_color, show_xaxis = FALSE) {
  if (nrow(ora_df) == 0) return(NULL)

  bars <- ora_df %>%
    slice_head(n = 3) %>%
    mutate(neg_log_padj = -log10(pmax(padj, 1e-20)),
           star = sig_stars(padj),
           y = rev(row_number()),
           bar_h = 0.7,
           text_size = case_when(
             nchar(pathway_label) > 35 ~ 2.8,
             nchar(pathway_label) > 25 ~ 3.4,
             TRUE ~ 4.0
           ))

  x_max <- max(bars$neg_log_padj, na.rm = TRUE)
  x_display_max <- x_max * 1.15

  ggplot(bars, aes(y = y)) +
    geom_rect(aes(xmin = 0, xmax = neg_log_padj,
                  ymin = y - bar_h / 2, ymax = y + bar_h / 2),
              fill = bar_color, color = NA) +
    geom_text(aes(x = neg_log_padj / 2, y = y, label = pathway_label),
              hjust = 0.5, size = bars$text_size, fontface = "bold",
              color = "white", lineheight = 0.85) +
    geom_text(aes(x = neg_log_padj + x_max * 0.03, label = star),
              hjust = 0, vjust = 0.5, size = 3.5, fontface = "bold",
              color = "black") +
    annotate("segment", x = 0, xend = x_display_max, y = -Inf, yend = -Inf,
             color = "grey40", linewidth = 0.3) +
    labs(title = set_label,
         x = if (show_xaxis) expression(-log[10](p[adj])) else NULL,
         y = NULL) +
    scale_x_continuous(limits = c(0, x_display_max),
                       breaks = scales::pretty_breaks(n = 3),
                       expand = expansion(mult = c(0, 0.02))) +
    scale_y_continuous(limits = c(0.5, nrow(bars) + 0.5), expand = c(0, 0)) +
    theme_minimal(base_size = 9) +
    theme(panel.grid    = element_blank(),
          axis.text.y   = element_blank(),
          axis.ticks.y  = element_blank(),
          axis.title.y  = element_blank(),
          axis.text.x   = if (show_xaxis) element_text(size = 8.5, face = "bold")
                          else element_blank(),
          axis.title.x  = if (show_xaxis) element_text(size = 9.5, face = "bold")
                          else element_blank(),
          axis.line.x   = element_blank(),
          axis.ticks.x  = if (show_xaxis) element_line(color = "grey40", linewidth = 0.3)
                          else element_blank(),
          plot.title    = element_text(face = "bold", size = 10, hjust = 0.5),
          plot.margin   = margin(2, 6, 2, 2, "mm"))
}

p_ora_up <- make_ora_bars(ora_leading_up, "Up DEPs", unname(DIR_COLORS["Up"]),
                           show_xaxis = FALSE)
p_ora_dn <- make_ora_bars(ora_leading_dn, "Down DEPs", unname(DIR_COLORS["Down"]),
                           show_xaxis = TRUE)

if (!is.null(p_ora_up) || !is.null(p_ora_dn)) {
  ora_panels <- Filter(Negate(is.null), list(p_ora_up, p_ora_dn))
  p_ora <- wrap_plots(ora_panels, ncol = 1) +
    plot_annotation(
      title = "Leading-Edge ORA: fry Driving Proteins",
      subtitle = "Hypergeometric ORA on concordant driving proteins | top 3 per set",
      theme = theme(plot.title    = element_text(size = 12, face = "bold", hjust = 0.5),
                    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey30")))

  ggsave(file.path(RPT, "supp", "panel_C_fry_ora_SUPP.pdf"), p_ora,
         width = 160, height = 100, units = "mm", device = pdf_device)
  ggsave(file.path(RPT, "supp", "panel_C_fry_ora_SUPP.png"), p_ora,
         width = 160, height = 100, units = "mm", dpi = 300)
  message("Panel C ORA (supp) done")
}

message("F06 CRvH Panel C (fry) done")
