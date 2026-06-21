# Reversal Panel D: fry Rotation Test (bidirectional)
#
# Test 1 (original): Cancer-sig gene sets tested on Training_CR contrast
#   "Do cancer-altered proteins reverse with training?"
#   Subsets to CR subjects only, blocks on Subject_ID
#
# Test 2 (reciprocal): Training-sig gene sets tested on Cancer_vs_Healthy contrast
#   "Do training-responsive proteins show enrichment in the cancer signature?"
#   Uses ALL samples (CR + Healthy), blocks on Subject_ID
#
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

library(tidyverse)
library(limma)
library(fgsea)
library(patchwork)

set.seed(42)

RPT_PNG     <- "04_Figures/Reversal/b_reports/main/png/panels"
RPT_PDF     <- "04_Figures/Reversal/b_reports/main/pdf/panels"
RPT_SUP_PNG <- "04_Figures/Reversal/b_reports/supp/png/panels"
RPT_SUP_PDF <- "04_Figures/Reversal/b_reports/supp/pdf/panels"
DAT         <- "04_Figures/Reversal/c_data"
PANEL_W     <- 178

dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUP_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUP_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_D_fry"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# ═══════════════════════════════════════════════════════════════════════════════
# SHARED: Load data
# ═══════════════════════════════════════════════════════════════════════════════
source("04_Figures/Reversal/a_script/reversal_inputs.R")   # dal (imp4p) + dep_df

mat  <- dal$data
meta <- as.data.frame(dal$metadata)
imp_ids <- rownames(mat)

# Circularity diagnostic (shared)
circ_r <- cor(dep_df$t_Cancer_vs_Healthy, dep_df$t_Training_CR,
              use = "complete.obs")
message(sprintf("Circularity: r(t_CvH, t_TR) = %.3f", circ_r))

# Shared helpers
running_es <- function(t_vals, in_set) {
  n <- length(t_vals); n_h <- sum(in_set)
  if (n_h == 0) return(rep(0, n))
  hit_w <- ifelse(in_set, abs(t_vals), 0)
  miss_w <- 1 / (n - n_h)
  cumsum(ifelse(in_set, hit_w / sum(hit_w), -miss_w))
}

run_fry_generic <- function(mat_sub, design_sub, contrast_vec, block_sub, cor_sub,
                             idx, set_name) {
  if (length(idx) < 3) return(tibble(set = set_name, n = length(idx),
                                      direction = NA_character_,
                                      PValue = NA_real_, PValue.Mixed = NA_real_))
  res <- fry(mat_sub, index = idx, design = design_sub,
             contrast = contrast_vec, block = block_sub,
             correlation = cor_sub)
  tibble(set = set_name, n = length(idx), direction = res$Direction[1],
         PValue = res$PValue[1], PValue.Mixed = res$PValue.Mixed[1])
}

make_barcode <- function(t_df, in_col, es_col, fry_row, title, color,
                          n_total, stat_corner = "topright") {
  marks <- t_df %>% filter(.data[[in_col]])
  is_sig <- !is.na(fry_row$PValue) && fry_row$PValue < 0.05
  line_color <- if (is_sig) color else scales::alpha(color, 0.4)

  p_label <- sprintf("fry %s, %s (n = %d)%s",
                      fry_row$direction, fmt_p(fry_row$PValue),
                      fry_row$n,
                      if (fry_row$consistent) "" else " \u2717")
  p_color <- if (fry_row$consistent) "grey20" else "#DC2626"

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
             size = 2.8, fontface = "bold", color = p_color) +
    labs(title = title, x = NULL, y = "ES") +
    scale_x_continuous(limits = c(1, n_total), expand = c(0.005, 0)) +
    FIG_THEME +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = FIG_AXIS_TEXT, face = "bold",
                                       margin = margin(r = 1)),
          plot.title    = element_text(size = 10, face = "bold", color = "grey15",
                                       margin = margin(b = 0.5, unit = "mm")),
          plot.margin   = margin(0, 1, 0, 0, "mm"))

  p_bc <- ggplot(marks, aes(x = rank, xend = rank, y = 0, yend = 1)) +
    geom_segment(color = line_color, linewidth = 0.3, alpha = 0.7) +
    scale_x_continuous(limits = c(1, n_total), expand = c(0.005, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    FIG_THEME +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey97"),
          plot.margin = margin(0, 1, 0, 2, "mm"))

  list(es = p_es, bc = p_bc)
}

shorten_ora_label <- function(x, max_chars = 22) {
  x <- gsub("Reference ", "", x)
  x <- gsub("Regulation Of ", "Reg. ", x)
  x <- gsub("Negative Regulation Of ", "Neg. Reg. ", x)
  x <- gsub("Positive Regulation Of ", "Pos. Reg. ", x)
  x <- gsub("Catabolic Process", "Catabolism", x)
  x <- gsub("Metabolic Process", "Metabolism", x)
  x <- gsub("Signalling Events", "Signaling", x)
  x <- gsub("Oxidative Phosphorylation", "OxPhos", x)
  x <- gsub("Extracellular Matrix", "ECM", x)
  x <- gsub("Organization", "Org.", x)
  x <- gsub("Aminoacylation", "Aminoacyl.", x)
  ifelse(nchar(x) > max_chars, stringr::str_trunc(x, max_chars), x)
}

make_flanking_ora <- function(ora_df, set_label, bar_color) {
  if (is.null(ora_df) || nrow(ora_df) == 0) {
    return(ggplot() + theme_void() +
             annotate("text", x = 0.5, y = 0.5, label = "No sig. pathways",
                      size = 3, color = "grey60"))
  }
  bars <- ora_df %>%
    slice_head(n = 5) %>%
    mutate(neg_log_padj = -log10(pmax(padj, 1e-20)),
           significant  = padj < 0.05,
           bar_fill     = ifelse(significant,
                                 scales::alpha(bar_color, 0.85),
                                 scales::alpha(bar_color, 0.30)),
           short_label  = shorten_ora_label(pathway_label),
           star = sig_stars(padj),
           y = rev(row_number()),
           bar_h = 0.85)

  x_max <- max(bars$neg_log_padj, na.rm = TRUE)
  x_display_max <- x_max * 1.15

  bars <- bars %>%
    mutate(label_inside = neg_log_padj >= x_max * 0.55,
           label_x      = ifelse(label_inside,
                                 neg_log_padj * 0.5,
                                 neg_log_padj + x_max * 0.03),
           label_hjust  = ifelse(label_inside, 0.5, 0),
           label_color  = ifelse(label_inside,
                                 ifelse(significant, "white", "grey15"),
                                 "grey20"),
           text_size    = 2.3)

  ggplot(bars, aes(y = y)) +
    geom_rect(aes(xmin = 0, xmax = neg_log_padj,
                  ymin = y - bar_h / 2, ymax = y + bar_h / 2),
              fill = bars$bar_fill, color = "black", linewidth = 0.3) +
    geom_text(aes(x = label_x, y = y, label = short_label),
              hjust = bars$label_hjust, size = bars$text_size, fontface = "bold",
              color = bars$label_color, lineheight = 0.85) +
    geom_text(aes(x = neg_log_padj + x_max * 0.01, label = star),
              hjust = 0, vjust = 0.5, size = 2.3, fontface = "bold",
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
          plot.title    = element_text(face = "bold", size = 9, hjust = 0.5),
          plot.margin   = margin(2, 1, 0, 0, "mm"))
}

# Shared pathway collection
pw_collection <- build_pathway_collection(min_size = 10, max_size = 500,
                                           include_goslim = TRUE,
                                           exclude_variants = TRUE)
all_genes <- unique(dep_df$gene[dep_df$uniprot_id %in% imp_ids])

# ═══════════════════════════════════════════════════════════════════════════════
# TEST 1: Cancer-sig sets → Training_CR contrast (CR subjects only)
# ═══════════════════════════════════════════════════════════════════════════════
message("\n--- Test 1: Cancer sets on Training_CR rank ---")

cr_idx  <- which(meta$Group_Time != "H_T1")
mat_cr  <- mat[, cr_idx]
meta_cr <- meta[cr_idx, ]
meta_cr$Group_Time <- factor(meta_cr$Group_Time,
                              levels = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2"))

design_cr <- model.matrix(~ 0 + Group_Time, data = meta_cr)
colnames(design_cr) <- gsub("^Group_Time", "", colnames(design_cr))

corfit_cr <- duplicateCorrelation(mat_cr, design_cr, block = meta_cr$Subject_ID)
cor_cr <- corfit_cr$consensus.correlation
message(sprintf("  Within-subject cor (CR): %.4f", cor_cr))

cm_tr <- makeContrasts(
  Training_CR = (CRE_T2 + PLA_T2) / 2 - (CRE_T1 + PLA_T1) / 2,
  levels = design_cr
)

# Cancer-sig gene sets
imp_ids_cr <- rownames(mat_cr)
sig_cancer <- dep_df %>%
  filter(pi_score_Cancer_vs_Healthy < 0.05, uniprot_id %in% imp_ids_cr)

sets_cancer <- list(
  up   = na.omit(match(sig_cancer$uniprot_id[sig_cancer$logFC_Cancer_vs_Healthy > 0], imp_ids_cr)),
  down = na.omit(match(sig_cancer$uniprot_id[sig_cancer$logFC_Cancer_vs_Healthy < 0], imp_ids_cr)),
  up_ids   = sig_cancer$uniprot_id[sig_cancer$logFC_Cancer_vs_Healthy > 0],
  down_ids = sig_cancer$uniprot_id[sig_cancer$logFC_Cancer_vs_Healthy < 0]
)

message(sprintf("  Cancer gene sets: up = %d, down = %d",
                length(sets_cancer$up), length(sets_cancer$down)))

fry_c_up <- run_fry_generic(mat_cr, design_cr, cm_tr[, "Training_CR"],
                              meta_cr$Subject_ID, cor_cr,
                              sets_cancer$up, "cancer_up") %>%
  mutate(expected = "Down", consistent = direction == expected, test = "cancer_on_training")
fry_c_dn <- run_fry_generic(mat_cr, design_cr, cm_tr[, "Training_CR"],
                              meta_cr$Subject_ID, cor_cr,
                              sets_cancer$down, "cancer_down") %>%
  mutate(expected = "Up", consistent = direction == expected, test = "cancer_on_training")

message(sprintf("  cancer_up  -> Training: fry %s, p = %.4f (exp Down) %s",
                fry_c_up$direction, fry_c_up$PValue,
                ifelse(fry_c_up$consistent, "OK", "UNEXPECTED")))
message(sprintf("  cancer_down -> Training: fry %s, p = %.4f (exp Up) %s",
                fry_c_dn$direction, fry_c_dn$PValue,
                ifelse(fry_c_dn$consistent, "OK", "UNEXPECTED")))

# Barcode rank data (Training_CR)
t_rank_tr <- dep_df %>%
  filter(uniprot_id %in% imp_ids_cr, !is.na(t_Training_CR)) %>%
  arrange(desc(t_Training_CR)) %>%
  mutate(rank = row_number(),
         in_up   = uniprot_id %in% sets_cancer$up_ids,
         in_down = uniprot_id %in% sets_cancer$down_ids)

t_rank_tr$es_up   <- running_es(t_rank_tr$t_Training_CR, t_rank_tr$in_up)
t_rank_tr$es_down <- running_es(t_rank_tr$t_Training_CR, t_rank_tr$in_down)
n_tr <- nrow(t_rank_tr)

# Driving proteins (cancer-up reversed by training)
driving_up_t1 <- dep_df %>%
  filter(uniprot_id %in% sets_cancer$up_ids, uniprot_id %in% imp_ids_cr,
         t_Training_CR < 0) %>%
  transmute(gene, uniprot_id, set = "cancer_up", test = "cancer_on_training",
            t_cancer = t_Cancer_vs_Healthy, t_training = t_Training_CR,
            logFC_Cancer_vs_Healthy, logFC_Training_CR, pi_score_Cancer_vs_Healthy)
driving_dn_t1 <- dep_df %>%
  filter(uniprot_id %in% sets_cancer$down_ids, uniprot_id %in% imp_ids_cr,
         t_Training_CR > 0) %>%
  transmute(gene, uniprot_id, set = "cancer_down", test = "cancer_on_training",
            t_cancer = t_Cancer_vs_Healthy, t_training = t_Training_CR,
            logFC_Cancer_vs_Healthy, logFC_Training_CR, pi_score_Cancer_vs_Healthy)

# Leading-edge ORA for Test 1
ora_t1_up <- if (nrow(driving_up_t1) >= 5) {
  tryCatch(
    run_ora_deduplicated(genes = unique(driving_up_t1$gene), universe = all_genes,
                          pathways = pw_collection, jaccard_cutoff = 0.5,
                          min_size = 10, max_size = 500, padj_cutoff = 0.25) %>%
      mutate(pathway_label = clean_pathway_name(pathway)) %>%
      slice_head(n = 5),
    error = function(e) tibble())
} else tibble()

ora_t1_dn <- if (nrow(driving_dn_t1) >= 5) {
  tryCatch(
    run_ora_deduplicated(genes = unique(driving_dn_t1$gene), universe = all_genes,
                          pathways = pw_collection, jaccard_cutoff = 0.5,
                          min_size = 10, max_size = 500, padj_cutoff = 0.25) %>%
      mutate(pathway_label = clean_pathway_name(pathway)) %>%
      slice_head(n = 5),
    error = function(e) tibble())
} else tibble()

# ═══════════════════════════════════════════════════════════════════════════════
# TEST 2: Training-sig sets → Cancer_vs_Healthy contrast (ALL samples)
# ═══════════════════════════════════════════════════════════════════════════════
message("\n--- Test 2: Training sets on Cancer_vs_Healthy rank ---")

meta_all <- meta
meta_all$Group_Time <- factor(meta_all$Group_Time,
                               levels = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1"))

design_all <- model.matrix(~ 0 + Group_Time, data = meta_all)
colnames(design_all) <- gsub("^Group_Time", "", colnames(design_all))

corfit_all <- duplicateCorrelation(mat, design_all, block = meta_all$Subject_ID)
cor_all <- corfit_all$consensus.correlation
message(sprintf("  Within-subject cor (all): %.4f", cor_all))

cm_cvh <- makeContrasts(
  Cancer_vs_Healthy = (CRE_T1 + PLA_T1) / 2 - H_T1,
  levels = design_all
)

# Training-sig gene sets
sig_training <- dep_df %>%
  filter(pi_score_Training_CR < 0.05, uniprot_id %in% imp_ids)

sets_training <- list(
  up   = na.omit(match(sig_training$uniprot_id[sig_training$logFC_Training_CR > 0], imp_ids)),
  down = na.omit(match(sig_training$uniprot_id[sig_training$logFC_Training_CR < 0], imp_ids)),
  up_ids   = sig_training$uniprot_id[sig_training$logFC_Training_CR > 0],
  down_ids = sig_training$uniprot_id[sig_training$logFC_Training_CR < 0]
)

message(sprintf("  Training gene sets: up = %d, down = %d",
                length(sets_training$up), length(sets_training$down)))

# Expected: training_up → cancer Down (cancer suppresses what training restores)
#           training_down → cancer Up (cancer elevates what training reduces)
fry_t_up <- run_fry_generic(mat, design_all, cm_cvh[, "Cancer_vs_Healthy"],
                              meta_all$Subject_ID, cor_all,
                              sets_training$up, "training_up") %>%
  mutate(expected = "Down", consistent = direction == expected, test = "training_on_cancer")
fry_t_dn <- run_fry_generic(mat, design_all, cm_cvh[, "Cancer_vs_Healthy"],
                              meta_all$Subject_ID, cor_all,
                              sets_training$down, "training_down") %>%
  mutate(expected = "Up", consistent = direction == expected, test = "training_on_cancer")

message(sprintf("  training_up  -> Cancer: fry %s, p = %.4f (exp Down) %s",
                fry_t_up$direction, fry_t_up$PValue,
                ifelse(fry_t_up$consistent, "OK", "UNEXPECTED")))
message(sprintf("  training_down -> Cancer: fry %s, p = %.4f (exp Up) %s",
                fry_t_dn$direction, fry_t_dn$PValue,
                ifelse(fry_t_dn$consistent, "OK", "UNEXPECTED")))

# Barcode rank data (Cancer_vs_Healthy)
t_rank_cvh <- dep_df %>%
  filter(uniprot_id %in% imp_ids, !is.na(t_Cancer_vs_Healthy)) %>%
  arrange(desc(t_Cancer_vs_Healthy)) %>%
  mutate(rank = row_number(),
         in_up   = uniprot_id %in% sets_training$up_ids,
         in_down = uniprot_id %in% sets_training$down_ids)

t_rank_cvh$es_up   <- running_es(t_rank_cvh$t_Cancer_vs_Healthy, t_rank_cvh$in_up)
t_rank_cvh$es_down <- running_es(t_rank_cvh$t_Cancer_vs_Healthy, t_rank_cvh$in_down)
n_cvh <- nrow(t_rank_cvh)

# Driving proteins (training-up appearing cancer-down = reversed)
driving_up_t2 <- dep_df %>%
  filter(uniprot_id %in% sets_training$up_ids, uniprot_id %in% imp_ids,
         t_Cancer_vs_Healthy < 0) %>%
  transmute(gene, uniprot_id, set = "training_up", test = "training_on_cancer",
            t_cancer = t_Cancer_vs_Healthy, t_training = t_Training_CR,
            logFC_Cancer_vs_Healthy, logFC_Training_CR, pi_score_Training_CR)
driving_dn_t2 <- dep_df %>%
  filter(uniprot_id %in% sets_training$down_ids, uniprot_id %in% imp_ids,
         t_Cancer_vs_Healthy > 0) %>%
  transmute(gene, uniprot_id, set = "training_down", test = "training_on_cancer",
            t_cancer = t_Cancer_vs_Healthy, t_training = t_Training_CR,
            logFC_Cancer_vs_Healthy, logFC_Training_CR, pi_score_Training_CR)

# --- Test 2: GO Slim composition bars (small n, ORA underpowered) ---
# Source go_slim ONLY if not already loaded (AnnotationDbi masking)
if (!exists("assign_go_slim_consolidated")) {
  source("04_Figures/shared/go_slim_categories.R")
}

make_goslim_bars <- function(gene_set, set_label, bar_color) {
  if (length(gene_set) < 3) {
    return(ggplot() + theme_void() +
             annotate("text", x = 0.5, y = 0.5, label = paste0("n = ", length(gene_set)),
                      size = 3, color = "grey60"))
  }
  slim <- assign_go_slim_consolidated(gene_set, all_genes)
  counts <- slim %>%
    count(consolidated, name = "n") %>%
    filter(consolidated != "Other", n > 0) %>%
    arrange(desc(n)) %>%
    slice_head(n = 5) %>%
    mutate(y = rev(row_number()), bar_h = 0.85)

  if (nrow(counts) == 0) {
    return(ggplot() + theme_void() +
             annotate("text", x = 0.5, y = 0.5, label = "No GO Slim",
                      size = 3, color = "grey60"))
  }

  ggplot(counts, aes(y = y)) +
    geom_rect(aes(xmin = 0, xmax = n,
                  ymin = y - bar_h / 2, ymax = y + bar_h / 2),
              fill = CONSOLIDATED_COLORS[counts$consolidated],
              color = "black", linewidth = 0.3) +
    geom_text(aes(x = n / 2, y = y,
                  label = stringr::str_trunc(as.character(consolidated), 20)),
              size = 2.2, fontface = "bold", color = "white", lineheight = 0.85) +
    geom_text(aes(x = n + max(n) * 0.05, y = y, label = n),
              hjust = 0, size = 2.2, fontface = "bold", color = "grey30") +
    labs(title = set_label, x = "Proteins", y = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    scale_y_continuous(limits = c(0.3, nrow(counts) + 0.7), expand = c(0, 0)) +
    theme_minimal(base_size = 9) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = FIG_AXIS_TEXT),
          axis.title.x = element_text(size = FIG_AXIS_TEXT),
          axis.line.x = element_line(color = "grey40", linewidth = 0.3),
          plot.title = element_text(face = "bold", size = 9, hjust = 0.5),
          plot.margin = margin(2, 1, 0, 0, "mm"))
}

# Build GO Slim bars for Test 2 driving proteins
p_ora_t2_up <- make_goslim_bars(unique(driving_up_t2$gene), "GO Slim: Tr-Up", COLOR_TRAINING)
p_ora_t2_dn <- make_goslim_bars(unique(driving_dn_t2$gene), "GO Slim: Tr-Down", COLOR_TRAINING)

# ═══════════════════════════════════════════════════════════════════════════════
# EFFECTIVE SIGNATURE SIZE (CAMERA formula + roastgsa cross-check)
# Wu & Smyth 2012 NAR (PMID 22638577): n_eff = n / (1 + (n-1) * r_bar)
# Caballé-Mestres et al. 2023 BMC Bioinform (PMID 37904108): roastgsa
#
# When proteins in a set are correlated, the "effective" number of
# independent observations is smaller than the nominal set size.
# This inflates rotation test power if unaccounted for.
# ═══════════════════════════════════════════════════════════════════════════════
compute_neff <- function(mat_sub, idx, label) {
  if (length(idx) < 3) return(list(n = length(idx), r_bar = NA, n_eff = length(idx)))
  set_mat <- mat_sub[idx, , drop = FALSE]
  # Mean pairwise Pearson correlation across samples
  cor_mat <- cor(t(set_mat), use = "pairwise.complete.obs")
  # Extract upper triangle (exclude diagonal)
  r_bar <- mean(cor_mat[upper.tri(cor_mat)], na.rm = TRUE)
  n <- length(idx)
  # CAMERA formula: variance inflation factor = 1 + (n-1)*r_bar
  vif <- 1 + (n - 1) * max(r_bar, 0)  # floor at 0 (negative r_bar → no inflation)
  n_eff <- n / vif
  message(sprintf("  %s: n = %d, mean r = %.3f, n_eff = %.1f (%.0f%% of nominal)",
                  label, n, r_bar, n_eff, 100 * n_eff / n))
  list(n = n, r_bar = round(r_bar, 4), n_eff = round(n_eff, 1))
}

message("\n--- Effective Signature Size ---")
neff_c_up <- compute_neff(mat_cr, sets_cancer$up, "Cancer-Up (Test 1)")
neff_c_dn <- compute_neff(mat_cr, sets_cancer$down, "Cancer-Down (Test 1)")
neff_t_up <- compute_neff(mat, sets_training$up, "Training-Up (Test 2)")
neff_t_dn <- compute_neff(mat, sets_training$down, "Training-Down (Test 2)")

neff_df <- tibble(
  set = c("cancer_up", "cancer_down", "training_up", "training_down"),
  n_nominal = c(neff_c_up$n, neff_c_dn$n, neff_t_up$n, neff_t_dn$n),
  mean_pairwise_r = c(neff_c_up$r_bar, neff_c_dn$r_bar, neff_t_up$r_bar, neff_t_dn$r_bar),
  n_effective = c(neff_c_up$n_eff, neff_c_dn$n_eff, neff_t_up$n_eff, neff_t_dn$n_eff),
  pct_of_nominal = round(100 * c(neff_c_up$n_eff, neff_c_dn$n_eff,
                                   neff_t_up$n_eff, neff_t_dn$n_eff) /
                           c(neff_c_up$n, neff_c_dn$n, neff_t_up$n, neff_t_dn$n), 1)
)
write_csv(neff_df, file.path(DAT, "panel_D_fry", "effective_signature_size.csv"))

# ═══════════════════════════════════════════════════════════════════════════════
# COMBINED EXPORTS
# ═══════════════════════════════════════════════════════════════════════════════
fry_all <- bind_rows(fry_c_up, fry_c_dn, fry_t_up, fry_t_dn) %>%
  mutate(cor_within_cr = cor_cr, cor_within_all = cor_all, circularity_r = circ_r,
         n_eff = c(neff_c_up$n_eff, neff_c_dn$n_eff, neff_t_up$n_eff, neff_t_dn$n_eff),
         mean_pairwise_r = c(neff_c_up$r_bar, neff_c_dn$r_bar,
                              neff_t_up$r_bar, neff_t_dn$r_bar))
write_csv(fry_all, file.path(DAT, "panel_D_fry", "fry_results_all.csv"))

driving_df <- bind_rows(driving_up_t1, driving_dn_t1, driving_up_t2, driving_dn_t2)
write_csv(driving_df, file.path(DAT, "panel_D_fry", "driving_proteins.csv"))

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALISATION: 2-column composite
# Left column:  Test 1 (cancer sets → Training_CR rank)
# Right column: Test 2 (training sets → Cancer_vs_Healthy rank)
# ═══════════════════════════════════════════════════════════════════════════════
COLOR_CANCER   <- unname(CONTRAST_COLORS["Cancer_vs_Healthy"])
COLOR_TRAINING <- unname(CONTRAST_COLORS["Training_CR"])

# --- Test 1 barcodes (left) ---
p1a <- make_barcode(t_rank_tr, "in_up", "es_up", fry_c_up,
                     sprintf("Cancer-Up (n = %d, n_eff = %.0f)",
                             length(sets_cancer$up), neff_c_up$n_eff),
                     COLOR_CANCER, n_tr, stat_corner = "bottomleft")
p1b <- make_barcode(t_rank_tr, "in_down", "es_down", fry_c_dn,
                     sprintf("Cancer-Down (n = %d, n_eff = %.0f)%s",
                             length(sets_cancer$down), neff_c_dn$n_eff,
                             if (fry_c_dn$PValue > 0.05) "  (n.s.)" else ""),
                     COLOR_CANCER, n_tr, stat_corner = "topright")

p_t1 <- ggplot(t_rank_tr, aes(x = rank, y = t_Training_CR)) +
  geom_area(fill = scales::alpha(COLOR_TRAINING, 0.20),
            color = COLOR_TRAINING, linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  labs(x = sprintf("Rank (Tr. CR t-stat, n = %d)", n_tr), y = NULL) +
  scale_x_continuous(limits = c(1, n_tr), expand = c(0.005, 0)) +
  FIG_THEME +
  theme(axis.title.x = element_text(size = 7, face = "bold"), axis.title.y = element_blank(),
        plot.margin = margin(0, 1, 1, 0, "mm"))

p_ora_t1_up <- make_flanking_ora(ora_t1_up, "Reversed Cancer-Up", DIR_COLORS["Up"])
p_ora_t1_dn <- make_flanking_ora(ora_t1_dn, "Reversed Cancer-Down", DIR_COLORS["Down"])

# --- Test 2 barcodes (right) ---
p2a <- make_barcode(t_rank_cvh, "in_up", "es_up", fry_t_up,
                     sprintf("Training-Up (n = %d)", length(sets_training$up)),
                     COLOR_TRAINING, n_cvh, stat_corner = "bottomleft")
p2b <- make_barcode(t_rank_cvh, "in_down", "es_down", fry_t_dn,
                     sprintf("Training-Down (n = %d)%s", length(sets_training$down),
                             if (fry_t_dn$PValue > 0.05) "  (n.s.)" else ""),
                     COLOR_TRAINING, n_cvh, stat_corner = "topright")

p_t2 <- ggplot(t_rank_cvh, aes(x = rank, y = t_Cancer_vs_Healthy)) +
  geom_area(fill = scales::alpha(COLOR_CANCER, 0.20),
            color = COLOR_CANCER, linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  labs(x = sprintf("Rank (Cancer t-stat, n = %d)", n_cvh), y = NULL) +
  scale_x_continuous(limits = c(1, n_cvh), expand = c(0.005, 0)) +
  FIG_THEME +
  theme(axis.title.x = element_text(size = 7, face = "bold"), axis.title.y = element_blank(),
        plot.margin = margin(0, 1, 1, 0, "mm"))

# p_ora_t2_up and p_ora_t2_dn already built above via make_goslim_bars()

# --- Build 2-column composite ---
# Each column: 2 ES curves + 2 barcodes + 1 t-stat + 2 ORA bars
fry_design_bidir <- c(
  # Column 1: Test 1
  area(1, 1, 1, 1),  area(2, 1, 2, 1),  # ES + barcode (cancer-up)
  area(3, 1, 3, 1),  area(4, 1, 4, 1),  # ES + barcode (cancer-down)
  area(5, 1, 5, 1),                       # t-stat rank (Training)
  area(1, 2, 2, 2),  area(3, 2, 4, 2),  # ORA bars (Test 1)
  # Column 2: Test 2
  area(1, 3, 1, 3),  area(2, 3, 2, 3),  # ES + barcode (training-up)
  area(3, 3, 3, 3),  area(4, 3, 4, 3),  # ES + barcode (training-down)
  area(5, 3, 5, 3),                       # t-stat rank (Cancer)
  area(1, 4, 2, 4),  area(3, 4, 4, 4)   # ORA bars (Test 2)
)

n_all <- max(n_tr, n_cvh)
cor_imp <- cor_cr  # backward compat for stitcher snapshot

pD_subtitle_full <- sprintf(
  "circ r = %.3f | dupCor(CR) = %.3f, dupCor(all) = %.3f | n = %d | n_eff: Up=%.0f, Dn=%.0f",
  circ_r, cor_cr, cor_all, n_all, neff_c_up$n_eff, neff_c_dn$n_eff)

pD_fry <- p1a$es + p1a$bc + p1b$es + p1b$bc + p_t1 +
  p_ora_t1_up + p_ora_t1_dn +
  p2a$es + p2a$bc + p2b$es + p2b$bc + p_t2 +
  p_ora_t2_up + p_ora_t2_dn +
  plot_layout(design = fry_design_bidir,
              heights = c(2.0, 0.25, 2.0, 0.25, 0.6),
              widths  = c(1.5, 1.3, 1.5, 1.3)) +
  plot_annotation(
    title = "fry: Cancer Recovery Reversal (bidirectional)",
    subtitle = pD_subtitle_full,
    theme = theme(
      plot.title = element_text(size = FIG_TITLE_SIZE, face = "bold", hjust = 0,
                                margin = margin(l = 6, unit = "mm")),
      plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, color = "grey30", hjust = 0,
                                   margin = margin(l = 6, unit = "mm")),
      plot.title.position = "panel"))

ggsave(file.path(RPT_PNG, "MAIN_panel_D_fry.png"), pD_fry,
       width = PANEL_W * 2, height = 175, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "MAIN_panel_D_fry.pdf"), pD_fry,
       width = PANEL_W * 2, height = 175, units = "mm", device = pdf_device)

# --- Export for composite (strip titles) ---
pD_title    <- "fry: Reversal (bidirectional)"
pD_subtitle <- pD_subtitle_full
pD_legend   <- NULL
pD_fry <- pD_fry &
  labs(title = NULL, subtitle = NULL, tag = NULL) &
  theme(legend.position = "none")
pD_fry <- pD_fry +
  plot_annotation(title = NULL, subtitle = NULL,
                  theme = theme(plot.title = element_blank(),
                                plot.subtitle = element_blank()))

# Backward compat
fry_up <- fry_c_up
fry_dn <- fry_c_dn

message("Reversal Panel D (fry bidirectional) done")
