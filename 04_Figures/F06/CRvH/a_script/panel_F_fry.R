# F06 CRvH Panel F_fry -- fry Rotation Test: CvH Concordance Barcode
# Tests whether Cancer_vs_Healthy-significant proteins respond concordantly
# in the Training_CR contrast, using limma's fry rotation framework.
#
# No circularity: Cancer_vs_Healthy uses CR_T1 vs H_T1 subjects;
# Training_CR uses CR_T2 vs CR_T1 on same CR subjects.
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
dir.create(file.path(DAT, "panel_F_fry"), recursive = TRUE, showWarnings = FALSE)

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
# CRvH model: Group_Time has levels for CR and Healthy subjects
# Subject extracted from Col_ID: sub("_T[12]$", "", Col_ID)

meta$subject <- sub("_T[12]$", "", meta$Col_ID)

# For fry we need the Training_CR contrast: CR_T2 - CR_T1
# We need a design with group-level coding
meta$Group_Time <- factor(paste0(meta$Group, "_", meta$Timepoint))
design <- model.matrix(~ 0 + Group_Time, data = meta)
colnames(design) <- gsub("^Group_Time", "", colnames(design))

block_id <- meta$subject

corfit_imp <- duplicateCorrelation(mat_imp, design, block = block_id)
cor_imp <- corfit_imp$consensus.correlation
message(sprintf("Within-subject cor: %.4f", cor_imp))

# Training_CR contrast: CR_T2 - CR_T1
# Need to identify which column names exist in the design
cr_t2_col <- grep("CR.*T2|CR_CRE_T2|CR_PLA_T2", colnames(design), value = TRUE)
cr_t1_col <- grep("CR.*T1|CR_CRE_T1|CR_PLA_T1", colnames(design), value = TRUE)

# Build contrast based on available columns
# The CRvH model groups: CR_CRE, CR_PLA (both are CR), PPS (Healthy)
# Timepoints: T1, T2 for CR subjects; T1 only for Healthy
available_cols <- colnames(design)
message(sprintf("Design columns: %s", paste(available_cols, collapse = ", ")))

# The Training_CR contrast averages T2-T1 across CR subjects
# Build contrast string dynamically from available columns
t2_cols <- available_cols[grepl("T2", available_cols)]
t1_cols <- available_cols[grepl("T1", available_cols) & !grepl("PPS|H_", available_cols)]

if (length(t2_cols) > 0 && length(t1_cols) > 0) {
  # Average over CR subgroups
  t2_str <- paste0("(", paste(t2_cols, collapse = " + "), ") / ", length(t2_cols))
  t1_str <- paste0("(", paste(t1_cols, collapse = " + "), ") / ", length(t1_cols))
  cm <- makeContrasts(
    Training_CR = list(t2_str, t1_str),
    levels = design
  )
  # If makeContrasts with list fails, try manual contrast vector
} else {
  stop("Could not identify T1/T2 columns for CR subjects")
}

# Fallback: build contrast vector manually
cm <- matrix(0, nrow = ncol(design), ncol = 1)
rownames(cm) <- colnames(design)
colnames(cm) <- "Training_CR"
for (col in t2_cols) cm[col, 1] <-  1 / length(t2_cols)
for (col in t1_cols) cm[col, 1] <- -1 / length(t1_cols)

# -- Step 4: Define Cancer_vs_Healthy gene sets --------------------------------

imp_ids <- rownames(mat_imp)

define_sets <- function(dep, ids, use_pi) {
  sig <- if (use_pi) dep %>% filter(pi_score_Cancer_vs_Healthy < 0.05)
         else dep %>% filter(P.Value_Cancer_vs_Healthy < 0.05)
  sig <- sig %>% filter(uniprot_id %in% ids)
  list(
    up   = match(sig$uniprot_id[sig$logFC_Cancer_vs_Healthy > 0], ids),
    down = match(sig$uniprot_id[sig$logFC_Cancer_vs_Healthy < 0], ids),
    up_ids   = sig$uniprot_id[sig$logFC_Cancer_vs_Healthy > 0],
    down_ids = sig$uniprot_id[sig$logFC_Cancer_vs_Healthy < 0]
  )
}

sets_pi <- define_sets(dep_df, imp_ids, TRUE)
sets_p  <- define_sets(dep_df, imp_ids, FALSE)

message(sprintf("Gene sets: Pi up=%d dn=%d | P up=%d dn=%d",
                length(sets_pi$up), length(sets_pi$down),
                length(sets_p$up),  length(sets_p$down)))

# -- Step 5: Pi subset of P verification --------------------------------------

pi_up <- dep_df$uniprot_id[!is.na(dep_df$pi_score_Cancer_vs_Healthy) &
                             dep_df$pi_score_Cancer_vs_Healthy < 0.05 &
                             dep_df$logFC_Cancer_vs_Healthy > 0]
p_up  <- dep_df$uniprot_id[!is.na(dep_df$P.Value_Cancer_vs_Healthy) &
                             dep_df$P.Value_Cancer_vs_Healthy < 0.05 &
                             dep_df$logFC_Cancer_vs_Healthy > 0]
pi_dn <- dep_df$uniprot_id[!is.na(dep_df$pi_score_Cancer_vs_Healthy) &
                             dep_df$pi_score_Cancer_vs_Healthy < 0.05 &
                             dep_df$logFC_Cancer_vs_Healthy < 0]
p_dn  <- dep_df$uniprot_id[!is.na(dep_df$P.Value_Cancer_vs_Healthy) &
                             dep_df$P.Value_Cancer_vs_Healthy < 0.05 &
                             dep_df$logFC_Cancer_vs_Healthy < 0]

overlap_df <- tibble(
  direction = c("up", "down"),
  n_pi = c(length(pi_up), length(pi_dn)),
  n_p  = c(length(p_up),  length(p_dn)),
  n_overlap = c(length(intersect(pi_up, p_up)), length(intersect(pi_dn, p_dn))),
  pi_subset_of_p = c(all(pi_up %in% p_up), all(pi_dn %in% p_dn))
)
write_csv(overlap_df, file.path(DAT, "panel_F_fry", "pi_p_overlap.csv"))

# -- Step 6: No circularity note ----------------------------------------------
# Cancer_vs_Healthy = CR_T1 - H_T1 (baseline comparison)
# Training_CR       = CR_T2 - CR_T1 (within-CR training effect)
# No shared terms -> no structural correlation bias

# -- Step 7: Run fry -- 2 configs x 2 sets ------------------------------------

run_fry <- function(mat, sets, design, cm, block, cor_val, config) {
  map_dfr(c("up", "down"), function(dir) {
    idx <- sets[[dir]]
    if (length(idx) < 3) return(tibble(config = config, set = paste0("cvh_", dir),
                                        n = length(idx), direction = NA_character_,
                                        PValue = NA_real_, PValue.Mixed = NA_real_))
    res <- fry(mat, index = idx, design = design,
               contrast = cm[, "Training_CR"], block = block, correlation = cor_val)
    tibble(config = config, set = paste0("cvh_", dir),
           n = length(idx), direction = res$Direction[1],
           PValue = res$PValue[1], PValue.Mixed = res$PValue.Mixed[1])
  })
}

fry_all <- bind_rows(
  run_fry(mat_imp, sets_pi, design, cm, block_id, cor_imp, "Imp_Pi"),
  run_fry(mat_imp, sets_p,  design, cm, block_id, cor_imp, "Imp_P")
) %>%
  mutate(
    expected = ifelse(set == "cvh_up", "Up", "Down"),
    consistent = direction == expected,
    cor_within = cor_imp
  )

write_csv(fry_all, file.path(DAT, "panel_F_fry", "fry_results_all.csv"))

# -- Step 8: Driving proteins -------------------------------------------------
# Driving = set members whose t_Training_CR is in the concordant direction

driving_df <- bind_rows(
  dep_df %>%
    filter(uniprot_id %in% sets_pi$up_ids, uniprot_id %in% imp_ids,
           t_Training_CR > 0) %>%
    transmute(gene, uniprot_id, set = "cvh_up",
              t_cancer_vs_healthy = t_Cancer_vs_Healthy, t_training_cr = t_Training_CR,
              logFC_Cancer_vs_Healthy, logFC_Training_CR, pi_score_Cancer_vs_Healthy),
  dep_df %>%
    filter(uniprot_id %in% sets_pi$down_ids, uniprot_id %in% imp_ids,
           t_Training_CR < 0) %>%
    transmute(gene, uniprot_id, set = "cvh_down",
              t_cancer_vs_healthy = t_Cancer_vs_Healthy, t_training_cr = t_Training_CR,
              logFC_Cancer_vs_Healthy, logFC_Training_CR, pi_score_Cancer_vs_Healthy)
)
write_csv(driving_df, file.path(DAT, "panel_F_fry", "driving_proteins.csv"))
message(sprintf("Driving proteins: %d", nrow(driving_df)))

# -- Step 9: Barcode data -----------------------------------------------------

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

# -- Step 10: Barcode visualization --------------------------------------------

fry_up <- fry_all %>% filter(config == "Imp_Pi", set == "cvh_up")
fry_dn <- fry_all %>% filter(config == "Imp_Pi", set == "cvh_down")

txt_s <- scale_text(BASE_STAT, PF_W)
n_all <- nrow(t_rank)

make_barcode <- function(t_df, in_col, es_col, fry_row, title, color) {
  marks <- t_df %>% filter(.data[[in_col]])

  p_es <- ggplot(t_df, aes(x = rank, y = .data[[es_col]])) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = color, alpha = 0.04) +
    geom_area(fill = scales::alpha(color, 0.2), color = NA) +
    geom_line(color = color, linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60",
               linewidth = 0.3) +
    annotate("text", x = n_all * 0.98, y = Inf,
             label = sprintf("fry %s, %s (n = %d)",
                              fry_row$direction, fmt_p(fry_row$PValue),
                              fry_row$n),
             hjust = 1, vjust = 1.3, size = txt_s * 1.15, fontface = "bold",
             color = ifelse(fry_row$consistent, "grey20", "#DC2626")) +
    labs(y = "ES", title = title) +
    scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    FIG_THEME +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(4, 4, 0, 4, "mm"),
          plot.title = element_text(size = 10, face = "bold"))

  p_bc <- ggplot(marks, aes(x = rank, xend = rank, y = 0, yend = 1)) +
    geom_segment(color = color, linewidth = 0.3, alpha = 0.7) +
    scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    FIG_THEME +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey97"),
          plot.margin = margin(0, 4, 0, 4, "mm"))

  list(es = p_es, bc = p_bc)
}

p1 <- make_barcode(t_rank, "in_up", "es_up", fry_up,
                    "CvH-Up DEPs (Pi < 0.05) \u2192 Tr.(CR) ranked t-statistics", "#D6604D")

dn_color <- ifelse(fry_dn$PValue > 0.05, "grey65", "#4393C3")
dn_title <- ifelse(fry_dn$PValue > 0.05,
                   "CvH-Down DEPs \u2192 Tr.(CR) ranked t-statistics (n.s.)",
                   "CvH-Down DEPs \u2192 Tr.(CR) ranked t-statistics")
p2 <- make_barcode(t_rank, "in_down", "es_down", fry_dn,
                    dn_title, dn_color)

p_t <- ggplot(t_rank, aes(x = rank, y = t_Training_CR)) +
  geom_area(fill = scales::alpha("#5DA5DA", 0.25), color = "#5DA5DA", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  labs(x = sprintf("Protein rank by t(Training CR)  [%d imputed proteins]", n_all),
       y = "t-stat") +
  scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
  FIG_THEME +
  theme(plot.margin = margin(2, 4, 4, 4, "mm"))

pF_fry <- p1$es / p1$bc / p_t +
  plot_layout(heights = c(3, 0.5, 1.5)) +
  plot_annotation(
    title = "fry Rotation Test: CvH DEP sets against Tr.(CR) ranked t-statistics",
    subtitle = sprintf("CvH Up set (n = %d) \u2192 Tr.(CR) | fry p = %.3f | missForest-imputed (n = %d) | dupCor = %.3f",
                        length(sets_pi$up), fry_up$PValue, n_imp, cor_imp),
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 9))
  )

ggsave(file.path(RPT, "panel_F_fry_SUPP.pdf"), pF_fry,
       width = PF_W, height = 130, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_F_fry_SUPP.png"), pF_fry,
       width = PF_W, height = 130, units = "mm", dpi = 300)

message("F06 CRvH Panel F (fry) done")
