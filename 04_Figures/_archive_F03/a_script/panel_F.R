# F03 Panel F — fry Rotation Test: Supplement Concordance Barcode
# Tests whether Training_CRE-significant proteins respond concordantly
# in PLA subjects, using limma's fry rotation framework.
#
# No circularity: Training_CRE (CRE_T2 - CRE_T1) and Training_PLA
# (PLA_T2 - PLA_T1) share no subjects.
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(limma)
  library(patchwork)
})

set.seed(42)

RPT <- "04_Figures/F03/b_reports"
DAT <- "04_Figures/F03/c_data"
dir.create(file.path(DAT, "panel_F_fry"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()
PF_W <- 220

# -- Step 1: Load data --------------------------------------------------------
dal      <- readRDS("02_Imputation/c_data/01_DAList_imputed.rds")
dep_df   <- read_csv("03_DEP/c_data/03_combined_results_CR.csv",
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

# Subset to CR subjects only (CRE + PLA, not H)
cr_meta <- meta %>% filter(Group %in% c("CR_CRE", "CR_PLA"))
cr_samples <- cr_meta$Col_ID
mat_cr <- mat_imp[, cr_samples]

n_imp <- nrow(mat_cr)
message(sprintf("CR matrix: %d proteins x %d samples", n_imp, ncol(mat_cr)))

# -- Step 3: Design matrix + duplicateCorrelation -----------------------------
cr_meta$Group_Time <- factor(cr_meta$Group_Time,
  levels = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2"))
design <- model.matrix(~ 0 + Group_Time, data = cr_meta)
colnames(design) <- gsub("^Group_Time", "", colnames(design))

block_id <- sub("_T[12]$", "", cr_meta$Col_ID)

corfit_imp <- duplicateCorrelation(mat_cr, design, block = block_id)
cor_imp <- corfit_imp$consensus.correlation
message(sprintf("Within-subject cor: %.4f", cor_imp))

cm <- makeContrasts(
  Training_PLA = PLA_T2 - PLA_T1,
  levels = design
)

# -- Step 4: Define Training_CRE gene sets -----------------------------------
imp_ids <- rownames(mat_cr)

define_sets <- function(dep, ids, use_pi) {
  sig <- if (use_pi) dep %>% filter(pi_score_Training_CRE < 0.05)
         else dep %>% filter(P.Value_Training_CRE < 0.05)
  sig <- sig %>% filter(uniprot_id %in% ids)
  list(
    up   = match(sig$uniprot_id[sig$logFC_Training_CRE > 0], ids),
    down = match(sig$uniprot_id[sig$logFC_Training_CRE < 0], ids),
    up_ids   = sig$uniprot_id[sig$logFC_Training_CRE > 0],
    down_ids = sig$uniprot_id[sig$logFC_Training_CRE < 0]
  )
}

sets_pi <- define_sets(dep_df, imp_ids, TRUE)
sets_p  <- define_sets(dep_df, imp_ids, FALSE)

message(sprintf("Gene sets: Pi up=%d dn=%d | P up=%d dn=%d",
                length(sets_pi$up), length(sets_pi$down),
                length(sets_p$up),  length(sets_p$down)))

# -- Step 5: Run fry -----------------------------------------------------------
run_fry <- function(mat, sets, design, cm, block, cor_val, config) {
  map_dfr(c("up", "down"), function(dir) {
    idx <- sets[[dir]]
    if (length(idx) < 3) return(tibble(config = config, set = paste0("tcre_", dir),
                                        n = length(idx), direction = NA_character_,
                                        PValue = NA_real_, PValue.Mixed = NA_real_))
    res <- fry(mat, index = idx, design = design,
               contrast = cm[, "Training_PLA"], block = block, correlation = cor_val)
    tibble(config = config, set = paste0("tcre_", dir),
           n = length(idx), direction = res$Direction[1],
           PValue = res$PValue[1], PValue.Mixed = res$PValue.Mixed[1])
  })
}

fry_all <- bind_rows(
  run_fry(mat_cr, sets_pi, design, cm, block_id, cor_imp, "Imp_Pi"),
  run_fry(mat_cr, sets_p,  design, cm, block_id, cor_imp, "Imp_P")
) %>%
  mutate(
    expected = ifelse(set == "tcre_up", "Up", "Down"),
    consistent = direction == expected,
    cor_within = cor_imp
  )

write_csv(fry_all, file.path(DAT, "panel_F_fry", "fry_results_all.csv"))

# -- Step 6: Driving proteins -------------------------------------------------
driving_df <- bind_rows(
  dep_df %>%
    filter(uniprot_id %in% sets_pi$up_ids, uniprot_id %in% imp_ids,
           t_Training_PLA > 0) %>%
    transmute(gene, uniprot_id, set = "tcre_up",
              t_training_cre = t_Training_CRE, t_training_pla = t_Training_PLA,
              logFC_Training_CRE, logFC_Training_PLA, pi_score_Training_CRE),
  dep_df %>%
    filter(uniprot_id %in% sets_pi$down_ids, uniprot_id %in% imp_ids,
           t_Training_PLA < 0) %>%
    transmute(gene, uniprot_id, set = "tcre_down",
              t_training_cre = t_Training_CRE, t_training_pla = t_Training_PLA,
              logFC_Training_CRE, logFC_Training_PLA, pi_score_Training_CRE)
)
write_csv(driving_df, file.path(DAT, "panel_F_fry", "driving_proteins.csv"))
message(sprintf("Driving proteins: %d", nrow(driving_df)))

# -- Step 7: Barcode data -----------------------------------------------------
t_rank <- dep_df %>%
  filter(uniprot_id %in% imp_ids, !is.na(t_Training_PLA)) %>%
  arrange(desc(t_Training_PLA)) %>%
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

t_rank$es_up   <- running_es(t_rank$t_Training_PLA, t_rank$in_up)
t_rank$es_down <- running_es(t_rank$t_Training_PLA, t_rank$in_down)

# -- Step 8: Barcode visualization ---------------------------------------------
fry_up <- fry_all %>% filter(config == "Imp_Pi", set == "tcre_up")
fry_dn <- fry_all %>% filter(config == "Imp_Pi", set == "tcre_down")

txt_s <- scale_text(BASE_STAT, PF_W)
n_all <- nrow(t_rank)

make_barcode <- function(t_df, in_col, es_col, fry_row, title, color) {
  marks <- t_df %>% filter(.data[[in_col]])

  p_es <- ggplot(t_df, aes(x = rank, y = .data[[es_col]])) +
    geom_line(color = color, linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60",
               linewidth = 0.3) +
    annotate("text", x = n_all * 0.98, y = Inf,
             label = sprintf("fry %s, %s (n = %d)",
                              fry_row$direction, fmt_p(fry_row$PValue),
                              fry_row$n),
             hjust = 1, vjust = 1.3, size = txt_s * 0.85, fontface = "bold",
             color = ifelse(fry_row$consistent, "grey20", "#DC2626")) +
    labs(y = "ES", title = title) +
    scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    FIG_THEME +
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(4, 4, 0, 4, "mm"),
          plot.title = element_text(size = 9, face = "bold"))

  p_bc <- ggplot(marks, aes(x = rank, xend = rank, y = 0, yend = 1)) +
    geom_segment(color = color, linewidth = 0.3, alpha = 0.7) +
    scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    FIG_THEME +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey97"),
          plot.margin = margin(0, 4, 0, 4, "mm"))

  p_es / p_bc + plot_layout(heights = c(3, 0.5))
}

p1 <- make_barcode(t_rank, "in_up", "es_up", fry_up,
                    "CRE-Up -> Training PLA  (expect Up)", "#D6604D")

p2 <- make_barcode(t_rank, "in_down", "es_down", fry_dn,
                    "CRE-Down -> Training PLA  (expect Down)", "#4393C3")

p_t <- ggplot(t_rank, aes(x = rank, y = t_Training_PLA)) +
  geom_area(fill = "grey85", color = "grey60", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  labs(x = sprintf("Protein rank by t(Training PLA)  [%d proteins]", n_all),
       y = "t-stat") +
  scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
  FIG_THEME +
  theme(plot.margin = margin(2, 4, 4, 4, "mm"))

pF_fry <- wrap_elements(p1) / wrap_elements(p2) / p_t +
  plot_layout(heights = c(3.5, 3.5, 1.5)) +
  plot_annotation(
    title = "fry Rotation Test: Supplement Concordance (CRE -> PLA)",
    subtitle = sprintf("missForest-imputed (n = %d) | Pi < 0.05 | dupCor = %.3f",
                        n_imp, cor_imp),
    theme = theme(plot.title = element_text(size = 11, face = "bold"),
                  plot.subtitle = element_text(size = 9))
  )

ggsave(file.path(RPT, "panel_F_fry.pdf"), pF_fry,
       width = PF_W, height = 250, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_F_fry.png"), pF_fry,
       width = PF_W, height = 250, units = "mm", dpi = 300)

message("F03 Panel F (fry) done")
