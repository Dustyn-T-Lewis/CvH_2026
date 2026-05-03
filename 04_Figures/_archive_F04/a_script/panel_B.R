# F04 Panel B: Reversal Concordance Scatter
# logFC_Cancer_vs_Healthy (x) vs logFC_Training_CR (y)
# INCLUDES: prepare_data subsections (melov, contingency, threshold, signed reversal)

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(boot)
})

set.seed(42)

PW <- 200; PH <- 200
RPT <- "04_Figures/F04/b_reports"
DAT <- "04_Figures/F04/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "reversal_tests"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_B"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# ============================================================================
# PREPARE DATA: Load DEP, imputation, and metadata
# ============================================================================

dep_df <- read_csv("03_DEP/c_data/03_combined_results_CRvH.csv", show_col_types = FALSE)
stopifnot(nrow(dep_df) > 2000)

# -- Load imputed matrix & metadata --
imp_data <- read_csv("02_Imputation/c_data/01_imputed.csv", show_col_types = FALSE)
imp_ann_cols  <- c("uniprot_id", "protein", "gene", "description")
imp_samp_cols <- setdiff(names(imp_data), imp_ann_cols)
imp_mat <- as.matrix(imp_data[, imp_samp_cols])
rownames(imp_mat) <- imp_data$uniprot_id

dal_meta <- as.data.frame(
  readRDS("02_Imputation/c_data/01_DAList_imputed.rds")$metadata)

meta <- tibble(
  sample_id = dal_meta$Col_ID,
  subject   = sub("_T[12]$", "", dal_meta$Col_ID),
  group     = dal_meta$Group,
  timepoint = dal_meta$Timepoint,
  group_time = dal_meta$Group_Time
)

# ============================================================================
# Section 1: Melov Reversal Permutation Test
# ============================================================================
# Cancer signature: proteins with nominal P < 0.05 for Cancer_vs_Healthy
cancer_sig <- dep_df %>%
  filter(!is.na(P.Value_Cancer_vs_Healthy) & P.Value_Cancer_vs_Healthy < 0.05) %>%
  pull(uniprot_id)
cancer_sig <- intersect(cancer_sig, rownames(imp_mat))
n_cancer <- length(cancer_sig)
message(sprintf("  Cancer signature: %d proteins (nominal P < 0.05)", n_cancer))

# Reference: Healthy at T1
h_t1_ids  <- meta$sample_id[meta$group == "PPS" & meta$timepoint == "T1"]
# Cancer Recovery at T1 (baseline) and T2 (post-training)
cr_t1_ids <- meta$sample_id[meta$group %in% c("CR_CRE", "CR_PLA") & meta$timepoint == "T1"]
cr_t2_ids <- meta$sample_id[meta$group %in% c("CR_CRE", "CR_PLA") & meta$timepoint == "T2"]

h_t1_mean  <- rowMeans(imp_mat[cancer_sig, intersect(h_t1_ids, colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
cr_t1_mean <- rowMeans(imp_mat[cancer_sig, intersect(cr_t1_ids, colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
cr_t2_mean <- rowMeans(imp_mat[cancer_sig, intersect(cr_t2_ids, colnames(imp_mat)), drop = FALSE], na.rm = TRUE)

d_pre  <- sqrt(sum((cr_t1_mean - h_t1_mean)^2))
d_post <- sqrt(sum((cr_t2_mean - h_t1_mean)^2))
reversal_pct <- (d_pre - d_post) / d_pre * 100

# Permutation: swap T1/T2 labels within CR subjects
set.seed(42)
n_perm <- 10000
perm_deltas <- numeric(n_perm)

cr_meta_all <- meta %>% filter(group %in% c("CR_CRE", "CR_PLA"))
cr_t1_meta  <- cr_meta_all %>% filter(timepoint == "T1")
cr_t2_meta  <- cr_meta_all %>% filter(timepoint == "T2")
cr_subjects <- intersect(cr_t1_meta$subject, cr_t2_meta$subject)

for (i in seq_len(n_perm)) {
  swap <- sample(c(TRUE, FALSE), length(cr_subjects), replace = TRUE)
  perm_t1_ids <- character(0)
  perm_t2_ids <- character(0)

  for (j in seq_along(cr_subjects)) {
    subj <- cr_subjects[j]
    t1_id <- cr_t1_meta$sample_id[cr_t1_meta$subject == subj]
    t2_id <- cr_t2_meta$sample_id[cr_t2_meta$subject == subj]
    if (swap[j]) {
      perm_t1_ids <- c(perm_t1_ids, t2_id)
      perm_t2_ids <- c(perm_t2_ids, t1_id)
    } else {
      perm_t1_ids <- c(perm_t1_ids, t1_id)
      perm_t2_ids <- c(perm_t2_ids, t2_id)
    }
  }

  perm_t1_mean <- rowMeans(imp_mat[cancer_sig, intersect(perm_t1_ids, colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
  perm_t2_mean <- rowMeans(imp_mat[cancer_sig, intersect(perm_t2_ids, colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
  d_pre_perm   <- sqrt(sum((perm_t1_mean - h_t1_mean)^2))
  d_post_perm  <- sqrt(sum((perm_t2_mean - h_t1_mean)^2))
  perm_deltas[i] <- d_pre_perm - d_post_perm
}

observed_delta <- d_pre - d_post
perm_pvalue <- mean(perm_deltas >= observed_delta)

n_exceed    <- as.integer(sum(perm_deltas >= observed_delta))
perm_pval_ci <- binom.test(n_exceed, as.integer(n_perm))$conf.int

# Bootstrap CI for reversal %
set.seed(42)
n_boot_rev <- 2000
boot_rev_pct <- replicate(n_boot_rev, {
  idx <- sample(seq_along(cancer_sig), replace = TRUE)
  boot_d_pre  <- sqrt(sum((cr_t1_mean[idx] - h_t1_mean[idx])^2))
  boot_d_post <- sqrt(sum((cr_t2_mean[idx] - h_t1_mean[idx])^2))
  (boot_d_pre - boot_d_post) / boot_d_pre * 100
})
rev_pct_ci <- quantile(boot_rev_pct, c(0.025, 0.975))

message(sprintf("  Permutation p-value: %.4f [%.4f, %.4f]",
                perm_pvalue, perm_pval_ci[1], perm_pval_ci[2]))
message(sprintf("  Reversal %%: %.1f%% [%.1f%%, %.1f%%]",
                reversal_pct, rev_pct_ci[1], rev_pct_ci[2]))

melov_df <- tibble(
  d_pre = d_pre, d_post = d_post,
  reversal_pct = round(reversal_pct, 2),
  reversal_pct_ci_lower = round(rev_pct_ci[1], 2),
  reversal_pct_ci_upper = round(rev_pct_ci[2], 2),
  observed_delta = observed_delta,
  p_value = perm_pvalue,
  p_value_ci_lower = round(perm_pval_ci[1], 6),
  p_value_ci_upper = round(perm_pval_ci[2], 6),
  n_cancer_proteins = n_cancer,
  n_permutations = n_perm,
  n_boot_reversal_pct = n_boot_rev
)
write_csv(melov_df, file.path(DAT, "reversal_tests", "melov_permutation.csv"))

# Permutation histogram
n_cr_subjects <- length(cr_subjects)
p_melov <- ggplot(tibble(delta = perm_deltas), aes(x = delta)) +
  geom_histogram(bins = 50, fill = "grey70", color = "grey50", linewidth = 0.2) +
  geom_vline(xintercept = observed_delta, color = "#D6604D",
             linewidth = 1, linetype = "solid") +
  annotate("text", x = observed_delta, y = Inf, vjust = 1.5, hjust = -0.1,
           label = sprintf("Observed = %.3f\np = %.4f [%.4f, %.4f]\nReversal = %.1f%% [%.1f, %.1f]",
                           observed_delta, perm_pvalue,
                           perm_pval_ci[1], perm_pval_ci[2],
                           reversal_pct, rev_pct_ci[1], rev_pct_ci[2]),
           size = 2.8, fontface = "bold", color = "#D6604D") +
  labs(title = "Melov Reversal Permutation Test",
       subtitle = sprintf("d(CR_T1, H_T1) - d(CR_T2, H_T1) | %d cancer-signature proteins | %d CR subjects",
                          n_cancer, n_cr_subjects),
       x = expression(d[pre] - d[post]), y = "Count (10,000 permutations)") +
  FIG_THEME +
  theme(plot.title = element_text(size = 10, face = "bold"))

ggsave(file.path(RPT, "melov_reversal_permutation.pdf"), p_melov,
       width = 180, height = 120, units = "mm", device = get_pdf_device())

# ============================================================================
# Section 2: Reversal Contingency (Fisher exact)
# ============================================================================
cancer_proteins <- dep_df %>%
  filter(!is.na(P.Value_Cancer_vs_Healthy) & P.Value_Cancer_vs_Healthy < 0.05 &
         !is.na(logFC_Cancer_vs_Healthy) & !is.na(logFC_Training_CR))

contingency <- cancer_proteins %>%
  mutate(
    cancer_dir   = ifelse(logFC_Cancer_vs_Healthy > 0, "Cancer_Up", "Cancer_Down"),
    training_dir = ifelse(logFC_Training_CR > 0, "Training_Up", "Training_Down"),
    pattern      = case_when(
      abs(logFC_Training_CR) <= 0.2 ~ "Negligible",
      sign(logFC_Cancer_vs_Healthy) != sign(logFC_Training_CR) ~ "Reversed",
      TRUE ~ "Exacerbated"
    )
  )

ct <- table(contingency$cancer_dir, contingency$training_dir)
fisher_res <- fisher.test(ct)
fisher_or_ci <- fisher_res$conf.int

n_rev_total <- sum(contingency$pattern == "Reversed")
rev_binom   <- binom.test(n_rev_total, nrow(contingency))
rev_pct_binom_ci <- rev_binom$conf.int * 100

contingency_summary <- tibble(
  cancer_up_training_down   = sum(contingency$cancer_dir == "Cancer_Up" &
                                   contingency$training_dir == "Training_Down"),
  cancer_up_training_up     = sum(contingency$cancer_dir == "Cancer_Up" &
                                   contingency$training_dir == "Training_Up"),
  cancer_down_training_up   = sum(contingency$cancer_dir == "Cancer_Down" &
                                   contingency$training_dir == "Training_Up"),
  cancer_down_training_down = sum(contingency$cancer_dir == "Cancer_Down" &
                                   contingency$training_dir == "Training_Down"),
  n_reversed    = n_rev_total,
  n_exacerbated = sum(contingency$pattern == "Exacerbated"),
  n_negligible  = sum(contingency$pattern == "Negligible"),
  pct_reversed  = round(mean(contingency$pattern == "Reversed") * 100, 1),
  pct_reversed_ci_lower = round(rev_pct_binom_ci[1], 1),
  pct_reversed_ci_upper = round(rev_pct_binom_ci[2], 1),
  fisher_or     = round(fisher_res$estimate, 3),
  fisher_or_ci_lower = round(fisher_or_ci[1], 3),
  fisher_or_ci_upper = round(fisher_or_ci[2], 3),
  fisher_p      = fisher_res$p.value,
  n_cancer_proteins = nrow(contingency)
)
write_csv(contingency_summary, file.path(DAT, "reversal_tests", "reversal_contingency.csv"))

# ============================================================================
# Section 3: Threshold sensitivity
# ============================================================================
sens_rows <- lapply(c(0.1, 0.2, 0.3), function(thresh) {
  patt <- case_when(
    abs(cancer_proteins$logFC_Training_CR) <= thresh ~ "Negligible",
    sign(cancer_proteins$logFC_Cancer_vs_Healthy) != sign(cancer_proteins$logFC_Training_CR) ~ "Reversed",
    TRUE ~ "Exacerbated"
  )
  tibble(threshold = thresh,
         n_reversed = sum(patt == "Reversed"),
         n_exacerbated = sum(patt == "Exacerbated"),
         n_negligible = sum(patt == "Negligible"),
         pct_reversed = round(mean(patt == "Reversed") * 100, 1))
})
write_csv(bind_rows(sens_rows), file.path(DAT, "reversal_tests", "threshold_sensitivity.csv"))

# ============================================================================
# Section 4: Signed reversal score (Pearson r)
# ============================================================================
reversal_cor <- dep_df %>%
  filter(!is.na(logFC_Cancer_vs_Healthy) & !is.na(logFC_Training_CR))

cor_res <- cor.test(reversal_cor$logFC_Cancer_vs_Healthy, reversal_cor$logFC_Training_CR,
                    method = "pearson")

signed_reversal <- tibble(
  r = round(cor_res$estimate, 4),
  ci_lower = round(cor_res$conf.int[1], 4),
  ci_upper = round(cor_res$conf.int[2], 4),
  p_value = cor_res$p.value,
  n_proteins = nrow(reversal_cor),
  interpretation = ifelse(cor_res$estimate < 0,
    "Negative: training opposes cancer signature globally",
    "Positive: training reinforces cancer direction")
)
write_csv(signed_reversal, file.path(DAT, "reversal_tests", "signed_reversal_score.csv"))

message("F04 prepare_data (sections 1-4) complete")

# ============================================================================
# PANEL B: Load imputation classification and prepare scatter data
# ============================================================================

imputation_df <- read_csv("02_Imputation/c_data/02_mar_mnar_classification.csv",
                           show_col_types = FALSE) %>%
  transmute(gene, imputed = classification != "Complete")

# melov_df already in scope from Section 1 above

scatter_df <- dep_df %>%
  transmute(gene,
            logFC_Cancer_vs_Healthy, logFC_Training_CR,
            pi_CvH = pi_score_Cancer_vs_Healthy,
            pi_TR  = pi_score_Training_CR) %>%
  filter(!is.na(logFC_Cancer_vs_Healthy), !is.na(logFC_Training_CR)) %>%
  left_join(imputation_df, by = "gene") %>%
  mutate(
    imputed      = replace_na(imputed, FALSE),
    significance = classify_proteins_f4(pi_CvH, pi_TR),
    quadrant = case_when(
      logFC_Cancer_vs_Healthy > 0 & logFC_Training_CR < 0 ~ "Reversed (Cancer Up / Training Down)",
      logFC_Cancer_vs_Healthy < 0 & logFC_Training_CR > 0 ~ "Reversed (Cancer Down / Training Up)",
      logFC_Cancer_vs_Healthy > 0 & logFC_Training_CR > 0 ~ "Exacerbated Up",
      TRUE ~ "Exacerbated Down"),
    border_col   = ifelse(imputed, "black", "grey75"),
    point_size   = ifelse(significance == "NS", 1.8, 2.3),
    point_stroke = ifelse(significance == "NS", 0.6, 0.9),
    bubble_alpha = case_when(
      significance == "NS"                  ~ 0.30,
      significance == "Sig Both"            ~ 0.75,
      TRUE                                  ~ 0.85)
  )

# -- Correlations --
cor_r   <- cor.test(scatter_df$logFC_Cancer_vs_Healthy, scatter_df$logFC_Training_CR,
                    method = "pearson", conf.level = 0.95)
cor_rho <- cor.test(scatter_df$logFC_Cancer_vs_Healthy, scatter_df$logFC_Training_CR,
                    method = "spearman", conf.level = 0.95)
n_obs   <- nrow(scatter_df)
rho_ci  <- tanh(atanh(cor_rho$estimate) + c(-1, 1) * qnorm(0.975) / sqrt(n_obs - 3))

reversal_pct <- mean(sign(scatter_df$logFC_Cancer_vs_Healthy) !=
                     sign(scatter_df$logFC_Training_CR)) * 100

set.seed(42)
boot_rev <- boot::boot(
  data = scatter_df,
  statistic = function(d, i)
    mean(sign(d$logFC_Cancer_vs_Healthy[i]) != sign(d$logFC_Training_CR[i])) * 100,
  R = 10000)
rev_ci <- tryCatch(
  boot::boot.ci(boot_rev, type = "bca", conf = 0.95)$bca[4:5],
  error = function(e) quantile(boot_rev$t, c(0.025, 0.975)))

# Sig-only correlations
sig_mask <- scatter_df$significance != "NS"
n_sig    <- sum(sig_mask)
cor_r_sig   <- cor.test(scatter_df$logFC_Cancer_vs_Healthy[sig_mask],
                        scatter_df$logFC_Training_CR[sig_mask],
                        method = "pearson", conf.level = 0.95)
cor_rho_sig <- cor.test(scatter_df$logFC_Cancer_vs_Healthy[sig_mask],
                        scatter_df$logFC_Training_CR[sig_mask],
                        method = "spearman", conf.level = 0.95)
set.seed(43)
boot_rho_sig <- boot::boot(
  data = scatter_df[sig_mask, ],
  statistic = function(d, i)
    cor(d$logFC_Cancer_vs_Healthy[i], d$logFC_Training_CR[i], method = "spearman"),
  R = 10000
)
rho_sig_ci <- tryCatch(
  boot::boot.ci(boot_rho_sig, type = "bca", conf = 0.95)$bca[4:5],
  error = function(e) quantile(boot_rho_sig$t, c(0.025, 0.975))
)
rev_sig <- mean(sign(scatter_df$logFC_Cancer_vs_Healthy[sig_mask]) !=
                sign(scatter_df$logFC_Training_CR[sig_mask])) * 100

message(sprintf("  Pearson r = %.3f [%.3f, %.3f], p = %.2g",
                cor_r$estimate, cor_r$conf.int[1], cor_r$conf.int[2], cor_r$p.value))
message(sprintf("  Spearman rho = %.3f [%.3f, %.3f]",
                cor_rho$estimate, rho_ci[1], rho_ci[2]))
message(sprintf("  Reversal %% = %.1f%% [%.1f, %.1f]",
                reversal_pct, rev_ci[1], rev_ci[2]))
message(sprintf("  Sig-only (n=%d): r = %.3f, rho = %.3f, reversal = %.1f%%",
                n_sig, cor_r_sig$estimate, cor_rho_sig$estimate, rev_sig))

txt_gene <- scale_text(BASE_GENE, PW)
txt_quad <- scale_text(BASE_QUADRANT, PW)
txt_stat <- scale_text(BASE_STAT, PW)

# -- Label top proteins --
label_df <- scatter_df %>%
  filter(significance != "NS") %>%
  group_by(significance) %>%
  arrange(desc(abs(logFC_Cancer_vs_Healthy) + abs(logFC_Training_CR))) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(label_fill     = SIG_LABEL_FILL_F4[as.character(significance)],
         label_text_col = SIG_LABEL_TEXT_F4[as.character(significance)])

# -- Quadrant counts --
q_df <- scatter_df %>%
  mutate(q = case_when(
    logFC_Cancer_vs_Healthy > 0 & logFC_Training_CR < 0 ~ "BR",
    logFC_Cancer_vs_Healthy < 0 & logFC_Training_CR > 0 ~ "TL",
    logFC_Cancer_vs_Healthy > 0 & logFC_Training_CR > 0 ~ "TR",
    TRUE ~ "BL"))
q_counts <- q_df %>% count(q) %>% deframe()
q_sig    <- q_df %>% filter(significance != "NS") %>% count(q) %>% deframe()
for (qq in c("BR", "TL", "TR", "BL")) if (is.na(q_sig[qq])) q_sig[qq] <- 0

ns_df  <- scatter_df %>% filter(significance == "NS")
sig_df <- scatter_df %>% filter(significance != "NS")

x_range <- range(scatter_df$logFC_Cancer_vs_Healthy, na.rm = TRUE)
y_range <- range(scatter_df$logFC_Training_CR, na.rm = TRUE)
x_pad   <- diff(x_range) * 0.10
y_pad   <- diff(y_range) * 0.10
xlim_range <- c(x_range[1] - x_pad, x_range[2] + x_pad)
ylim_range <- c(y_range[1] - y_pad, y_range[2] + y_pad * 1.5)

melov_rev_pct <- melov_df$reversal_pct
melov_p       <- melov_df$p_value
melov_n       <- melov_df$n_cancer_proteins

sub_txt <- sprintf(
  "All (n = %s): r = %.2f [%.2f, %.2f], \u03c1 = %.2f [%.2f, %.2f] | reversal = %.0f%%\nSig. (n = %d): r = %.2f [%.2f, %.2f], \u03c1 = %.2f [%.2f, %.2f] | reversal = %.0f%%",
  format(n_obs, big.mark = ","),
  cor_r$estimate, cor_r$conf.int[1], cor_r$conf.int[2],
  cor_rho$estimate, rho_ci[1], rho_ci[2], reversal_pct,
  n_sig,
  cor_r_sig$estimate, cor_r_sig$conf.int[1], cor_r_sig$conf.int[2],
  cor_rho_sig$estimate, rho_sig_ci[1], rho_sig_ci[2], rev_sig)

melov_txt <- sprintf("Melov magnitude reversal: %.1f%%, p = %.2f (%d cancer-sig. proteins)",
                     melov_rev_pct, melov_p, melov_n)

pB <- ggplot(mapping = aes(x = logFC_Cancer_vs_Healthy, y = logFC_Training_CR)) +
  # Reversal quadrants (blue shading)
  annotate("rect", xmin = 0, xmax = Inf,  ymin = -Inf, ymax = 0,
           fill = "#DCEEFF", alpha = 0.55) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#DCEEFF", alpha = 0.55) +
  # Exacerbation quadrants (red shading)
  annotate("rect", xmin = 0, xmax = Inf,  ymin = 0, ymax = Inf,
           fill = "#FFE0E0", alpha = 0.55) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#FFE0E0", alpha = 0.55) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  # Anti-diagonal: perfect reversal line
  geom_abline(slope = -1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  geom_point(data = ns_df, color = "grey70", size = 0.8, alpha = 0.25) +
  geom_point(data = sig_df, aes(fill = significance), shape = 21,
             size = sig_df$point_size, color = sig_df$border_col,
             alpha = sig_df$bubble_alpha, stroke = sig_df$point_stroke) +
  scale_fill_manual(values = SIG_COLORS_F4, name = "Significance") +
  geom_label_repel(data = label_df, aes(label = gene),
                   fill = label_df$label_fill, color = label_df$label_text_col,
                   size = txt_gene, fontface = "italic", max.overlaps = 40,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.6, force = 3, force_pull = 0.5,
                   label.padding = unit(1.5, "pt"), label.r = unit(1, "pt"),
                   label.size = 0.15, seed = 42) +
  # Quadrant labels
  annotate("label", x = Inf, y = -Inf,
           label = sprintf("Reversed\u2002n = %s/%s", q_sig["BR"], q_counts["BR"]),
           hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = Inf,
           label = sprintf("Reversed\u2002n = %s/%s", q_sig["TL"], q_counts["TL"]),
           hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("Exacerbated\u2002n = %s/%s", q_sig["TR"], q_counts["TR"]),
           hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = -Inf,
           label = sprintf("Exacerbated\u2002n = %s/%s", q_sig["BL"], q_counts["BL"]),
           hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("text", x = -Inf, y = -Inf, label = melov_txt,
           hjust = -0.02, vjust = -1.5,
           size = txt_stat, color = "grey30", fontface = "bold") +
  coord_fixed(ratio = 1, xlim = xlim_range, ylim = ylim_range, expand = FALSE) +
  labs(title = "Protein-Level Reversal  |  limma + dupCor, missForest-imputed",
       subtitle = sub_txt,
       x = expression(log[2]*FC ~ "(Cancer vs Healthy)"),
       y = expression(log[2]*FC ~ "(Training CR)")) +
  FIG_THEME +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 8, face = "bold"),
    legend.text     = element_text(size = 7),
    legend.key.size = unit(3, "mm"),
    legend.margin   = margin(0, 0, 0, 0)
  ) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 3, alpha = 0.8)))

ggsave(file.path(RPT, "panel_B_reversal_scatter.pdf"), pB,
       width = PW, height = PH, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_B_reversal_scatter.png"), pB,
       width = PW, height = PH, units = "mm", dpi = 300)

# -- Data export --
scatter_df %>%
  transmute(gene,
            logFC_Cancer_vs_Healthy = round(logFC_Cancer_vs_Healthy, 4),
            logFC_Training_CR       = round(logFC_Training_CR, 4),
            pi_score_Cancer_vs_Healthy = round(pi_CvH, 6),
            pi_score_Training_CR       = round(pi_TR, 6),
            significance               = as.character(significance),
            quadrant, imputed) %>%
  arrange(significance, desc(abs(logFC_Cancer_vs_Healthy) + abs(logFC_Training_CR))) %>%
  write_csv(file.path(DAT, "panel_B", "reversal_scatter.csv"))

tibble(
  metric   = c("Pearson_r", "Spearman_rho", "Reversal_pct",
               "Pearson_r_sig", "Spearman_rho_sig", "Reversal_pct_sig",
               "Melov_magnitude_reversal_pct"),
  estimate = c(cor_r$estimate, cor_rho$estimate, reversal_pct,
               cor_r_sig$estimate, cor_rho_sig$estimate, rev_sig,
               melov_rev_pct),
  ci_lower = c(cor_r$conf.int[1], rho_ci[1], rev_ci[1],
               cor_r_sig$conf.int[1], rho_sig_ci[1], NA_real_,
               melov_df$reversal_pct_ci_lower),
  ci_upper = c(cor_r$conf.int[2], rho_ci[2], rev_ci[2],
               cor_r_sig$conf.int[2], rho_sig_ci[2], NA_real_,
               melov_df$reversal_pct_ci_upper),
  p_value  = c(cor_r$p.value, cor_rho$p.value, NA_real_,
               cor_r_sig$p.value, cor_rho_sig$p.value, NA_real_,
               melov_p),
  n        = c(n_obs, n_obs, n_obs, n_sig, n_sig, n_sig, melov_n),
  note     = c("95% CI from cor.test()",
               "95% CI via Fisher z-transformation",
               "95% BCa bootstrap CI (10000 replicates, all proteins)",
               "Sig proteins only -- 95% CI from cor.test()",
               "Sig proteins only -- 95% BCa bootstrap CI (10000 replicates)",
               "Sig proteins only -- no CI",
               sprintf("Melov permutation test (%d perms)", melov_df$n_permutations))
) %>%
  write_csv(file.path(DAT, "panel_B", "reversal_scatter_stats.csv"))

cat("F04 Panel B done\n")
