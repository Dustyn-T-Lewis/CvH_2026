# F05/CR Melov-style concordance tests -- Cosine Similarity + Permutation
# Cosine similarity of CRE vs PLA training response vectors.
# Standalone exploratory script -- run from project root.

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

set.seed(42)

# -- Load data -----------------------------------------------------------------
dep_df <- read_csv("03_DEP/c_data/03_combined_results_CR.csv", show_col_types = FALSE)

imp_data <- read_csv("02_Imputation/c_data/01_imputed.csv", show_col_types = FALSE)
imp_ann_cols <- c("uniprot_id", "protein", "gene", "description")
imp_samp_cols <- setdiff(names(imp_data), imp_ann_cols)
imp_mat <- as.matrix(imp_data[, imp_samp_cols])
rownames(imp_mat) <- imp_data$uniprot_id

dal_meta <- as.data.frame(
  readRDS("02_Imputation/c_data/01_DAList_imputed.rds")$metadata)
meta <- tibble(
  sample_id = dal_meta$Col_ID,
  subject   = sub("_T[12]$", "", dal_meta$Col_ID),
  group     = dal_meta$Group,
  time      = dal_meta$Timepoint,
  group_time = dal_meta$Group_Time
)

# -- Sample IDs ----------------------------------------------------------------
cre_t1_ids <- meta$sample_id[meta$group == "CR_CRE" & meta$time == "T1"]
cre_t2_ids <- meta$sample_id[meta$group == "CR_CRE" & meta$time == "T2"]
pla_t1_ids <- meta$sample_id[meta$group == "CR_PLA" & meta$time == "T1"]
pla_t2_ids <- meta$sample_id[meta$group == "CR_PLA" & meta$time == "T2"]

cat(sprintf("Subjects: %d CRE, %d PLA\n",
            length(unique(meta$subject[meta$group == "CR_CRE"])),
            length(unique(meta$subject[meta$group == "CR_PLA"]))))

# -- Define training-responsive signature (CRE P < 0.05) ----------------------
training_sig <- dep_df %>%
  filter(!is.na(P.Value_Training_CRE) & P.Value_Training_CRE < 0.05) %>%
  pull(uniprot_id)
training_sig <- intersect(training_sig, rownames(imp_mat))
n_sig <- length(training_sig)
cat(sprintf("Training signature: %d proteins (Training_CRE nominal P < 0.05)\n\n", n_sig))

# -- Compute group centroids on training signature -----------------------------
cre_t1_mean <- rowMeans(imp_mat[training_sig, intersect(cre_t1_ids, colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
cre_t2_mean <- rowMeans(imp_mat[training_sig, intersect(cre_t2_ids, colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
pla_t1_mean <- rowMeans(imp_mat[training_sig, intersect(pla_t1_ids, colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
pla_t2_mean <- rowMeans(imp_mat[training_sig, intersect(pla_t2_ids, colnames(imp_mat)), drop = FALSE], na.rm = TRUE)

# Training response vectors
cre_delta <- cre_t2_mean - cre_t1_mean
pla_delta <- pla_t2_mean - pla_t1_mean

# ==============================================================================
# OPTION A: Cosine Similarity Permutation Test
# Question: Is the PLA training response directionally aligned with CRE's?
# ==============================================================================
cosine_sim <- function(a, b) sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))

observed_cosine <- cosine_sim(cre_delta, pla_delta)
cat("=== OPTION A: Cosine Similarity Test ===\n")
cat(sprintf("Observed cosine similarity: %.4f\n", observed_cosine))

# Permutation: shuffle T1/T2 within PLA subjects
n_perm <- 10000
set.seed(42)

pla_subjects  <- unique(meta$subject[meta$group == "CR_PLA"])
pla_t1_meta   <- meta %>% filter(group == "CR_PLA", time == "T1")
pla_t2_meta   <- meta %>% filter(group == "CR_PLA", time == "T2")

perm_cosines <- numeric(n_perm)
for (i in seq_len(n_perm)) {
  swap <- sample(c(TRUE, FALSE), length(pla_subjects), replace = TRUE)
  perm_t1_ids <- character(0)
  perm_t2_ids <- character(0)

  for (j in seq_along(pla_subjects)) {
    subj <- pla_subjects[j]
    t1_id <- pla_t1_meta$sample_id[pla_t1_meta$subject == subj]
    t2_id <- pla_t2_meta$sample_id[pla_t2_meta$subject == subj]
    if (swap[j]) {
      perm_t1_ids <- c(perm_t1_ids, t2_id)
      perm_t2_ids <- c(perm_t2_ids, t1_id)
    } else {
      perm_t1_ids <- c(perm_t1_ids, t1_id)
      perm_t2_ids <- c(perm_t2_ids, t2_id)
    }
  }

  perm_t1_mean <- rowMeans(imp_mat[training_sig, intersect(perm_t1_ids, colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
  perm_t2_mean <- rowMeans(imp_mat[training_sig, intersect(perm_t2_ids, colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
  perm_pla_delta <- perm_t2_mean - perm_t1_mean
  perm_cosines[i] <- cosine_sim(cre_delta, perm_pla_delta)
}

cosine_p <- mean(perm_cosines >= observed_cosine)
cosine_p_ci <- binom.test(sum(perm_cosines >= observed_cosine), n_perm)$conf.int

cat(sprintf("Permutation p-value: %.4f [%.4f, %.4f]\n", cosine_p, cosine_p_ci[1], cosine_p_ci[2]))
cat(sprintf("Null cosine mean: %.4f, SD: %.4f\n", mean(perm_cosines), sd(perm_cosines)))
cat(sprintf("Observed is %.1f SDs above null mean\n\n",
            (observed_cosine - mean(perm_cosines)) / sd(perm_cosines)))

# ==============================================================================
# OPTION B: Distance-to-Reference Test
# Question: Does training in PLA move the proteome toward the CRE_T1 state?
# (Melov-style adaptation for supplement concordance)
# ==============================================================================
cat("=== OPTION B: Distance-to-Reference Test ===\n")

d_pla_t1 <- sqrt(sum((pla_t1_mean - cre_t1_mean)^2))
d_pla_t2 <- sqrt(sum((pla_t2_mean - cre_t1_mean)^2))
observed_delta_B <- d_pla_t1 - d_pla_t2
reversal_pct_B <- (d_pla_t1 - d_pla_t2) / d_pla_t1 * 100

cat(sprintf("d(PLA_T1, CRE_T1) = %.4f\n", d_pla_t1))
cat(sprintf("d(PLA_T2, CRE_T1) = %.4f\n", d_pla_t2))
cat(sprintf("Observed delta (d_t1 - d_t2) = %.4f\n", observed_delta_B))
cat(sprintf("Magnitude reversal: %.1f%%\n", reversal_pct_B))

set.seed(42)
perm_deltas_B <- numeric(n_perm)

for (i in seq_len(n_perm)) {
  swap <- sample(c(TRUE, FALSE), length(pla_subjects), replace = TRUE)
  perm_t1_ids <- character(0)
  perm_t2_ids <- character(0)

  for (j in seq_along(pla_subjects)) {
    subj <- pla_subjects[j]
    t1_id <- pla_t1_meta$sample_id[pla_t1_meta$subject == subj]
    t2_id <- pla_t2_meta$sample_id[pla_t2_meta$subject == subj]
    if (swap[j]) {
      perm_t1_ids <- c(perm_t1_ids, t2_id)
      perm_t2_ids <- c(perm_t2_ids, t1_id)
    } else {
      perm_t1_ids <- c(perm_t1_ids, t1_id)
      perm_t2_ids <- c(perm_t2_ids, t2_id)
    }
  }

  perm_t1_mean <- rowMeans(imp_mat[training_sig, intersect(perm_t1_ids, colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
  perm_t2_mean <- rowMeans(imp_mat[training_sig, intersect(perm_t2_ids, colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
  d_t1_perm <- sqrt(sum((perm_t1_mean - cre_t1_mean)^2))
  d_t2_perm <- sqrt(sum((perm_t2_mean - cre_t1_mean)^2))
  perm_deltas_B[i] <- d_t1_perm - d_t2_perm
}

melov_p_B <- mean(perm_deltas_B >= observed_delta_B)
melov_p_ci_B <- binom.test(sum(perm_deltas_B >= observed_delta_B), n_perm)$conf.int

# Bootstrap CI on reversal %
set.seed(42)
boot_rev_B <- replicate(2000, {
  idx <- sample(seq_along(training_sig), replace = TRUE)
  b_d_t1 <- sqrt(sum((pla_t1_mean[idx] - cre_t1_mean[idx])^2))
  b_d_t2 <- sqrt(sum((pla_t2_mean[idx] - cre_t1_mean[idx])^2))
  (b_d_t1 - b_d_t2) / b_d_t1 * 100
})
rev_ci_B <- quantile(boot_rev_B, c(0.025, 0.975))

cat(sprintf("Permutation p-value: %.4f [%.4f, %.4f]\n", melov_p_B, melov_p_ci_B[1], melov_p_ci_B[2]))
cat(sprintf("Reversal %% = %.1f%% [%.1f, %.1f]\n\n", reversal_pct_B, rev_ci_B[1], rev_ci_B[2]))

# ==============================================================================
# BONUS: Cosine on Baseline_Supplement signature
# ==============================================================================
cat("=== BONUS: Cosine on Baseline_Supplement signature ===\n")

bl_sig <- dep_df %>%
  filter(!is.na(P.Value_Baseline_Supplement) & P.Value_Baseline_Supplement < 0.05) %>%
  pull(uniprot_id)
bl_sig <- intersect(bl_sig, rownames(imp_mat))
cat(sprintf("Baseline_Supplement signature: %d proteins\n", length(bl_sig)))

if (length(bl_sig) > 5) {
  cre_t1_bl <- rowMeans(imp_mat[bl_sig, intersect(cre_t1_ids, colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
  cre_t2_bl <- rowMeans(imp_mat[bl_sig, intersect(cre_t2_ids, colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
  pla_t1_bl <- rowMeans(imp_mat[bl_sig, intersect(pla_t1_ids, colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
  pla_t2_bl <- rowMeans(imp_mat[bl_sig, intersect(pla_t2_ids, colnames(imp_mat)), drop = FALSE], na.rm = TRUE)

  cre_training_bl <- cre_t2_bl - cre_t1_bl
  pla_training_bl <- pla_t2_bl - pla_t1_bl
  cos_bl <- cosine_sim(cre_training_bl, pla_training_bl)
  cat(sprintf("Cosine(Training_CRE, Training_PLA) on BL sig: %.4f\n", cos_bl))
  cat("  (positive = concordant training response on BL-affected proteins)\n\n")
} else {
  cat("  Too few BL signature proteins -- skipping\n\n")
}

# ==============================================================================
# Summary table
# ==============================================================================
cat("=== SUMMARY ===\n")
summary_df <- tibble(
  test = c("Option A: Cosine concordance (training sig)",
           "Option B: Distance-to-CRE_T1 (training sig)"),
  signature = c(sprintf("%d proteins (CRE P<0.05)", n_sig),
                sprintf("%d proteins (CRE P<0.05)", n_sig)),
  statistic = c(sprintf("cos = %.4f", observed_cosine),
                sprintf("delta = %.4f (%.1f%%)", observed_delta_B, reversal_pct_B)),
  p_value = c(sprintf("%.4f", cosine_p),
              sprintf("%.4f", melov_p_B)),
  interpretation = c(
    ifelse(cosine_p < 0.05, "SIGNIFICANT -- PLA training aligns with CRE", "NS -- directional alignment not distinguishable from noise"),
    ifelse(melov_p_B < 0.05, "SIGNIFICANT -- PLA moves toward CRE_T1", "NS -- distance reduction not distinguishable from noise")
  )
)

for (r in seq_len(nrow(summary_df))) {
  cat(sprintf("\n%s\n  Signature: %s\n  Statistic: %s\n  p = %s\n  -> %s\n",
              summary_df$test[r], summary_df$signature[r],
              summary_df$statistic[r], summary_df$p_value[r],
              summary_df$interpretation[r]))
}

# -- Save results --------------------------------------------------------------
DAT <- "04_Figures/F05/CR/c_data"
dir.create(file.path(DAT, "melov"), recursive = TRUE, showWarnings = FALSE)

write_csv(summary_df, file.path(DAT, "melov", "concordance_test_summary.csv"))

cat("\nDone.\n")
