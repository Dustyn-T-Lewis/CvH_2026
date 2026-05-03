# F05 CRvH Melov-style concordance test -- Cosine Similarity Permutation
# Cosine similarity between Cancer_vs_Healthy and Training_CR logFC vectors.
# Permutation test: shuffle T1/T2 labels within CR subjects (10,000 reps).
# Standalone exploratory script -- run from project root.

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

set.seed(42)

# -- Load data --
dep_df <- read_csv("03_DEP/c_data/03_combined_results_CRvH.csv", show_col_types = FALSE)

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
  time      = dal_meta$Timepoint
)

# -- Sample IDs --
# CRvH model: CR subjects have T1+T2, Healthy (PPS) subjects have T1 only
# Groups in metadata: CR_CRE, CR_PLA (both are Cancer Recovery), PPS (Healthy)
cr_subjects <- unique(meta$subject[meta$group %in% c("CR_CRE", "CR_PLA")])
h_subjects  <- unique(meta$subject[meta$group == "PPS"])

cr_t1_ids <- meta$sample_id[meta$subject %in% cr_subjects & meta$time == "T1"]
cr_t2_ids <- meta$sample_id[meta$subject %in% cr_subjects & meta$time == "T2"]
h_t1_ids  <- meta$sample_id[meta$subject %in% h_subjects  & meta$time == "T1"]

cat(sprintf("Subjects: %d CR, %d Healthy\n",
            length(cr_subjects), length(h_subjects)))

# -- Define cancer-different signature (Cancer_vs_Healthy P < 0.05) --
cancer_sig <- dep_df %>%
  filter(!is.na(P.Value_Cancer_vs_Healthy) & P.Value_Cancer_vs_Healthy < 0.05) %>%
  pull(uniprot_id)
cancer_sig <- intersect(cancer_sig, rownames(imp_mat))
n_sig <- length(cancer_sig)
cat(sprintf("Cancer signature: %d proteins (Cancer_vs_Healthy nominal P < 0.05)\n\n", n_sig))

# -- Compute group centroids on cancer signature --
cr_t1_mean <- rowMeans(imp_mat[cancer_sig, intersect(cr_t1_ids, colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
cr_t2_mean <- rowMeans(imp_mat[cancer_sig, intersect(cr_t2_ids, colnames(imp_mat)), drop = FALSE], na.rm = TRUE)
h_t1_mean  <- rowMeans(imp_mat[cancer_sig, intersect(h_t1_ids,  colnames(imp_mat)), drop = FALSE], na.rm = TRUE)

# Cancer difference vector (CRvH baseline)
cancer_delta <- cr_t1_mean - h_t1_mean

# Training response vector (CR training effect)
training_delta <- cr_t2_mean - cr_t1_mean

# ======================================================================
# OPTION A: Cosine Similarity Permutation Test
# Question: Is the CR training response directionally aligned with the
# cancer difference? (Positive = training moves in same direction as cancer)
# ======================================================================
cosine_sim <- function(a, b) sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2)))

observed_cosine <- cosine_sim(cancer_delta, training_delta)
cat("=== OPTION A: Cosine Similarity Test ===\n")
cat(sprintf("Observed cosine similarity: %.4f\n", observed_cosine))

# Permutation: shuffle T1/T2 within CR subjects
n_perm <- 10000
set.seed(42)

cr_t1_meta <- meta %>% filter(subject %in% cr_subjects, time == "T1")
cr_t2_meta <- meta %>% filter(subject %in% cr_subjects, time == "T2")

perm_cosines <- numeric(n_perm)
for (i in seq_len(n_perm)) {
  swap <- sample(c(TRUE, FALSE), length(cr_subjects), replace = TRUE)
  perm_t1_ids <- character(0)
  perm_t2_ids <- character(0)

  for (j in seq_along(cr_subjects)) {
    subj <- cr_subjects[j]
    t1_id <- cr_t1_meta$sample_id[cr_t1_meta$subject == subj]
    t2_id <- cr_t2_meta$sample_id[cr_t2_meta$subject == subj]
    if (length(t1_id) == 0 || length(t2_id) == 0) next
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
  perm_training_delta <- perm_t2_mean - perm_t1_mean
  perm_cosines[i] <- cosine_sim(cancer_delta, perm_training_delta)
}

cosine_p <- mean(perm_cosines >= observed_cosine)
cosine_p_ci <- binom.test(sum(perm_cosines >= observed_cosine), n_perm)$conf.int

cat(sprintf("Permutation p-value: %.4f [%.4f, %.4f]\n", cosine_p, cosine_p_ci[1], cosine_p_ci[2]))
cat(sprintf("Null cosine mean: %.4f, SD: %.4f\n", mean(perm_cosines), sd(perm_cosines)))
cat(sprintf("Observed is %.1f SDs above null mean\n\n",
            (observed_cosine - mean(perm_cosines)) / sd(perm_cosines)))

# ======================================================================
# OPTION B: Distance to Healthy Reference
# Question: Does training in CR move the proteome toward the Healthy state?
# ======================================================================
cat("=== OPTION B: Distance-to-Healthy Test ===\n")

d_cr_t1 <- sqrt(sum((cr_t1_mean - h_t1_mean)^2))
d_cr_t2 <- sqrt(sum((cr_t2_mean - h_t1_mean)^2))
observed_delta_B <- d_cr_t1 - d_cr_t2
reversal_pct_B <- (d_cr_t1 - d_cr_t2) / d_cr_t1 * 100

cat(sprintf("d(CR_T1, H_T1) = %.4f\n", d_cr_t1))
cat(sprintf("d(CR_T2, H_T1) = %.4f\n", d_cr_t2))
cat(sprintf("Observed delta (d_T1 - d_T2) = %.4f\n", observed_delta_B))
cat(sprintf("Magnitude reversal: %.1f%%\n", reversal_pct_B))

set.seed(42)
perm_deltas_B <- numeric(n_perm)

for (i in seq_len(n_perm)) {
  swap <- sample(c(TRUE, FALSE), length(cr_subjects), replace = TRUE)
  perm_t1_ids <- character(0)
  perm_t2_ids <- character(0)

  for (j in seq_along(cr_subjects)) {
    subj <- cr_subjects[j]
    t1_id <- cr_t1_meta$sample_id[cr_t1_meta$subject == subj]
    t2_id <- cr_t2_meta$sample_id[cr_t2_meta$subject == subj]
    if (length(t1_id) == 0 || length(t2_id) == 0) next
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
  d_t1_perm <- sqrt(sum((perm_t1_mean - h_t1_mean)^2))
  d_t2_perm <- sqrt(sum((perm_t2_mean - h_t1_mean)^2))
  perm_deltas_B[i] <- d_t1_perm - d_t2_perm
}

melov_p_B <- mean(perm_deltas_B >= observed_delta_B)
melov_p_ci_B <- binom.test(sum(perm_deltas_B >= observed_delta_B), n_perm)$conf.int

# Bootstrap CI on reversal %
set.seed(42)
boot_rev_B <- replicate(2000, {
  idx <- sample(seq_along(cancer_sig), replace = TRUE)
  b_d_t1 <- sqrt(sum((cr_t1_mean[idx] - h_t1_mean[idx])^2))
  b_d_t2 <- sqrt(sum((cr_t2_mean[idx] - h_t1_mean[idx])^2))
  (b_d_t1 - b_d_t2) / b_d_t1 * 100
})
rev_ci_B <- quantile(boot_rev_B, c(0.025, 0.975))

cat(sprintf("Permutation p-value: %.4f [%.4f, %.4f]\n", melov_p_B, melov_p_ci_B[1], melov_p_ci_B[2]))
cat(sprintf("Reversal %% = %.1f%% [%.1f, %.1f]\n\n", reversal_pct_B, rev_ci_B[1], rev_ci_B[2]))

# ======================================================================
# Summary table
# ======================================================================
cat("=== SUMMARY ===\n")
summary_df <- tibble(
  test = c("Option A: Cosine concordance (cancer sig)",
           "Option B: Distance-to-Healthy (cancer sig)"),
  signature = c(sprintf("%d proteins (CvH P<0.05)", n_sig),
                sprintf("%d proteins (CvH P<0.05)", n_sig)),
  statistic = c(sprintf("cos = %.4f", observed_cosine),
                sprintf("delta = %.4f (%.1f%%)", observed_delta_B, reversal_pct_B)),
  p_value = c(sprintf("%.4f", cosine_p),
              sprintf("%.4f", melov_p_B)),
  interpretation = c(
    ifelse(cosine_p < 0.05, "SIGNIFICANT -- CR training aligns with cancer difference", "NS -- directional alignment not distinguishable from noise"),
    ifelse(melov_p_B < 0.05, "SIGNIFICANT -- Training moves CR toward Healthy", "NS -- distance reduction not distinguishable from noise")
  )
)

for (r in seq_len(nrow(summary_df))) {
  cat(sprintf("\n%s\n  Signature: %s\n  Statistic: %s\n  p = %s\n  -> %s\n",
              summary_df$test[r], summary_df$signature[r],
              summary_df$statistic[r], summary_df$p_value[r],
              summary_df$interpretation[r]))
}

# Export summary
DAT <- "04_Figures/F05/CRvH/c_data"
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)
write_csv(summary_df, file.path(DAT, "melov_concordance_summary.csv"))

cat("\nDone.\n")
