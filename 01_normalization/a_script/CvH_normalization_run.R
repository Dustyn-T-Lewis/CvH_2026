#!/usr/bin/env Rscript
# ============================================================================
# CvH Proteomics — Normalization Pipeline (Standalone Runner)
# Cancer Survivor (CR) vs. Healthy (PPS) Skeletal Muscle DIA-MS Proteomics
#
# 1:1 reproduction of CvH_normalization.qmd as a clean R script.
# Run from: 01_normalization/a_script/
# Usage:    source("CvH_normalization_run.R")  OR  Rscript CvH_normalization_run.R
# ============================================================================

cat("=== CvH Normalization Pipeline ===\n\n")

# --- 0: Setup ----------------------------------------------------------------
cat(">> 0 — Setup\n")

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  proteoDA,
  readxl, readr, dplyr, tidyr, stringr,
  httr, xml2,
  ggplot2, patchwork, knitr
)

# Paths
base_dir   <- normalizePath(file.path(dirname(getwd()), ".."), mustWork = TRUE)
input_dir  <- file.path(base_dir, "00_input")
report_dir <- file.path(base_dir, "01_normalization", "b_reports")
data_dir   <- file.path(base_dir, "01_normalization", "c_data")

dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir,   recursive = TRUE, showWarnings = FALSE)

# Color palette
pal_group <- c(CR_CRE = "#2166AC", CR_PLA = "#67A9CF", PPS = "#B2182B")
pal_group_time <- c(
  CR_T1 = "#92C5DE", CR_T2 = "#2166AC",
  H_T1  = "#B2182B"
)
pal_supplement <- c(CRE = "#2166AC", PLA = "#67A9CF", `NA` = "#999999")
shape_tp <- c(T1 = 16, T2 = 17)

cat("   Base directory:", base_dir, "\n")

# --- 1: Load & Validate Data -------------------------------------------------
cat("\n>> 1 — Load & Validate Data\n")

# Raw intensity data
raw <- read_excel(file.path(input_dir, "CvH_raw.xlsx"))

# Separate annotation from intensity
annot_cols <- c("uniprot_id", "protein", "gene", "description", "n_seq")
annotation <- raw[, annot_cols]
intensity  <- raw[, setdiff(names(raw), annot_cols)]

cat(sprintf("   Raw data: %d proteins x %d samples\n", nrow(raw), ncol(intensity)))

# Metadata
metadata <- read.csv(file.path(input_dir, "CvH_meta.csv"), stringsAsFactors = FALSE)
rownames(metadata) <- metadata$Col_ID

cat(sprintf("   Metadata: %d samples\n", nrow(metadata)))

# Validate alignment
data_samples <- colnames(intensity)
meta_samples <- metadata$Col_ID

stopifnot(
  "Sample mismatch between data and metadata" =
    setequal(data_samples, meta_samples)
)

# Reorder intensity columns to match metadata row order
intensity <- intensity[, meta_samples]

cat("   Validation passed: all", length(meta_samples), "samples aligned.\n")

# Sample summary
cat("\n   Sample distribution:\n")
print(metadata %>% count(Group, Timepoint, Group_Time, Supplement))

# --- 2: CRAPome Contaminant Annotation ---------------------------------------
cat("\n>> 2 — CRAPome Contaminant Annotation\n")

muscle_protect <- c(
  "MYH1", "MYH2", "MYH7", "MYH4", "MYL1", "MYL2", "MYL3",
  "MYBPC1", "MYBPC2",
  "ACTA1", "ACTC1", "ACTN2", "ACTN3", "ACTB", "ACTG1",
  "TPM1", "TPM2", "TPM3", "TNNC1", "TNNC2", "TNNI1", "TNNI2",
  "TNNT1", "TNNT3",
  "TTN", "NEB",
  "CKM", "CKB", "CKMT2", "MB", "PYGM", "ALDOA", "ENO3",
  "GAPDH", "PKM", "LDHA", "LDHB", "ENO1", "PFKM", "GPI",
  "ATP5F1A", "ATP5F1B", "COX4I1", "COX5A", "UQCRC1", "UQCRC2",
  "SDHA", "SDHB", "CS", "MDH2", "IDH2", "OGDH",
  "DES", "VIM", "FLNC", "LDB3", "MYOZ1", "MYOZ2", "SYNPO2L",
  "COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "COL6A3",
  "HSPA1A", "HSPA8", "HSPB1", "HSPB6", "HSPD1", "HSPE1",
  "HSP90AA1", "HSP90AB1",
  "HBA1", "HBA2", "HBB",
  "ATP2A1", "ATP2A2", "CASQ1", "CALM1", "CALM2",
  "TPI1", "PGAM2", "PGK1", "MDH1"
)

# CRAPome query (with caching)
crapome_cache <- file.path(data_dir, "00_crapome_annotations.csv")

if (file.exists(crapome_cache)) {
  cat("   Loading cached CRAPome annotations...\n")
  crapome <- read_csv(crapome_cache, show_col_types = FALSE)
} else {
  cat("   Querying CRAPome API for", length(unique(annotation$gene)), "genes...\n")
  cat("   This will take ~10-15 minutes.\n\n")

  genes <- unique(annotation$gene)
  crapome <- tibble(gene = genes, crapome_n_expt = NA_integer_)

  for (i in seq_along(genes)) {
    g <- genes[i]
    if (i %% 100 == 0) cat(sprintf("     %d / %d ...\n", i, length(genes)))

    tryCatch({
      url <- sprintf(
        "https://reprint-apms.org/?q=ws/proteindetail/%s/human/singleStep",
        URLencode(g, reserved = TRUE)
      )
      resp <- GET(url, timeout(15))

      if (status_code(resp) == 200) {
        xml_text <- content(resp, as = "text", encoding = "UTF-8")
        if (grepl("^<\\?xml", xml_text)) {
          doc <- read_xml(xml_text)
          expts <- xml_attr(xml_find_all(doc, ".//protein"), "expt")
          crapome$crapome_n_expt[i] <- length(unique(expts))
        } else {
          crapome$crapome_n_expt[i] <- 0L
        }
      } else {
        crapome$crapome_n_expt[i] <- 0L
      }
    }, error = function(e) {
      crapome$crapome_n_expt[i] <<- 0L
    })

    Sys.sleep(0.2)
  }

  max_expt <- max(crapome$crapome_n_expt, na.rm = TRUE)
  crapome$crapome_pct <- round(crapome$crapome_n_expt / max_expt * 100, 1)

  write_csv(crapome, crapome_cache)
  cat(sprintf("   CRAPome annotation complete. Max experiments: %d\n", max_expt))
}

# Ensure pct is computed even from cache
if (!"crapome_pct" %in% names(crapome)) {
  max_expt <- max(crapome$crapome_n_expt, na.rm = TRUE)
  crapome$crapome_pct <- round(crapome$crapome_n_expt / max_expt * 100, 1)
}

# Flag contaminants
crapome <- crapome %>%
  mutate(
    is_muscle_protected = gene %in% muscle_protect,
    is_contaminant = (crapome_pct > 20) & !is_muscle_protected
  )

# Join to annotation
annotation <- annotation %>%
  left_join(crapome, by = "gene") %>%
  mutate(
    crapome_n_expt = replace_na(crapome_n_expt, 0L),
    crapome_pct    = replace_na(crapome_pct, 0),
    is_contaminant = replace_na(is_contaminant, FALSE)
  )

n_flagged   <- sum(annotation$is_contaminant)
n_protected <- sum(annotation$is_muscle_protected & annotation$crapome_pct > 20,
                   na.rm = TRUE)

cat(sprintf("   Proteins flagged as contaminants: %d\n", n_flagged))
cat(sprintf("   Muscle proteins protected (CRAPome >20%% but retained): %d\n", n_protected))

# Print contaminant tables
cat("\n   Flagged contaminants:\n")
print(annotation %>%
  filter(is_contaminant) %>%
  select(gene, uniprot_id, description, crapome_n_expt, crapome_pct) %>%
  arrange(desc(crapome_pct)))

cat("\n   Muscle-protected proteins retained:\n")
print(annotation %>%
  filter(is_muscle_protected, crapome_pct > 20) %>%
  select(gene, uniprot_id, description, crapome_n_expt, crapome_pct) %>%
  arrange(desc(crapome_pct)))

# --- 3: Assemble DAList -------------------------------------------------------
cat("\n>> 3 — Assemble DAList\n")

# Check for duplicate uniprot_ids
dup_ids <- annotation$uniprot_id[duplicated(annotation$uniprot_id)]

if (length(dup_ids) > 0) {
  cat(sprintf("   Found %d duplicate uniprot_ids — deduplicating by highest mean intensity.\n",
              length(dup_ids)))

  intensity_num <- as.data.frame(lapply(intensity, as.numeric))
  annotation$row_mean <- rowMeans(intensity_num, na.rm = TRUE)

  keep_idx <- annotation %>%
    mutate(row_idx = row_number()) %>%
    group_by(uniprot_id) %>%
    slice_max(row_mean, n = 1, with_ties = FALSE) %>%
    pull(row_idx)

  annotation <- annotation[keep_idx, ]
  intensity  <- intensity[keep_idx, ]
  annotation$row_mean <- NULL

  cat(sprintf("   After deduplication: %d unique proteins\n", nrow(annotation)))
} else {
  cat("   No duplicate uniprot_ids found.\n")
}

# Convert intensity to numeric matrix
intensity_mat <- as.data.frame(lapply(intensity, as.numeric))
rownames(intensity_mat) <- annotation$uniprot_id

# Annotation df
annot_df <- as.data.frame(annotation)
rownames(annot_df) <- annotation$uniprot_id

# Metadata df
meta_df <- as.data.frame(metadata)
rownames(meta_df) <- metadata$Col_ID

# Assemble
dal <- DAList(
  data       = intensity_mat,
  annotation = annot_df,
  metadata   = meta_df
)

cat(sprintf("   DAList assembled: %d proteins x %d samples\n",
            nrow(dal$data), ncol(dal$data)))

# --- 4: Quality Filtering ----------------------------------------------------
cat("\n>> 4 — Quality Filtering\n")

filter_log <- tibble(
  step      = character(),
  n_before  = integer(),
  n_after   = integer(),
  n_removed = integer()
)

n_start <- nrow(dal$data)

# Convert 0s to NA
dal <- zero_to_missing(dal)
cat(sprintf("   Converted zeros to NA. Total proteins: %d\n", nrow(dal$data)))

# CRAPome contaminant removal
n_before <- nrow(dal$data)
dal <- filter_proteins_by_annotation(dal, !is_contaminant)
n_after <- nrow(dal$data)

filter_log <- bind_rows(filter_log, tibble(
  step = "CRAPome contaminant removal",
  n_before = n_before, n_after = n_after, n_removed = n_before - n_after
))

cat(sprintf("   CRAPome filtering: %d -> %d proteins (%d removed)\n",
            n_before, n_after, n_before - n_after))

# Missingness filter
n_before <- nrow(dal$data)
dal <- filter_proteins_by_proportion(
  dal,
  min_prop = 0.66,
  grouping_column = "Group_Time"
)
n_after <- nrow(dal$data)

filter_log <- bind_rows(filter_log, tibble(
  step = "Missingness filter (66% per Group_Time)",
  n_before = n_before, n_after = n_after, n_removed = n_before - n_after
))

cat(sprintf("   Missingness filtering: %d -> %d proteins (%d removed)\n",
            n_before, n_after, n_before - n_after))

# Filtering outputs
filter_log <- bind_rows(
  tibble(step = "Raw input", n_before = NA_integer_,
         n_after = n_start, n_removed = NA_integer_),
  filter_log
) %>%
  mutate(pct_retained = round(n_after / n_start * 100, 1))

write_csv(filter_log, file.path(data_dir, "03_filtering_effects.csv"))

all_uniprots <- annot_df$uniprot_id
kept_uniprots <- rownames(dal$data)
removed_uniprots <- setdiff(all_uniprots, kept_uniprots)

filtered_proteins <- annot_df %>%
  filter(uniprot_id %in% removed_uniprots) %>%
  select(uniprot_id, gene, description, crapome_pct, is_contaminant) %>%
  mutate(reason = if_else(is_contaminant, "CRAPome contaminant", "Missingness"))

write_csv(filtered_proteins, file.path(data_dir, "03_filtered_proteins.csv"))

cat("\n   Filtering summary:\n")
print(filter_log)

# Filtering barplot
filter_plot_data <- filter_log %>%
  mutate(step = factor(step, levels = step))

p_filter <- ggplot(filter_plot_data, aes(x = step, y = n_after)) +
  geom_col(fill = "#2166AC", width = 0.6) +
  geom_text(aes(label = n_after), vjust = -0.3, size = 4) +
  labs(x = NULL, y = "Number of proteins",
       title = "Protein retention through filtering pipeline") +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

ggsave(file.path(report_dir, "03_filtering_effects.pdf"),
       p_filter, width = 8, height = 5)
cat("   Saved: 03_filtering_effects.pdf\n")

# --- 5: Outlier Detection & Removal ------------------------------------------
cat("\n>> 5 — Outlier Detection & Removal\n")

# Per-sample percent missing
pct_missing <- colMeans(is.na(dal$data)) * 100

# Paired missingness (delta) for CR subjects with both T1 and T2
paired_subjects <- dal$metadata %>%
  filter(Group != "PPS") %>%
  count(Subject_ID) %>%
  filter(n == 2) %>%
  pull(Subject_ID)

delta_missing <- tibble(
  Col_ID = character(),
  Subject_ID = character(),
  pct_missing = numeric(),
  delta_missing = numeric()
)

for (subj in paired_subjects) {
  rows <- dal$metadata %>% filter(Subject_ID == subj)
  t1_id <- rows$Col_ID[rows$Timepoint == "T1"]
  t2_id <- rows$Col_ID[rows$Timepoint == "T2"]
  if (length(t1_id) == 1 && length(t2_id) == 1) {
    m1 <- pct_missing[t1_id]
    m2 <- pct_missing[t2_id]
    delta_missing <- bind_rows(delta_missing, tibble(
      Col_ID = c(t1_id, t2_id),
      Subject_ID = subj,
      pct_missing = c(m1, m2),
      delta_missing = m2 - m1
    ))
  }
}

# For unpaired samples (PPS), add with delta = NA
unpaired_ids <- dal$metadata %>%
  filter(Group == "PPS") %>%
  pull(Col_ID)

delta_missing <- bind_rows(delta_missing, tibble(
  Col_ID = unpaired_ids,
  Subject_ID = dal$metadata$Subject_ID[match(unpaired_ids, dal$metadata$Col_ID)],
  pct_missing = pct_missing[unpaired_ids],
  delta_missing = NA_real_
))

# Also add unpaired CR subjects
unpaired_cr <- setdiff(dal$metadata$Col_ID, delta_missing$Col_ID)
if (length(unpaired_cr) > 0) {
  delta_missing <- bind_rows(delta_missing, tibble(
    Col_ID = unpaired_cr,
    Subject_ID = dal$metadata$Subject_ID[match(unpaired_cr, dal$metadata$Col_ID)],
    pct_missing = pct_missing[unpaired_cr],
    delta_missing = NA_real_
  ))
}

# IQR-based flags
miss_q3  <- quantile(pct_missing, 0.75)
miss_iqr <- IQR(pct_missing)
miss_threshold <- miss_q3 + 1.5 * miss_iqr

delta_vals <- delta_missing$delta_missing[!is.na(delta_missing$delta_missing)]
if (length(delta_vals) > 2) {
  delta_q3  <- quantile(delta_vals, 0.75)
  delta_iqr <- IQR(delta_vals)
  delta_threshold <- delta_q3 + 1.5 * delta_iqr
  delta_lower     <- quantile(delta_vals, 0.25) - 1.5 * delta_iqr
} else {
  delta_threshold <- Inf
  delta_lower <- -Inf
}

delta_missing <- delta_missing %>%
  mutate(
    miss_flag = pct_missing > miss_threshold |
      (!is.na(delta_missing) & (delta_missing > delta_threshold | delta_missing < delta_lower))
  )

# PCA outlier detection
data_for_pca <- dal$data
for (j in seq_len(ncol(data_for_pca))) {
  nas <- is.na(data_for_pca[, j])
  if (any(nas)) {
    data_for_pca[nas, j] <- median(data_for_pca[, j], na.rm = TRUE)
  }
}
data_log2 <- log2(data_for_pca + 1)

pca_res <- prcomp(t(data_log2), center = TRUE, scale. = TRUE)
pc_scores <- pca_res$x[, 1:3]

center <- colMeans(pc_scores)
cov_mat <- cov(pc_scores)
mahal_dist <- mahalanobis(pc_scores, center, cov_mat)

mahal_threshold <- qchisq(0.99, df = 3)

pca_flags <- tibble(
  Col_ID = colnames(dal$data),
  mahal_dist = mahal_dist,
  pca_flag = mahal_dist > mahal_threshold
)

# MAD-based intensity outlier detection
sample_medians <- apply(log2(dal$data + 1), 2, median, na.rm = TRUE)
global_median  <- median(sample_medians)
mad_val        <- mad(sample_medians)

mad_flags <- tibble(
  Col_ID = names(sample_medians),
  sample_median = sample_medians,
  mad_deviation = abs(sample_medians - global_median),
  mad_flag = abs(sample_medians - global_median) > 3 * mad_val
)

# Consensus
outlier_diag <- delta_missing %>%
  select(Col_ID, pct_missing, delta_missing, miss_flag) %>%
  left_join(pca_flags, by = "Col_ID") %>%
  left_join(mad_flags, by = "Col_ID") %>%
  mutate(
    n_flags = miss_flag + pca_flag + mad_flag,
    consensus_outlier = n_flags >= 2
  )

write_csv(outlier_diag, file.path(data_dir, "04_outlier_diagnostics.csv"))

n_outliers <- sum(outlier_diag$consensus_outlier)
cat(sprintf("   Outlier consensus: %d sample(s) flagged by >= 2 methods\n", n_outliers))

if (n_outliers > 0) {
  cat("\n   Consensus outliers:\n")
  print(outlier_diag %>%
    filter(consensus_outlier) %>%
    select(Col_ID, pct_missing, delta_missing, mahal_dist, miss_flag, pca_flag, mad_flag))
}

# Outlier diagnostic plots
p1 <- ggplot(outlier_diag %>% filter(!is.na(delta_missing)),
             aes(x = pct_missing, y = delta_missing)) +
  geom_point(aes(color = consensus_outlier), size = 3) +
  geom_hline(yintercept = c(delta_lower, delta_threshold),
             linetype = "dashed", color = "red", alpha = 0.5) +
  geom_vline(xintercept = miss_threshold,
             linetype = "dashed", color = "red", alpha = 0.5) +
  scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
  labs(x = "% Missing", y = "Delta Missing (T2 - T1)",
       title = "A: Paired Missingness") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

pc_df <- as.data.frame(pca_res$x[, 1:2])
pc_df$Col_ID <- rownames(pc_df)
pc_df <- left_join(pc_df, outlier_diag %>% select(Col_ID, consensus_outlier), by = "Col_ID")
var_explained <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

p2 <- ggplot(pc_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = consensus_outlier), size = 3) +
  scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
  labs(x = sprintf("PC1 (%.1f%%)", var_explained[1]),
       y = sprintf("PC2 (%.1f%%)", var_explained[2]),
       title = "B: PCA Outliers (Mahalanobis)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

p3 <- ggplot(outlier_diag, aes(x = reorder(Col_ID, sample_median), y = sample_median)) +
  geom_point(aes(color = consensus_outlier), size = 2.5) +
  geom_hline(yintercept = global_median, color = "black") +
  geom_hline(yintercept = global_median + c(-3, 3) * mad_val,
             linetype = "dashed", color = "red", alpha = 0.5) +
  scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
  labs(x = "Sample", y = "Median log2 intensity",
       title = "C: MAD Intensity Outliers") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        legend.position = "none")

p_outlier <- p1 + p2 + p3 + plot_layout(ncol = 1)

ggsave(file.path(report_dir, "04_outlier_diagnostics.pdf"),
       p_outlier, width = 10, height = 12)
cat("   Saved: 04_outlier_diagnostics.pdf\n")

# Remove outliers
if (n_outliers > 0) {
  outlier_ids <- outlier_diag %>%
    filter(consensus_outlier) %>%
    pull(Col_ID)

  cat(sprintf("   Removing %d outlier sample(s): %s\n",
              length(outlier_ids), paste(outlier_ids, collapse = ", ")))

  dal <- filter_samples(dal, !(Col_ID %in% outlier_ids))

  cat(sprintf("   After outlier removal: %d samples remain\n", ncol(dal$data)))
} else {
  cat("   No outlier samples removed.\n")
}

# --- 6: Normalization ---------------------------------------------------------
cat("\n>> 6 — Normalization\n")

write_norm_report(
  dal,
  grouping_column = "Group_Time",
  output_dir = report_dir,
  filename = "01_normalization_report.pdf",
  overwrite = TRUE
)
cat("   Saved: 01_normalization_report.pdf\n")

dal <- normalize_data(dal, norm_method = "cycloess")

cat(sprintf("   Normalization complete (cycloess): %d proteins x %d samples\n",
            nrow(dal$data), ncol(dal$data)))

# --- 7: Post-Normalization QC -------------------------------------------------
cat("\n>> 7 — Post-Normalization QC\n")

write_qc_report(
  dal,
  color_column = "Group_Time",
  label_column = "Col_ID",
  output_dir = report_dir,
  filename = "02_qc_report.pdf",
  overwrite = TRUE
)
cat("   Saved: 02_qc_report.pdf\n")

# Custom PCA
norm_mat <- dal$data
for (j in seq_len(ncol(norm_mat))) {
  nas <- is.na(norm_mat[, j])
  if (any(nas)) {
    norm_mat[nas, j] <- median(norm_mat[, j], na.rm = TRUE)
  }
}

pca_norm <- prcomp(t(norm_mat), center = TRUE, scale. = TRUE)
pc_df <- as.data.frame(pca_norm$x[, 1:3])
pc_df$Col_ID <- rownames(pc_df)
pc_df <- left_join(pc_df, dal$metadata, by = "Col_ID")

var_exp <- round(summary(pca_norm)$importance[2, 1:3] * 100, 1)

p_pca1 <- ggplot(pc_df, aes(x = PC1, y = PC2, color = Group, shape = Timepoint)) +
  geom_point(size = 3.5, alpha = 0.85) +
  stat_ellipse(aes(group = Group), type = "norm", level = 0.68,
               linetype = "solid", linewidth = 0.7) +
  scale_color_manual(values = pal_group) +
  scale_shape_manual(values = shape_tp) +
  labs(x = sprintf("PC1 (%.1f%%)", var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp[2]),
       title = "A: By Group") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

p_pca2 <- ggplot(pc_df, aes(x = PC1, y = PC2, color = Group_Time, shape = Timepoint)) +
  geom_point(size = 3.5, alpha = 0.85) +
  stat_ellipse(aes(group = Group_Time), type = "norm", level = 0.68,
               linetype = "solid", linewidth = 0.7) +
  scale_color_manual(values = pal_group_time) +
  scale_shape_manual(values = shape_tp) +
  labs(x = sprintf("PC1 (%.1f%%)", var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp[2]),
       title = "B: By Group x Timepoint") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

pc_df$Supplement_label <- replace_na(pc_df$Supplement, "NA")

p_pca3 <- ggplot(pc_df, aes(x = PC1, y = PC2,
                              color = Supplement_label, shape = Timepoint)) +
  geom_point(size = 3.5, alpha = 0.85) +
  stat_ellipse(aes(group = Supplement_label), type = "norm", level = 0.68,
               linetype = "solid", linewidth = 0.7) +
  scale_color_manual(values = pal_supplement) +
  scale_shape_manual(values = shape_tp) +
  labs(x = sprintf("PC1 (%.1f%%)", var_exp[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp[2]),
       title = "C: By Supplement",
       color = "Supplement") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

p_pca_combined <- p_pca1 / p_pca2 / p_pca3

ggsave(file.path(report_dir, "05_post_norm_pca.pdf"),
       p_pca_combined, width = 9, height = 16)
cat("   Saved: 05_post_norm_pca.pdf\n")

# --- 8: Export ----------------------------------------------------------------
cat("\n>> 8 — Export\n")

export_df <- bind_cols(
  as_tibble(dal$annotation) %>%
    select(uniprot_id, protein, gene, description),
  as_tibble(dal$data)
)

write_csv(export_df, file.path(data_dir, "01_normalized.csv"))

saveRDS(dal, file.path(data_dir, "01_DAList_normalized.rds"))

cat(sprintf("   Exported: %d proteins x %d samples\n",
            nrow(dal$data), ncol(dal$data)))
cat(sprintf("   CSV: %s\n", file.path(data_dir, "01_normalized.csv")))
cat(sprintf("   RDS: %s\n", file.path(data_dir, "01_DAList_normalized.rds")))

# --- 9: Session Info ----------------------------------------------------------
cat("\n>> 9 — Session Info\n")
sessionInfo()

cat("\n=== CvH Normalization Pipeline Complete ===\n")
