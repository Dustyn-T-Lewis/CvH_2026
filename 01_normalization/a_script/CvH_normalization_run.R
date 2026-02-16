#!/usr/bin/env Rscript
# CvH_normalization_run.R — Normalization pipeline (proteoDA)
# INPUT:  00_input/CvH_raw.xlsx, CvH_meta.csv, HPA_skeletal_muscle_annotations.tsv
# OUTPUT: c_data/01_normalized.csv, 01_DAList_normalized.rds, 00_DAList_pre_norm.rds
#         c_data/03_filtering_effects.csv, 03_filtered_proteins.csv
#         c_data/04_outlier_diagnostics.csv, 05_norm_quality_scores.csv
#         b_reports/01_normalization_report.pdf
#         b_reports/02_qc_report_pre_norm.pdf, 03_qc_report_post_norm.pdf

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(proteoDA, readxl, readr, dplyr, tidyr, stringr, ggplot2)

# --- Paths ---
# Run from project root
input_dir  <- "00_input"
report_dir <- "01_normalization/b_reports"
data_dir   <- "01_normalization/c_data"
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir,   recursive = TRUE, showWarnings = FALSE)

# --- 1. Load & validate ---
cat("Step 1: Load & validate\n")

raw      <- read_excel(file.path(input_dir, "CvH_raw.xlsx"))
metadata <- read.csv(file.path(input_dir, "CvH_meta.csv"), stringsAsFactors = FALSE)
rownames(metadata) <- metadata$Col_ID

annot_cols <- c("uniprot_id", "protein", "gene", "description", "n_seq")
annotation <- raw[, annot_cols]
intensity  <- raw[, setdiff(names(raw), annot_cols)]

stopifnot("Sample mismatch" = setequal(colnames(intensity), metadata$Col_ID))
intensity <- intensity[, metadata$Col_ID]

cat(sprintf("  %d proteins x %d samples\n", nrow(raw), ncol(intensity)))

stopifnot("Intensity columns must equal metadata rows" =
            ncol(intensity) == nrow(metadata))

# --- 2. Annotate: HPA skeletal muscle ---
cat("Step 2: HPA annotation\n")

hpa_file <- file.path(input_dir, "HPA_skeletal_muscle_annotations.tsv")
stopifnot("HPA annotations file not found" = file.exists(hpa_file))
hpa_genes <- unique(read_tsv(hpa_file, show_col_types = FALSE)$Gene)

annotation <- annotation %>%
  mutate(in_hpa_muscle = gene %in% hpa_genes)

# --- 3. DAList + filter cascade ---
cat("Step 3: DAList + filter cascade\n")

if (any(duplicated(annotation$uniprot_id))) {
  annotation$rm <- rowMeans(as.data.frame(lapply(intensity, as.numeric)), na.rm = TRUE)
  keep <- annotation %>% mutate(i = row_number()) %>% group_by(uniprot_id) %>%
    slice_max(rm, n = 1, with_ties = FALSE) %>% pull(i)
  annotation <- annotation[keep, ]; intensity <- intensity[keep, ]
  annotation$rm <- NULL
  cat(sprintf("  Deduplicated: %d unique proteins\n", nrow(annotation)))
}

int_mat <- as.data.frame(lapply(intensity, as.numeric))
rownames(int_mat) <- annotation$uniprot_id
annot_df <- as.data.frame(annotation); rownames(annot_df) <- annotation$uniprot_id
dal <- DAList(data = int_mat, annotation = annot_df, metadata = metadata)

counts <- c(Raw = nrow(dal$data))

dal <- zero_to_missing(dal)

dal <- filter_proteins_by_annotation(dal, in_hpa_muscle)
counts["HPA"] <- nrow(dal$data)
cat(sprintf("  HPA filter: %d -> %d\n", counts["Raw"], counts["HPA"]))

dal <- filter_proteins_by_proportion(dal, min_prop = 0.50, grouping_column = "Group_Time")
counts["Miss50"] <- nrow(dal$data)
cat(sprintf("  50%% filter: %d -> %d\n", counts["HPA"], counts["Miss50"]))

stopifnot("Filter cascade must only remove proteins" = all(diff(counts) <= 0))

n_start <- counts["Raw"]
filter_log <- bind_rows(
  tibble(step = "Raw input", n_before = NA_integer_,
         n_after = as.integer(counts["Raw"]), n_removed = NA_integer_),
  tibble(step = "HPA muscle tissue filter",
         n_before = as.integer(counts["Raw"]),
         n_after  = as.integer(counts["HPA"]),
         n_removed = as.integer(counts["Raw"] - counts["HPA"])),
  tibble(step = "Missingness filter (50% per Group_Time)",
         n_before = as.integer(counts["HPA"]),
         n_after  = as.integer(counts["Miss50"]),
         n_removed = as.integer(counts["HPA"] - counts["Miss50"]))
) %>%
  mutate(pct_retained = round(n_after / n_start * 100, 1))
write_csv(filter_log, file.path(data_dir, "03_filtering_effects.csv"))

removed <- annot_df %>%
  filter(!uniprot_id %in% rownames(dal$data)) %>%
  dplyr::select(uniprot_id, gene, description, in_hpa_muscle) %>%
  mutate(reason = if_else(!in_hpa_muscle, "Not in HPA muscle", "Missingness"))
write_csv(removed, file.path(data_dir, "03_filtered_proteins.csv"))

# --- 4. Outlier detection ---
cat("Step 4: Outlier detection\n")

pm <- colMeans(is.na(dal$data)) * 100
f_miss <- pm > (quantile(pm, 0.75) + 1.5 * IQR(pm))

# Compute paired abs_delta for CR subjects (like QMD does)
paired_subjects <- dal$metadata %>%
  filter(Group != "PPS") %>%
  count(Subject_ID) %>%
  filter(n == 2) %>%
  pull(Subject_ID)

abs_delta_map <- setNames(rep(NA_real_, ncol(dal$data)), colnames(dal$data))
for (subj in paired_subjects) {
  rows <- dal$metadata %>% filter(Subject_ID == subj)
  t1_id <- rows$Col_ID[rows$Timepoint == "T1"]
  t2_id <- rows$Col_ID[rows$Timepoint == "T2"]
  if (length(t1_id) == 1 && length(t2_id) == 1) {
    d <- abs(pm[t2_id] - pm[t1_id])
    abs_delta_map[t1_id] <- d
    abs_delta_map[t2_id] <- d
  }
}

tmp <- log2(dal$data)
for (j in 1:ncol(tmp)) tmp[is.na(tmp[, j]), j] <- median(tmp[, j], na.rm = TRUE)
pc <- prcomp(t(tmp), center = TRUE, scale. = TRUE)
md <- mahalanobis(pc$x[, 1:3], colMeans(pc$x[, 1:3]), cov(pc$x[, 1:3]))
f_pca <- md > qchisq(0.99, df = 3)

sm <- apply(log2(dal$data), 2, median, na.rm = TRUE)
mad_dev <- abs(sm - median(sm))
f_mad <- mad_dev > 3 * mad(sm)

nf <- f_miss + f_pca + f_mad
outlier_diag <- tibble(
  Col_ID = colnames(dal$data), pct_missing = pm, abs_delta = abs_delta_map[colnames(dal$data)],
  miss_flag = f_miss,
  mahal_dist = md, pca_flag = f_pca, sample_median = sm, mad_deviation = mad_dev,
  mad_flag = f_mad,
  n_flags = nf,
  tier = case_when(
    nf >= 3 & pct_missing > 50 ~ "catastrophic",
    nf >= 2                    ~ "soft_outlier",
    TRUE                       ~ "normal"
  ))
write_csv(outlier_diag, file.path(data_dir, "04_outlier_diagnostics.csv"))

hard_remove  <- outlier_diag$Col_ID[outlier_diag$tier == "catastrophic"]
soft_outlier <- outlier_diag$Col_ID[outlier_diag$tier == "soft_outlier"]

cat(sprintf("  Flagged: %d catastrophic, %d soft outlier\n",
            length(hard_remove), length(soft_outlier)))
if (length(hard_remove) + length(soft_outlier) > 0) {
  print(filter(outlier_diag, tier != "normal") %>%
        dplyr::select(Col_ID, pct_missing, abs_delta, mahal_dist, n_flags, tier))
}
if (length(hard_remove) > 0) {
  dal <- filter_samples(dal, !(Col_ID %in% hard_remove))
  cat(sprintf("  Removed: %s\n", paste(hard_remove, collapse = ", ")))
}

cat(sprintf("  After outlier removal: %d proteins x %d samples\n",
            nrow(dal$data), ncol(dal$data)))

# --- 5. Normalize ---
cat("Step 5: Cycloess normalization\n")

saveRDS(dal, file.path(data_dir, "00_DAList_pre_norm.rds"))

write_norm_report(dal, grouping_column = "Group_Time", output_dir = report_dir,
                  filename = "01_normalization_report.pdf", overwrite = TRUE)

write_qc_report(dal, color_column = "Group_Time", label_column = "Col_ID",
                output_dir = report_dir, filename = "02_qc_report_pre_norm.pdf", overwrite = TRUE)

dal <- normalize_data(dal, norm_method = "cycloess")

write_qc_report(dal, color_column = "Group_Time", label_column = "Col_ID",
                output_dir = report_dir, filename = "03_qc_report_post_norm.pdf", overwrite = TRUE)

cat(sprintf("  Final: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

# --- 5b. Norm quality scoring ---
cat("Step 5b: Norm quality scoring\n")
compute_pcv <- function(mat, groups) {
  cvs <- c()
  for (grp in unique(groups)) {
    sub <- mat[, groups == grp, drop = FALSE]
    cv <- apply(sub, 1, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2 || mean(x) == 0) return(NA_real_)
      sd(x) / abs(mean(x))
    })
    cvs <- c(cvs, cv)
  }
  median(cvs, na.rm = TRUE)
}

compute_pmad <- function(mat, groups) {
  mads <- c()
  for (grp in unique(groups)) {
    sub <- mat[, groups == grp, drop = FALSE]
    m <- apply(sub, 1, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) return(NA_real_)
      mad(x, constant = 1)
    })
    mads <- c(mads, m)
  }
  median(mads, na.rm = TRUE)
}

compute_cor <- function(mat, groups) {
  cors <- c()
  for (grp in unique(groups)) {
    idx <- which(groups == grp)
    if (length(idx) < 2) next
    cm <- cor(mat[, idx], use = "pairwise.complete.obs")
    cors <- c(cors, cm[lower.tri(cm)])
  }
  mean(cors, na.rm = TRUE)
}

dal_pre <- readRDS(file.path(data_dir, "00_DAList_pre_norm.rds"))
dal_pre$metadata$group <- factor(dal_pre$metadata$Group_Time,
                                  levels = c("CRE_T1","CRE_T2","PLA_T1","PLA_T2","H_T1"))

NORMS <- c("log2", "median", "mean", "vsn", "quantile", "cycloess", "rlr", "gi")

norm_scores <- tibble()
for (nm in NORMS) {
  dal_n <- tryCatch(normalize_data(dal_pre, norm_method = nm), error = function(e) NULL)
  if (is.null(dal_n)) { cat(sprintf("  %s: FAILED\n", nm)); next }

  mat_n <- as.matrix(dal_n$data)
  grps  <- dal_n$metadata$group

  pcv     <- compute_pcv(mat_n, grps)
  pmad    <- compute_pmad(mat_n, grps)
  cor_val <- compute_cor(mat_n, grps)

  norm_scores <- bind_rows(norm_scores, tibble(
    norm = nm, PCV = round(pcv, 4), PMAD = round(pmad, 4), COR = round(cor_val, 4)
  ))
  cat(sprintf("  %s: PCV=%.4f  PMAD=%.4f  COR=%.4f\n", nm, pcv, pmad, cor_val))
}

norm_scores <- norm_scores %>%
  mutate(PCV_rank  = rank(PCV),
         PMAD_rank = rank(PMAD),
         COR_rank  = rank(-COR),
         norm_composite = (PCV_rank + PMAD_rank + COR_rank) / 3) %>%
  arrange(norm_composite)

cat("\n  Normalization rankings:\n")
print(norm_scores %>% dplyr::select(norm, PCV, PMAD, COR, norm_composite))

write_csv(norm_scores, file.path(data_dir, "05_norm_quality_scores.csv"))

# --- 6. Export ---
cat("Step 6: Export\n")

export_df <- bind_cols(
  dal$annotation %>% dplyr::select(uniprot_id, protein, gene, description, n_seq),
  as_tibble(dal$data)
)

stopifnot("Export row count must match DAList" = nrow(export_df) == nrow(dal$data))
stopifnot("Export columns must match DAList"   = ncol(export_df) - 5 == ncol(dal$data))
stopifnot("No duplicate uniprot_ids in export" = !any(duplicated(export_df$uniprot_id)))

write_csv(export_df, file.path(data_dir, "01_normalized.csv"))
saveRDS(dal, file.path(data_dir, "01_DAList_normalized.rds"))

cat(sprintf("  Exported: %d proteins x %d samples\n", nrow(export_df), ncol(export_df) - 5))
