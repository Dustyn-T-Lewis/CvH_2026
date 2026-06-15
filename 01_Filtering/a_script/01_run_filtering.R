#!/usr/bin/env Rscript
# =============================================================================
# 01_run_filtering.R  --  CvH Stage 01: Filtering (no normalization)
#
# Chain: load -> HPA tissue filter -> HPA-derived blood removal -> UniProt dedup
#        -> missingness filter (proteoDA, by group_time) -> 4-method outlier
#        consensus (>=3/4). Outputs a FILTERED, outlier-removed, UN-normalized
#        DAList for Stage 02 (Normalization).
#
# Design: 3-group pooled scheme  group_time in {H_pre, CR_pre, CR_post}
#         (supplement pooled; CRE/PLA retained as a column for DEP sensitivity).
# =============================================================================

suppressPackageStartupMessages({
  library(proteoDA); library(here); library(readxl); library(readr)
  library(dplyr); library(tidyr)
})
set.seed(42)
source(here("R", "cvh_design.R"))
source(here("R", "blood_filter.R"))

cfg <- list(
  raw_file   = here("00_input", "CvH_raw.xlsx"),
  meta_file  = here("00_input", "CvH_meta.csv"),
  pheno_file = here("00_input", "CRm_meta.csv"),
  hpa_file   = here("00_input", "HPA_annotations.tsv"),  # single HPA export (presence + blood)
  data_dir   = here("01_Filtering", "c_data"),
  min_reps = 5L, min_groups = 1L, outlier_k = 3, mad_k = 3, mahal_p = 0.01
)
dir.create(cfg$data_dir, recursive = TRUE, showWarnings = FALSE)

run_pca <- function(mat) {                       # median-impute for PCA only
  for (j in seq_len(ncol(mat))) mat[is.na(mat[, j]), j] <- median(mat[, j], na.rm = TRUE)
  prcomp(t(log2(mat)), center = TRUE, scale. = TRUE)
}

# --- 1. Load -----------------------------------------------------------------
raw <- read_excel(cfg$raw_file)
annot_cols <- c("uniprot_id", "protein", "gene", "description", "n_seq")
annotation <- raw[, annot_cols]
intensity  <- raw[, setdiff(names(raw), annot_cols)]

metadata <- as.data.frame(load_cvh_analysis_metadata(
  meta_file = cfg$meta_file, pheno_file = cfg$pheno_file,
  raw_sample_ids = colnames(intensity)))
# 3-group pooled scheme (supplement collapsed)
metadata$group_time <- dplyr::case_when(
  metadata$Group_Time %in% c("CRE_T1", "PLA_T1") ~ "CR_pre",
  metadata$Group_Time %in% c("CRE_T2", "PLA_T2") ~ "CR_post",
  metadata$Group_Time == "H_T1"                  ~ "H_pre"
)
stopifnot(!anyNA(metadata$group_time))
rownames(metadata) <- metadata$Col_ID
intensity <- intensity[, metadata$Col_ID]

n_raw <- nrow(annotation)
flog <- tibble(step = "Raw input", n_after = n_raw, n_removed = NA_integer_)
log_step <- function(flog, step, before, after)
  bind_rows(flog, tibble(step = step, n_after = after, n_removed = before - after))
cat(sprintf("Raw: %d proteins x %d samples\n", n_raw, ncol(intensity)))

# --- 2. HPA tissue-context filter --------------------------------------------
hpa <- read_tsv(cfg$hpa_file, show_col_types = FALSE) |>
  select(Gene, Protein_class = `Protein class`, Secretome = `Secretome location`) |>
  distinct(Gene, .keep_all = TRUE)
n0 <- nrow(annotation); keep <- annotation$gene %in% hpa$Gene
intensity <- intensity[keep, ]; annotation <- annotation[keep, ] |> left_join(hpa, by = c("gene" = "Gene"))
flog <- log_step(flog, "HPA tissue filter", n0, nrow(annotation))
cat(sprintf("HPA: %d -> %d\n", n0, nrow(annotation)))

# --- 2b. HPA-derived blood-contaminant removal -------------------------------
blood_genes <- blood_contaminant_genes(cfg$hpa_file)
n0 <- nrow(annotation); keep <- !annotation$gene %in% blood_genes
intensity <- intensity[keep, ]; annotation <- annotation[keep, ]
flog <- log_step(flog, "Blood contaminant removal", n0, nrow(annotation))
cat(sprintf("Blood: %d -> %d (-%d)\n", n0, nrow(annotation), n0 - nrow(annotation)))

# --- 3. Deduplicate by UniProt ID --------------------------------------------
if (any(duplicated(annotation$uniprot_id))) {
  n0 <- nrow(annotation)
  annotation$row_mean <- rowMeans(data.matrix(intensity), na.rm = TRUE)
  keep_idx <- annotation |> mutate(row_idx = row_number()) |>
    group_by(uniprot_id) |> slice_max(row_mean, n = 1, with_ties = FALSE) |> pull(row_idx)
  annotation <- annotation[keep_idx, ]; intensity <- intensity[keep_idx, ]; annotation$row_mean <- NULL
  flog <- log_step(flog, "Deduplication", n0, nrow(annotation))
}

# --- 4. Assemble DAList + missingness filter (by group_time) ------------------
int_mat <- as.data.frame(data.matrix(intensity)); rownames(int_mat) <- annotation$uniprot_id
annot_df <- as.data.frame(annotation); rownames(annot_df) <- annotation$uniprot_id
meta_df  <- as.data.frame(metadata);   rownames(meta_df)  <- metadata$Col_ID

dal <- zero_to_missing(DAList(data = int_mat, annotation = annot_df, metadata = meta_df))
n0 <- nrow(dal$data)
dal <- filter_proteins_by_group(dal, min_reps = cfg$min_reps,
                                min_groups = cfg$min_groups, grouping_column = "group_time")
flog <- log_step(flog, sprintf("Missingness (>=%d in >=%d group)", cfg$min_reps, cfg$min_groups),
                 n0, nrow(dal$data))
flog <- flog |> mutate(pct_of_raw = round(n_after / n_raw * 100, 1))
cat(sprintf("Missingness: %d -> %d\n", n0, nrow(dal$data)))

# --- 5. Outlier detection (4-method consensus, >=3/4) ------------------------
m <- dal$metadata
pct_missing <- colMeans(is.na(dal$data)) * 100
# (1) sample missingness, pooled + paired delta within CR subjects
paired <- m |> filter(Group != "PPS") |> count(Subject_ID) |> filter(n == 2) |> pull(Subject_ID)
delta <- setNames(rep(NA_real_, nrow(m)), m$Col_ID)
for (s in paired) {
  r <- m[m$Subject_ID == s, ]; t1 <- r$Col_ID[r$Timepoint == "T1"]; t2 <- r$Col_ID[r$Timepoint == "T2"]
  if (length(t1) == 1 && length(t2) == 1) delta[c(t1, t2)] <- abs(pct_missing[t2] - pct_missing[t1])
}
miss_thr  <- quantile(pct_missing, 0.75) + 1.5 * IQR(pct_missing)
dvals <- delta[!is.na(delta)]; delta_thr <- if (length(dvals) > 2) quantile(dvals, 0.75) + 1.5 * IQR(dvals) else Inf
miss_flag <- pct_missing > miss_thr | (!is.na(delta) & delta > delta_thr)
# (2) MAD median intensity
sm <- apply(log2(dal$data), 2, median, na.rm = TRUE)
mad_flag <- abs(sm - median(sm)) > cfg$mad_k * mad(sm)
# (3) PCA Mahalanobis on complete cases
pc3 <- run_pca(dal$data[rowSums(is.na(dal$data)) == 0, ])$x[, 1:3]
pca_flag <- mahalanobis(pc3, colMeans(pc3), cov(pc3)) > qchisq(1 - cfg$mahal_p, df = 3)
# (4) inter-sample correlation
cm <- cor(log2(dal$data), use = "pairwise.complete.obs")
mc <- apply(cm, 2, function(x) median(x[x < 1], na.rm = TRUE))
cor_flag <- mc < median(mc) - cfg$mad_k * mad(mc)

outlier_diag <- tibble(
  Col_ID = colnames(dal$data),
  miss_flag = miss_flag[Col_ID], mad_flag = mad_flag[Col_ID],
  pca_flag = pca_flag[Col_ID], cor_flag = cor_flag[Col_ID]) |>
  mutate(n_flags = miss_flag + mad_flag + pca_flag + cor_flag,
         consensus_outlier = n_flags >= cfg$outlier_k)
outlier_ids <- outlier_diag$Col_ID[outlier_diag$consensus_outlier]
cat(sprintf("Outliers: %d flagged (>=%d/4)%s\n", length(outlier_ids), cfg$outlier_k,
            if (length(outlier_ids)) paste0(": ", paste(outlier_ids, collapse = ", ")) else ""))
if (length(outlier_ids)) dal <- filter_samples(dal, !(Col_ID %in% outlier_ids))

# --- 6. Export filtered (un-normalized) artifacts ----------------------------
saveRDS(dal, file.path(cfg$data_dir, "01_DAList_filtered.rds"))
write_csv(bind_cols(as_tibble(dal$annotation) |> select(uniprot_id, protein, gene, description),
                    as_tibble(dal$data)),
          file.path(cfg$data_dir, "02_filtered_matrix.csv"))
write_csv(flog, file.path(cfg$data_dir, "03_filter_log.csv"))
saveRDS(list(filter_log = flog, outlier_diag = outlier_diag, outlier_ids = outlier_ids,
             blood_genes = blood_genes, n_raw = n_raw),
        file.path(cfg$data_dir, "00_filter_intermediates.rds"))
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
cat(sprintf("Done: %d proteins x %d samples -> %s/\n", nrow(dal$data), ncol(dal$data), cfg$data_dir))
print(flog)
