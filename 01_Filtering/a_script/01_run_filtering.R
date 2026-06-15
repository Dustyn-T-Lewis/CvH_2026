#!/usr/bin/env Rscript
# =============================================================================
# 01_run_filtering.R  --  CvH Stage 01: Filtering (proteoDA-native)
#
# load -> dedup -> build DAList -> annotate (HPA presence + blood verdict) ->
# native filters: filter_proteins_by_annotation (HPA, blood) ->
# filter_proteins_by_group (missingness) -> 4-method outlier consensus ->
# filter_samples. Outputs a FILTERED, un-normalized DAList for Stage 02.
#
# Native proteoDA throughout; only the contaminant classification (R/blood_filter.R)
# and the 4-method outlier consensus are custom (no native equivalent).
# Design: 3-group pooled scheme group_time in {H_pre, CR_pre, CR_post}.
# =============================================================================

suppressPackageStartupMessages({
  library(proteoDA); library(here); library(readxl); library(readr); library(dplyr); library(tidyr)
})
set.seed(42)
source(here("R", "cvh_design.R"))
source(here("R", "blood_filter.R"))

cfg <- list(
  raw_file   = here("00_input", "CvH_raw.xlsx"),
  meta_file  = here("00_input", "CvH_meta.csv"),
  pheno_file = here("00_input", "CRm_meta.csv"),
  hpa_file   = here("00_input", "HPA_annotations.tsv"),
  data_dir   = here("01_Filtering", "c_data"),
  min_reps = 5L, min_groups = 1L, outlier_k = 3, mad_k = 3, mahal_p = 0.01
)
dir.create(cfg$data_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load + metadata + 3-group scheme + dedup -----------------------------
raw <- read_excel(cfg$raw_file)
annot_cols <- c("uniprot_id", "protein", "gene", "description", "n_seq")
annotation <- raw[, annot_cols]
intensity  <- raw[, setdiff(names(raw), annot_cols)]

metadata <- as.data.frame(load_cvh_analysis_metadata(
  meta_file = cfg$meta_file, pheno_file = cfg$pheno_file, raw_sample_ids = colnames(intensity)))
metadata$group_time <- dplyr::case_when(
  metadata$Group_Time %in% c("CRE_T1", "PLA_T1") ~ "CR_pre",
  metadata$Group_Time %in% c("CRE_T2", "PLA_T2") ~ "CR_post",
  metadata$Group_Time == "H_T1"                  ~ "H_pre")
stopifnot(!anyNA(metadata$group_time))
rownames(metadata) <- metadata$Col_ID
intensity <- intensity[, metadata$Col_ID]
n_raw <- nrow(annotation)

# dedup by UniProt (keep highest mean intensity; no native equivalent)
if (any(duplicated(annotation$uniprot_id))) {
  rm_mean <- rowMeans(data.matrix(intensity), na.rm = TRUE)
  keep_idx <- tibble(i = seq_along(rm_mean), id = annotation$uniprot_id, m = rm_mean) |>
    group_by(id) |> slice_max(m, n = 1, with_ties = FALSE) |> pull(i)
  annotation <- annotation[keep_idx, ]; intensity <- intensity[keep_idx, ]
}

# --- 2. Contaminant + HPA-presence removal -----------------------------------
# NB: proteoDA::filter_proteins_by_annotation() is broken under R >= 4.4 (its
# internal length-2 class() check errors), so the annotation-based protein
# removal is done by direct subsetting here; the group/sample filters below stay
# native (they use inherits() and work fine).
bl <- classify_blood_reference(cfg$hpa_file)             # full HPA table + verdict
hpa_genes    <- unique(bl$gene)
remove_genes <- bl$gene[bl$verdict == "remove"]

flog <- tibble(step = "Raw input", n_after = nrow(annotation), n_removed = NA_integer_)
keep <- annotation$gene %in% hpa_genes
flog <- bind_rows(flog, tibble(step = "HPA presence", n_after = sum(keep), n_removed = sum(!keep)))
annotation <- annotation[keep, ]; intensity <- intensity[keep, ]
keep <- !(annotation$gene %in% remove_genes)
flog <- bind_rows(flog, tibble(step = "Blood contaminant removal", n_after = sum(keep), n_removed = sum(!keep)))
annotation <- annotation[keep, ]; intensity <- intensity[keep, ]

# --- 3. Build DAList + native missingness filter -----------------------------
int_mat <- as.data.frame(data.matrix(intensity)); rownames(int_mat) <- annotation$uniprot_id
annot_df <- as.data.frame(annotation); rownames(annot_df) <- annotation$uniprot_id
meta_df  <- as.data.frame(metadata);   rownames(meta_df)  <- metadata$Col_ID
dal <- zero_to_missing(DAList(data = int_mat, annotation = annot_df, metadata = meta_df))

n0 <- nrow(dal$data)
dal <- filter_proteins_by_group(dal, min_reps = cfg$min_reps,
                                min_groups = cfg$min_groups, grouping_column = "group_time")
flog <- bind_rows(flog, tibble(step = sprintf("Missingness (>=%d in >=%d group)", cfg$min_reps, cfg$min_groups),
                               n_after = nrow(dal$data), n_removed = n0 - nrow(dal$data)))
flog <- flog |> mutate(pct_of_raw = round(n_after / n_raw * 100, 1))

# --- 4. Outlier consensus (4-method, >=3/4) -> native filter_samples ---------
m <- dal$metadata; pct_missing <- colMeans(is.na(dal$data)) * 100
paired <- m |> filter(Group != "PPS") |> count(Subject_ID) |> filter(n == 2) |> pull(Subject_ID)
delta <- setNames(rep(NA_real_, nrow(m)), m$Col_ID)
for (s in paired) {
  r <- m[m$Subject_ID == s, ]; t1 <- r$Col_ID[r$Timepoint == "T1"]; t2 <- r$Col_ID[r$Timepoint == "T2"]
  if (length(t1) == 1 && length(t2) == 1) delta[c(t1, t2)] <- abs(pct_missing[t2] - pct_missing[t1])
}
miss_thr <- quantile(pct_missing, 0.75) + 1.5 * IQR(pct_missing)
dv <- delta[!is.na(delta)]; delta_thr <- if (length(dv) > 2) quantile(dv, 0.75) + 1.5 * IQR(dv) else Inf
miss_flag <- pct_missing > miss_thr | (!is.na(delta) & delta > delta_thr)
sm <- apply(log2(dal$data), 2, median, na.rm = TRUE); mad_flag <- abs(sm - median(sm)) > cfg$mad_k * mad(sm)
pc3 <- prcomp(t(log2(dal$data[rowSums(is.na(dal$data)) == 0, ])), center = TRUE, scale. = TRUE)$x[, 1:3]
pca_flag <- mahalanobis(pc3, colMeans(pc3), cov(pc3)) > qchisq(1 - cfg$mahal_p, df = 3)
cm <- cor(log2(dal$data), use = "pairwise.complete.obs"); mc <- apply(cm, 2, function(x) median(x[x < 1], na.rm = TRUE))
cor_flag <- mc < median(mc) - cfg$mad_k * mad(mc)
outlier_diag <- tibble(Col_ID = colnames(dal$data), miss_flag = miss_flag[Col_ID], mad_flag = mad_flag[Col_ID],
                       pca_flag = pca_flag[Col_ID], cor_flag = cor_flag[Col_ID]) |>
  mutate(n_flags = miss_flag + mad_flag + pca_flag + cor_flag, consensus_outlier = n_flags >= cfg$outlier_k)
outlier_ids <- outlier_diag$Col_ID[outlier_diag$consensus_outlier]
cat(sprintf("Outliers (>=%d/4): %s\n", cfg$outlier_k, if (length(outlier_ids)) paste(outlier_ids, collapse = ", ") else "none"))
if (length(outlier_ids)) dal <- filter_samples(dal, !(Col_ID %in% outlier_ids))

# --- 5. Export ---------------------------------------------------------------
saveRDS(dal, file.path(cfg$data_dir, "01_DAList_filtered.rds"))
write_csv(bind_cols(as_tibble(dal$annotation) |> select(uniprot_id, protein, gene, description),
                    as_tibble(dal$data)), file.path(cfg$data_dir, "02_filtered_matrix.csv"))
write_csv(flog, file.path(cfg$data_dir, "03_filter_log.csv"))
saveRDS(list(filter_log = flog, outlier_diag = outlier_diag, outlier_ids = outlier_ids, n_raw = n_raw),
        file.path(cfg$data_dir, "00_filter_intermediates.rds"))
if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
cat(sprintf("Done: %d proteins x %d samples -> %s/\n", nrow(dal$data), ncol(dal$data), cfg$data_dir))
print(as.data.frame(flog))
