#!/usr/bin/env Rscript
# apply_missforest.R -- YvO-style MAR/MNAR classification + missForest imputation
#
# Outputs:
#   c_data/01_imputed.csv
#   c_data/01_DAList_imputed.rds
#   c_data/02_imputation.xlsx
#   c_data/02_mar_mnar_classification.csv
#   c_data/07_imputation_mask.csv
#   c_data/08_mnar_imputation_audit.csv
#   c_data/09_imputation_summary.txt
#   c_data/00_report_intermediates.rds

library(missForest)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(openxlsx)

set.seed(42)
setwd(rprojroot::find_rstudio_root_file())
source("R/cvh_design.R")

cfg <- list(
  norm_csv = "01_normalization/c_data/02_normalized.csv",
  norm_rds = "01_normalization/c_data/03_DAList_normalized.rds",
  data_dir = "02_Imputation/c_data",
  miss_unreliable = 50
)
dir.create(cfg$data_dir, showWarnings = FALSE, recursive = TRUE)

PAL_GT <- c(
  CRE_T1 = "#2166AC", CRE_T2 = "#67A9CF",
  PLA_T1 = "#D6604D", PLA_T2 = "#F4A582",
  H_T1   = "#4DAF4A"
)
PAL_MAR <- c(MAR = "#4393C3", MNAR = "#D6604D")
PAL_CLASS <- c(Complete = "#4DAF4A", MAR = "#4393C3", MNAR = "#D6604D")

write_sheet <- function(wb, name, data) {
  addWorksheet(wb, name)
  writeData(
    wb, name, data,
    headerStyle = createStyle(textDecoration = "bold", fgFill = "#DCE6F1")
  )
  freezePane(wb, name, firstRow = TRUE)
  setColWidths(wb, name, cols = seq_len(ncol(data)), widths = "auto")
}

df <- readr::read_csv(cfg$norm_csv, show_col_types = FALSE)
ann <- df |>
  dplyr::select(uniprot_id, gene, protein, description)
mat <- as.matrix(df[, setdiff(names(df), names(ann))])
rownames(mat) <- ann$uniprot_id

dal_norm <- readRDS(cfg$norm_rds)
meta <- tibble::as_tibble(dal_norm$metadata) |>
  dplyr::select(Col_ID, Subject_ID, Group, Supplement, Timepoint, Group_Time)
assert_cvh_design_rules(
  meta,
  context = "normalized CvH metadata",
  allow_t2_only_singletons = TRUE
)

stopifnot(setequal(meta$Col_ID, colnames(mat)))
cat(sprintf("Loaded: %d proteins x %d samples\n", nrow(mat), ncol(mat)))

prot_miss <- rowSums(is.na(mat))
prot_pct <- prot_miss / ncol(mat) * 100
obs_means <- rowMeans(mat, na.rm = TRUE)
pct_miss <- round(sum(is.na(mat)) / length(mat) * 100, 2)
cat(sprintf("Missing: %d / %d (%.2f%%)\n", sum(is.na(mat)), length(mat), pct_miss))

miss_by_group <- sapply(unique(meta$Group_Time), function(gt) {
  cols <- meta$Col_ID[meta$Group_Time == gt]
  rowSums(is.na(mat[, cols, drop = FALSE])) / length(cols) * 100
})

has_na <- which(prot_miss > 0 & prot_miss < ncol(mat))
inc_mean <- obs_means[has_na]
inc_pct <- prot_pct[has_na]

set.seed(42)
km <- kmeans(scale(cbind(inc_mean, inc_pct)), centers = 2, nstart = 25)
km_mnar <- km$cluster == which.min(tapply(inc_mean, km$cluster, mean))

lr_df <- data.frame(
  is_miss = as.integer(is.na(as.vector(mat))),
  intensity = rep(obs_means, ncol(mat))
)
lr_fit <- glm(is_miss ~ intensity, data = lr_df, family = binomial)
lr_pred <- predict(
  lr_fit,
  newdata = data.frame(intensity = inc_mean),
  type = "response"
)
lr_mnar <- lr_pred > median(lr_pred)

global_q25 <- quantile(mat, 0.25, na.rm = TRUE)
tail_frac <- vapply(
  has_na,
  function(i) mean(mat[i, !is.na(mat[i, ])] < global_q25),
  numeric(1)
)
lt_mnar <- (tail_frac * inc_pct / 100) > median(tail_frac * inc_pct / 100)

votes <- as.integer(km_mnar) + as.integer(lr_mnar) + as.integer(lt_mnar)
classification_method <- "3-method consensus (kmeans + logistic + left-tail)"

miss_class <- tibble(
  uniprot_id = rownames(mat),
  gene = ann$gene,
  n_miss = prot_miss,
  pct_miss = prot_pct,
  mean_intensity = obs_means,
  vote_kmeans = NA_integer_,
  vote_logistic = NA_integer_,
  vote_lefttail = NA_integer_,
  n_mnar_votes = NA_integer_
)
miss_class$vote_kmeans[has_na] <- as.integer(km_mnar)
miss_class$vote_logistic[has_na] <- as.integer(lr_mnar)
miss_class$vote_lefttail[has_na] <- as.integer(lt_mnar)
miss_class$n_mnar_votes[has_na] <- votes

miss_class <- miss_class |>
  dplyr::mutate(
    classification = dplyr::case_when(
      n_miss == 0 ~ "Complete",
      n_miss >= ncol(mat) ~ "MNAR",
      n_mnar_votes >= 2 ~ "MNAR",
      TRUE ~ "MAR"
    ),
    imputation_reliable = classification == "Complete" | pct_miss < cfg$miss_unreliable
  )

mnar_ids <- miss_class$uniprot_id[miss_class$classification == "MNAR"]

group_miss_pval <- setNames(rep(NA_real_, nrow(miss_class)), miss_class$uniprot_id)
for (uid in mnar_ids) {
  ct <- sapply(unique(meta$Group_Time), function(gt) {
    cols <- meta$Col_ID[meta$Group_Time == gt]
    c(
      missing = sum(is.na(mat[uid, cols])),
      observed = sum(!is.na(mat[uid, cols]))
    )
  })
  group_miss_pval[uid] <- tryCatch(
    fisher.test(ct, simulate.p.value = TRUE, B = 2000)$p.value,
    error = function(e) NA_real_
  )
}
miss_class$group_miss_pval <- group_miss_pval[miss_class$uniprot_id]

n_mar_prots <- sum(miss_class$classification == "MAR")
n_mnar_prots <- length(mnar_ids)
n_comp_prots <- sum(miss_class$classification == "Complete")
mar_miss_vals <- sum(miss_class$n_miss[miss_class$classification == "MAR"])
mnar_miss_vals <- sum(miss_class$n_miss[miss_class$classification == "MNAR"])
total_miss_vals <- mar_miss_vals + mnar_miss_vals
cat(sprintf(
  "Classification: MAR %d | MNAR %d | Complete %d\n",
  n_mar_prots, n_mnar_prots, n_comp_prots
))

cat("Imputing with missForest...\n")
gene_order <- order(rownames(mat))
mat <- mat[gene_order, , drop = FALSE]
ann <- ann[gene_order, , drop = FALSE]

set.seed(42)
mf_result <- missForest::missForest(t(mat), maxiter = 10, ntree = 100, verbose = TRUE)
mat_imp <- t(mf_result$ximp)
rownames(mat_imp) <- rownames(mat)
colnames(mat_imp) <- colnames(mat)
stopifnot(sum(is.na(mat_imp)) == 0)

oob_error <- as.numeric(mf_result$OOBerror[1])
cat(sprintf("missForest OOB error: %.4f\n", oob_error))

mnar_audit <- miss_class |>
  dplyr::filter(uniprot_id %in% mnar_ids) |>
  dplyr::transmute(
    uniprot_id,
    gene,
    pre_mean = rowMeans(mat[uniprot_id, , drop = FALSE], na.rm = TRUE),
    post_mean = rowMeans(mat_imp[uniprot_id, , drop = FALSE]),
    pre_sd = apply(mat[uniprot_id, , drop = FALSE], 1, sd, na.rm = TRUE),
    pct_miss,
    shift = post_mean - pre_mean,
    effect_d = (post_mean - pre_mean) / pre_sd,
    imputation_reliable = pct_miss < cfg$miss_unreliable
  )

imp_df <- bind_cols(ann, as_tibble(mat_imp))
mask_df <- bind_cols(tibble(uniprot_id = rownames(mat)), as_tibble(is.na(mat) + 0L))
summary_df <- tibble(
  metric = c(
    "n_proteins", "n_samples", "pct_missing", "n_complete",
    "n_mar_proteins", "n_mnar_proteins", "n_mar_values", "n_mnar_values",
    "classification_method", "method", "oob_error"
  ),
  value = c(
    nrow(mat), ncol(mat), pct_miss, n_comp_prots,
    n_mar_prots, n_mnar_prots, mar_miss_vals, mnar_miss_vals,
    classification_method, "missForest", round(oob_error, 4)
  )
)

wb <- createWorkbook()
write_sheet(wb, "imputed_matrix", imp_df)
write_sheet(wb, "mar_mnar_classification", as.data.frame(miss_class))
write_sheet(wb, "imputation_mask", mask_df)
write_sheet(wb, "mnar_audit", as.data.frame(mnar_audit))
write_sheet(wb, "imputation_summary", summary_df)

benchmark_path <- file.path(cfg$data_dir, "benchmark", "04_composite_ranking.csv")
if (file.exists(benchmark_path)) {
  benchmark_df <- readr::read_csv(benchmark_path, show_col_types = FALSE)
  write_sheet(wb, "benchmark_ranking", as.data.frame(benchmark_df))
}
saveWorkbook(wb, file.path(cfg$data_dir, "02_imputation.xlsx"), overwrite = TRUE)

write_csv(imp_df, file.path(cfg$data_dir, "01_imputed.csv"))
write_csv(miss_class, file.path(cfg$data_dir, "02_mar_mnar_classification.csv"))
write_csv(mask_df, file.path(cfg$data_dir, "07_imputation_mask.csv"))
write_csv(mnar_audit, file.path(cfg$data_dir, "08_mnar_imputation_audit.csv"))

summary_lines <- c(
  paste("n_proteins =", nrow(mat)),
  paste("n_samples =", ncol(mat)),
  paste("pct_missing =", pct_miss),
  paste("n_complete =", n_comp_prots),
  paste("n_mar_proteins =", n_mar_prots),
  paste("n_mnar_proteins =", n_mnar_prots),
  paste("n_mar_values =", mar_miss_vals),
  paste("n_mnar_values =", mnar_miss_vals),
  paste(
    "pct_mar_values =",
    if (total_miss_vals > 0) round(mar_miss_vals / total_miss_vals * 100, 1) else NA_real_
  ),
  paste("classification_method =", classification_method),
  paste("n_unreliable =", sum(!miss_class$imputation_reliable)),
  "best_method = missForest",
  paste("oob_error =", round(oob_error, 4))
)
writeLines(summary_lines, file.path(cfg$data_dir, "09_imputation_summary.txt"))

dal <- dal_norm
dal$data <- mat_imp
n_ann <- nrow(dal$annotation)
dal$annotation <- merge(
  dal$annotation,
  miss_class |>
    dplyr::select(
      uniprot_id, n_miss, pct_miss,
      miss_classification = classification,
      imputation_reliable
    ),
  by = "uniprot_id",
  all.x = TRUE,
  sort = FALSE
)
stopifnot(nrow(dal$annotation) == n_ann)
# Re-align $annotation rows to $data row order. mat_imp was reordered by
# gene_order for missForest determinism; merge() preserves left-frame order.
# Without this match() step the saved DAList has the same set of proteins
# in $data and $annotation but at different row positions.
dal$annotation <- dal$annotation[
  match(rownames(dal$data), dal$annotation$uniprot_id), , drop = FALSE]
rownames(dal$annotation) <- dal$annotation$uniprot_id
stopifnot(identical(rownames(dal$data), dal$annotation$uniprot_id))
saveRDS(dal, file.path(cfg$data_dir, "01_DAList_imputed.rds"))

saveRDS(
  list(
    mat = mat,
    mat_imp = mat_imp,
    was_na = is.na(mat),
    ann = ann,
    meta = meta,
    miss_class = miss_class,
    miss_by_group = miss_by_group,
    prot_pct = prot_pct,
    pct_miss = pct_miss,
    mnar_ids = mnar_ids,
    mnar_audit = mnar_audit,
    best = "missForest",
    n_mar_prots = n_mar_prots,
    n_mnar_prots = n_mnar_prots,
    mar_miss_vals = mar_miss_vals,
    mnar_miss_vals = mnar_miss_vals,
    total_miss_vals = total_miss_vals,
    classification_method = classification_method,
    oob_error = oob_error,
    PAL_GT = PAL_GT,
    PAL_MAR = PAL_MAR,
    PAL_CLASS = PAL_CLASS
  ),
  file.path(cfg$data_dir, "00_report_intermediates.rds")
)

writeLines(capture.output(sessionInfo()), file.path(cfg$data_dir, "sessionInfo.txt"))
cat(sprintf(
  "Done: missForest | %d proteins x %d samples | OOB=%.4f\n",
  nrow(mat_imp), ncol(mat_imp), oob_error
))
