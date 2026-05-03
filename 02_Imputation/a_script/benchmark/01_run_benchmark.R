# CvH Imputation — MAR/MNAR classification, 12-method benchmark, best-method apply
# Refs: Lazar 2016 (hybrid), Hediyeh-zadeh 2023 (EBM), Wei 2018 (Procrustes)

library(MsCoreUtils)
library(msImpute)
library(pcaMethods)
library(missMDA)
library(vegan)
library(dplyr)
library(tidyr)
library(readr)

setwd(rprojroot::find_rstudio_root_file())

cfg <- list(
  norm_csv        = "01_normalization/c_data/02_normalized.csv",
  norm_rds        = "01_normalization/c_data/03_DAList_normalized.rds",
  data_dir        = "02_Imputation/c_data/benchmark",
  n_iter          = 20L,
  mask_frac       = 0.10,
  miss_unreliable = 50,
  methods = list(
    MinProb = list(method = "MinProb"), MinDet = list(method = "MinDet"),
    QRILC   = list(method = "QRILC"),   zero   = list(method = "zero"),
    knn     = list(method = "knn"),      bpca   = list(method = "bpca"),
    RF      = list(method = "RF"),       SVD    = list(method = "SVD"),
    imputePCA      = list(method = "imputePCA"),
    mix_bpca_QRILC = list(method = "mixed", mar = "bpca", mnar = "QRILC"),
    mix_knn_QRILC  = list(method = "mixed", mar = "knn",  mnar = "QRILC"),
    mix_RF_QRILC   = list(method = "mixed", mar = "RF",   mnar = "QRILC")
  )
)

METHOD_TYPE <- tibble(method = names(cfg$methods)) |>
  mutate(type = case_when(
    method %in% c("MinProb", "MinDet", "QRILC", "zero") ~ "MNAR",
    startsWith(method, "mix_") ~ "Hybrid",
    TRUE ~ "MAR"))

dir.create(cfg$data_dir, showWarnings = FALSE, recursive = TRUE)

# --- Helpers -----------------------------------------------------------------
run_impute <- function(m, mat, randna) {
  if (m$method == "mixed")
    MsCoreUtils::impute_matrix(mat, method = "mixed", randna = randna,
                               mar = m$mar, mnar = m$mnar)
  else if (m$method == "imputePCA")
    missMDA::imputePCA(mat, ncp = 2, method = "Regularized")$completeObs
  else if (m$method == "SVD")
    pcaMethods::pca(mat, method = "svdImpute", nPcs = 2, verbose = FALSE)@completeObs
  else
    MsCoreUtils::impute_matrix(mat, method = m$method)
}

nrmse <- function(true_vals, imp_vals) {
  sqrt(mean((true_vals - imp_vals)^2)) / sd(true_vals)
}

pss_metric <- function(original, imputed, n_pc = 5) {
  complete_prots <- complete.cases(original)
  if (sum(complete_prots) < n_pc + 1) return(NA_real_)
  n_pc_use <- min(n_pc, ncol(original) - 1)
  pca_orig <- prcomp(t(original[complete_prots, ]), center = TRUE, scale. = TRUE)
  pca_imp  <- prcomp(t(imputed[complete_prots, ]),  center = TRUE, scale. = TRUE)
  vegan::procrustes(pca_orig$x[, seq_len(n_pc_use)],
                    pca_imp$x[, seq_len(n_pc_use)])$ss
}

# --- 1. Load data ------------------------------------------------------------
df  <- read_csv(cfg$norm_csv, show_col_types = FALSE)
ann <- df |> select(uniprot_id, gene, protein, description)
mat <- as.matrix(df[, -(1:4)])
rownames(mat) <- ann$uniprot_id
cat(sprintf("%d proteins x %d samples\n", nrow(mat), ncol(mat)))

# Metadata from normalized DAList (not regex-derived)
dal_norm <- readRDS(cfg$norm_rds)
meta <- tibble(
  Col_ID     = dal_norm$metadata$Col_ID,
  Subject_ID = dal_norm$metadata$Subject_ID,
  Group      = dal_norm$metadata$Group,
  Timepoint  = dal_norm$metadata$Timepoint,
  Group_Time = dal_norm$metadata$Group_Time
)
stopifnot(setequal(meta$Col_ID, colnames(mat)))

# --- 2. Missingness profiling ------------------------------------------------
prot_miss <- rowSums(is.na(mat))
prot_pct  <- prot_miss / ncol(mat) * 100
obs_means <- rowMeans(mat, na.rm = TRUE)
pct_miss  <- round(sum(is.na(mat)) / length(mat) * 100, 2)

cat(sprintf("Missing: %d / %d (%.2f%%) | Complete: %d\n", sum(is.na(mat)), length(mat), pct_miss, sum(prot_miss == 0)))

miss_by_group <- sapply(unique(meta$Group_Time), function(g) {
  cols <- meta$Col_ID[meta$Group_Time == g]
  rowSums(is.na(mat[, cols, drop = FALSE])) / length(cols) * 100
})

# --- 3. MAR/MNAR classification ----------------------------------------------
has_na <- which(prot_miss > 0 & prot_miss < ncol(mat))
miss_class <- tibble(uniprot_id = rownames(mat), gene = ann$gene,
                     n_miss = prot_miss, pct_miss = prot_pct,
                     mean_intensity = obs_means)

mar_result <- tryCatch({
  feat <- msImpute::selectFeatures(mat[has_na, ], method = "ebm",
                                   group = meta$Group_Time)
  mar_names <- feat$name[feat$msImpute_feature]
  n_incomplete <- nrow(feat)
  cat(sprintf("EBM result: %d MAR / %d MNAR (of %d incomplete proteins)\n",
              length(mar_names), n_incomplete - length(mar_names), n_incomplete))
  # Guard against degenerate classification
  if (length(mar_names) < 0.05 * n_incomplete) {
    cat("  EBM degenerate (<5% MAR) -- falling back to k-means\n")
    NULL
  } else {
    list(mar_ids = mar_names, method = "msImpute_ebm")
  }
}, error = function(e) {
  cat(sprintf("msImpute EBM failed (%s) -- falling back to k-means\n",
              conditionMessage(e)))
  NULL
})

if (is.null(mar_result)) {
  mc_sub <- miss_class |> filter(n_miss > 0, n_miss < ncol(mat))
  km <- kmeans(scale(cbind(mc_sub$mean_intensity, mc_sub$pct_miss)),
               centers = 2, nstart = 25)
  cl_means <- tapply(mc_sub$mean_intensity, km$cluster, mean)
  mnar_cl <- which.min(cl_means)
  mar_result <- list(
    mar_ids = mc_sub$uniprot_id[km$cluster != mnar_cl],
    method = sprintf("k-means (cluster means: %.1f vs %.1f)",
                     cl_means[mnar_cl], cl_means[-mnar_cl]))
}

cat(sprintf("Classification: %s\n", mar_result$method))

miss_class <- miss_class |>
  mutate(
    classification = case_when(
      n_miss == 0 ~ "Complete",
      uniprot_id %in% mar_result$mar_ids ~ "MAR",
      TRUE ~ "MNAR"),
    imputation_reliable = classification == "Complete" | pct_miss < cfg$miss_unreliable
  )

# Group-stratified missingness (Fisher test per MNAR protein)
mnar_ids <- miss_class$uniprot_id[miss_class$classification == "MNAR"]
group_miss_pval <- setNames(rep(NA_real_, nrow(miss_class)), miss_class$uniprot_id)

for (g in mnar_ids) {
  ct <- sapply(unique(meta$Group_Time), function(gt) {
    cols <- meta$Col_ID[meta$Group_Time == gt]
    c(missing = sum(is.na(mat[g, cols])), observed = sum(!is.na(mat[g, cols])))
  })
  group_miss_pval[g] <- tryCatch(
    fisher.test(ct, simulate.p.value = TRUE, B = 2000)$p.value,
    error = function(e) NA_real_)
}
miss_class$group_miss_pval <- group_miss_pval[miss_class$uniprot_id]

n_sig <- sum(group_miss_pval[mnar_ids] < 0.05, na.rm = TRUE)

n_mar_prots  <- sum(miss_class$classification == "MAR")
n_mnar_prots <- length(mnar_ids)
n_comp_prots <- sum(miss_class$classification == "Complete")
mar_miss_vals  <- sum(miss_class$n_miss[miss_class$classification == "MAR"])
mnar_miss_vals <- sum(miss_class$n_miss[miss_class$classification == "MNAR"])
total_miss_vals <- mar_miss_vals + mnar_miss_vals

cat(sprintf("Proteins: MAR %d | MNAR %d (%d group-biased) | Complete %d | Values: MAR %d (%.0f%%) MNAR %d (%.0f%%)\n",
            n_mar_prots, n_mnar_prots, n_sig, n_comp_prots,
            mar_miss_vals, mar_miss_vals / total_miss_vals * 100,
            mnar_miss_vals, mnar_miss_vals / total_miss_vals * 100))

# --- 4. Benchmark ------------------------------------------------------------
set.seed(42)
randna <- setNames(miss_class$classification != "MNAR", miss_class$uniprot_id)
bench_ids   <- miss_class$uniprot_id[randna[miss_class$uniprot_id]]
bench_means <- rowMeans(mat[bench_ids, , drop = FALSE], na.rm = TRUE)
bin_breaks  <- unique(quantile(bench_means, probs = c(0, 1/3, 2/3, 1)))
bin_labels  <- c("low", "mid", "high")[seq_len(length(bin_breaks) - 1)]
prot_bins   <- cut(bench_means, breaks = bin_breaks, labels = bin_labels,
                   include.lowest = TRUE)
names(prot_bins) <- bench_ids

cat(sprintf("Benchmarking: %d methods x %d iterations\n",
            length(cfg$methods), cfg$n_iter))
bench_raw <- vector("list", length(cfg$methods) * cfg$n_iter)
bin_raw   <- vector("list", length(cfg$methods) * cfg$n_iter * length(bin_labels))
k <- 0L; kb <- 0L

for (iter in seq_len(cfg$n_iter)) {
  if (iter %% 5 == 0) cat(sprintf("  Iter %d/%d\n", iter, cfg$n_iter))

  mar_obs_idx <- which(!is.na(mat) & randna[row(mat)])
  mask_idx    <- sample(mar_obs_idx, round(length(mar_obs_idx) * cfg$mask_frac))
  true_v      <- mat[mask_idx]
  mm <- mat; mm[mask_idx] <- NA
  mask_rows <- rownames(mat)[((mask_idx - 1) %% nrow(mat)) + 1]
  mask_bin  <- prot_bins[mask_rows]

  for (nm in names(cfg$methods)) {
    imp <- tryCatch(run_impute(cfg$methods[[nm]], mm, randna),
                    error = function(e) NULL)
    if (is.null(imp)) next

    k <- k + 1L
    nrmse_val <- nrmse(true_v, imp[mask_idx])
    pss_val   <- tryCatch(pss_metric(mat, imp), error = function(e) NA_real_)
    bench_raw[[k]] <- tibble(method = nm, iter = iter,
                             nrmse = nrmse_val, pss = pss_val)

    for (b in bin_labels) {
      sel <- which(mask_bin == b & !is.na(mask_bin))
      if (length(sel) < 5) next
      kb <- kb + 1L
      bin_raw[[kb]] <- tibble(method = nm, iter = iter, bin = b,
                              nrmse = nrmse(true_v[sel], imp[mask_idx[sel]]))
    }
  }
}

bench_df  <- bind_rows(bench_raw)
bench_sum <- bench_df |>
  group_by(method) |>
  summarise(mean_nrmse   = mean(nrmse),
            sd_nrmse     = sd(nrmse),
            median_nrmse = median(nrmse),
            mean_pss     = mean(pss, na.rm = TRUE),
            sd_pss       = sd(pss, na.rm = TRUE),
            .groups = "drop") |>
  arrange(mean_nrmse)

bin_df  <- bind_rows(bin_raw)
bin_sum <- bin_df |>
  group_by(method, bin) |>
  summarise(mean_nrmse = mean(nrmse), sd_nrmse = sd(nrmse), .groups = "drop")

top5 <- bench_sum$method[1:min(5, nrow(bench_sum))]
ext_sum <- bench_sum |>
  filter(method %in% top5) |>
  select(method, nrmse = mean_nrmse, pss = mean_pss) |>
  left_join(
    bin_sum |> filter(method %in% top5) |>
      select(method, bin, mean_nrmse) |>
      pivot_wider(names_from = bin, values_from = mean_nrmse, names_prefix = "nrmse_"),
    by = "method")

best <- bench_sum$method[1]
cat(sprintf("Best: %s (NRMSE %.4f, PSS %.4f)\n", best, bench_sum$mean_nrmse[1], bench_sum$mean_pss[1]))

# --- 5. Apply best method ----------------------------------------------------
set.seed(42)
mat_imp <- run_impute(cfg$methods[[best]], mat, randna)
stopifnot(sum(is.na(mat_imp)) == 0)

# --- 6. Export data ----------------------------------------------------------
was_na <- is.na(mat)
stopifnot(identical(ann$uniprot_id, rownames(mat_imp)))

mnar_audit <- tibble(
  uniprot_id = mnar_ids,
  gene       = ann$gene[match(mnar_ids, ann$uniprot_id)],
  pre_mean   = rowMeans(mat[mnar_ids, ], na.rm = TRUE),
  post_mean  = rowMeans(mat_imp[mnar_ids, ]),
  pre_sd     = apply(mat[mnar_ids, ], 1, sd, na.rm = TRUE),
  pct_miss   = prot_pct[mnar_ids],
  shift      = post_mean - pre_mean,
  effect_d   = shift / pre_sd,
  imputation_reliable = prot_pct[mnar_ids] < cfg$miss_unreliable
)

write_csv(bind_cols(ann, as_tibble(mat_imp)),
          file.path(cfg$data_dir, "01_imputed.csv"))

dal <- readRDS(cfg$norm_rds)
dal$data <- mat_imp
dal$annotation <- merge(
  dal$annotation,
  miss_class |> select(uniprot_id, n_miss, pct_miss,
                       miss_classification = classification,
                       imputation_reliable),
  by = "uniprot_id", all.x = TRUE, sort = FALSE)
saveRDS(dal, file.path(cfg$data_dir, "01_DAList_imputed.rds"))

write_csv(miss_class, file.path(cfg$data_dir, "02_mar_mnar_classification.csv"))
write_csv(bench_sum,  file.path(cfg$data_dir, "03_benchmark_summary.csv"))
write_csv(bench_df,   file.path(cfg$data_dir, "04_benchmark_raw_iterations.csv"))
write_csv(ext_sum,    file.path(cfg$data_dir, "05_benchmark_extended.csv"))
write_csv(bin_sum,    file.path(cfg$data_dir, "06_benchmark_per_intensity.csv"))
write_csv(bind_cols(tibble(uniprot_id = rownames(was_na)), as_tibble(was_na)),
          file.path(cfg$data_dir, "07_imputation_mask.csv"))
write_csv(mnar_audit, file.path(cfg$data_dir, "08_mnar_imputation_audit.csv"))

info <- list(
  n_proteins   = nrow(mat), n_samples = ncol(mat), pct_missing = pct_miss,
  n_complete   = n_comp_prots, n_mar_proteins = n_mar_prots,
  n_mnar_proteins = n_mnar_prots, n_mar_values = mar_miss_vals,
  n_mnar_values = mnar_miss_vals,
  pct_mar_values = round(mar_miss_vals / total_miss_vals * 100, 1),
  classification_method = mar_result$method,
  n_unreliable = sum(!miss_class$imputation_reliable),
  best_method  = best,
  best_nrmse   = round(bench_sum$mean_nrmse[1], 4),
  best_pss     = round(bench_sum$mean_pss[1], 4))
writeLines(paste(names(info), info, sep = " = "),
           file.path(cfg$data_dir, "09_imputation_summary.txt"))

# --- 7. Per-Sample LOOCV (best method, 20 iter x 10% mask) ------------------
cat(sprintf("Computing per-sample LOOCV (%s, 20 iter x 10%% mask)...\n", best))

samp_miss <- tibble(
  Col_ID   = colnames(mat),
  pct_miss = colSums(is.na(mat)) / nrow(mat) * 100) |>
  left_join(meta |> select(Col_ID, Group_Time), by = "Col_ID")

loocv_samples <- bind_rows(
  samp_miss |> slice_max(pct_miss, n = 5, with_ties = FALSE) |> mutate(group = "high_miss"),
  samp_miss |> slice_min(pct_miss, n = 5, with_ties = FALSE) |> mutate(group = "low_miss"))

set.seed(42)
loocv_res <- do.call(rbind, lapply(seq_len(nrow(loocv_samples)), function(i) {
  sid     <- loocv_samples$Col_ID[i]
  col_idx <- which(colnames(mat) == sid)
  obs_rows <- which(!is.na(mat[, col_idx]))

  nrmse_vals <- vapply(seq_len(20L), function(iter) {
    n_mask   <- max(1L, round(length(obs_rows) * 0.10))
    mask_idx <- sample(obs_rows, n_mask)
    mat_masked <- mat
    mat_masked[mask_idx, col_idx] <- NA
    mat_re <- run_impute(cfg$methods[[best]], mat_masked, randna)
    truth   <- mat[mask_idx, col_idx]
    imputed <- mat_re[mask_idx, col_idx]
    sqrt(mean((truth - imputed)^2)) / sd(truth)
  }, numeric(1))

  tibble(Col_ID = sid, Group_Time = loocv_samples$Group_Time[i],
         pct_miss = loocv_samples$pct_miss[i],
         mean_nrmse = mean(nrmse_vals), sd_nrmse = sd(nrmse_vals),
         group = loocv_samples$group[i])
}))

write_csv(loocv_res, file.path(cfg$data_dir, "11_per_sample_loocv.csv"))
cat(sprintf("  LOOCV saved: %s/11_per_sample_loocv.csv\n", cfg$data_dir))

# --- 8. Save intermediates for report script ---------------------------------
saveRDS(list(
  mat          = mat,
  mat_imp      = mat_imp,
  was_na       = was_na,
  ann          = ann,
  meta         = meta,
  miss_class   = miss_class,
  miss_by_group = miss_by_group,
  prot_pct     = prot_pct,
  pct_miss     = pct_miss,
  mnar_ids     = mnar_ids,
  mnar_audit   = mnar_audit,
  bench_sum    = bench_sum,
  bench_df     = bench_df,
  bin_sum      = bin_sum,
  ext_sum      = ext_sum,
  top5         = top5,
  best         = best,
  METHOD_TYPE  = METHOD_TYPE,
  n_iter       = cfg$n_iter,
  mask_frac    = cfg$mask_frac,
  n_mar_prots  = n_mar_prots,
  n_mnar_prots = n_mnar_prots,
  mar_miss_vals  = mar_miss_vals,
  mnar_miss_vals = mnar_miss_vals,
  total_miss_vals = total_miss_vals,
  loocv_res    = loocv_res
), file.path(cfg$data_dir, "00_report_intermediates.rds"))

cat(sprintf("Done: %s (NRMSE %.4f) | %d unreliable proteins\n",
            best, bench_sum$mean_nrmse[1], sum(!miss_class$imputation_reliable)))
