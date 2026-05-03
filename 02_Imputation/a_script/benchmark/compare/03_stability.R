# compare/03_stability.R — Q1/Q4 ratio, KS statistic, jackknife retention
# Reads imp_list from parent environment (or CACHE_RDS)
# Writes: 02_Imputation/c_data/benchmark/03_stability.csv

if (!exists("imp_list")) {
  source("02_Imputation/a_script/benchmark/_common.R")
  imp_list <- readRDS(CACHE_RDS)
}

suppressPackageStartupMessages({
  library(limma)
})
select <- dplyr::select

# --- Helper: run limma Cancer_vs_Healthy contrast ---
run_limma_benchmark_quick <- function(mat, meta_df) {
  meta_df$group <- factor(meta_df$Group_Time,
    levels = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1"))

  design <- model.matrix(~ 0 + group, data = meta_df)
  colnames(design) <- gsub("^group", "", colnames(design))

  dupcor <- duplicateCorrelation(mat, design, block = meta_df$Subject_ID)
  fit <- lmFit(mat, design, block = meta_df$Subject_ID,
               correlation = dupcor$consensus.correlation)
  cm <- makeContrasts(
    Cancer_vs_Healthy = (CRE_T1 + PLA_T1)/2 - H_T1,
    levels = design)
  fit2 <- contrasts.fit(fit, cm)
  fit2 <- eBayes(fit2)
  topTable(fit2, coef = "Cancer_vs_Healthy", number = Inf, sort.by = "none")
}

# Truly fast methods eligible for jackknife (N = ncol iterations each)
# Only instant methods; BPCA/KNN/imputePCA/msImpute take too long per iteration
FAST_METHODS <- c("MinProb", "Half_minimum", "Non_imputed",
                   "KNN_QRILC_km", "KNN_QRILC_logistic")

meta_df <- as.data.frame(meta)
rownames(meta_df) <- meta_df$Col_ID

# --- Compute metrics per method ---
stability_rows <- list()

for (mname in names(imp_list)) {
  cat(sprintf("  Stability: %s... ", mname))
  imp_mat <- imp_list[[mname]]

  # --- Q1/Q4 ratio ---
  # Classify proteins by median observed intensity quartile
  obs_median <- apply(norm_mat, 1, median, na.rm = TRUE)
  quartiles <- cut(obs_median, breaks = quantile(obs_median, na.rm = TRUE),
                   labels = c("Q1", "Q2", "Q3", "Q4"), include.lowest = TRUE)
  names(quartiles) <- rownames(norm_mat)

  # Run limma for DEP identification
  tt <- tryCatch(run_limma_benchmark_quick(imp_mat, meta_df),
                  error = function(e) NULL)
  q1q4 <- NA_real_
  if (!is.null(tt)) {
    deps <- rownames(tt)[tt$adj.P.Val < 0.05]
    n_q1 <- sum(quartiles[deps] == "Q1", na.rm = TRUE)
    n_q4 <- sum(quartiles[deps] == "Q4", na.rm = TRUE)
    q1q4 <- if (n_q4 > 0) n_q1 / n_q4 else NA_real_
  }

  # --- KS statistic ---
  # Per protein: KS test of imputed values vs observed values
  ks_stats <- numeric(0)
  has_na_rows <- which(rowSums(is.na(norm_mat)) > 0 & rowSums(is.na(norm_mat)) < ncol(norm_mat))

  if (mname != "Non_imputed" && length(has_na_rows) > 0) {
    for (i in has_na_rows) {
      obs <- norm_mat[i, !is.na(norm_mat[i, ])]
      imp_vals <- imp_mat[i, is.na(norm_mat[i, ])]
      if (length(obs) >= 3 && length(imp_vals) >= 2) {
        ks_stats <- c(ks_stats, ks.test(imp_vals, obs)$statistic)
      }
    }
  }
  ks_median <- if (length(ks_stats) > 0) median(ks_stats) else NA_real_

  # --- Jackknife retention (fast methods only) ---
  jackknife_rho <- NA_real_
  if (mname %in% FAST_METHODS && !is.null(tt)) {
    full_logfc <- setNames(tt$logFC, rownames(tt))

    info <- find_method_fn(mname)
    if (!is.null(info)) {
      fn <- get(info$fn_name)
      is_mnar_use <- CLASSIFIERS[[info$classifier]]

      loo_rhos <- numeric(ncol(norm_mat))
      for (s in seq_len(ncol(norm_mat))) {
        loo_mat <- norm_mat[, -s, drop = FALSE]
        loo_meta <- meta_df[-s, , drop = FALSE]

        loo_imp <- tryCatch(
          fn(mat = loo_mat, meta = loo_meta, is_mnar = is_mnar_use),
          error = function(e) NULL)
        if (is.null(loo_imp)) { loo_rhos[s] <- NA; next }

        loo_tt <- tryCatch(run_limma_benchmark_quick(loo_imp, loo_meta),
                            error = function(e) NULL)
        if (is.null(loo_tt)) { loo_rhos[s] <- NA; next }

        loo_logfc <- setNames(loo_tt$logFC, rownames(loo_tt))
        shared <- intersect(names(full_logfc), names(loo_logfc))
        loo_rhos[s] <- cor(full_logfc[shared], loo_logfc[shared],
                           method = "spearman", use = "complete.obs")
      }
      jackknife_rho <- median(loo_rhos, na.rm = TRUE)
    }
  }

  stability_rows[[mname]] <- data.frame(
    method = mname,
    q1_q4_ratio = q1q4,
    ks_median = ks_median,
    jackknife_rho = jackknife_rho,
    stringsAsFactors = FALSE
  )
  cat(sprintf("Q1/Q4=%.2f KS=%.3f Jack=%.3f\n",
              ifelse(is.na(q1q4), 0, q1q4),
              ifelse(is.na(ks_median), 0, ks_median),
              ifelse(is.na(jackknife_rho), 0, jackknife_rho)))
}

stability_df <- do.call(rbind, stability_rows)
rownames(stability_df) <- NULL
write.csv(stability_df, file.path(BENCH_DIR, "03_stability.csv"), row.names = FALSE)
cat(sprintf("Wrote %d rows to 03_stability.csv\n", nrow(stability_df)))
