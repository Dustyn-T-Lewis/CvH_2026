# compare/01_reconstruct.R — NRMSE + Procrustes (artificial masking)
# Reads imp_list from parent environment (or CACHE_RDS)
# Writes: 02_Imputation/c_data/benchmark/01_reconstruction.csv

if (!exists("imp_list")) {
  source("02_Imputation/a_script/benchmark/_common.R")
  imp_list <- readRDS(CACHE_RDS)
}

library(vegan)
select <- dplyr::select

N_REPS <- 5  # reduced from 20 for practical runtime; each rep requires full re-imputation
MASK_FRAC <- 0.05
set.seed(42)

nrmse <- function(true_vals, imp_vals) {
  sqrt(mean((true_vals - imp_vals)^2)) / sd(true_vals)
}

# Methods to evaluate (exclude Non_imputed — can't reconstruct without a method)
eval_methods <- setdiff(names(imp_list), "Non_imputed")

# Pre-compute observed indices and Q1 threshold for MNAR masking
obs_idx <- which(!is.na(norm_mat))
obs_vals_all <- norm_mat[obs_idx]
q1_threshold <- quantile(obs_vals_all, 0.25)
q1_idx <- obs_idx[obs_vals_all <= q1_threshold]

results <- vector("list", length(eval_methods))
names(results) <- eval_methods

for (mname in eval_methods) {
  cat(sprintf("  Reconstruct: %s ", mname))

  info <- find_method_fn(mname)
  if (is.null(info)) {
    cat("SKIP (no method file)\n")
    next
  }
  fn <- get(info$fn_name)
  is_mnar_use <- CLASSIFIERS[[info$classifier]]

  nrmse_mcar <- numeric(N_REPS)
  nrmse_mnar <- numeric(N_REPS)
  procrustes_m2 <- numeric(N_REPS)

  for (rep in seq_len(N_REPS)) {
    # MCAR mask: uniform random 5% of observed values
    mcar_mask <- sample(obs_idx, round(length(obs_idx) * MASK_FRAC))
    mcar_true <- norm_mat[mcar_mask]
    mm_mcar <- norm_mat
    mm_mcar[mcar_mask] <- NA

    # MNAR mask: 5% from Q1 (bottom 25% intensity)
    n_mnar_mask <- round(length(obs_idx) * MASK_FRAC)
    mnar_mask <- sample(q1_idx, min(n_mnar_mask, length(q1_idx)))
    mnar_true <- norm_mat[mnar_mask]
    mm_mnar <- norm_mat
    mm_mnar[mnar_mask] <- NA

    # Re-impute masked data
    imp_mcar <- tryCatch(
      fn(mat = mm_mcar, meta = meta, is_mnar = is_mnar_use),
      error = function(e) NULL)
    imp_mnar <- tryCatch(
      fn(mat = mm_mnar, meta = meta, is_mnar = is_mnar_use),
      error = function(e) NULL)

    if (!is.null(imp_mcar)) {
      nrmse_mcar[rep] <- nrmse(mcar_true, imp_mcar[mcar_mask])
    } else {
      nrmse_mcar[rep] <- NA
    }

    if (!is.null(imp_mnar)) {
      nrmse_mnar[rep] <- nrmse(mnar_true, imp_mnar[mnar_mask])
    } else {
      nrmse_mnar[rep] <- NA
    }

    # Procrustes M² (on MCAR-imputed matrix vs original complete proteins)
    # Use procrustes() directly (no permutation test for speed)
    if (!is.null(imp_mcar)) {
      complete_prots <- complete.cases(norm_mat)
      if (sum(complete_prots) > 10) {
        n_pc <- min(5, ncol(norm_mat) - 1)
        pca_true <- prcomp(t(norm_mat[complete_prots, ]), center = TRUE, scale. = TRUE)
        pca_imp  <- prcomp(t(imp_mcar[complete_prots, ]), center = TRUE, scale. = TRUE)
        proc <- vegan::procrustes(pca_true$x[, 1:n_pc], pca_imp$x[, 1:n_pc])
        procrustes_m2[rep] <- 1 - proc$ss  # M² correlation (higher = better)
      } else {
        procrustes_m2[rep] <- NA
      }
    } else {
      procrustes_m2[rep] <- NA
    }
  }
  cat(sprintf("NRMSE-MCAR=%.3f NRMSE-MNAR=%.3f Procrustes=%.3f\n",
              median(nrmse_mcar, na.rm = TRUE),
              median(nrmse_mnar, na.rm = TRUE),
              median(procrustes_m2, na.rm = TRUE)))

  results[[mname]] <- data.frame(
    method = mname,
    nrmse_mcar = median(nrmse_mcar, na.rm = TRUE),
    nrmse_mcar_sd = sd(nrmse_mcar, na.rm = TRUE),
    nrmse_mnar = median(nrmse_mnar, na.rm = TRUE),
    nrmse_mnar_sd = sd(nrmse_mnar, na.rm = TRUE),
    procrustes_m2 = median(procrustes_m2, na.rm = TRUE),
    procrustes_m2_sd = sd(procrustes_m2, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

recon_df <- do.call(rbind, results)
rownames(recon_df) <- NULL
write.csv(recon_df, file.path(BENCH_DIR, "01_reconstruction.csv"), row.names = FALSE)
cat(sprintf("Wrote %d rows to 01_reconstruction.csv\n", nrow(recon_df)))
