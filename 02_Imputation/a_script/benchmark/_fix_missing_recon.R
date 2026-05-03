# _fix_missing_recon.R — Reconstruct the 3 methods that were skipped
source("02_Imputation/a_script/benchmark/_common.R")
imp_list <- readRDS(CACHE_RDS)

suppressPackageStartupMessages({
  library(MsCoreUtils)
  library(vegan)
})
select <- dplyr::select

N_REPS <- 5
MASK_FRAC <- 0.05
set.seed(42)

nrmse <- function(true_vals, imp_vals) {
  sqrt(mean((true_vals - imp_vals)^2)) / sd(true_vals)
}

obs_idx <- which(!is.na(norm_mat))
obs_vals_all <- norm_mat[obs_idx]
q1_threshold <- quantile(obs_vals_all, 0.25)
q1_idx <- obs_idx[obs_vals_all <= q1_threshold]

missing_methods <- c("knn_standalone", "bpca_standalone", "RF_MsCoreUtils")
recon <- read.csv(file.path(BENCH_DIR, "01_reconstruction.csv"), stringsAsFactors = FALSE)

for (mname in missing_methods) {
  cat(sprintf("Reconstructing %s... ", mname))
  info <- find_method_fn(mname)
  if (is.null(info)) {
    cat("SKIP (fn not found)\n")
    next
  }
  fn <- get(info$fn_name)
  is_mnar_use <- CLASSIFIERS[[info$classifier]]

  nrmse_mcar <- nrmse_mnar <- numeric(N_REPS)
  for (rep in 1:N_REPS) {
    mcar_mask <- sample(obs_idx, round(length(obs_idx) * MASK_FRAC))
    mm <- norm_mat; mm[mcar_mask] <- NA
    imp_mcar <- tryCatch(fn(mat = mm, meta = meta, is_mnar = is_mnar_use),
                          error = function(e) NULL)
    if (!is.null(imp_mcar)) {
      nrmse_mcar[rep] <- nrmse(norm_mat[mcar_mask], imp_mcar[mcar_mask])
    } else {
      nrmse_mcar[rep] <- NA
    }

    n_mask <- round(length(obs_idx) * MASK_FRAC)
    mnar_mask <- sample(q1_idx, min(n_mask, length(q1_idx)))
    mm2 <- norm_mat; mm2[mnar_mask] <- NA
    imp_mnar <- tryCatch(fn(mat = mm2, meta = meta, is_mnar = is_mnar_use),
                          error = function(e) NULL)
    if (!is.null(imp_mnar)) {
      nrmse_mnar[rep] <- nrmse(norm_mat[mnar_mask], imp_mnar[mnar_mask])
    } else {
      nrmse_mnar[rep] <- NA
    }
  }

  med_mcar <- median(nrmse_mcar, na.rm = TRUE)
  med_mnar <- median(nrmse_mnar, na.rm = TRUE)
  cat(sprintf("NRMSE-MCAR=%.3f NRMSE-MNAR=%.3f\n", med_mcar, med_mnar))

  # Update or append to recon
  if (mname %in% recon$method) {
    idx <- which(recon$method == mname)
    recon$nrmse_mcar[idx] <- med_mcar
    recon$nrmse_mcar_sd[idx] <- sd(nrmse_mcar, na.rm = TRUE)
    recon$nrmse_mnar[idx] <- med_mnar
    recon$nrmse_mnar_sd[idx] <- sd(nrmse_mnar, na.rm = TRUE)
  } else {
    recon <- rbind(recon, data.frame(
      method = mname,
      nrmse_mcar = med_mcar, nrmse_mcar_sd = sd(nrmse_mcar, na.rm = TRUE),
      nrmse_mnar = med_mnar, nrmse_mnar_sd = sd(nrmse_mnar, na.rm = TRUE),
      procrustes_m2 = NA, procrustes_m2_sd = NA,
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(recon, file.path(BENCH_DIR, "01_reconstruction.csv"), row.names = FALSE)
cat(sprintf("Updated 01_reconstruction.csv: %d rows\n", nrow(recon)))

# Now re-run composite
cat("\n--- Regenerating composite ---\n")
source("02_Imputation/a_script/benchmark/compare/04_composite.R")
