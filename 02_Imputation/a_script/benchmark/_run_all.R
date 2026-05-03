# _run_all.R — Run full benchmark
# Usage: Rscript 02_Imputation/a_script/benchmark/_run_all.R
# Or with cache: SKIP_IMPUTE=1 Rscript 02_Imputation/a_script/benchmark/_run_all.R

source("02_Imputation/a_script/benchmark/_common.R")

SKIP_IMPUTE <- nzchar(Sys.getenv("SKIP_IMPUTE")) && file.exists(CACHE_RDS)

if (SKIP_IMPUTE) {
  cat("Loading cached imputed matrices...\n")
  imp_list <- readRDS(CACHE_RDS)
} else {
  # Run imputation phase
  source("02_Imputation/a_script/benchmark/_run_impute_only.R")
  imp_list <- readRDS(CACHE_RDS)
}

cat(sprintf("\n=== Imputation complete: %d methods ===\n", length(imp_list)))
for (nm in names(imp_list)) {
  m <- imp_list[[nm]]
  cat(sprintf("  %-25s  %dx%d  NAs=%d\n", nm, nrow(m), ncol(m), sum(is.na(m))))
}

# Run comparisons
cat("\n--- Reconstruction metrics ---\n")
source("02_Imputation/a_script/benchmark/compare/01_reconstruct.R")

cat("\n--- Downstream metrics ---\n")
source("02_Imputation/a_script/benchmark/compare/02_downstream.R")

cat("\n--- Stability metrics ---\n")
source("02_Imputation/a_script/benchmark/compare/03_stability.R")

cat("\n--- Composite ranking ---\n")
source("02_Imputation/a_script/benchmark/compare/04_composite.R")

cat("\nBenchmark complete. Results in: ", BENCH_DIR, "\n")
