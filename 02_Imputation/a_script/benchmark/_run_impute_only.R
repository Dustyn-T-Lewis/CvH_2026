# _run_impute_only.R — Run only the imputation phase (no comparison scripts)
# Usage: Rscript 02_Imputation/a_script/benchmark/_run_impute_only.R

source("02_Imputation/a_script/benchmark/_common.R")

suppressPackageStartupMessages({
  library(MsCoreUtils)
  library(impute)
  library(missForest)
  library(missMDA)
  library(msImpute)
})
select <- dplyr::select

method_files <- sort(list.files("02_Imputation/a_script/benchmark/methods",
                                 pattern = "\\.R$", full.names = TRUE))
imp_list <- list()

for (f in method_files) {
  # Source the file and find the impute_* function it defines
  env_before <- ls(envir = .GlobalEnv)
  source(f, local = FALSE)
  env_after <- ls(envir = .GlobalEnv)
  new_fns <- setdiff(env_after, env_before)
  fn_name <- grep("^impute_", new_fns, value = TRUE)[1]
  if (is.na(fn_name)) {
    cat(sprintf("SKIP %s: no impute_* function found\n", basename(f)))
    next
  }

  # Derive canon_name from function name: impute_BPCA_QRILC -> BPCA_QRILC
  canon_name <- sub("^impute_", "", fn_name)
  base_lower <- tolower(gsub("_", "", canon_name))
  is_hybrid <- any(sapply(HYBRID_METHODS, function(h) tolower(gsub("_", "", h)) == base_lower))

  if (is_hybrid) {
    for (cl_name in names(CLASSIFIERS)) {
      key <- paste0(canon_name, "_", cl_name)
      cat(sprintf("[%s] Running %s... ", format(Sys.time(), "%H:%M:%S"), key))
      t0 <- proc.time()["elapsed"]
      result <- tryCatch(
        get(fn_name)(mat = norm_mat, meta = meta, is_mnar = CLASSIFIERS[[cl_name]]),
        error = function(e) { cat(sprintf("FAILED: %s\n", e$message)); NULL }
      )
      dt <- round(proc.time()["elapsed"] - t0, 1)
      if (!is.null(result)) {
        imp_list[[key]] <- result
        cat(sprintf("OK (%.1fs, NAs=%d)\n", dt, sum(is.na(result))))
      }
    }
  } else {
    cat(sprintf("[%s] Running %s... ", format(Sys.time(), "%H:%M:%S"), canon_name))
    t0 <- proc.time()["elapsed"]
    result <- tryCatch(
      get(fn_name)(mat = norm_mat, meta = meta, is_mnar = CLASSIFIERS$km),
      error = function(e) { cat(sprintf("FAILED: %s\n", e$message)); NULL }
    )
    dt <- round(proc.time()["elapsed"] - t0, 1)
    if (!is.null(result)) {
      imp_list[[canon_name]] <- result
      cat(sprintf("OK (%.1fs, NAs=%d)\n", dt, sum(is.na(result))))
    }
  }
}

saveRDS(imp_list, CACHE_RDS)
cat(sprintf("\nCached %d entries to %s\n", length(imp_list), CACHE_RDS))
for (nm in names(imp_list)) {
  m <- imp_list[[nm]]
  cat(sprintf("  %-30s %dx%d NAs=%d\n", nm, nrow(m), ncol(m), sum(is.na(m))))
}
