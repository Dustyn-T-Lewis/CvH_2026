# _run_new_methods.R — Add new methods 17-23 to existing cache
# Usage: Rscript 02_Imputation/a_script/benchmark/_run_new_methods.R

source("02_Imputation/a_script/benchmark/_common.R")
suppressPackageStartupMessages({
  library(MsCoreUtils)
  library(pcaMethods)
})
select <- dplyr::select

imp_list <- readRDS(CACHE_RDS)
cat(sprintf("Loaded %d existing methods\n", length(imp_list)))

new_files <- list.files("02_Imputation/a_script/benchmark/methods",
                        pattern = "^(1[789]|2[0123])_.*\\.R$", full.names = TRUE)
new_files <- sort(new_files)
cat(sprintf("New method files: %s\n", paste(basename(new_files), collapse = ", ")))

for (f in new_files) {
  env_before <- ls(envir = .GlobalEnv)
  source(f, local = FALSE)
  env_after <- ls(envir = .GlobalEnv)
  new_fns <- grep("^impute_", setdiff(env_after, env_before), value = TRUE)
  if (length(new_fns) == 0) {
    raw_name <- sub("^\\d+_", "", tools::file_path_sans_ext(basename(f)))
    all_fns <- grep("^impute_", ls(envir = .GlobalEnv), value = TRUE)
    new_fns <- all_fns[tolower(gsub("_", "", sub("^impute_", "", all_fns))) ==
                        tolower(gsub("_", "", raw_name))]
  }
  fn_name <- new_fns[1]
  canon_name <- sub("^impute_", "", fn_name)

  if (canon_name %in% names(imp_list)) {
    cat(sprintf("  SKIP %s (already cached)\n", canon_name))
    next
  }

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

saveRDS(imp_list, CACHE_RDS)
cat(sprintf("\nUpdated cache: %d entries\n", length(imp_list)))
for (nm in names(imp_list)) {
  m <- imp_list[[nm]]
  cat(sprintf("  %-30s %dx%d NAs=%d\n", nm, nrow(m), ncol(m), sum(is.na(m))))
}
