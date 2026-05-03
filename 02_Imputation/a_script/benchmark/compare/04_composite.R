# compare/04_composite.R — composite ranking from all metrics
# Reads: 01_reconstruction.csv, 02_downstream.csv, 03_stability.csv
# Writes: 04_composite_ranking.csv, 04_full_report.txt

if (!exists("BENCH_DIR")) {
  source("02_Imputation/a_script/benchmark/_common.R")
}

# --- Load all metric tables ---
recon <- read.csv(file.path(BENCH_DIR, "01_reconstruction.csv"), stringsAsFactors = FALSE)
down  <- read.csv(file.path(BENCH_DIR, "02_downstream.csv"), stringsAsFactors = FALSE)
stab  <- read.csv(file.path(BENCH_DIR, "03_stability.csv"), stringsAsFactors = FALSE)

# --- Merge ---
df <- merge(down, stab, by = "method", all = TRUE)
df <- merge(df, recon, by = "method", all = TRUE)

# --- Classifier sensitivity (hybrid methods only) ---
# For each hybrid base method, compute |FC_rho_km - FC_rho_logistic|
hybrid_bases <- c("BPCA_QRILC", "KNN_QRILC", "imp4p_mixed", "msImpute_v2_mnar")
df$class_sens_fc <- NA_real_
df$class_sens_dep <- NA_real_

for (hb in hybrid_bases) {
  km_row <- which(df$method == paste0(hb, "_km"))
  log_row <- which(df$method == paste0(hb, "_logistic"))
  if (length(km_row) == 1 && length(log_row) == 1) {
    fc_diff <- abs(df$fc_rho[km_row] - df$fc_rho[log_row])
    dep_diff <- abs(df$dep_count[km_row] - df$dep_count[log_row])
    df$class_sens_fc[km_row] <- fc_diff
    df$class_sens_fc[log_row] <- fc_diff
    df$class_sens_dep[km_row] <- dep_diff
    df$class_sens_dep[log_row] <- dep_diff
  }
}

# --- Min-max normalize each metric to [0,1] ---
# Higher = better for all; invert where lower is better
minmax <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (r[2] == r[1]) return(rep(0.5, length(x)))
  (x - r[1]) / (r[2] - r[1])
}

# Invert: NRMSE-MCAR, NRMSE-MNAR, Q1/Q4, KS (lower is better -> invert)
df$norm_nrmse_mcar <- 1 - minmax(ifelse(is.na(df$nrmse_mcar), max(df$nrmse_mcar, na.rm = TRUE), df$nrmse_mcar))
df$norm_nrmse_mnar <- 1 - minmax(ifelse(is.na(df$nrmse_mnar), max(df$nrmse_mnar, na.rm = TRUE), df$nrmse_mnar))
df$norm_procrustes <- minmax(ifelse(is.na(df$procrustes_m2), min(df$procrustes_m2, na.rm = TRUE), df$procrustes_m2))
df$norm_fc_rho     <- minmax(df$fc_rho)
df$norm_nes_rho    <- minmax(df$nes_rho)

# Q1/Q4: ideal is 1.0 (no enrichment), penalize deviations
df$norm_q1q4 <- 1 - minmax(abs(ifelse(is.na(df$q1_q4_ratio), 1, df$q1_q4_ratio) - 1))
df$norm_ks   <- 1 - minmax(ifelse(is.na(df$ks_median), 0, df$ks_median))

# Jackknife (higher = better, NA -> 0)
df$norm_jackknife <- minmax(ifelse(is.na(df$jackknife_rho), min(df$jackknife_rho, na.rm = TRUE), df$jackknife_rho))

# DEP count: distance from Non-imputed median (closer = better)
non_imp_dep <- df$dep_count[df$method == "Non_imputed"]
if (length(non_imp_dep) == 1) {
  df$norm_dep <- 1 - minmax(abs(df$dep_count - non_imp_dep))
} else {
  df$norm_dep <- minmax(df$dep_count)
}

# --- Composite score ---
# Weights: reconstruction 30% (NRMSE only; Procrustes broken), downstream 25%,
# artifact 20%, discovery 15%, stability/jackknife 10% (only 5/20 have values).
# Procrustes weight = 0 (raw SS, not M²; unbounded negative values).
# Jackknife weight = 0 (13/20 methods NA; biases toward the 7 that have values).
df$composite <- (
  0.20 * df$norm_nrmse_mcar +
  0.20 * df$norm_nrmse_mnar +
  0.00 * df$norm_procrustes +
  0.15 * df$norm_fc_rho +
  0.15 * df$norm_nes_rho +
  0.10 * df$norm_q1q4 +
  0.05 * df$norm_ks +
  0.00 * df$norm_jackknife +
  0.15 * df$norm_dep
)

df$rank <- rank(-df$composite, ties.method = "first")
df <- df[order(df$rank), ]

# --- Add method class ---
df$class <- NA_character_
for (i in seq_len(nrow(df))) {
  m <- df$method[i]
  # Strip classifier suffix for class lookup
  base <- sub("_(km|logistic)$", "", m)
  if (base %in% names(METHOD_CLASS)) {
    df$class[i] <- METHOD_CLASS[base]
  } else {
    df$class[i] <- "Unknown"
  }
}

# --- Write ranking CSV ---
out_cols <- c("rank", "method", "class", "composite",
              "nrmse_mcar", "nrmse_mnar", "procrustes_m2",
              "fc_rho", "nes_rho", "dep_count",
              "q1_q4_ratio", "ks_median", "jackknife_rho",
              "class_sens_fc", "class_sens_dep")
out_cols <- intersect(out_cols, names(df))
write.csv(df[, out_cols], file.path(BENCH_DIR, "04_composite_ranking.csv"), row.names = FALSE)

# --- Full report to stdout and file ---
report_lines <- character(0)
add <- function(...) report_lines <<- c(report_lines, sprintf(...))

add("=== PER-METHOD PROFILES ===\n")
for (i in seq_len(nrow(df))) {
  r <- df[i, ]
  add("%s [%s]", r$method, r$class)
  add("  Reconstruction:  NRMSE-MCAR=%.3f  NRMSE-MNAR=%.3f  Procrustes=%.3f",
      ifelse(is.na(r$nrmse_mcar), NA, r$nrmse_mcar),
      ifelse(is.na(r$nrmse_mnar), NA, r$nrmse_mnar),
      ifelse(is.na(r$procrustes_m2), NA, r$procrustes_m2))
  add("  Downstream:      FC_rho=%.3f  NES_rho=%.3f  DEP_count=%d",
      r$fc_rho, r$nes_rho, r$dep_count)
  add("  Artifacts:       Q1_Q4=%.2f  KS=%.3f",
      ifelse(is.na(r$q1_q4_ratio), NA, r$q1_q4_ratio),
      ifelse(is.na(r$ks_median), NA, r$ks_median))
  add("  Stability:       Jackknife=%.3f",
      ifelse(is.na(r$jackknife_rho), NA, r$jackknife_rho))
  if (!is.na(r$class_sens_fc)) {
    add("  Classifier:      FC_rho_diff=%.3f  DEP_diff=%d  [%s]",
        r$class_sens_fc, r$class_sens_dep,
        ifelse(r$class_sens_fc < 0.05, "STABLE", "SENSITIVE"))
  }
  add("  Composite:       %.3f (rank #%d)", r$composite, r$rank)

  # Flags
  flags <- character(0)
  if (!is.na(r$q1_q4_ratio) && r$q1_q4_ratio > 1.5) flags <- c(flags, "Q1/Q4>1.5")
  if (!is.na(r$nrmse_mcar) && r$nrmse_mcar > 0.8) flags <- c(flags, "NRMSE-MCAR>0.8")
  if (!is.na(r$nrmse_mnar) && r$nrmse_mnar > 0.8) flags <- c(flags, "NRMSE-MNAR>0.8")
  if (!is.na(r$class_sens_fc) && r$class_sens_fc > 0.05) flags <- c(flags, "Classifier-sensitive")
  if (length(flags) > 0) add("  Flags:           %s", paste(flags, collapse = ", "))
  add("")
}

# Unbiased metrics only
recon_sorted <- df[!is.na(df$nrmse_mcar), ]
recon_sorted <- recon_sorted[order(recon_sorted$nrmse_mcar), ]
add("=== UNBIASED METRICS ONLY (reconstruction) ===")
for (i in seq_len(nrow(recon_sorted))) {
  r <- recon_sorted[i, ]
  add("  #%d  %-25s  NRMSE-MCAR=%.3f  NRMSE-MNAR=%.3f  Procrustes=%.3f",
      i, r$method, r$nrmse_mcar, r$nrmse_mnar, r$procrustes_m2)
}

add("\n=== DOWNSTREAM AGREEMENT (reference-dependent, interpret with caution) ===")
down_sorted <- df[order(-df$fc_rho), ]
for (i in seq_len(nrow(down_sorted))) {
  r <- down_sorted[i, ]
  add("  #%d  %-25s  FC_rho=%.3f  NES_rho=%.3f", i, r$method, r$fc_rho, r$nes_rho)
}

add("\n=== CLASSIFIER SENSITIVITY (hybrid methods only) ===")
for (hb in hybrid_bases) {
  km_row <- df[df$method == paste0(hb, "_km"), ]
  if (nrow(km_row) == 1 && !is.na(km_row$class_sens_fc)) {
    add("  %-25s  FC_rho_diff=%.3f  DEP_diff=%d  [%s]",
        hb, km_row$class_sens_fc, km_row$class_sens_dep,
        ifelse(km_row$class_sens_fc < 0.05, "STABLE", "SENSITIVE"))
  }
}

add("\n=== COMPOSITE RANKING ===")
for (i in seq_len(nrow(df))) {
  add("  #%d  %-25s  score=%.3f", df$rank[i], df$method[i], df$composite[i])
}

report_text <- paste(report_lines, collapse = "\n")
cat(report_text, "\n")
writeLines(report_text, file.path(BENCH_DIR, "04_full_report.txt"))

cat(sprintf("\nWrote %d rows to 04_composite_ranking.csv\n", nrow(df)))
cat(sprintf("Wrote full report to 04_full_report.txt\n"))
