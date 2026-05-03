# CvH DEP — Robustness: response differential, bootstrap CIs, power, imputation sensitivity
# Refs: Conover 1999 (KS, Fligner), Romano et al. 2006 (Cliff's delta),
#       Efron 1993 (bootstrap), Cohen 1988 (power),
#       Karpievitch 2012 BMC Bioinform 13(S16):S5 (imputation sensitivity)

library(dplyr)
library(readr)
library(proteoDA)
library(openxlsx)
library(boot)
library(pwr)

set.seed(42)
setwd(rprojroot::find_rstudio_root_file())

`%||%` <- function(x, y) if (!is.null(x)) x else y

add_sheet <- function(wb, name, df, title = NULL, notes = NULL) {
  openxlsx::addWorksheet(wb, name)
  start_row <- 1L
  if (!is.null(title)) {
    openxlsx::writeData(wb, name, title, startRow = 1)
    openxlsx::addStyle(wb, name,
      openxlsx::createStyle(textDecoration = "bold", fontSize = 12),
      rows = 1, cols = 1)
    start_row <- start_row + 1L
  }
  if (!is.null(notes)) {
    for (i in seq_along(notes)) {
      openxlsx::writeData(wb, name, notes[i], startRow = start_row)
      openxlsx::addStyle(wb, name,
        openxlsx::createStyle(fontSize = 10, fontColour = "#555555", wrapText = TRUE),
        rows = start_row, cols = 1)
      start_row <- start_row + 1L
    }
    start_row <- start_row + 1L
  }
  hs <- openxlsx::createStyle(textDecoration = "bold", border = "Bottom",
                               fgFill = "#DCE6F1")
  openxlsx::writeData(wb, name, df, startRow = start_row, headerStyle = hs)
  openxlsx::freezePane(wb, name, firstActiveRow = start_row + 1L,
                        firstActiveCol = 2)
  openxlsx::setColWidths(wb, name, cols = seq_len(ncol(df)), widths = "auto")
}

cfg <- list(
  data_dir    = "03_DEP/c_data",
  per_dir     = "03_DEP/c_data/04_per_contrast_results",
  imp_path    = "02_Imputation/c_data/01_imputed.csv",
  pi_thresh   = 0.05,
  pval_thresh = 0.10
)

dal_crvh <- readRDS(file.path(cfg$data_dir, "01_limma_DAList_CRvH.rds"))
dal_cr   <- readRDS(file.path(cfg$data_dir, "01_limma_DAList_CR.rds"))

contrast_names_m1 <- names(dal_crvh$results)
contrast_names_m2 <- names(dal_cr$results)
all_contrast_names <- c(contrast_names_m1, contrast_names_m2)

required_contrasts <- c(
  "Cancer_vs_Healthy", "Training_CR",
  "Baseline_Supplement", "Training_CRE", "Training_PLA", "Supplement_Interaction"
)
missing_contrasts <- setdiff(required_contrasts, all_contrast_names)
if (length(missing_contrasts) > 0) {
  stop("Robustness inputs are missing required contrasts: ",
       paste(missing_contrasts, collapse = ", "))
}

results_list <- setNames(
  lapply(all_contrast_names, function(cname)
    read_csv(file.path(cfg$per_dir, paste0(cname, ".csv")), show_col_types = FALSE)),
  all_contrast_names)

meta_m1 <- as.data.frame(dal_crvh$metadata)
meta_m2 <- as.data.frame(dal_cr$metadata)

# --- 1. RESPONSE DIFFERENTIAL DIAGNOSTICS ---
# Tests whether |logFC| distributions differ between Training_CRE and
# Training_PLA, providing evidence for supplement-dependent response modulation.
# Non-parametric complement to the Supplement_Interaction contrast.
# Refs: Conover 1999 (KS, Fligner), Romano et al. 2006 (Cliff's delta)

diff_df <- tibble(
  gene        = results_list[["Training_CRE"]]$gene,
  abs_lfc_cre = abs(results_list[["Training_CRE"]]$logFC),
  abs_lfc_pla = abs(results_list[["Training_PLA"]]$logFC)
) |> filter(!is.na(abs_lfc_cre) & !is.na(abs_lfc_pla))

ks_res <- ks.test(diff_df$abs_lfc_cre, diff_df$abs_lfc_pla)

fk_data <- data.frame(
  abs_lfc  = c(diff_df$abs_lfc_cre, diff_df$abs_lfc_pla),
  contrast = factor(rep(c("Training_CRE", "Training_PLA"), each = nrow(diff_df)))
)
fk_res <- fligner.test(abs_lfc ~ contrast, data = fk_data)

wx_res <- wilcox.test(diff_df$abs_lfc_cre, diff_df$abs_lfc_pla, paired = TRUE)

diffs <- diff_df$abs_lfc_cre - diff_df$abs_lfc_pla
cliff_delta <- (sum(diffs > 0) - sum(diffs < 0)) / length(diffs)
cliff_mag <- dplyr::case_when(
  abs(cliff_delta) < 0.147 ~ "negligible",
  abs(cliff_delta) < 0.33  ~ "small",
  abs(cliff_delta) < 0.474 ~ "medium",
  TRUE                     ~ "large"
)

quants <- c(0.25, 0.50, 0.75, 0.90)
quant_df <- tibble(
  quantile     = paste0("Q", quants * 100),
  Training_CRE = quantile(diff_df$abs_lfc_cre, quants),
  Training_PLA = quantile(diff_df$abs_lfc_pla, quants),
  ratio_CRE_PLA = round(Training_CRE / Training_PLA, 3)
)

resp_diag <- tibble(
  test = c("Kolmogorov-Smirnov", "Fligner-Killeen",
           "Wilcoxon signed-rank", "Cliff's delta"),
  statistic = c(ks_res$statistic, fk_res$statistic,
                wx_res$statistic, cliff_delta),
  p_value = c(ks_res$p.value, fk_res$p.value, wx_res$p.value, NA),
  interpretation = c(
    ifelse(ks_res$p.value < 0.05,
           "Distributions differ significantly (shape/location)",
           "No significant distributional difference"),
    ifelse(fk_res$p.value < 0.05,
           "Variance differs (consistent with differential response)",
           "No significant variance difference"),
    ifelse(wx_res$p.value < 0.05,
           "Paired |logFC| shift significant (CRE vs PLA response differs)",
           "No significant paired shift"),
    sprintf("%s effect (delta = %.3f; >0 = CRE responds more)",
            tools::toTitleCase(cliff_mag), cliff_delta)
  )
)
write_csv(resp_diag, file.path(cfg$data_dir, "06_response_differential.csv"))

cat(sprintf("  Response differential: KS p=%.2g, Wilcoxon p=%.2g, Cliff d=%.3f (%s)\n",
            ks_res$p.value, wx_res$p.value, cliff_delta, cliff_mag))

# --- 2. BOOTSTRAP CI (Effect Sizes) ---
median_fn <- function(d, i) median(d[i], na.rm = TRUE)

boot_df <- lapply(all_contrast_names, function(cname) {
  vals <- abs(results_list[[cname]]$logFC)
  vals <- vals[!is.na(vals)]
  b  <- boot(vals, median_fn, R = 10000)
  ci <- tryCatch(boot.ci(b, type = "bca"),
                 error = function(e) boot.ci(b, type = "perc"))
  ci_lo <- if (!is.null(ci$bca)) ci$bca[4] else ci$percent[4]
  ci_hi <- if (!is.null(ci$bca)) ci$bca[5] else ci$percent[5]
  tibble(contrast = cname, median_absLFC = median(vals),
         ci_lower = ci_lo, ci_upper = ci_hi,
         boot_se = sd(b$t), n_proteins = length(vals))
}) |> bind_rows()
write_csv(boot_df, file.path(cfg$data_dir, "07_effect_size_bootstrap.csv"))

print(boot_df)

# --- 3. POWER ANALYSIS ---

fit_m1 <- dal_crvh$eBayes_fit
fit_m2 <- dal_cr$eBayes_fit
cor_m1 <- dal_crvh$eBayes_fit$correlation %||%
  dal_crvh$tags$duplicate_correlation %||% NA_real_
cor_m2 <- dal_cr$eBayes_fit$correlation %||%
  dal_cr$tags$duplicate_correlation %||% NA_real_
sigma_m1 <- sqrt(mean(fit_m1$sigma^2, na.rm = TRUE))
sigma_m2 <- sqrt(mean(fit_m2$sigma^2, na.rm = TRUE))

n_cr_t1     <- sum(meta_m1$group %in% c("CRE_T1", "PLA_T1"))
n_h_t1      <- sum(meta_m1$group == "H_T1")
n_cr_paired <- length(unique(meta_m2$subject[meta_m2$timepoint == "T2"]))
n_cre       <- sum(meta_m2$supplement == "CRE" & meta_m2$timepoint == "T1")
n_pla       <- sum(meta_m2$supplement == "PLA" & meta_m2$timepoint == "T1")

power_specs <- list(
  list(contrast = "Cancer_vs_Healthy", n = min(n_cr_t1, n_h_t1),
       sigma = sigma_m1, cor = NA, paired = FALSE, model = "CRvH"),
  list(contrast = "Training_CR", n = n_cr_paired,
       sigma = sigma_m1, cor = cor_m1, paired = TRUE, model = "CRvH"),
  list(contrast = "Training_CRE", n = n_cre,
       sigma = sigma_m2, cor = cor_m2, paired = TRUE, model = "CR"),
  list(contrast = "Training_PLA", n = n_pla,
       sigma = sigma_m2, cor = cor_m2, paired = TRUE, model = "CR"),
  list(contrast = "Supplement_Interaction", n = min(n_cre, n_pla),
       sigma = sigma_m2, cor = cor_m2, paired = TRUE, model = "CR")
)

power_df <- lapply(power_specs, function(spec) {
  eff_sigma <- if (spec$paired && !is.na(spec$cor))
    spec$sigma * sqrt(2 * (1 - spec$cor)) else
    spec$sigma * sqrt(2)
  pw <- pwr.t.test(n = spec$n, d = NULL, sig.level = cfg$pval_thresh,
                   power = 0.80,
                   type = if (spec$paired) "paired" else "two.sample")
  tibble(model = spec$model, contrast = spec$contrast,
         n_subjects = spec$n,
         within_cor = ifelse(spec$paired, spec$cor, NA_real_),
         effective_sigma = round(eff_sigma, 4),
         min_detectable_d = round(pw$d, 4),
         min_detectable_logFC = round(pw$d * eff_sigma, 4),
         power = 0.80, alpha = cfg$pval_thresh)
}) |> bind_rows()
write_csv(power_df, file.path(cfg$data_dir, "08_power_analysis.csv"))
print(as.data.frame(power_df))

# --- 4. IMPUTATION SENSITIVITY ---

sens_df <- NULL

if (file.exists(cfg$imp_path)) {
  ann_cols <- c("uniprot_id", "protein", "gene", "description")
  imp_data <- read_csv(cfg$imp_path, show_col_types = FALSE)
  imp_samp <- setdiff(names(imp_data), ann_cols)
  imp_mat  <- as.matrix(imp_data[, imp_samp])
  rownames(imp_mat) <- imp_data$uniprot_id

  # Align imp_mat row order to each model's annotation BEFORE DAList().
  # Imputation reorders proteins by gene_order for missForest determinism;
  # DAList() silently re-labels rownames(data) to match annotation while
  # leaving values in place, scrambling protein labels relative to intensities
  # if data and annotation row orders differ.
  align_imp_mat <- function(imp_mat, ann_uid) {
    ord <- match(ann_uid, rownames(imp_mat))
    if (any(is.na(ord))) {
      stop("Imputed CSV is missing proteins from the DEP annotation.")
    }
    out <- imp_mat[ord, , drop = FALSE]
    stopifnot(identical(rownames(out), ann_uid))
    out
  }

  # Model 1 imputed
  ann_crvh <- as.data.frame(dal_crvh$annotation)
  imp_mat_m1 <- align_imp_mat(imp_mat, ann_crvh$uniprot_id)
  shared1 <- intersect(colnames(dal_crvh$data), colnames(imp_mat_m1))
  dal_imp1 <- DAList(
    data       = imp_mat_m1[, shared1],
    annotation = ann_crvh,
    metadata   = as.data.frame(dal_crvh$metadata)[
      as.data.frame(dal_crvh$metadata)$sample_id %in% shared1, ],
    tags       = list(norm_method = "cycloess_imputed")
  )
  dal_imp1$design <- dal_crvh$design
  dal_imp1 <- fit_limma_model(dal_imp1)
  dal_imp1 <- extract_DA_results(dal_imp1, pval_thresh = cfg$pval_thresh,
                                  lfc_thresh = 0, adj_method = "BH")

  # Model 2 imputed
  ann_cr <- as.data.frame(dal_cr$annotation)
  imp_mat_m2 <- align_imp_mat(imp_mat, ann_cr$uniprot_id)
  shared2 <- intersect(colnames(dal_cr$data), colnames(imp_mat_m2))
  dal_imp2 <- DAList(
    data       = imp_mat_m2[, shared2],
    annotation = ann_cr,
    metadata   = as.data.frame(dal_cr$metadata)[
      as.data.frame(dal_cr$metadata)$sample_id %in% shared2, ],
    tags       = list(norm_method = "cycloess_imputed")
  )
  dal_imp2$design <- dal_cr$design
  dal_imp2 <- fit_limma_model(dal_imp2)
  dal_imp2 <- extract_DA_results(dal_imp2, pval_thresh = cfg$pval_thresh,
                                  lfc_thresh = 0, adj_method = "BH")

  sens_rows <- list()
  for (cname in all_contrast_names) {
    comb_file <- if (cname %in% contrast_names_m1)
      "03_combined_results_CRvH.csv" else "03_combined_results_CR.csv"
    comb <- read_csv(file.path(cfg$data_dir, comb_file), show_col_types = FALSE)

    t_col   <- paste0("t_", cname)
    adj_col <- paste0("adj.P.Val_", cname)

    dal_imp <- if (cname %in% contrast_names_m1) dal_imp1 else dal_imp2
    if (!(t_col %in% names(comb)) || !(cname %in% names(dal_imp$results))) next

    imp_df <- dal_imp$results[[cname]] |>
      tibble::rownames_to_column("uniprot_id") |>
      select(uniprot_id, t_imp = t, padj_imp = adj.P.Val)

    merged <- inner_join(
      comb |> select(uniprot_id, t_nonimp = all_of(t_col),
                      padj_nonimp = all_of(adj_col)),
      imp_df,
      by = "uniprot_id"
    ) |> filter(!is.na(t_nonimp) & !is.na(t_imp))

    sp <- cor.test(merged$t_nonimp, merged$t_imp, method = "spearman")

    sens_rows[[cname]] <- tibble(contrast = cname,
           spearman_rho = round(sp$estimate, 4),
           p_value = sp$p.value, n_proteins = nrow(merged))
  }
  sens_df <- bind_rows(sens_rows)
  write_csv(sens_df, file.path(cfg$data_dir, "09_imputation_sensitivity.csv"))

  print(as.data.frame(sens_df))
} else {
  cat("  Imputed data not found — skipping\n")
  write_csv(tibble(contrast = character(), spearman_rho = numeric(),
                   p_value = numeric(), n_proteins = integer()),
            file.path(cfg$data_dir, "09_imputation_sensitivity.csv"))
}

# --- 5. SAVE INTERMEDIATES ---

saveRDS(list(
  resp_diag  = resp_diag,
  quant_df   = quant_df,
  boot_df    = boot_df,
  power_df   = power_df,
  sens_df    = sens_df,
  all_contrast_names = all_contrast_names,
  contrast_names_m1 = contrast_names_m1,
  contrast_names_m2 = contrast_names_m2,
  results_list = results_list
), file.path(cfg$data_dir, "00_report_intermediates_robustness.rds"))

# --- 6. SUPPLEMENTARY EXCEL ---
wb_supp <- createWorkbook()

cor_m1_val <- if (!is.na(cor_m1)) sprintf("%.3f", cor_m1) else "N/A"
cor_m2_val <- if (!is.na(cor_m2)) sprintf("%.3f", cor_m2) else "N/A"

overview <- tibble(
  Parameter = c(
    "Model 1 (CRvH)", "Model 1 Design", "Model 1 Contrasts",
    "Model 1 Within-subject cor", "Model 1 N proteins", "Model 1 N samples",
    "Model 2 (CR)", "Model 2 Design", "Model 2 Contrasts",
    "Model 2 Within-subject cor", "Model 2 N proteins", "Model 2 N samples",
    "FDR method", "FDR threshold", "Pi-score threshold",
    "Normalization", "Imputation strategy"),
  Value = c(
    "All subjects (CRE, PLA, H)",
    "~ 0 + group + (1 | Subject_ID)",
    paste(contrast_names_m1, collapse = "; "),
    cor_m1_val, as.character(nrow(dal_crvh$data)), as.character(ncol(dal_crvh$data)),
    "CR subjects only (CRE, PLA)",
    "~ 0 + group + (1 | Subject_ID)",
    paste(contrast_names_m2, collapse = "; "),
    cor_m2_val, as.character(nrow(dal_cr$data)), as.character(ncol(dal_cr$data)),
    "Benjamini-Hochberg",
    "0.10 (exploratory; pi-score provides secondary filter)",
    "Capital Pi < 0.05 (= pi-value > 1.3; Xiao et al. 2014, Bioinformatics 30:801)",
    "Cycloess (selected by PRONE benchmarking)",
    "Non-imputed; limma handles NAs per-protein (Karpievitch et al. 2012)")
)
add_sheet(wb_supp, "Methods_Overview", overview,
  title = "CvH DEP Analysis \u2014 Methods Overview",
  notes = "Parameters and settings for the differential expression analysis.")

da_sum_m1 <- read_csv(file.path(cfg$data_dir, "02_DA_summary_CRvH.csv"),
                       show_col_types = FALSE) |> mutate(model = "CRvH", .before = 1)
da_sum_m2 <- read_csv(file.path(cfg$data_dir, "02_DA_summary_CR.csv"),
                       show_col_types = FALSE) |> mutate(model = "CR", .before = 1)
da_summary <- bind_rows(da_sum_m1, da_sum_m2)

add_sheet(wb_supp, "Significance_Summary", da_summary,
  title = "Significance counts by model, contrast, direction, and criterion.",
  notes = c("sig.PVal: nominal P < pval_thresh | sig.FDR: adj.P < pval_thresh (BH)",
            "sig.Pi: Pi-score < 0.05 | sig.FDR.05/10: FDR at 0.05 and 0.10 thresholds."))

keep_cols <- c("uniprot_id", "gene", "protein", "description",
               "logFC", "CI.L", "CI.R", "t", "P.Value", "adj.P.Val",
               "pi_score", "sig_pi")

contrast_labels <- c(
  Cancer_vs_Healthy       = "(CRE_T1 + PLA_T1)/2 - H_T1",
  Training_CR             = "(CRE_T2 + PLA_T2)/2 - (CRE_T1 + PLA_T1)/2",
  Baseline_Supplement     = "CRE_T1 - PLA_T1",
  Training_CRE            = "CRE_T2 - CRE_T1",
  Training_PLA            = "PLA_T2 - PLA_T1",
  Supplement_Interaction   = "(CRE_T2 - CRE_T1) - (PLA_T2 - PLA_T1)"
)

for (cname in all_contrast_names) {
  res <- results_list[[cname]]
  cols_present <- intersect(keep_cols, names(res))
  n_sig_fdr <- sum(res$adj.P.Val < 0.10, na.rm = TRUE)
  n_sig_pi  <- sum(res$sig_pi != 0, na.rm = TRUE)
  model_label <- if (cname %in% contrast_names_m1) "CRvH" else "CR"
  add_sheet(wb_supp, cname, res[order(res$pi_score), cols_present],
    title = sprintf("[%s] %s: %s", model_label, cname, contrast_labels[[cname]]),
    notes = c(sprintf("Proteins tested: %d | FDR<0.10: %d | Pi<0.05: %d",
                      nrow(res), n_sig_fdr, n_sig_pi),
              "Sorted by Pi-score (most significant first)."))
}

# Response Differential sheet (tests + quantiles, matching YvO Blunting_Diagnostics layout)
addWorksheet(wb_supp, "Response_Differential")
rd_header <- c(
  "Supplement-dependent modulation of training response.",
  "Tests compare |logFC| distributions: Training_CRE vs Training_PLA.",
  "KS: shape/location | Fligner-Killeen: variance | Wilcoxon: paired shift",
  "Cliff's delta: non-parametric effect size (>0 = CRE responds more)."
)
for (i in seq_along(rd_header))
  writeData(wb_supp, "Response_Differential", rd_header[i], startRow = i, startCol = 1)
addStyle(wb_supp, "Response_Differential",
         createStyle(textDecoration = "italic", fontColour = "#555555", wrapText = TRUE),
         rows = seq_along(rd_header), cols = 1, stack = TRUE)
writeData(wb_supp, "Response_Differential", resp_diag,
          startRow = length(rd_header) + 2,
          headerStyle = createStyle(textDecoration = "bold"))
quant_start <- length(rd_header) + 2 + nrow(resp_diag) + 2
writeData(wb_supp, "Response_Differential", "Quantile comparison (|logFC|):",
          startRow = quant_start)
writeData(wb_supp, "Response_Differential", quant_df,
          startRow = quant_start + 1,
          headerStyle = createStyle(textDecoration = "bold"))
setColWidths(wb_supp, "Response_Differential",
             cols = 1:max(ncol(resp_diag), ncol(quant_df)), widths = "auto")

add_sheet(wb_supp, "Effect_Size_CI", boot_df,
  title = "Median |logFC| with 95% BCa bootstrap confidence intervals (10,000 resamples).",
  notes = "Summarises global effect magnitude per contrast (Efron & Tibshirani 1993).")

add_sheet(wb_supp, "Power_Analysis", power_df,
  title = "Minimum detectable logFC at 80% power (alpha = 0.10).",
  notes = c("Conservative: uses pwr.t.test (standard t); limma's moderated t has higher power.",
            "Paired contrasts adjust for within-subject correlation; between-group do not.",
            "Interpret as approximate lower bounds on detectable effect sizes."))

if (!is.null(sens_df) && nrow(sens_df) > 0) {
  add_sheet(wb_supp, "Imputation_Sensitivity", sens_df,
    title = "Spearman correlation of t-statistics: non-imputed (main) vs imputed limma.",
    notes = c("High rho across all contrasts indicates robustness to imputation choice.",
              "Reference: Karpievitch et al. 2012, BMC Bioinform 13(S16):S5."))
}

saveWorkbook(wb_supp, file.path(cfg$data_dir, "10_DEP_supplementary.xlsx"),
             overwrite = TRUE)
cat("Done: 03_dep_robustness.R\n")
