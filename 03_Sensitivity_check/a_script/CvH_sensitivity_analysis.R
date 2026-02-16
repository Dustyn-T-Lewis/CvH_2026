#!/usr/bin/env Rscript
# CvH_sensitivity_analysis.R — Top-3 norms x (top-3 imps + none) = 12 combos
# INPUT:  01_normalization/c_data/00_DAList_pre_norm.rds, 05_norm_quality_scores.csv
#         02_Imputation/c_data/benchmark_summary.csv
#         01_normalization/c_data/01_DAList_normalized.rds
# OUTPUT: c_data/sensitivity_comparison.csv, pairwise_spearman.csv
#         c_data/cat_curves.csv, dep_counts.csv
#         b_reports/sensitivity_comparison.pdf

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(proteoDA, MsCoreUtils, msImpute, pcaMethods, imputeLCMD, missForest,
               limma, readr, dplyr, tidyr, tibble, stringr,
               ggplot2, patchwork, scales, here)

# --- Paths ---
base_dir   <- here::here()
report_dir <- file.path(base_dir, "03_Sensitivity_check", "b_reports")
data_dir   <- file.path(base_dir, "03_Sensitivity_check", "c_data")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir,   recursive = TRUE, showWarnings = FALSE)

norm_scores_file <- file.path(base_dir, "01_normalization", "c_data", "05_norm_quality_scores.csv")
bench_file       <- file.path(base_dir, "02_Imputation", "c_data", "benchmark_summary.csv")
pre_norm_file    <- file.path(base_dir, "01_normalization", "c_data", "00_DAList_pre_norm.rds")
dal_norm_file    <- file.path(base_dir, "01_normalization", "c_data", "01_DAList_normalized.rds")

stopifnot("Norm quality scores not found — run 01_normalization first" = file.exists(norm_scores_file))
stopifnot("Benchmark summary not found — run 02_Imputation first"     = file.exists(bench_file))
stopifnot("Pre-norm DAList not found"    = file.exists(pre_norm_file))
stopifnot("Normalized DAList not found"  = file.exists(dal_norm_file))

# --- 1. Read upstream rankings ---
cat("Step 1: Read upstream rankings\n")

norm_scores <- read_csv(norm_scores_file, show_col_types = FALSE)
top_norms   <- norm_scores$norm[1:3]
cat(sprintf("  Top 3 norms: %s\n", paste(top_norms, collapse = ", ")))

bench <- read_csv(bench_file, show_col_types = FALSE) %>% arrange(mean)
top_imps <- bench$method[1:3]
cat(sprintf("  Top 3 imps:  %s\n", paste(top_imps, collapse = ", ")))

IMP_SPECS <- list(none = NULL)
for (imp_name in top_imps) {
  if (grepl("^mix_", imp_name)) {
    parts <- strsplit(sub("^mix_", "", imp_name), "_")[[1]]
    IMP_SPECS[[imp_name]] <- list(method = "mixed", mar = parts[1], mnar = parts[2])
  } else {
    IMP_SPECS[[imp_name]] <- list(method = imp_name)
  }
}

# --- 2. Load data & classify MAR/MNAR ---
cat("Step 2: Load data\n")

dal_raw  <- readRDS(pre_norm_file)
dal_norm <- readRDS(dal_norm_file)

dal_raw$metadata$group <- factor(dal_raw$metadata$Group_Time,
                                  levels = c("CRE_T1","CRE_T2","PLA_T1","PLA_T2","H_T1"))
meta <- dal_raw$metadata

cat(sprintf("  %d proteins x %d samples\n", nrow(dal_raw$data), ncol(dal_raw$data)))

ref_mat  <- as.matrix(dal_norm$data)
rownames(ref_mat) <- dal_norm$annotation$gene
prot_pct <- rowSums(is.na(ref_mat)) / ncol(ref_mat) * 100
obs_mean <- rowMeans(ref_mat, na.rm = TRUE)

randna <- tryCatch({
  has_na <- which(prot_pct > 0 & prot_pct < 100)
  sf <- msImpute::selectFeatures(ref_mat[has_na, ], method = "ebm",
                                  group = dal_norm$metadata$Group_Time)
  mnar_genes <- rownames(ref_mat[has_na, ])[sf$msImpute_feature]
  is_mnar <- (rownames(ref_mat) %in% mnar_genes)
  setNames(is_mnar, rownames(ref_mat))
}, error = function(e) {
  cat(sprintf("  selectFeatures failed: %s -> intensity heuristic\n", e$message))
  med <- median(obs_mean, na.rm = TRUE)
  is_mnar <- (prot_pct > 30 & obs_mean < med) | prot_pct > 50
  setNames(is_mnar, rownames(ref_mat))
})

CONTRASTS <- c(
  Training_CR           = "(CRE_T2 + PLA_T2)/2 - (CRE_T1 + PLA_T1)/2",
  Cancer_vs_Healthy     = "(CRE_T1 + PLA_T1)/2 - H_T1",
  Supplement_Interaction = "(CRE_T2 - CRE_T1) - (PLA_T2 - PLA_T1)"
)

# --- 3. Run 12-combo grid ---
cat("Step 3: Run 12-combo grid\n")

run_pipeline <- function(dal, norm_method, imp_spec, randna_vec) {
  dal_n <- tryCatch(normalize_data(dal, norm_method = norm_method),
                    error = function(e) NULL)
  if (is.null(dal_n)) return(NULL)

  mat <- as.matrix(dal_n$data)

  if (!is.null(imp_spec)) {
    rn <- randna_vec[rownames(mat)]
    rn[is.na(rn)] <- FALSE

    mat <- tryCatch({
      if (imp_spec$method == "mixed")
        impute_matrix(mat, method = "mixed", randna = rn,
                      mar = imp_spec$mar, mnar = imp_spec$mnar)
      else
        impute_matrix(mat, method = imp_spec$method)
    }, error = function(e) NULL)
    if (is.null(mat)) return(NULL)
  }

  meta   <- dal_n$metadata
  design <- model.matrix(~ 0 + group, data = meta)
  colnames(design) <- gsub("^group", "", colnames(design))

  contrast_mat <- makeContrasts(
    Training_CR            = (CRE_T2 + PLA_T2)/2 - (CRE_T1 + PLA_T1)/2,
    Cancer_vs_Healthy      = (CRE_T1 + PLA_T1)/2 - H_T1,
    Supplement_Interaction = (CRE_T2 - CRE_T1) - (PLA_T2 - PLA_T1),
    levels = design)

  dupcor <- tryCatch(
    duplicateCorrelation(mat, design, block = meta$Subject_ID),
    error = function(e) {
      warning(sprintf("[%s_%s] duplicateCorrelation failed: %s\n  -> Falling back to correlation = 0.",
                      norm_method, deparse(imp_spec), e$message), call. = FALSE)
      list(consensus.correlation = 0)
    })

  fit <- lmFit(mat, design, block = meta$Subject_ID,
               correlation = dupcor$consensus.correlation) |>
    contrasts.fit(contrast_mat) |> eBayes()

  dep_results <- list()
  for (i in seq_along(CONTRASTS)) {
    cname <- names(CONTRASTS)[i]
    tt <- topTable(fit, coef = i, number = Inf, sort.by = "none")
    tt$gene <- rownames(mat)
    dep_results[[cname]] <- tt
  }

  dep_results
}

combos <- expand.grid(norm = top_norms, imp = names(IMP_SPECS), stringsAsFactors = FALSE)
combos$label <- paste0(combos$norm, "_", combos$imp)
imp_clustered_order <- combos %>%
  arrange(factor(imp, levels = names(IMP_SPECS)), norm) %>%
  pull(label)

all_results <- list()
for (i in seq_len(nrow(combos))) {
  lab <- combos$label[i]
  cat(sprintf("  [%d/%d] %s\n", i, nrow(combos), lab))

  res <- tryCatch(
    run_pipeline(dal_raw, combos$norm[i], IMP_SPECS[[combos$imp[i]]], randna),
    error = function(e) { cat(sprintf("  FAILED: %s\n", e$message)); NULL })

  if (!is.null(res)) all_results[[lab]] <- res
}

cat(sprintf("  Completed: %d / %d combinations\n", length(all_results), nrow(combos)))
method_names <- names(all_results)
n_methods    <- length(method_names)

# --- 4. Pairwise Spearman rho ---
cat("Step 4: Pairwise Spearman rho\n")

all_genes <- unique(unlist(lapply(all_results, function(res) res[[1]]$gene)))

spearman_by_contrast <- list()
for (cname in names(CONTRASTS)) {
  t_mat <- sapply(method_names, function(lab) {
    tt <- all_results[[lab]][[cname]]
    t_vals <- setNames(rep(NA_real_, length(all_genes)), all_genes)
    t_vals[tt$gene] <- tt$t
    t_vals
  })

  spearman_mat <- cor(t_mat, method = "spearman", use = "pairwise.complete.obs")
  spearman_by_contrast[[cname]] <- spearman_mat

  mean_rho <- mean(spearman_mat[lower.tri(spearman_mat)], na.rm = TRUE)
  min_rho  <- min(spearman_mat[lower.tri(spearman_mat)], na.rm = TRUE)
  cat(sprintf("  %s: mean rho = %.4f, min rho = %.4f\n", cname, mean_rho, min_rho))
}

spearman_long <- list()
for (cname in names(spearman_by_contrast)) {
  sm <- spearman_by_contrast[[cname]]
  pairs <- which(lower.tri(sm), arr.ind = TRUE)
  spearman_long[[cname]] <- tibble(
    contrast = cname,
    method_a = rownames(sm)[pairs[, 1]],
    method_b = colnames(sm)[pairs[, 2]],
    spearman_rho = sm[pairs]
  )
}
spearman_df <- bind_rows(spearman_long)
write_csv(spearman_df, file.path(data_dir, "pairwise_spearman.csv"))

# --- 5. CAT curves ---
cat("Step 5: CAT curves\n")

compute_cat <- function(ranks_a, ranks_b, max_k = 300) {
  shared <- intersect(names(ranks_a), names(ranks_b))
  if (length(shared) < 2) return(rep(NA_real_, max_k))
  ord_a <- names(sort(ranks_a[shared]))
  ord_b <- names(sort(ranks_b[shared]))
  ks <- seq_len(min(max_k, length(shared)))
  sapply(ks, function(k) length(intersect(ord_a[1:k], ord_b[1:k])) / k)
}

cat_data_all <- list()
for (cname in names(CONTRASTS)) {
  consensus_ranks <- setNames(rep(NA_real_, length(all_genes)), all_genes)
  for (g in all_genes) {
    ranks_g <- sapply(all_results, function(res) {
      tt <- res[[cname]]
      rank(-abs(tt$t))[match(g, tt$gene)]
    })
    consensus_ranks[g] <- mean(ranks_g, na.rm = TRUE)
  }
  consensus_ranks <- consensus_ranks[!is.na(consensus_ranks)]

  cat_data <- lapply(method_names, function(lab) {
    tt <- all_results[[lab]][[cname]]
    method_ranks <- setNames(rank(-abs(tt$t)), tt$gene)
    cat_vals <- compute_cat(method_ranks, consensus_ranks, max_k = 300)
    tibble(k = seq_along(cat_vals), concordance = cat_vals, method = lab, contrast = cname)
  }) %>% bind_rows()

  cat_data_all[[cname]] <- cat_data
}

write_csv(bind_rows(cat_data_all), file.path(data_dir, "cat_curves.csv"))

# --- 6. DEP counts ---
cat("Step 6: DEP counts\n")
dep_summary <- list()
for (lab in method_names) {
  for (cname in names(CONTRASTS)) {
    tt <- all_results[[lab]][[cname]]

    tt$pi_score <- -log10(tt$P.Value) * abs(tt$logFC)
    tt$pi_pval  <- rank(-tt$pi_score) / nrow(tt)

    criteria <- list(
      nominal_p05  = tt %>% filter(!is.na(P.Value)   & P.Value   < 0.05),
      FDR_05       = tt %>% filter(!is.na(adj.P.Val) & adj.P.Val < 0.05),
      FDR_10       = tt %>% filter(!is.na(adj.P.Val) & adj.P.Val < 0.10),
      pi_score_05  = tt %>% filter(pi_pval < 0.05)
    )

    for (crit_name in names(criteria)) {
      sig <- criteria[[crit_name]]
      dep_summary <- c(dep_summary, list(tibble(
        method    = lab,
        contrast  = cname,
        criterion = crit_name,
        n_up      = sum(sig$logFC > 0),
        n_down    = sum(sig$logFC < 0),
        n_total   = nrow(sig)
      )))
    }
  }
}
dep_df <- bind_rows(dep_summary)

write_csv(dep_df, file.path(data_dir, "dep_counts.csv"))

# --- 7. Summary table ---
cat("Step 7: Summary table\n")

summary_rows <- list()
for (lab in method_names) {
  mean_rhos <- sapply(names(spearman_by_contrast), function(cname) {
    sm <- spearman_by_contrast[[cname]]
    mean(sm[lab, setdiff(method_names, lab)], na.rm = TRUE)
  })

  dep_fdr05 <- dep_df %>% filter(method == lab, criterion == "FDR_05") %>%
    summarise(
      n_DEP_Training    = n_total[contrast == "Training_CR"],
      n_DEP_CvH         = n_total[contrast == "Cancer_vs_Healthy"],
      n_DEP_Supplement   = n_total[contrast == "Supplement_Interaction"]
    )

  dep_extra <- dep_df %>%
    filter(method == lab, contrast == "Cancer_vs_Healthy") %>%
    summarise(
      n_CvH_nominal_p05 = n_total[criterion == "nominal_p05"],
      n_CvH_FDR_10      = n_total[criterion == "FDR_10"],
      n_CvH_pi_score    = n_total[criterion == "pi_score_05"]
    )

  parts <- strsplit(lab, "_", fixed = TRUE)[[1]]

  summary_rows <- c(summary_rows, list(tibble(
    method           = lab,
    norm             = parts[1],
    imputation       = paste(parts[-1], collapse = "_"),
    mean_spearman    = round(mean(mean_rhos), 4),
    rho_Training     = round(mean_rhos["Training_CR"], 4),
    rho_CvH          = round(mean_rhos["Cancer_vs_Healthy"], 4),
    rho_Supplement   = round(mean_rhos["Supplement_Interaction"], 4),
    n_DEP_Training   = dep_fdr05$n_DEP_Training,
    n_DEP_CvH        = dep_fdr05$n_DEP_CvH,
    n_DEP_Supplement = dep_fdr05$n_DEP_Supplement,
    n_CvH_nominal_p05 = dep_extra$n_CvH_nominal_p05,
    n_CvH_FDR_10      = dep_extra$n_CvH_FDR_10,
    n_CvH_pi_score    = dep_extra$n_CvH_pi_score
  )))
}
summary_df <- bind_rows(summary_rows) %>% arrange(desc(mean_spearman))

print(summary_df)

write_csv(summary_df, file.path(data_dir, "sensitivity_comparison.csv"))

# --- 8. Visuals ---
cat("Step 8: Visuals\n")

pdf(file.path(report_dir, "sensitivity_comparison.pdf"), width = 14, height = 10)

heat_data <- list()
for (cname in names(spearman_by_contrast)) {
  sm <- spearman_by_contrast[[cname]]
  heat_data[[cname]] <- expand.grid(
    method_a = factor(rownames(sm), levels = imp_clustered_order),
    method_b = factor(colnames(sm), levels = imp_clustered_order),
    stringsAsFactors = FALSE
  ) %>%
    mutate(rho = as.vector(sm), contrast = cname)
}
heat_df <- bind_rows(heat_data)

p_heat <- ggplot(heat_df, aes(x = method_a, y = method_b, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.3f", rho)), size = 2) +
  facet_wrap(~ contrast, ncol = 2) +
  scale_fill_gradient2(low = "#B2182B", mid = "#F7F7F7", high = "#2166AC",
                       midpoint = 0.95, limits = c(
                         min(heat_df$rho, na.rm = TRUE),
                         max(heat_df$rho, na.rm = TRUE)),
                       name = "Spearman\nrho") +
  labs(title = "Pairwise Spearman rho of t-statistics (Karpievitch et al. 2012)",
       subtitle = "Higher = methods rank proteins more similarly",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold"))
print(p_heat)

cat_df <- bind_rows(cat_data_all)
cat_plots <- list()
for (cname in names(CONTRASTS)) {
  cat_plots[[cname]] <- ggplot(
    cat_df %>% filter(contrast == cname),
    aes(x = k, y = concordance, color = factor(method, levels = imp_clustered_order))) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = 0.7, linetype = "dashed", alpha = 0.4) +
    labs(x = "Top k proteins", y = "Concordance", title = cname, color = "Pipeline") +
    theme_minimal(base_size = 9) +
    scale_color_brewer(palette = "Paired") +
    ylim(0, 1)
}
p_cat <- wrap_plots(cat_plots, ncol = 1) +
  plot_annotation(
    title = "Concordance At the Top (CAT) — each method vs consensus ranking",
    subtitle = "Irizarry et al. (2005) | Dashed line = 0.70 concordance",
    theme = theme(plot.title = element_text(face = "bold")))
print(p_cat)

dep_fdr05_long <- dep_df %>%
  filter(criterion == "FDR_05") %>%
  mutate(method = factor(method, levels = imp_clustered_order)) %>%
  pivot_longer(c(n_up, n_down), names_to = "direction", values_to = "count") %>%
  mutate(direction = ifelse(direction == "n_up", "Up", "Down"),
         count = ifelse(direction == "Down", -count, count))

p_dep <- ggplot(dep_fdr05_long, aes(x = method, y = count, fill = direction)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, color = "grey30") +
  facet_wrap(~ contrast, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c(Up = "#D6604D", Down = "#4393C3")) +
  labs(x = NULL, y = "DEP count (up / down)",
       title = "DEP counts per method per contrast (FDR < 0.05)",
       fill = "Direction") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))
print(p_dep)

dep_criteria <- dep_df %>%
  mutate(method = factor(method, levels = imp_clustered_order),
         criterion = factor(criterion,
    levels = c("nominal_p05", "FDR_10", "FDR_05", "pi_score_05"),
    labels = c("Nominal p<0.05", "FDR<0.10", "FDR<0.05", "Pi-score<0.05")))

p_criteria <- ggplot(dep_criteria, aes(x = method, y = n_total, fill = criterion)) +
  geom_col(position = "dodge", width = 0.7) +
  facet_wrap(~ contrast, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("Nominal p<0.05" = "#92C5DE",
                                "FDR<0.10"       = "#F4A582",
                                "FDR<0.05"       = "#D6604D",
                                "Pi-score<0.05"  = "#5AAE61")) +
  labs(x = NULL, y = "Total DEPs",
       title = "DEP counts across significance criteria",
       subtitle = "Pi-score = -log10(p) x |logFC| (Xiao et al. 2014); empirical p-value < 0.05",
       fill = "Criterion") +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"))
print(p_criteria)

dev.off()

# --- 9. Interpretation ---
overall_rho <- mean(spearman_df$spearman_rho, na.rm = TRUE)
cat(sprintf("Overall mean Spearman rho: %.4f\n", overall_rho))
if (overall_rho > 0.95) {
  cat("Methods agree closely (rho > 0.95). Choice matters little.\n")
} else if (overall_rho > 0.85) {
  cat("Methods mostly agree (rho > 0.85). Some variation in extremes.\n")
} else {
  cat("Methods diverge meaningfully (rho < 0.85). Inspect contrasts individually.\n")
}

cat(sprintf("Top pipeline: %s (rho = %.4f)\n",
            summary_df$method[1], summary_df$mean_spearman[1]))
cat(sprintf("Done: %d pipelines x %d contrasts\n", n_methods, length(CONTRASTS)))
