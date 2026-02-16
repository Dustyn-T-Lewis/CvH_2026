#!/usr/bin/env Rscript
# CvH_imputation_run.R — MAR/MNAR classification + benchmark + imputation
# INPUT:  01_normalization/c_data/01_normalized.csv, 00_input/CvH_meta.csv
# OUTPUT: c_data/01_imputed.csv, 01_DAList_imputed.rds, benchmark_summary.csv
#         c_data/mar_mnar_classification.csv, imputation_effect_audit.csv
#         c_data/imputation_summary.txt
#         b_reports/01_missingness_report.pdf, 02_imputation_report.pdf

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(MsCoreUtils, msImpute, pcaMethods, imputeLCMD, missForest,
               readr, dplyr, tidyr, tibble, stringr,
               ggplot2, patchwork, scales)

# --- Paths ---
# Run from project root
input_csv  <- "01_normalization/c_data/01_normalized.csv"
meta_file  <- "00_input/CvH_meta.csv"
report_dir <- "02_Imputation/b_reports"
data_dir   <- "02_Imputation/c_data"
dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_dir,   showWarnings = FALSE, recursive = TRUE)

# --- Aesthetics ---
pal_mar <- c(MAR = "#4393C3", MNAR = "#D6604D")
pal_typ <- c(MNAR = "#D6604D", MAR = "#4393C3", Hybrid = "#5AAE61")
thm     <- theme_minimal(base_size = 11)

run_impute <- function(spec, mat, randna) {
  if (spec$method == "mixed")
    impute_matrix(mat, method = "mixed", randna = randna, mar = spec$mar, mnar = spec$mnar)
  else
    impute_matrix(mat, method = spec$method)
}

# --- 1. Load ---
cat("Step 1: Load\n")
stopifnot("Normalized CSV not found" = file.exists(input_csv))
df  <- read_csv(input_csv, show_col_types = FALSE)
ann_cols <- c("uniprot_id", "protein", "gene", "description", "n_seq")
ann <- df[, intersect(ann_cols, names(df))]
mat <- as.matrix(df[, setdiff(names(df), ann_cols)])
rownames(mat) <- df$uniprot_id

stopifnot("Metadata file not found" = file.exists(meta_file))
meta <- read_csv(meta_file, show_col_types = FALSE) %>% filter(Col_ID %in% colnames(mat))
mat  <- mat[, meta$Col_ID]

cat(sprintf("  %d proteins x %d samples\n", nrow(mat), ncol(mat)))

stopifnot("Annotation-matrix row mismatch" = nrow(ann) == nrow(mat))
stopifnot("UniProt IDs must match" = identical(ann$uniprot_id, rownames(mat)))

# --- 2. Missingness characterization ---
cat("Step 2: Missingness characterization\n")
pct_all  <- sum(is.na(mat)) / length(mat) * 100
prot_pct <- rowSums(is.na(mat)) / ncol(mat) * 100
obs_mean <- rowMeans(mat, na.rm = TRUE)

cat(sprintf("  Overall: %.1f%% missing\n", pct_all))

# --- 3. MAR/MNAR classification ---
cat("Step 3: MAR/MNAR classification\n")
has_na <- which(prot_pct > 0 & prot_pct < 100)

mar_result <- tryCatch({
  sf <- msImpute::selectFeatures(mat[has_na, ], method = "ebm",
                                  group = meta$Group_Time)
  mar_ids <- rownames(mat[has_na, ])[!sf$msImpute_feature]
  list(ok = TRUE, mar_ids = mar_ids)
}, error = function(e) {
  cat(sprintf("  selectFeatures failed: %s -> intensity heuristic\n", e$message))
  list(ok = FALSE)
})

miss_class <- tibble(uniprot_id = rownames(mat), gene = ann$gene,
                     n_miss = rowSums(is.na(mat)),
                     pct_miss = prot_pct, mean_int = obs_mean)

if (mar_result$ok) {
  miss_class <- miss_class %>%
    mutate(class = case_when(n_miss == 0 ~ "Complete",
                             uniprot_id %in% mar_result$mar_ids ~ "MAR",
                             TRUE ~ "MNAR"))
} else {
  med_int <- median(obs_mean, na.rm = TRUE)
  miss_class <- miss_class %>%
    mutate(class = case_when(
      n_miss == 0 ~ "Complete",
      pct_miss > 30 & mean_int < med_int ~ "MNAR",
      pct_miss > 50 ~ "MNAR",
      TRUE ~ "MAR"))
}

print(count(miss_class, class))

randna <- setNames(miss_class$class == "MNAR", miss_class$uniprot_id)

# --- Missingness report ---
cat("Step 3b: Missingness report\n")
pdf(file.path(report_dir, "01_missingness_report.pdf"), width = 12, height = 10)

n_total   <- length(mat)
n_missing <- sum(is.na(mat))
n_obs     <- n_total - n_missing
frac_df <- tibble(
  Status = factor(c("Observed", "Missing"), levels = c("Observed", "Missing")),
  n      = c(n_obs, n_missing),
  pct    = round(c(n_obs, n_missing) / n_total * 100, 1),
  label  = sprintf("%s\n%s (%.1f%%)", c("Observed", "Missing"),
                   format(c(n_obs, n_missing), big.mark = ","),
                   round(c(n_obs, n_missing) / n_total * 100, 1))
)

pFrac <- ggplot(frac_df, aes(x = "", y = n, fill = Status)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 4) +
  scale_fill_manual(values = c(Observed = "#2166AC", Missing = "#D6604D")) +
  coord_flip() +
  labs(x = NULL, y = "Total values (proteins x samples)",
       title = sprintf("A: Overall missingness — %d proteins x %d samples",
                       nrow(mat), ncol(mat))) +
  thm + theme(legend.position = "none", axis.text.y = element_blank(),
              axis.ticks.y = element_blank())

pA <- tibble(x = obs_mean, y = prot_pct) %>% ggplot(aes(x, y)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "loess", color = "#B2182B") +
  labs(x = "Mean log2 intensity", y = "% missing",
       title = "B: Protein missingness vs abundance") + thm

samp_miss <- tibble(
  Col_ID = colnames(mat),
  Missing  = colSums(is.na(mat)),
  Observed = nrow(mat) - colSums(is.na(mat))) %>%
  left_join(meta %>% dplyr::select(Col_ID, Group_Time), by = "Col_ID") %>%
  pivot_longer(c(Observed, Missing), names_to = "Status", values_to = "n") %>%
  mutate(Status = factor(Status, levels = c("Missing", "Observed")))

pB <- ggplot(samp_miss, aes(x = reorder(Col_ID, -n), y = n, fill = Status)) +
  geom_col(position = "stack", width = 0.8) +
  scale_fill_manual(values = c(Observed = "#2166AC", Missing = "#D6604D")) +
  labs(x = NULL, y = "Protein count",
       title = "C: Sample-level missingness (observed vs to-be-imputed)") +
  thm + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))

print(pFrac / pA / pB)

donut_df <- count(miss_class, class) %>%
  mutate(frac = n / sum(n),
         ymax = cumsum(frac), ymin = c(0, head(ymax, -1)),
         mid  = (ymin + ymax) / 2,
         label = sprintf("%s\n%d (%.0f%%)", class, n, frac * 100))
pal_donut <- c(Complete = "#999999", MAR = "#4393C3", MNAR = "#D6604D")
pC <- ggplot(donut_df, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 2.5, fill = class)) +
  geom_rect(alpha = 0.85) +
  geom_text(aes(x = 3.25, y = mid, label = label), size = 3.5) +
  scale_fill_manual(values = pal_donut) +
  coord_polar(theta = "y") + xlim(c(1, 4)) +
  labs(title = "C: Missingness classification") +
  theme_void(base_size = 11) + theme(legend.position = "none")

pD <- ggplot(filter(miss_class, class != "Complete"),
             aes(x = class, y = mean_int, fill = class)) +
  geom_boxplot(alpha = 0.6, width = 0.5) +
  scale_fill_manual(values = pal_mar) +
  labs(x = NULL, y = "Mean log2 intensity",
       title = "D: Intensity by missingness class") +
  thm + theme(legend.position = "none")

pE <- ggplot(filter(miss_class, class != "Complete"),
             aes(mean_int, fill = class)) +
  geom_density(alpha = 0.5) + scale_fill_manual(values = pal_mar) +
  labs(x = "Mean log2 intensity", y = "Density",
       title = "E: Intensity density by class") + thm

print((pC | pD) / pE)
dev.off()

# --- 4. Benchmark ---
cat("Step 4: Benchmark (NRMSE, 10%% masked, 20 iters)\n")
set.seed(42); N_ITER <- 20

METHODS <- list(
  MinProb = list(method = "MinProb"),
  MinDet  = list(method = "MinDet"),
  QRILC   = list(method = "QRILC"),
  knn     = list(method = "knn"),
  bpca    = list(method = "bpca"),
  MLE     = list(method = "MLE"),
  RF      = list(method = "RF"),
  mix_knn_MinProb  = list(method = "mixed", mar = "knn",  mnar = "MinProb"),
  mix_knn_MinDet   = list(method = "mixed", mar = "knn",  mnar = "MinDet"),
  mix_knn_QRILC    = list(method = "mixed", mar = "knn",  mnar = "QRILC"),
  mix_bpca_MinProb = list(method = "mixed", mar = "bpca", mnar = "MinProb"),
  mix_bpca_QRILC   = list(method = "mixed", mar = "bpca", mnar = "QRILC"),
  mix_RF_MinProb   = list(method = "mixed", mar = "RF",   mnar = "MinProb"),
  mix_RF_QRILC     = list(method = "mixed", mar = "RF",   mnar = "QRILC")
)

nrmse_fn <- function(true, imp) sqrt(mean((true - imp)^2)) / sd(true)

res <- list(); k <- 0L
for (it in seq_len(N_ITER)) {
  obs   <- which(!is.na(mat))
  mask  <- sample(obs, round(length(obs) * 0.10))
  truth <- mat[mask]
  mm    <- mat; mm[mask] <- NA

  masked_rows <- unique((mask - 1) %% nrow(mat) + 1)
  rn_iter     <- randna; rn_iter[masked_rows] <- FALSE

  for (nm in names(METHODS)) {
    imp <- tryCatch(run_impute(METHODS[[nm]], mm, rn_iter), error = function(e) NULL)
    if (is.null(imp)) next
    k <- k + 1L
    res[[k]] <- tibble(method = nm, iter = it, nrmse = nrmse_fn(truth, imp[mask]))
  }
  if (it %% 5 == 0) cat(sprintf("  iter %d/%d\n", it, N_ITER))
}

bench <- bind_rows(res) %>% group_by(method) %>%
  summarise(mean = mean(nrmse), sd = sd(nrmse), median = median(nrmse), .groups = "drop") %>%
  arrange(mean)
print(bench)

best <- bench$method[1]
cat(sprintf("  Best: %s (NRMSE = %.4f)\n", best, bench$mean[1]))

mtype <- tibble(method = names(METHODS)) %>%
  mutate(type = case_when(
    method %in% c("MinProb","MinDet","QRILC") ~ "MNAR",
    method %in% c("knn","bpca","MLE","RF")    ~ "MAR",
    TRUE                                       ~ "Hybrid"))

# --- 5. Apply best method ---
cat("Step 5: Impute\n")
set.seed(42)
mat_imp <- run_impute(METHODS[[best]], mat, randna)
cat(sprintf("  Remaining NAs: %d\n", sum(is.na(mat_imp))))

stopifnot("Imputed matrix dims changed" = identical(dim(mat), dim(mat_imp)))
stopifnot("Row names altered" = identical(rownames(mat), rownames(mat_imp)))

# --- 6. Imputation report ---
cat("Step 6: Imputation report\n")
was_na <- is.na(mat)

meta_pca <- meta %>%
  mutate(BioGroup = case_when(
    Group_Time %in% c("CRE_T1", "PLA_T1") ~ "CR_T1",
    Group_Time %in% c("CRE_T2", "PLA_T2") ~ "CR_T2",
    TRUE ~ "PPS_T1"))
pal_bio <- c(CR_T1 = "#4393C3", CR_T2 = "#2166AC", PPS_T1 = "#B2182B")

pdf(file.path(report_dir, "02_imputation_report.pdf"), width = 12, height = 10)

bench_df <- bind_rows(res) %>% left_join(mtype, by = "method")

pBench <- ggplot(bench_df, aes(x = nrmse, y = reorder(method, nrmse, median), fill = type)) +
  geom_boxplot(alpha = 0.6) +
  scale_fill_manual(values = pal_typ) +
  scale_x_log10() +
  annotation_logticks(sides = "b") +
  labs(y = NULL, x = "NRMSE (log scale)", title = sprintf("A: Benchmark (best: %s)", best)) +
  thm + theme(legend.position = "bottom")

top3 <- bench$method[1:3]
pTop3 <- bench_df %>% filter(method %in% top3) %>%
  ggplot(aes(x = nrmse, y = reorder(method, nrmse, median), fill = type)) +
  geom_boxplot(alpha = 0.6, width = 0.5) +
  geom_vline(xintercept = bench$mean[1], linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = pal_typ) +
  labs(y = NULL, x = "NRMSE",
       title = sprintf("B: Top 3 — %s %.4f | %s %.4f | %s %.4f",
                        top3[1], bench$mean[bench$method == top3[1]],
                        top3[2], bench$mean[bench$method == top3[2]],
                        top3[3], bench$mean[bench$method == top3[3]])) +
  thm + theme(legend.position = "none")

print(pBench); print(pTop3)

obs_vals <- as.vector(mat[!was_na])
imp_vals <- mat_imp[was_na]

brks <- seq(floor(min(obs_vals, na.rm = TRUE)),
            ceiling(max(c(obs_vals, imp_vals), na.rm = TRUE)), by = 0.25)
hist_df <- bind_rows(
  tibble(log2 = obs_vals, Source = "Observed"),
  tibble(log2 = imp_vals, Source = "Imputed"))

pHist <- ggplot(hist_df, aes(x = log2, fill = Source)) +
  geom_histogram(breaks = brks, position = "identity", alpha = 0.7) +
  scale_fill_manual(values = c(Observed = "#2166AC", Imputed = "#D6604D")) +
  labs(x = "log2 intensity", y = "Count",
       title = sprintf("C: Intensity distribution — observed (n=%s) vs imputed (n=%s)",
                        format(length(obs_vals), big.mark = ","),
                        format(length(imp_vals), big.mark = ","))) +
  thm + theme(legend.position = "bottom")
print(pHist)

med_fill <- mat
for (j in 1:ncol(med_fill)) med_fill[is.na(med_fill[,j]),j] <- median(med_fill[,j], na.rm=TRUE)

make_pca <- function(m, label) {
  pca_out <- prcomp(t(m), center = TRUE, scale. = TRUE)
  ve <- round(summary(pca_out)$importance[2, 1:2] * 100, 1)
  as.data.frame(pca_out$x[, 1:2]) %>%
    mutate(Col_ID = rownames(.)) %>%
    left_join(meta_pca, by = "Col_ID") %>%
    ggplot(aes(PC1, PC2, color = BioGroup, fill = BioGroup)) +
    geom_point(size = 3) +
    stat_ellipse(geom = "polygon", alpha = 0.1, level = 0.68, linewidth = 0.6) +
    scale_color_manual(values = pal_bio) +
    scale_fill_manual(values = pal_bio) +
    labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]),
         title = label, color = "Group", fill = "Group") +
    thm + theme(legend.position = "bottom")
}

pPCA_pre  <- make_pca(med_fill, "D: PCA pre-imputation (normalized, outlier-removed)")
pPCA_post <- make_pca(mat_imp,  sprintf("E: PCA post-imputation (%s)", best))
print(pPCA_pre | pPCA_post)

samp_pre  <- tibble(Col_ID = rep(colnames(mat), each = nrow(mat)),
                    log2 = as.vector(mat), stage = "Pre")
samp_post <- tibble(Col_ID = rep(colnames(mat_imp), each = nrow(mat_imp)),
                    log2 = as.vector(mat_imp), stage = "Post")
box_df <- bind_rows(samp_pre, samp_post) %>%
  left_join(meta_pca %>% dplyr::select(Col_ID, BioGroup), by = "Col_ID") %>%
  mutate(stage = factor(stage, levels = c("Pre", "Post")))

pBox <- ggplot(box_df, aes(x = Col_ID, y = log2, fill = stage)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0.3, width = 0.7) +
  scale_fill_manual(values = c(Pre = "#92C5DE", Post = "#2166AC")) +
  labs(x = NULL, y = "log2 intensity",
       title = "F: Per-sample log2 distributions (pre vs post imputation)") +
  thm + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))
print(pBox)

inc_ids <- miss_class$uniprot_id[miss_class$n_miss > 0]
audit <- tibble(
  uniprot_id = inc_ids,
  gene = ann$gene[match(inc_ids, ann$uniprot_id)],
  pre_mean  = rowMeans(mat[inc_ids, ], na.rm = TRUE),
  post_mean = rowMeans(mat_imp[inc_ids, ]),
  pre_sd    = apply(mat[inc_ids, ], 1, sd, na.rm = TRUE),
  pct_miss  = prot_pct[inc_ids]) %>%
  mutate(shift    = post_mean - pre_mean,
         effect_d = shift / pre_sd)

pEff <- ggplot(audit, aes(pct_miss, effect_d)) +
  geom_point(alpha = 0.5, size = 1, color = "#4393C3") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = "dotted", color = "grey50") +
  annotate("text", x = max(audit$pct_miss) * 0.9, y = 0.22,
           label = "|d| = 0.2 (small)", size = 3, color = "grey40") +
  labs(x = "% missing", y = "Effect size (shift / SD)",
       title = "G: Imputation effect on protein means (Cohen's d)") + thm

samp_means <- tibble(
  Col_ID = colnames(mat),
  Pre  = colMeans(mat, na.rm = TRUE),
  Post = colMeans(mat_imp)) %>%
  left_join(meta_pca %>% dplyr::select(Col_ID, BioGroup), by = "Col_ID")

pSampMean <- ggplot(samp_means, aes(Pre, Post, color = BioGroup)) +
  geom_point(size = 3) + geom_abline(linetype = "dashed") +
  scale_color_manual(values = pal_bio) +
  labs(x = "Pre-imputation sample mean (log2)",
       y = "Post-imputation sample mean (log2)",
       title = "H: Sample-level mean shift (identity line = no change)") +
  thm + theme(legend.position = "bottom")

print(pEff | pSampMean)
dev.off()

write_csv(audit, file.path(data_dir, "imputation_effect_audit.csv"))
write_csv(bench, file.path(data_dir, "benchmark_summary.csv"))

# --- 7. Export ---
cat("Step 7: Export\n")

stopifnot("UniProt ID order mismatch" = identical(ann$uniprot_id, rownames(mat_imp)))
stopifnot("Observed values altered" = all.equal(mat[!was_na], mat_imp[!was_na], tolerance = 1e-10))

write_csv(bind_cols(ann, as_tibble(mat_imp)), file.path(data_dir, "01_imputed.csv"))
write_csv(miss_class, file.path(data_dir, "mar_mnar_classification.csv"))

dal_path <- "01_normalization/c_data/01_DAList_normalized.rds"
if (file.exists(dal_path)) {
  dal <- readRDS(dal_path)
  rownames(mat_imp) <- ann$uniprot_id
  dal$data <- mat_imp
  saveRDS(dal, file.path(data_dir, "01_DAList_imputed.rds"))
}

info <- c(proteins = nrow(mat), samples = ncol(mat), pct_missing = round(pct_all, 2),
          n_MAR = sum(miss_class$class == "MAR"),
          n_MNAR = sum(miss_class$class == "MNAR"),
          best_method = best, best_NRMSE = round(bench$mean[1], 4))
writeLines(paste(names(info), info, sep = " = "), file.path(data_dir, "imputation_summary.txt"))

cat(sprintf("  Exported: %s | NRMSE %.4f\n", best, bench$mean[1]))
