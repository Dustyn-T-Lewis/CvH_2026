#!/usr/bin/env Rscript
# =============================================================================
# 01_impute.R  --  CvH Stage 03: Imputation (figures/WGCNA only)
#
# MAR/MNAR consensus (reporting) -> missForest on the full matrix. The imputed
# matrix feeds figures/WGCNA ONLY; Stage 04 DEP uses the NON-imputed normalized
# matrix (limma handles per-protein NAs). missForest: Stekhoven 2012 (PMID 22039212).
# =============================================================================

suppressPackageStartupMessages({
  library(missForest); library(here); library(readr); library(dplyr); library(tibble)
})
set.seed(42)
source(here("R", "mar_mnar.R"))

norm_csv <- here("02_Normalization", "c_data", "02_normalized.csv")
norm_rds <- here("02_Normalization", "c_data", "03_DAList_normalized.rds")
data_dir <- here("03_Imputation", "c_data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

df  <- read_csv(norm_csv, show_col_types = FALSE)
ann <- df |> select(uniprot_id, gene, protein, description)
mat <- as.matrix(df[, setdiff(names(df), names(ann))]); rownames(mat) <- ann$uniprot_id
cat(sprintf("Loaded: %d proteins x %d samples | %.2f%% missing\n",
            nrow(mat), ncol(mat), sum(is.na(mat)) / length(mat) * 100))

# --- MAR/MNAR classification (reporting) -------------------------------------
miss_class <- classify_mar_mnar(mat) |> mutate(gene = ann$gene, .after = uniprot_id)
cat(sprintf("Classification: Complete %d | MAR %d | MNAR %d\n",
            sum(miss_class$classification == "Complete"),
            sum(miss_class$classification == "MAR"),
            sum(miss_class$classification == "MNAR")))

# --- missForest on the full matrix (reorder by uniprot for determinism) ------
ord <- order(rownames(mat)); mat <- mat[ord, , drop = FALSE]; ann <- ann[ord, ]
set.seed(42)
mf  <- missForest(t(mat), maxiter = 10, ntree = 100, verbose = FALSE)
mat_imp <- t(mf$ximp); dimnames(mat_imp) <- dimnames(mat)
stopifnot(sum(is.na(mat_imp)) == 0)
oob <- as.numeric(mf$OOBerror[1])
cat(sprintf("missForest OOB error: %.4f\n", oob))

# --- Build imputed DAList (realign $annotation to $data row order) -----------
dal <- readRDS(norm_rds)
dal$data <- mat_imp
dal$annotation <- merge(dal$annotation,
                        miss_class |> select(uniprot_id, n_miss, pct_miss,
                                             miss_classification = classification, imputation_reliable),
                        by = "uniprot_id", all.x = TRUE, sort = FALSE)
dal$annotation <- dal$annotation[match(rownames(dal$data), dal$annotation$uniprot_id), , drop = FALSE]
rownames(dal$annotation) <- dal$annotation$uniprot_id
stopifnot(identical(rownames(dal$data), dal$annotation$uniprot_id))

# --- Export ------------------------------------------------------------------
write_csv(bind_cols(ann, as_tibble(mat_imp)), file.path(data_dir, "01_imputed.csv"))
write_csv(miss_class, file.path(data_dir, "02_mar_mnar_classification.csv"))
saveRDS(dal, file.path(data_dir, "01_DAList_imputed.rds"))
saveRDS(list(mat = mat, mat_imp = mat_imp, was_na = is.na(mat), miss_class = miss_class, oob = oob),
        file.path(data_dir, "00_report_intermediates.rds"))
writeLines(capture.output(sessionInfo()), file.path(data_dir, "sessionInfo.txt"))
cat(sprintf("Done: %d proteins x %d samples | OOB=%.4f -> %s/\n",
            nrow(mat_imp), ncol(mat_imp), oob, data_dir))
