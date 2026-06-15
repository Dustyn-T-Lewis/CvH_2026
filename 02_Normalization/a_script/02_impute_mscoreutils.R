#!/usr/bin/env Rscript
# =============================================================================
# 02_impute_mscoreutils.R  --  CvH Stage 02: hybrid imputation (MsCoreUtils)
#
# Mechanism-aware HYBRID imputation (Lazar et al. 2016, J Proteome Res):
#   MsCoreUtils::impute_matrix(method = "mixed", randna = <MAR features>,
#                              mar = "knn", mnar = "QRILC")
# Per the MsCoreUtils docs, `mixed` applies `mar` to the MCAR/MAR feature subset
# and `mnar` to the MNAR subset, split by the logical `randna` (length = nrow,
# TRUE = missing-at-random). We define `randna` from our 3-method MAR/MNAR
# consensus (R/mar_mnar.R): a protein is MAR unless classified MNAR.
#   - mar  = "knn"   : k-nearest-neighbour averaging (MAR; correlated features)
#   - mnar = "QRILC" : quantile regression imputation of left-censored data
#                      (imputeLCMD::impute.QRILC; the canonical MNAR method)
# Imputed matrix feeds 03_DEP/b_imputed and the figures; DEP/a_non_imputed uses
# the non-imputed matrix (imputing before testing is not done there).
# =============================================================================

suppressPackageStartupMessages({
  library(proteoDA); library(here); library(MsCoreUtils); library(dplyr); library(tibble); library(readr)
})
set.seed(42)
source(here("R", "mar_mnar.R"))
data_dir <- here("02_Normalization", "c_data")

dal <- readRDS(file.path(data_dir, "03_DAList_normalized.rds"))
mat <- as.matrix(dal$data)                              # proteins x samples, log2, with NAs
cat(sprintf("Loaded normalized matrix: %d x %d | %.1f%% missing\n",
            nrow(mat), ncol(mat), mean(is.na(mat)) * 100))

# MAR/MNAR consensus -> randna (TRUE = MAR feature). Complete features (no NA) are
# treated as MAR (they are untouched by imputation anyway).
cls <- classify_mar_mnar(mat)
randna <- cls$classification != "MNAR"
cat(sprintf("MAR/MNAR split: %d MAR(+Complete) / %d MNAR features\n", sum(randna), sum(!randna)))

imp <- impute_matrix(mat, method = "mixed", randna = randna, mar = "knn", mnar = "QRILC")
stopifnot(sum(is.na(imp)) == 0, identical(dim(imp), dim(mat)))

dal$data <- imp
dal$imputation <- list(method = "MsCoreUtils mixed", mar = "knn", mnar = "QRILC",
                       n_mar = sum(randna), n_mnar = sum(!randna))
saveRDS(dal, file.path(data_dir, "04_DAList_imputed_mscoreutils.rds"))
write_csv(bind_cols(as_tibble(dal$annotation) |> select(uniprot_id, protein, gene, description),
                    as_tibble(imp)), file.path(data_dir, "04_imputed_mscoreutils.csv"))
cat(sprintf("Done: hybrid (knn/QRILC) imputed %d cells -> 04_DAList_imputed_mscoreutils.rds\n", sum(is.na(mat))))
