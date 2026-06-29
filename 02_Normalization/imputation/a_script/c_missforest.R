#!/usr/bin/env Rscript
# STANDALONE imputation option C: missForest (exploratory arm; not the canonical DEP input).
#
# missForest (Stekhoven & Buhlmann 2012) is a non-parametric random-forest imputer: each
# protein's missing values are predicted from the multivariate structure of the other samples,
# with no left-censoring model. Stage-02 missingness here is largely left-censored (low-abundance
# dropout), which missForest does not model, yet it tracks the non-imputed effect sizes closely
# and barely adds DEPs. The CvH figures/WGCNA read the imp4p arm; missForest is kept alongside
# imp4p and the MsCoreUtils hybrid as a comparison, with the DEP input staying non-imputed.

pacman::p_load(proteoDA, here, missForest)
set.seed(42) # missForest is stochastic (random forests)
norm_dir <- here("02_Normalization", "c_data") # read stage-02 normalized matrix
data_dir <- here("02_Normalization", "imputation", "c_data") # write imputed DAList here
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

#### Load normalized DAList ####

dal <- readRDS(file.path(norm_dir, "DAList_normalized.rds"))
mat <- as.matrix(dal$data)
cat(sprintf("[missforest] %d x %d | %.1f%% missing\n", nrow(mat), ncol(mat), mean(is.na(mat)) * 100))

#### Impute ####

mf <- missForest(mat, maxiter = 10, ntree = 100, verbose = FALSE)
imp <- mf$ximp
dimnames(imp) <- dimnames(mat)
stopifnot(sum(is.na(imp)) == 0, identical(dim(imp), dim(mat)))

#### Export ####

dal$data <- imp
dal$imputation <- list(
  method = "missForest::missForest", maxiter = 10, ntree = 100,
  OOB_NRMSE = unname(mf$OOBerror["NRMSE"])
)
saveRDS(dal, file.path(data_dir, "DAList_imputed_missforest.rds"))
cat(sprintf(
  "[missforest] done: imputed %d cells (OOB NRMSE = %.4f) -> DAList_imputed_missforest.rds\n",
  sum(is.na(mat)), unname(mf$OOBerror["NRMSE"])
))
