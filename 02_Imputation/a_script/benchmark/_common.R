# _common.R — shared constants and data loading for benchmark
# Source this at the top of every method and comparison script

library(dplyr)
select <- dplyr::select  # prevent AnnotationDbi masking

setwd(rprojroot::find_rstudio_root_file())

# --- Paths ---
NORM_CSV  <- "01_normalization/c_data/02_normalized.csv"
DAL_RDS   <- "01_normalization/c_data/03_DAList_normalized.rds"
MNAR_CSV  <- "02_Imputation/c_data/02_mar_mnar_classification.csv"
MASK_CSV  <- "02_Imputation/c_data/07_imputation_mask.csv"
BENCH_DIR <- "02_Imputation/c_data/benchmark"
CACHE_RDS <- file.path(BENCH_DIR, "imputed_matrices.rds")

# --- Load data ---
dal  <- readRDS(DAL_RDS)
meta <- dal$metadata

# IMPORTANT: 02_normalized.csv has annotation columns (uniprot_id, protein, gene, description).
# Must drop them to get a numeric matrix.
raw <- read.csv(NORM_CSV, check.names = FALSE)
annot_cols <- c("uniprot_id", "protein", "gene", "description")
num_cols <- setdiff(names(raw), annot_cols)
norm_mat <- apply(raw[, num_cols], 2, as.numeric)
rownames(norm_mat) <- raw$uniprot_id
norm_mat <- norm_mat[, meta$Col_ID]  # ensure column order matches metadata

cat(sprintf("Loaded: %d proteins x %d samples, %.1f%% missing\n",
            nrow(norm_mat), ncol(norm_mat),
            100 * sum(is.na(norm_mat)) / length(norm_mat)))

# --- MAR/MNAR classifiers ---
# Classifier 1: K-means (current pipeline default)
mnar_df <- read.csv(MNAR_CSV)
is_mnar_km <- setNames(mnar_df$classification == "MNAR", mnar_df$uniprot_id)
is_mnar_km <- is_mnar_km[rownames(norm_mat)]

# Classifier 2: Logistic regression P(missing | intensity)
obs_mean <- rowMeans(norm_mat, na.rm = TRUE)
miss_vec <- as.vector(is.na(norm_mat))
int_vec  <- rep(obs_mean, ncol(norm_mat))
logistic_fit <- glm(miss_vec ~ int_vec, family = binomial)
pred_prob <- predict(logistic_fit, newdata = data.frame(int_vec = obs_mean), type = "response")
incomplete <- rowSums(is.na(norm_mat)) > 0
median_threshold <- median(pred_prob[incomplete])
is_mnar_logistic <- pred_prob > median_threshold & incomplete
is_mnar_logistic[!incomplete] <- FALSE
names(is_mnar_logistic) <- rownames(norm_mat)

CLASSIFIERS <- list(km = is_mnar_km, logistic = is_mnar_logistic)

# Report classifier agreement
both_incomplete <- incomplete & !is.na(is_mnar_km)
agree <- mean(is_mnar_km[both_incomplete] == is_mnar_logistic[both_incomplete])
cat(sprintf("Classifiers: km (%d MNAR) vs logistic (%d MNAR), %.1f%% agreement\n",
            sum(is_mnar_km, na.rm = TRUE), sum(is_mnar_logistic),
            agree * 100))

# --- Method registry ---
HYBRID_METHODS <- c("BPCA_QRILC", "KNN_QRILC", "imp4p_mixed", "msImpute_v2_mnar")

BASE_METHODS <- c(
  "BPCA_QRILC", "KNN_QRILC", "MinProb", "MICE", "ProJect",
  "missForest", "imp4p_RF", "imp4p_mixed", "GSimp", "DreamAI",
  "msImpute_v2", "msImpute_v2_mnar", "MAI", "Half_minimum",
  "imputePCA", "Non_imputed",
  "MinDet", "QRILC", "zero", "knn_standalone", "bpca_standalone",
  "RF_MsCoreUtils", "SVD"
)

METHOD_CLASS <- c(
  BPCA_QRILC = "Hybrid", KNN_QRILC = "Hybrid", MinProb = "MNAR",
  MICE = "MAR", ProJect = "Model-free", missForest = "MAR",
  imp4p_RF = "MAR", imp4p_mixed = "Hybrid", GSimp = "MNAR",
  DreamAI = "Ensemble", msImpute_v2 = "MAR", msImpute_v2_mnar = "Hybrid",
  MAI = "Hybrid", Half_minimum = "MNAR", imputePCA = "MAR",
  Non_imputed = "Reference",
  MinDet = "MNAR", QRILC = "MNAR", zero = "Baseline",
  knn_standalone = "MAR", bpca_standalone = "MAR",
  RF_MsCoreUtils = "MAR", SVD = "MAR"
)

# --- Ensure output directory ---
dir.create(BENCH_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BENCH_DIR, "figures"), showWarnings = FALSE)

# --- Helper: find and source method function ---
# Given a method name like "BPCA_QRILC_km", find the script, source it,
# and return the impute_* function name.
find_method_fn <- function(mname) {
  method_files <- sort(list.files("02_Imputation/a_script/benchmark/methods",
                                   pattern = "\\.R$", full.names = TRUE))
  # Split off classifier suffix
  parts <- strsplit(mname, "_(?=(km|logistic)$)", perl = TRUE)[[1]]
  if (length(parts) == 2) {
    base_method <- parts[1]
    classifier <- parts[2]
  } else {
    base_method <- mname
    classifier <- "km"
  }

  # Try exact function name first: impute_<base_method>
  target_fn <- paste0("impute_", base_method)

  # Source all method files if needed, then look for the function
  for (f in method_files) {
    source(f, local = FALSE)
  }

  # Direct match
  if (exists(target_fn, envir = .GlobalEnv)) {
    return(list(fn_name = target_fn, classifier = classifier))
  }

  # Fuzzy match: strip underscores and compare lowercase
  base_lower <- tolower(gsub("_", "", base_method))
  all_fns <- grep("^impute_", ls(envir = .GlobalEnv), value = TRUE)
  for (fn in all_fns) {
    fn_base <- sub("^impute_", "", fn)
    if (tolower(gsub("_", "", fn_base)) == base_lower) {
      return(list(fn_name = fn, classifier = classifier))
    }
  }

  NULL
}
