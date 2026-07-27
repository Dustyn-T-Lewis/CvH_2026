# F05 stage 00: build the signed co-abundance network from the imputed DAList.
# Inputs:  02_Normalization/imputation/c_data/DAList_imputed_missforest.rds
# Outputs: 04_Figures/F05_WGCNA/c_data/wgcna_network.rds (consolidated list)
#          plus panel-ready rds/csv in c_data/
#
# Signed network and TOM, minModuleSize 30, mergeCutHeight 0.25 (Cahill 2018 /
# Langfelder & Horvath 2008). Correlation is biweight midcorrelation
# (maxPOutliers 0.05), robust to the outliers and heavy tails proteomics
# abundances carry. Soft power = first power with signed R^2 > 0.90, which also
# satisfies the WGCNA sample-size heuristic (>= 14 for a signed network at
# 30-40 samples) and keeps mean connectivity low. Module-level inference
# (eigengene contrasts, fry, NES, phenotype) lives in 01_module_stats.R.

pacman::p_load(WGCNA, tidyverse)

allowWGCNAThreads()
set.seed(42)
setwd(here::here())

DALIST_RDS <- "02_Normalization/imputation/c_data/DAList_imputed_missforest.rds"
PANEL_DIR <- "04_Figures/F05_WGCNA/c_data"
DATA_DIR <- file.path(PANEL_DIR, "wgcna")
REPORT_DIR <- "04_Figures/F05_WGCNA/b_reports/supp/01_QC"

dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(DALIST_RDS))

dal <- readRDS(DALIST_RDS)
ann_cols <- c("uniprot_id", "protein", "gene", "description")
ann <- as.data.frame(dal$annotation)[, ann_cols]
mat <- as.matrix(dal$data)
rownames(mat) <- ann$uniprot_id

datExpr <- t(mat) # samples as rows, proteins as columns

dal_meta <- as.data.frame(dal$metadata)

meta <- tibble(
  sample_id = dal_meta$Col_ID,
  subject = dal_meta$Subject_ID, # CvH IDs are cohort-unique
  group = dal_meta$Group,
  timepoint = dal_meta$Timepoint,
  time = dal_meta$Timepoint, # alias for downstream supp panels
  group_time = dal_meta$Group_Time,
  supplement = if ("Supplement" %in% names(dal_meta)) dal_meta$Supplement else NA_character_,
  cancer = dal_meta$cancer, # CTL / SURV
  cancer_time = dal_meta$cancer_time # CTL_T1 / SURV_T1 / SURV_T2
)

gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  ann <- ann |> filter(uniprot_id %in% colnames(datExpr))
  message(sprintf(
    "After goodSamplesGenes: %d samples x %d proteins",
    nrow(datExpr), ncol(datExpr)
  ))
}

cor <- WGCNA::cor

RSQ_CUT <- 0.90

# WGCNA's published minimum soft power for a signed network, by sample count. At small n
# the scale-free criterion alone settles too low, so the heuristic acts as a floor.
signed_power_floor <- function(n) {
  if (n < 20) {
    18L
  } else if (n < 30) {
    16L
  } else if (n < 40) {
    14L
  } else {
    12L
  }
}

powers <- 1:20
sft <- pickSoftThreshold(datExpr,
  powerVector = powers, networkType = "signed",
  corFnc = bicor, corOptions = list(maxPOutliers = 0.05), verbose = 2
)
saveRDS(sft$fitIndices, file.path(DATA_DIR, "sft_fitIndices.rds"))

r2_values <- -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq
power_idx <- which(r2_values > RSQ_CUT)[1]
r2_power <- if (!is.na(power_idx)) powers[power_idx] else NA_integer_
soft_power <- max(r2_power, signed_power_floor(nrow(datExpr)), na.rm = TRUE)
message(sprintf(
  "Soft power: %d (R^2 criterion gave %s, signed-network floor for n=%d is %d)",
  soft_power, r2_power, nrow(datExpr), signed_power_floor(nrow(datExpr))
))

png(file.path(REPORT_DIR, "SUPP_soft_threshold.png"),
  width = 3000, height = 1500, res = 300
)
par(mfrow = c(1, 2))
plot(sft$fitIndices$Power, r2_values,
  xlab = "Soft Threshold (power)",
  ylab = expression(paste("Scale Free Topology Model Fit (", R^2, ")")),
  main = "Scale independence", type = "n"
)
text(sft$fitIndices$Power, r2_values, labels = powers, cex = 0.9, col = "red")
abline(h = RSQ_CUT, col = "red", lty = 2)
plot(sft$fitIndices$Power, sft$fitIndices$mean.k.,
  xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
  main = "Mean connectivity", type = "n"
)
text(sft$fitIndices$Power, sft$fitIndices$mean.k.,
  labels = powers,
  cex = 0.9, col = "red"
)
dev.off()

net <- blockwiseModules(
  datExpr,
  power             = soft_power,
  networkType       = "signed",
  TOMType           = "signed",
  corType           = "bicor",
  maxPOutliers      = 0.05,
  minModuleSize     = 30,
  mergeCutHeight    = 0.25,
  numericLabels     = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs          = FALSE,
  verbose           = 3
)

module_colors <- labels2colors(net$colors)
n_modules <- length(unique(net$colors)) - (0 %in% net$colors)
message(sprintf("Modules detected: %d (+ grey/unassigned)", n_modules))

MEs <- moduleEigengenes(datExpr, colors = module_colors)$eigengenes
MEs <- orderMEs(MEs)
kME <- signedKME(datExpr, MEs, corFnc = "bicor", corOptions = "maxPOutliers = 0.05")

module_df <- tibble(
  uniprot_id   = colnames(datExpr),
  module_color = module_colors,
  module_num   = net$colors
) |>
  left_join(ann |> dplyr::select(uniprot_id, gene), by = "uniprot_id")

mod_sizes <- sort(table(module_colors[module_colors != "grey"]),
  decreasing = TRUE
)
mod_bio_labels <- tibble(
  module_color = names(mod_sizes),
  module_id = paste0("M", seq_along(mod_sizes)),
  bio_label = NA_character_,
  n_proteins = as.integer(mod_sizes),
  display_label = paste0(
    "M", seq_along(mod_sizes), ": ",
    names(mod_sizes), " (n=", as.integer(mod_sizes), ")"
  )
)

cor <- stats::cor

sft_summary <- tibble(
  selected_power    = soft_power,
  R_squared         = r2_values[soft_power],
  mean_connectivity = sft$fitIndices$mean.k.[soft_power],
  n_proteins        = ncol(datExpr),
  n_samples         = nrow(datExpr)
)


imp_mat <- t(datExpr)
saveRDS(MEs, file.path(PANEL_DIR, "MEs.rds"))
saveRDS(datExpr, file.path(PANEL_DIR, "datExpr.rds"))
saveRDS(imp_mat, file.path(PANEL_DIR, "imp_mat.rds"))
saveRDS(module_colors, file.path(PANEL_DIR, "module_colors.rds"))
write_csv(module_df, file.path(PANEL_DIR, "wgcna_module_assignments.csv"))
write_csv(mod_bio_labels, file.path(PANEL_DIR, "mod_bio_labels.csv"))
write_csv(meta, file.path(PANEL_DIR, "meta.csv"))
write_csv(ann, file.path(PANEL_DIR, "imp_annotations.csv"))

# Trajectory table: Ctl + pooled CR (T1, T2) group-mean eigengenes for the card.
traj_summary <- as.data.frame(MEs) |>
  tibble::rownames_to_column("sample_id") |>
  pivot_longer(-sample_id, names_to = "module", values_to = "eigengene") |>
  left_join(meta[, c("sample_id", "cancer_time")], by = "sample_id") |>
  mutate(traj = recode(cancer_time,
    CTL_T1 = "Ctl", SURV_T1 = "CR_T1", SURV_T2 = "CR_T2"
  )) |>
  group_by(module, traj) |>
  summarise(
    mean_eig = mean(eigengene),
    se       = sd(eigengene) / sqrt(dplyr::n()),
    n        = dplyr::n(),
    .groups  = "drop"
  ) |>
  mutate(module_color = sub("^ME", "", module))
write_csv(traj_summary, file.path(PANEL_DIR, "trajectory_eigengenes.csv"))


# Consolidated network object mirroring MITO's wgcna_network.rds list.
w <- list(
  net = net, module_colors = module_colors, MEs = MEs, kME = kME,
  mod_bio_labels = mod_bio_labels, module_df = module_df,
  traj_summary = traj_summary, sft_summary = sft_summary,
  chosen_power = soft_power, meta = meta, ann = ann
)
saveRDS(w, file.path(PANEL_DIR, "wgcna_network.rds"))

message(sprintf(
  "Done: %d modules x %d proteins x %d samples (power = %d)",
  n_modules, ncol(datExpr), nrow(datExpr), soft_power
))
