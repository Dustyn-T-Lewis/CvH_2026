# CvH WGCNA Runner — produces network and module artifacts for F06
# Inputs:  02_Normalization/imputation/c_data/DAList_imputed_imp4p.rds
# Outputs: 04_Figures/F06/c_data/wgcna/  (network, sft_summary)
#          04_Figures/F06/c_data/        (panel-ready: MEs, kME, module_colors, datExpr,
#                                         module_assignments, mod_bio_labels)
#
# Mirrors YvO_WGCNA_run.R parameters (Cahill 2018 / Langfelder & Horvath 2008):
#   networkType = "signed", TOMType = "signed",
#   minModuleSize = 30, mergeCutHeight = 0.25,
#   pickSoftThreshold over powers 1:20, R^2 > 0.87 cutoff (default 6 fallback).

suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
})

allowWGCNAThreads()
set.seed(42)
setwd(here::here())

DALIST_RDS <- "02_Normalization/imputation/c_data/DAList_imputed_imp4p.rds"
PANEL_DIR  <- "04_Figures/F06/c_data"
DATA_DIR   <- file.path(PANEL_DIR, "wgcna")
REPORT_DIR <- "04_Figures/F06/b_reports/supp/01_QC"

dir.create(PANEL_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR,   recursive = TRUE, showWarnings = FALSE)
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(DALIST_RDS))

dal        <- readRDS(DALIST_RDS)
ann_cols   <- c("uniprot_id", "protein", "gene", "description")
ann        <- as.data.frame(dal$annotation)[, ann_cols]
mat        <- as.matrix(dal$data)
rownames(mat) <- ann$uniprot_id

datExpr <- t(mat)  # samples as rows, proteins as columns

dal_meta <- as.data.frame(dal$metadata)

meta <- tibble(
  sample_id  = dal_meta$Col_ID,
  subject    = dal_meta$Subject_ID,         # CvH IDs are cohort-unique
  group      = dal_meta$Group,
  timepoint  = dal_meta$Timepoint,
  time       = dal_meta$Timepoint,          # alias for downstream supp panels
  group_time = dal_meta$Group_Time,
  supplement = if ("Supplement" %in% names(dal_meta)) dal_meta$Supplement else NA_character_
)

gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  ann <- ann |> filter(uniprot_id %in% colnames(datExpr))
  message(sprintf("After goodSamplesGenes: %d samples x %d proteins",
                  nrow(datExpr), ncol(datExpr)))
}

cor <- WGCNA::cor

powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                         networkType = "signed", verbose = 2)
saveRDS(sft$fitIndices, file.path(DATA_DIR, "sft_fitIndices.rds"))

# Soft power: first power with signed R^2 > 0.87 (Langfelder & Horvath 2008
# small-n guidance; same as YvO). Default 6 fallback if no power qualifies.
r2_values <- -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq
power_idx <- which(r2_values > 0.87)[1]
soft_power <- if (!is.na(power_idx)) powers[power_idx] else 6L
message(sprintf("Soft power: %d (R^2 = %.3f)", soft_power, r2_values[soft_power]))

png(file.path(REPORT_DIR, "SUPP_soft_threshold.png"),
    width = 3000, height = 1500, res = 300)
par(mfrow = c(1, 2))
plot(sft$fitIndices$Power, r2_values,
     xlab = "Soft Threshold (power)",
     ylab = expression(paste("Scale Free Topology Model Fit (", R^2, ")")),
     main = "Scale independence", type = "n")
text(sft$fitIndices$Power, r2_values, labels = powers, cex = 0.9, col = "red")
abline(h = 0.85, col = "red", lty = 2)
plot(sft$fitIndices$Power, sft$fitIndices$mean.k.,
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     main = "Mean connectivity", type = "n")
text(sft$fitIndices$Power, sft$fitIndices$mean.k., labels = powers,
     cex = 0.9, col = "red")
dev.off()

net <- blockwiseModules(
  datExpr,
  power             = soft_power,
  networkType       = "signed",
  TOMType           = "signed",
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
kME <- signedKME(datExpr, MEs)

module_df <- tibble(
  uniprot_id   = colnames(datExpr),
  module_color = module_colors,
  module_num   = net$colors
) |>
  left_join(ann |> dplyr::select(uniprot_id, gene), by = "uniprot_id")

mod_sizes <- sort(table(module_colors[module_colors != "grey"]),
                  decreasing = TRUE)
mod_bio_labels <- tibble(
  module_color  = names(mod_sizes),
  module_id     = paste0("M", seq_along(mod_sizes)),
  bio_label     = NA_character_,
  n_proteins    = as.integer(mod_sizes),
  display_label = paste0("M", seq_along(mod_sizes), ": ",
                          names(mod_sizes), " (n=", as.integer(mod_sizes), ")")
)

# Restore stats::cor before any non-WGCNA correlation calls downstream
cor <- stats::cor

# --- LMM contrasts: per-model audit consumed by F08/panel_A.R ---
# CvH fits two models (per 03_DEP): CRvH (all 5 group_time levels) and CR
# (cancer-only 4 levels). Each module's eigengene is fit ME ~ group_time +
# (1|subject); contrasts mirror the DEP design exactly.
fit_lmm_contrasts <- function(model_label, gt_levels, contrasts_named) {
  lmm_meta <- meta |>
    filter(group_time %in% gt_levels) |>
    mutate(group_time = factor(group_time, levels = gt_levels))
  if (nrow(lmm_meta) == 0) return(tibble())
  rows <- list()
  for (mod in colnames(MEs)) {
    lmm_meta[[mod]] <- MEs[lmm_meta$sample_id, mod]
    fit <- tryCatch(suppressWarnings(
      lme4::lmer(as.formula(paste0("`", mod, "` ~ group_time + (1 | subject)")),
                 data = lmm_meta, REML = TRUE)
    ), error = function(e) NULL)
    if (is.null(fit)) next
    sing <- lme4::isSingular(fit)
    emm  <- emmeans::emmeans(fit, ~ group_time)
    for (cname in names(contrasts_named)) {
      cv  <- contrasts_named[[cname]]
      ctr <- emmeans::contrast(emm, list(ctr = cv))
      s   <- summary(ctr, ddf = "Kenward-Roger")
      r_eq <- sign(s$estimate) * sqrt(s$t.ratio^2 / (s$t.ratio^2 + s$df))
      rows <- c(rows, list(tibble(
        module   = mod,
        model    = model_label,
        contrast = cname,
        estimate = round(s$estimate, 5),
        SE       = round(s$SE, 5),
        df       = round(s$df, 2),
        t_ratio  = round(s$t.ratio, 4),
        p_raw    = s$p.value,
        r_equiv  = round(r_eq, 4),
        singular = sing
      )))
    }
  }
  bind_rows(rows)
}

crvh_levels <- c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1")
crvh_levels <- intersect(crvh_levels, unique(meta$group_time))
crvh_contrasts <- list(
  Cancer_vs_Healthy = c(0.5,    0, 0.5,   0, -1),
  Training_CR       = c(-0.5, 0.5, -0.5, 0.5,  0)
)
# Trim contrast vectors to match available levels (in case H_T1 is missing)
crvh_contrasts <- lapply(crvh_contrasts, function(v) v[seq_along(crvh_levels)])

lmm_df <- if (length(crvh_levels) >= 4)
  fit_lmm_contrasts("CRvH", crvh_levels, crvh_contrasts) else tibble()
if (nrow(lmm_df) > 0) {
  lmm_df$p_bh <- p.adjust(lmm_df$p_raw, method = "BH")
  write_csv(lmm_df, file.path(DATA_DIR, "wgcna_lmm_contrast_audit.csv"))
  write_csv(lmm_df, file.path(PANEL_DIR, "wgcna_lmm_contrast_audit.csv"))
  message(sprintf("LMM audit: %d tests across %d modules",
                  nrow(lmm_df), length(colnames(MEs))))
}

sft_summary <- tibble(
  selected_power    = soft_power,
  R_squared         = r2_values[soft_power],
  mean_connectivity = sft$fitIndices$mean.k.[soft_power],
  n_proteins        = ncol(datExpr),
  n_samples         = nrow(datExpr)
)

# WGCNA network artifacts (DATA_DIR)
saveRDS(net,         file.path(DATA_DIR, "wgcna_network.rds"))
write_csv(sft_summary, file.path(DATA_DIR, "wgcna_sft_summary.csv"))
write_csv(module_df,   file.path(DATA_DIR, "wgcna_module_assignments.csv"))

# Panel-ready artifacts (PANEL_DIR)
imp_mat <- t(datExpr)
saveRDS(MEs,           file.path(PANEL_DIR, "MEs.rds"))
saveRDS(kME,           file.path(PANEL_DIR, "kME_all.rds"))
saveRDS(datExpr,       file.path(PANEL_DIR, "datExpr.rds"))
saveRDS(imp_mat,       file.path(PANEL_DIR, "imp_mat.rds"))
saveRDS(module_colors, file.path(PANEL_DIR, "module_colors.rds"))
write_csv(module_df,   file.path(PANEL_DIR, "wgcna_module_assignments.csv"))
write_csv(mod_bio_labels, file.path(PANEL_DIR, "mod_bio_labels.csv"))
write_csv(meta,        file.path(PANEL_DIR, "meta.csv"))
write_csv(ann,         file.path(PANEL_DIR, "imp_annotations.csv"))

# Top 5 modules by size as default key set
key_modules <- head(mod_bio_labels$module_color, 5)
writeLines(key_modules, file.path(PANEL_DIR, "key_modules.txt"))

message(sprintf("Done: %d modules x %d proteins x %d samples (power = %d)",
                n_modules, ncol(datExpr), nrow(datExpr), soft_power))
