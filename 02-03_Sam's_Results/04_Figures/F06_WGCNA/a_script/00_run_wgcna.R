# Sam-parallel WGCNA runner.
# Reads Sam's N=35 normalized DAList, runs blockwiseModules (signed Pearson),
# computes module eigengenes, module-trait Pearson correlations, and LMM
# contrasts (eigengene ~ group + (1|subject), Kenward-Roger df).
# Saves all downstream artifacts consumed by 01_main_panels.R and 02_supp_panels.R.
#
# Run from A_CvH_2026/ root:
#   Rscript "02-03_Sam's_Results/04_Figures/F06_WGCNA/a_script/00_run_wgcna.R"

setwd(rprojroot::find_rstudio_root_file())

suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
  library(lme4)
  library(emmeans)
})

source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

allowWGCNAThreads()
set.seed(42)

# ── Paths ────────────────────────────────────────────────────────────────────

SAM_RDS   <- "02-03_Sam's_Results/00_input/01_normalized_DAList_SURV_stringent_muscle.RDS"
DEP_CRvH  <- "02-03_Sam's_Results/03_DEP/c_data/01_limma_DAList_CRvH.rds"
DEP_CR    <- "02-03_Sam's_Results/03_DEP/c_data/01_limma_DAList_CR.rds"

BASE   <- "02-03_Sam's_Results/04_Figures/F06_WGCNA"
WGCNA_DIR  <- file.path(BASE, "c_data", "wgcna")
PANEL_DIR  <- file.path(BASE, "c_data")
SUPP_DIR   <- file.path(BASE, "b_reports", "supp", "png", "panels")
SUPP_PDF   <- file.path(BASE, "b_reports", "supp", "pdf", "panels")

for (d in c(WGCNA_DIR, PANEL_DIR, SUPP_DIR, SUPP_PDF)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

stopifnot(file.exists(SAM_RDS))

# ── Load Sam's DAList ─────────────────────────────────────────────────────────

sam   <- readRDS(SAM_RDS)
mat   <- sam$data                           # 1944 proteins x 35 samples
ann   <- as.data.frame(sam$annotation)[, c("uniprot_id", "protein", "gene", "description")]
md    <- as.data.frame(sam$metadata)        # 35 x 45

message(sprintf("Sam's data: %d proteins x %d samples", nrow(mat), ncol(mat)))

# ── Build metadata table ──────────────────────────────────────────────────────

meta <- tibble(
  sample_id  = md$sample_id,
  pid        = md$pid,
  cancer     = md$cancer,          # SURV / CTL
  supp       = md$supp,            # CRE / PLA / NA
  timepoint  = md$timepoint,       # T1 / T2
  group      = dplyr::case_when(
    md$cancer == "SURV" & md$supp == "CRE" & md$timepoint == "T1" ~ "CRE_T1",
    md$cancer == "SURV" & md$supp == "CRE" & md$timepoint == "T2" ~ "CRE_T2",
    md$cancer == "SURV" & md$supp == "PLA" & md$timepoint == "T1" ~ "PLA_T1",
    md$cancer == "SURV" & md$supp == "PLA" & md$timepoint == "T2" ~ "PLA_T2",
    md$cancer == "CTL"                                             ~ "H_T1",
    TRUE ~ NA_character_
  ),
  age        = md$age,
  sex        = md$sex,
  B_M_ratio  = md$B_M_ratio
)

# Per-sample phenotype columns (pre_*/post_*) already in metadata wide format.
# We carry them for module-trait correlation.
pheno_cols_pre <- c("pre_BMI", "pre_LBM_kg", "pre_ALM_kg",
                    "pre_sts_max_pwr", "pre_grip_lbs",
                    "pre_chest_press_lbs", "pre_leg_ext_lbs")
pheno_cols_post <- c("post_BMI", "post_LBM_kg", "post_ALM_kg",
                     "post_sts_max_pwr", "post_grip_lbs",
                     "post_chest_press_lbs", "post_leg_ext_lbs")
for (pc in c(pheno_cols_pre, pheno_cols_post, "B_M_ratio")) {
  if (pc %in% names(md)) meta[[pc]] <- md[[pc]]
}

# ── WGCNA expression matrix (samples x proteins) ─────────────────────────────

datExpr <- t(mat)                           # 35 x 1944
rownames(datExpr) <- colnames(mat)

# goodSamplesGenes check
cor <- WGCNA::cor
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  ann <- ann |> dplyr::filter(uniprot_id %in% colnames(datExpr))
  message(sprintf("After goodSamplesGenes: %d samples x %d proteins",
                  nrow(datExpr), ncol(datExpr)))
}

# ── Soft-power selection ──────────────────────────────────────────────────────

powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers,
                         networkType = "signed", verbose = 2)
saveRDS(sft$fitIndices, file.path(WGCNA_DIR, "sft_fitIndices.rds"))

r2_values <- -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq
power_idx <- which(r2_values > 0.87)[1]
soft_power <- if (!is.na(power_idx)) powers[power_idx] else {
  # elbow fallback: first power where R^2 > 0.80
  p80 <- which(r2_values > 0.80)[1]
  if (!is.na(p80)) powers[p80] else 6L
}

sft_slope <- sft$fitIndices$slope[soft_power]
r2_note <- paste0(
  ifelse(is.na(sft_slope) || sft_slope > -1,
         " [NOTE: slope > -1, weak scale-free fit]", ""),
  ifelse(r2_values[soft_power] < 0.90,
         " [below 0.90 convention; acceptable for small-n, see Langfelder & Horvath 2008]", ""))
message(sprintf("Soft power: %d  (R^2 = %.3f, slope = %.2f%s)",
                soft_power, r2_values[soft_power],
                if (is.na(sft_slope)) NaN else sft_slope, r2_note))

# Soft-power plot saved to supp
png(file.path(SUPP_DIR, "SUPP_soft_threshold.png"),
    width = 3000, height = 1500, res = 300)
par(mfrow = c(1, 2))
plot(sft$fitIndices$Power, r2_values,
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit (R^2)",
     main = "Scale independence", type = "n")
text(sft$fitIndices$Power, r2_values, labels = powers, cex = 0.9, col = "red")
abline(h = 0.85, col = "red", lty = 2)
abline(v = soft_power, col = "blue", lty = 3)
plot(sft$fitIndices$Power, sft$fitIndices$mean.k.,
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     main = "Mean connectivity", type = "n")
text(sft$fitIndices$Power, sft$fitIndices$mean.k., labels = powers, cex = 0.9, col = "red")
abline(v = soft_power, col = "blue", lty = 3)
dev.off()

# ── blockwiseModules ──────────────────────────────────────────────────────────

net <- blockwiseModules(
  datExpr,
  power             = soft_power,
  networkType       = "signed",
  TOMType           = "signed",
  corType           = "pearson",
  minModuleSize     = 30,
  mergeCutHeight    = 0.25,
  deepSplit         = 2,
  numericLabels     = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs          = FALSE,
  verbose           = 3
)

module_colors <- labels2colors(net$colors)
n_modules <- length(unique(net$colors[net$colors != 0]))
message(sprintf("Modules detected: %d (+grey/unassigned)", n_modules))

# Module sizes (excluding grey)
mod_sizes <- sort(table(module_colors[module_colors != "grey"]), decreasing = TRUE)
message("Top 5 modules:")
print(head(mod_sizes, 5))

# ── Dendrogram plot ───────────────────────────────────────────────────────────

png(file.path(SUPP_DIR, "SUPP_dendrogram.png"),
    width = 3000, height = 1800, res = 300)
plotDendroAndColors(net$dendrograms[[1]],
                    module_colors[net$blockGenes[[1]]],
                    "Module Colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = sprintf("Sam CvH — WGCNA Dendrogram (n=35, power=%d)", soft_power))
dev.off()

# ── Module eigengenes ─────────────────────────────────────────────────────────

MEs <- moduleEigengenes(datExpr, colors = module_colors)$eigengenes
MEs <- orderMEs(MEs)

# ── Simple module-trait Pearson correlations ──────────────────────────────────
# Traits: design (cancer, supp, timepoint encoded numeric) + phenotype deltas

# Encode design traits numerically
meta <- meta |>
  mutate(
    cancer_num    = if_else(cancer == "SURV", 1L, 0L),
    supp_num      = case_when(
      supp == "CRE" ~ 1L,
      supp == "PLA" ~ 0L,
      TRUE          ~ NA_integer_),
    timepoint_num = if_else(timepoint == "T2", 1L, 0L)
  )

# Compute pre/post deltas per subject (only for SURV who have T1+T2)
# subject_key = strip _T# suffix from sample_id
meta <- meta |>
  mutate(subject_key = sub("_T[12]$", "", sample_id))

delta_pheno_vars <- c("BMI", "LBM_kg", "ALM_kg", "sts_max_pwr",
                      "grip_lbs", "chest_press_lbs", "leg_ext_lbs")

delta_df <- tibble(subject_key = character())
for (v in delta_pheno_vars) {
  pre_col  <- paste0("pre_", v)
  post_col <- paste0("post_", v)
  if (pre_col %in% names(meta) && post_col %in% names(meta)) {
    tmp <- meta |>
      dplyr::select(subject_key, pid, all_of(pre_col), all_of(post_col)) |>
      distinct() |>
      mutate(delta = .data[[post_col]] - .data[[pre_col]]) |>
      dplyr::select(subject_key, delta) |>
      dplyr::rename(!!paste0("delta_", v) := delta)
    delta_df <- if (nrow(delta_df) == 0) tmp else
      left_join(delta_df, tmp, by = "subject_key")
  }
}

# Merge deltas back to per-sample meta
if (ncol(delta_df) > 1) {
  meta <- left_join(meta, delta_df, by = "subject_key")
}

# Build trait matrix (one row per sample, aligned to datExpr rows)
trait_cols_design <- c("cancer_num", "supp_num", "timepoint_num", "B_M_ratio", "age")
trait_cols_pheno  <- paste0("delta_", delta_pheno_vars)
trait_cols_pheno  <- intersect(trait_cols_pheno, names(meta))
all_trait_cols    <- c(trait_cols_design, trait_cols_pheno)
all_trait_cols    <- intersect(all_trait_cols, names(meta))

traits_mat <- meta |>
  dplyr::select(sample_id, all_of(all_trait_cols)) |>
  column_to_rownames("sample_id") |>
  mutate(across(everything(), as.numeric))
traits_mat <- traits_mat[rownames(datExpr), , drop = FALSE]

n_per_trait <- colSums(!is.na(traits_mat))
module_trait_cor <- cor(MEs, traits_mat, use = "pairwise.complete.obs")

module_trait_pval <- module_trait_cor
for (j in seq_len(ncol(module_trait_cor))) {
  module_trait_pval[, j] <- corPvalueStudent(module_trait_cor[, j], n_per_trait[j])
}
pval_vec    <- as.vector(module_trait_pval)
pval_bh_vec <- p.adjust(pval_vec, method = "BH")
module_trait_pval_bh <- matrix(pval_bh_vec,
                                nrow = nrow(module_trait_pval),
                                ncol = ncol(module_trait_pval),
                                dimnames = dimnames(module_trait_pval))

# Save module-trait correlations
write_csv(
  as.data.frame(module_trait_cor) |> rownames_to_column("module"),
  file.path(WGCNA_DIR, "wgcna_module_trait_correlations.csv"))
write_csv(
  as.data.frame(module_trait_pval_bh) |> rownames_to_column("module"),
  file.path(WGCNA_DIR, "wgcna_module_trait_pvalues_bh.csv"))

# Report top 3 module-trait correlations
cor_long <- as.data.frame(module_trait_cor) |>
  rownames_to_column("module") |>
  pivot_longer(-module, names_to = "trait", values_to = "r") |>
  left_join(
    as.data.frame(module_trait_pval_bh) |>
      rownames_to_column("module") |>
      pivot_longer(-module, names_to = "trait", values_to = "p_bh"),
    by = c("module", "trait")) |>
  arrange(p_bh)
message("Top 3 module-trait correlations:")
print(head(cor_long, 3))

# Quick heatmap to SUPP (base R labeledHeatmap, mirrors YvO)
star_matrix <- ifelse(module_trait_pval_bh < 0.001, "***",
               ifelse(module_trait_pval_bh < 0.01,  "**",
               ifelse(module_trait_pval_bh < 0.05,  "*", "")))
text_matrix <- paste(signif(module_trait_cor, 2), star_matrix, sep = "\n")
dim(text_matrix) <- dim(module_trait_cor)

png(file.path(SUPP_DIR, "SUPP_module_trait_heatmap_simple.png"),
    width = 3200, height = 3000, res = 300)
par(mar = c(6, 12, 3, 3))
labeledHeatmap(
  Matrix     = module_trait_cor,
  xLabels    = colnames(traits_mat),
  yLabels    = colnames(MEs),
  ySymbols   = colnames(MEs),
  colorLabels = FALSE,
  colors     = blueWhiteRed(50),
  textMatrix = text_matrix,
  setStdMargins = FALSE,
  cex.text   = 0.5,
  zlim       = c(-1, 1),
  main       = "Sam CvH — Module-trait correlations (* BH < 0.05)"
)
dev.off()

# ── Hub genes (kME) ───────────────────────────────────────────────────────────

kME <- signedKME(datExpr, MEs)

module_df <- tibble(
  uniprot_id   = colnames(datExpr),
  module_color = module_colors,
  module_num   = net$colors
) |> left_join(ann |> dplyr::select(uniprot_id, gene), by = "uniprot_id")

unique_modules <- setdiff(unique(module_colors), "grey")

hub_rows <- list()
for (mod in unique_modules) {
  mod_proteins <- module_df$uniprot_id[module_df$module_color == mod]
  kme_col <- paste0("kME", mod)
  if (!(kme_col %in% colnames(kME))) next
  mod_kme <- tibble(
    uniprot_id = rownames(kME),
    kME        = kME[, kme_col]
  ) |>
    dplyr::filter(uniprot_id %in% mod_proteins) |>
    arrange(desc(abs(kME))) |>
    head(15) |>
    mutate(module = mod) |>
    left_join(ann |> dplyr::select(uniprot_id, gene), by = "uniprot_id")
  hub_rows <- c(hub_rows, list(mod_kme))
}
hub_df <- bind_rows(hub_rows)

# ── ORA per module ────────────────────────────────────────────────────────────

cor <- stats::cor    # restore base cor after WGCNA block

bg_genes <- ann$gene[ann$uniprot_id %in% colnames(datExpr)]
bg_genes <- unique(bg_genes[!is.na(bg_genes) & bg_genes != ""])

pw_collection <- build_pathway_collection(min_size = 15, max_size = 500,
                                          include_goslim = FALSE)

ora_results_list <- list()
for (mod in unique_modules) {
  mod_genes <- module_df$gene[module_df$module_color == mod]
  mod_genes <- unique(mod_genes[!is.na(mod_genes) & mod_genes != ""])
  if (length(mod_genes) < 5) next
  ora_res <- tryCatch(
    run_ora_deduplicated(
      genes          = mod_genes,
      universe       = bg_genes,
      pathways       = pw_collection,
      jaccard_cutoff = 0.5,
      min_size       = 15,
      max_size       = 500,
      padj_cutoff    = 0.10
    ),
    error = function(e) {
      warning(sprintf("ORA failed for '%s': %s", mod, e$message)); NULL
    }
  )
  if (!is.null(ora_res) && nrow(ora_res) > 0) {
    ora_res$module      <- mod
    ora_res$Description <- clean_pathway_name(ora_res$pathway)
    ora_res$geneID      <- vapply(ora_res$overlapGenes,
                                   function(g) paste(g, collapse = "/"), character(1))
    ora_res$Count       <- ora_res$overlap
    ora_res$p.adjust    <- ora_res$padj
    ora_res$ID          <- ora_res$pathway
    ora_results_list    <- c(ora_results_list, list(ora_res))
  }
}
enrich_df <- bind_rows(ora_results_list)
if ("overlapGenes" %in% names(enrich_df)) enrich_df$overlapGenes <- NULL

# ── Module bio labels ─────────────────────────────────────────────────────────

mod_bio_labels <- tibble(
  module_color = names(mod_sizes),
  module_id    = paste0("M", seq_along(mod_sizes)),
  n_proteins   = as.integer(mod_sizes)
) |>
  mutate(
    top_pathway = map_chr(module_color, function(mc) {
      sub_e <- enrich_df |> dplyr::filter(module == mc) |> arrange(padj)
      if (nrow(sub_e) > 0) sub_e$Description[1] else "Uncharacterised"
    }),
    display_label = paste0(module_id, ": ", top_pathway,
                           " (n=", n_proteins, ")")
  )

# ── LMM contrasts (eigengene ~ group + (1|subject)) ──────────────────────────

# All samples with repeated measures. CTL (H_T1) are T1-only; they contribute
# to the intercept but have NA subject random effects (no post measurement).
# Use the same group factor levels as our main CvH pipeline.
lmm_data <- meta |>
  mutate(
    subject = pid,
    group   = factor(group, levels = c("CRE_T1", "CRE_T2",
                                        "PLA_T1", "PLA_T2",
                                        "H_T1"))
  )

# Two-model contrast sets (matching main CvH DEP contrasts)
lmm_contrast_list <- list(
  Cancer_vs_Healthy   = c(-0.5, -0.5, -0.5, -0.5, 2) / 2,
  # (CRE_T1 + PLA_T1)/2 - H_T1 normalised to sum-to-zero equivalent
  Training_CR         = c(-1, 1, -1, 1, 0) / 2
)

# CR-only contrasts (subset to SURV samples)
cr_contrasts <- list(
  Baseline_Supplement    = c(-1, -1,  1,  1,  0) / 2,  # CRE_T1 - PLA_T1
  Training_CRE           = c(-1,  1,  0,  0,  0),       # CRE_T2 - CRE_T1
  Training_PLA           = c( 0,  0, -1,  1,  0),       # PLA_T2 - PLA_T1
  Supplement_Interaction = c( 1, -1, -1,  1,  0) / 2   # (CRE_T2-CRE_T1) - (PLA_T2-PLA_T1)
)

all_mods <- colnames(MEs)

lmm_rows <- list()
for (mod in all_mods) {
  lmm_data[[mod]] <- MEs[lmm_data$sample_id, mod]

  # Full model (all 5 groups)
  fit_full <- tryCatch(
    suppressWarnings(
      lmer(as.formula(paste0("`", mod, "` ~ group + (1 | subject)")),
           data = lmm_data, REML = TRUE)
    ),
    error = function(e) { warning(sprintf("LMM failed for %s: %s", mod, e$message)); NULL }
  )
  if (!is.null(fit_full)) {
    singular_full <- isSingular(fit_full)
    emm_full <- emmeans(fit_full, ~ group)
    for (cname in names(lmm_contrast_list)) {
      ctr <- contrast(emm_full, list(ctr = lmm_contrast_list[[cname]]))
      s   <- summary(ctr, ddf = "Kenward-Roger")
      t_v <- s$t.ratio; df_v <- s$df
      r_e <- sign(s$estimate) * sqrt(t_v^2 / (t_v^2 + df_v))
      lmm_rows <- c(lmm_rows, list(tibble(
        module = mod, contrast = cname, model = "CRvH",
        estimate = round(s$estimate, 5), SE = round(s$SE, 5),
        df = round(df_v, 2), t_ratio = round(t_v, 4),
        p_raw = s$p.value, r_equiv = round(r_e, 4), singular = singular_full
      )))
    }
  }

  # CR-only model (SURV samples only)
  cr_data <- lmm_data |>
    dplyr::filter(cancer == "SURV") |>
    mutate(group = droplevels(group))

  fit_cr <- tryCatch(
    suppressWarnings(
      lmer(as.formula(paste0("`", mod, "` ~ group + (1 | subject)")),
           data = cr_data, REML = TRUE)
    ),
    error = function(e) NULL
  )
  if (!is.null(fit_cr)) {
    singular_cr <- isSingular(fit_cr)
    emm_cr <- emmeans(fit_cr, ~ group)
    cr_levels <- levels(cr_data$group)
    for (cname in names(cr_contrasts)) {
      ctr_vec_full <- cr_contrasts[[cname]]
      # Subset to the CR group levels (drop H_T1)
      idx_in_cr <- which(c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1") %in% cr_levels)
      ctr_vec   <- ctr_vec_full[idx_in_cr]
      if (length(ctr_vec) != length(cr_levels)) next
      ctr <- contrast(emm_cr, list(ctr = ctr_vec))
      s   <- summary(ctr, ddf = "Kenward-Roger")
      t_v <- s$t.ratio; df_v <- s$df
      r_e <- sign(s$estimate) * sqrt(t_v^2 / (t_v^2 + df_v))
      lmm_rows <- c(lmm_rows, list(tibble(
        module = mod, contrast = cname, model = "CR",
        estimate = round(s$estimate, 5), SE = round(s$SE, 5),
        df = round(df_v, 2), t_ratio = round(t_v, 4),
        p_raw = s$p.value, r_equiv = round(r_e, 4), singular = singular_cr
      )))
    }
  }
}

lmm_df <- bind_rows(lmm_rows)
lmm_df$p_bh <- p.adjust(lmm_df$p_raw, method = "BH")
message(sprintf("LMM contrasts: %d tests", nrow(lmm_df)))
print(lmm_df |> arrange(p_raw) |> head(5))

# ── SFT summary ───────────────────────────────────────────────────────────────

sft_summary <- tibble(
  selected_power    = soft_power,
  R_squared         = r2_values[soft_power],
  mean_connectivity = sft$fitIndices$mean.k.[soft_power],
  n_proteins        = ncol(datExpr),
  n_samples         = nrow(datExpr)
)
write_csv(sft_summary, file.path(WGCNA_DIR, "wgcna_sft_summary.csv"))

# ── Save WGCNA result objects ─────────────────────────────────────────────────

saveRDS(net,        file.path(WGCNA_DIR, "wgcna_network.rds"))
write_csv(module_df, file.path(WGCNA_DIR, "wgcna_module_assignments.csv"))
write_csv(hub_df,    file.path(WGCNA_DIR, "wgcna_hub_proteins.csv"))
write_csv(enrich_df, file.path(WGCNA_DIR, "wgcna_module_enrichment.csv"))
write_csv(lmm_df,    file.path(WGCNA_DIR, "wgcna_lmm_contrasts.csv"))
write_csv(mod_bio_labels, file.path(WGCNA_DIR, "wgcna_mod_bio_labels.csv"))

# Key modules (top by enrichment path count)
key_mod_counts <- enrich_df |>
  dplyr::count(module, sort = TRUE) |>
  head(5) |>
  dplyr::pull(module)
if (length(key_mod_counts) == 0) key_mod_counts <- names(mod_sizes)[1:min(3, length(mod_sizes))]
writeLines(key_mod_counts, file.path(WGCNA_DIR, "key_modules.txt"))

# ── Panel-level objects ───────────────────────────────────────────────────────

# Pre/post eigengene matrices for repeated-measures panels
surv_meta <- meta |> dplyr::filter(cancer == "SURV")
pre_meta   <- surv_meta |> dplyr::filter(timepoint == "T1")
post_meta  <- surv_meta |> dplyr::filter(timepoint == "T2")

pre_subj   <- pre_meta$subject_key
post_subj  <- post_meta$subject_key
common_subj <- intersect(pre_subj, post_subj)

me_pre  <- MEs[pre_meta$sample_id[match(common_subj,  pre_meta$subject_key)], , drop = FALSE]
me_post <- MEs[post_meta$sample_id[match(common_subj, post_meta$subject_key)], , drop = FALSE]
rownames(me_pre)  <- common_subj
rownames(me_post) <- common_subj
delta_me <- me_post - me_pre

# z-scored expression per protein (for module triptych supp panel)
expr_g <- t(datExpr)
rownames(expr_g) <- ann$gene[match(rownames(expr_g), ann$uniprot_id)]
expr_g <- expr_g[!is.na(rownames(expr_g)) & rownames(expr_g) != "", ]
expr_g <- expr_g[!duplicated(rownames(expr_g)), ]
z_g    <- t(scale(t(expr_g)))

group_levels <- c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1")
group_z <- vapply(
  group_levels,
  function(g) {
    s <- meta$sample_id[meta$group == g]
    s <- intersect(s, colnames(z_g))
    if (length(s) == 0) return(rep(NA_real_, nrow(z_g)))
    rowMeans(z_g[, s, drop = FALSE], na.rm = TRUE)
  },
  numeric(nrow(z_g))
)
colnames(group_z) <- group_levels

saveRDS(MEs,          file.path(PANEL_DIR, "MEs.rds"))
saveRDS(kME,          file.path(PANEL_DIR, "kME_all.rds"))
saveRDS(datExpr,      file.path(PANEL_DIR, "datExpr.rds"))
saveRDS(module_colors, file.path(PANEL_DIR, "module_colors.rds"))
saveRDS(me_pre,        file.path(PANEL_DIR, "me_pre.rds"))
saveRDS(me_post,       file.path(PANEL_DIR, "me_post.rds"))
saveRDS(delta_me,      file.path(PANEL_DIR, "delta_me.rds"))
saveRDS(group_z,       file.path(PANEL_DIR, "group_z.rds"))

write_csv(meta,           file.path(PANEL_DIR, "meta.csv"))
write_csv(ann,            file.path(PANEL_DIR, "imp_annotations.csv"))
write_csv(mod_bio_labels, file.path(PANEL_DIR, "mod_bio_labels.csv"))
if (length(common_subj) > 0) {
  write_csv(tibble(subject_key = common_subj), file.path(PANEL_DIR, "common_subj.csv"))
}

saveRDS(list(
  common_subj = common_subj,
  pre_subj    = pre_subj,
  post_subj   = post_subj,
  n_modules   = n_modules,
  soft_power  = soft_power,
  mod_sizes   = as.list(mod_sizes),
  key_modules = key_mod_counts
), file.path(PANEL_DIR, "shared_objects.rds"))

message(sprintf("Done: %d modules, %d hub proteins, %d enriched pathways, %d paired subjects",
                n_modules, nrow(hub_df), nrow(enrich_df), length(common_subj)))
