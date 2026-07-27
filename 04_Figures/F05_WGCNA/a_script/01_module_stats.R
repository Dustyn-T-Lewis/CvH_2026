# F05 stage 01: module-level inference on the pooled CR-vs-Ctl design.
#   - eigengene LMM contrasts (r_equiv, feeds the trajectory brackets)
#   - fry (self-contained gate) + camera (reported) per contrast
#   - eigengene omnibus moderated-F (the main-figure gate)
#   - per-module fGSEA NES: module as gene set, genes ranked by DE moderated-t
#   - module-eigengene association with clinical outcomes (three arms)
# Same subject block (duplicateCorrelation / random intercept) throughout.

setwd(here::here())
pacman::p_load(limma, lme4, emmeans, dplyr, tidyr, tibble, readr, purrr)
source("04_Figures/shared/wgcna_stats.R")

DAT <- "04_Figures/F05_WGCNA/c_data"
w <- readRDS(file.path(DAT, "wgcna_network.rds"))
MEs <- w$MEs
meta <- w$meta
expr <- readRDS(file.path(DAT, "imp_mat.rds"))
mc <- readRDS(file.path(DAT, "module_colors.rds"))
names(mc) <- rownames(expr)
module_df <- read_csv(file.path(DAT, "wgcna_module_assignments.csv"), show_col_types = FALSE)

gt_levels <- intersect(
  c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1"), unique(meta$group_time)
)
# Weights mirror 03_DEP exactly: the NES row of each tile is read from the DE table by
# contrast name, so the fry row has to test the same thing or the tile compares two
# different contrasts under one label.
crvh_contrasts <- list(
  CRvH_Baseline = c(CRE_T1 = 0.5, CRE_T2 = 0, PLA_T1 = 0.5, PLA_T2 = 0, H_T1 = -1),
  CR_Training = c(CRE_T1 = -0.5, CRE_T2 = 0.5, PLA_T1 = -0.5, PLA_T2 = 0.5, H_T1 = 0)
)

# Eigengene LMM contrasts: ME ~ group_time + (1 | subject)
fit_lmm_contrasts <- function(gt_levels, contrasts_named) {
  lmm_meta <- meta |>
    filter(group_time %in% gt_levels) |>
    mutate(group_time = factor(group_time, levels = gt_levels))
  cv_named <- lapply(contrasts_named, function(v) unname(v[gt_levels]))
  rows <- list()
  for (mod in colnames(MEs)) {
    lmm_meta[[mod]] <- MEs[lmm_meta$sample_id, mod]
    fit <- tryCatch(
      suppressWarnings(lme4::lmer(
        as.formula(paste0("`", mod, "` ~ group_time + (1 | subject)")),
        data = lmm_meta, REML = TRUE
      )),
      error = function(e) NULL
    )
    if (is.null(fit)) next
    sing <- lme4::isSingular(fit)
    emm <- emmeans::emmeans(fit, ~group_time)
    for (cname in names(cv_named)) {
      ctr <- emmeans::contrast(emm, list(ctr = cv_named[[cname]]))
      s <- summary(ctr, ddf = "Kenward-Roger")
      r_eq <- sign(s$estimate) * sqrt(s$t.ratio^2 / (s$t.ratio^2 + s$df))
      rows <- c(rows, list(tibble(
        module = mod, model = "CRvCtl", contrast = cname,
        estimate = round(s$estimate, 5), SE = round(s$SE, 5),
        df = round(s$df, 2), t_ratio = round(s$t.ratio, 4),
        p_raw = s$p.value, r_equiv = round(r_eq, 4), singular = sing
      )))
    }
  }
  bind_rows(rows)
}

lmm_df <- fit_lmm_contrasts(gt_levels, crvh_contrasts)
lmm_df$p_bh <- p.adjust(lmm_df$p_raw, method = "BH")
write_csv(lmm_df, file.path(DAT, "wgcna_lmm_contrast_audit.csv"))

# fry + camera + omnibus F
samp <- colnames(expr)
grp <- factor(meta$group_time[match(samp, meta$sample_id)], levels = gt_levels)
block <- meta$subject[match(samp, meta$sample_id)]
stopifnot(!anyNA(grp))

design <- model.matrix(~ 0 + grp)
colnames(design) <- levels(grp)
cm <- vapply(crvh_contrasts, function(v) v[colnames(design)], numeric(ncol(design)))
rownames(cm) <- colnames(design)

non_grey <- setdiff(unique(mc), "grey")
idx <- split(seq_len(nrow(expr)), mc[rownames(expr)])[non_grey]
idx <- idx[lengths(idx) >= 3]
set_rho <- max(duplicateCorrelation(expr, design, block = block)$consensus, 0)
settests <- bind_rows(lapply(colnames(cm), function(cn) {
  fr <- fry(expr, idx, design, contrast = cm[, cn], block = block, correlation = set_rho)
  ca <- camera(expr, idx, design, contrast = cm[, cn])
  tibble(
    module_color = rownames(fr), contrast = cn,
    direction = fr$Direction, fry_p = fr$PValue, fry_fdr = fr$FDR,
    camera_p = ca[rownames(fr), "PValue"], camera_fdr = ca[rownames(fr), "FDR"]
  )
}))
write_csv(settests, file.path(DAT, "module_set_tests.csv"))

me <- t(as.matrix(MEs[samp, , drop = FALSE]))
me <- me[sub("^ME", "", rownames(me)) != "grey", , drop = FALSE]
eig_rho <- max(duplicateCorrelation(me, design, block = block)$consensus, 0)
omni_cm <- makeContrasts(
  CRE_T1 - H_T1, CRE_T2 - H_T1, PLA_T1 - H_T1, PLA_T2 - H_T1,
  levels = design
)
omni_fit <- eBayes(contrasts.fit(
  lmFit(me, design, block = block, correlation = eig_rho), omni_cm
))
omnibus <- tibble(
  module_color = sub("^ME", "", rownames(me)),
  F = omni_fit$F, p = omni_fit$F.p.value, fdr = p.adjust(omni_fit$F.p.value, "BH")
)
write_csv(omnibus, file.path(DAT, "module_omnibus_F.csv"))

# per-module fGSEA NES: modules as gene sets, ranked by DE moderated-t
rank_wide <- read_csv(
  here::here("03_DEP/a_non_imputed/c_data/combined_results_pi.csv"),
  show_col_types = FALSE
) |>
  filter(!is.na(gene), gene != "", contrast %in% names(crvh_contrasts)) |>
  select(gene, contrast, t) |>
  pivot_wider(names_from = contrast, values_from = t, names_prefix = "t_", values_fn = mean)
module_genes <- split(module_df$gene, module_df$module_color)
module_genes <- lapply(
  module_genes[setdiff(names(module_genes), "grey")], \(g) unique(na.omit(g))
)
nes <- run_module_fgsea(rank_wide, module_genes, names(crvh_contrasts))
write_csv(nes, file.path(DAT, "module_fgsea_nes.csv"))

# supplement arm: the CR-only 2x2, with creatine and placebo kept apart
# The main card pools the two arms because only 6 subjects per arm are paired. This
# block keeps them separate so the supplement can show what pooling hides.
cr_levels <- intersect(c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2"), gt_levels)
cr_idx <- which(grp %in% cr_levels)
grp_cr <- droplevels(factor(grp[cr_idx], levels = cr_levels))
design_cr <- model.matrix(~ 0 + grp_cr)
colnames(design_cr) <- levels(grp_cr)
block_cr <- block[cr_idx]
expr_cr <- expr[, cr_idx, drop = FALSE]

supp_contrasts <- makeContrasts(
  Baseline_Supplement = CRE_T1 - PLA_T1,
  Training_CRE = CRE_T2 - CRE_T1,
  Training_PLA = PLA_T2 - PLA_T1,
  Supplement_Interaction = (CRE_T2 - CRE_T1) - (PLA_T2 - PLA_T1),
  levels = design_cr
)
rho_cr <- max(duplicateCorrelation(expr_cr, design_cr, block = block_cr)$consensus, 0)
settests_supp <- bind_rows(lapply(colnames(supp_contrasts), function(cn) {
  fr <- fry(expr_cr, idx, design_cr,
    contrast = supp_contrasts[, cn], block = block_cr, correlation = rho_cr
  )
  tibble(
    module_color = rownames(fr), contrast = cn,
    direction = fr$Direction, fry_p = fr$PValue, fry_fdr = fr$FDR
  )
}))
write_csv(settests_supp, file.path(DAT, "module_set_tests_supp.csv"))

nes_supp <- run_module_fgsea(
  read_csv(here::here("03_DEP/a_non_imputed/c_data/combined_results_pi.csv"),
    show_col_types = FALSE
  ) |>
    filter(!is.na(gene), gene != "", contrast %in% colnames(supp_contrasts)) |>
    select(gene, contrast, t) |>
    pivot_wider(
      names_from = contrast, values_from = t,
      names_prefix = "t_", values_fn = mean
    ),
  module_genes, colnames(supp_contrasts)
)
write_csv(nes_supp, file.path(DAT, "module_fgsea_nes_supp.csv"))

# Five-cell eigengene means: healthy is shared by both arms, so each arm's line starts
# from the same control point.
traj_supp <- as.data.frame(MEs) |>
  rownames_to_column("sample_id") |>
  pivot_longer(-sample_id, names_to = "module", values_to = "eigengene") |>
  left_join(select(meta, sample_id, group_time, supplement, cancer), by = "sample_id") |>
  filter(module != "MEgrey") |>
  mutate(
    arm = if_else(cancer == "CTL", "Ctl", supplement),
    stage = case_when(
      cancer == "CTL" ~ "Ctl",
      grepl("_T1$", group_time) ~ "pre",
      TRUE ~ "post"
    )
  ) |>
  group_by(module, arm, stage) |>
  summarise(
    mean_eig = mean(eigengene), se = sd(eigengene) / sqrt(dplyr::n()),
    n = dplyr::n(), .groups = "drop"
  ) |>
  mutate(module_color = sub("^ME", "", module))
write_csv(traj_supp, file.path(DAT, "trajectory_eigengenes_supp.csv"))

# module-eigengene x clinical outcome association
eig <- as.matrix(MEs)
colnames(eig) <- sub("^ME", "", colnames(eig))
eig <- eig[, colnames(eig) != "grey", drop = FALSE]

cvh <- read_csv(here::here("00_input/CvH_meta.csv"), show_col_types = FALSE)
cvh <- cvh[match(meta$sample_id, cvh$Col_ID), ]
outcome_bases <- c(
  "ALM_kg", "sts_max_pwr", "LBM_kg", "chest_press_lbs", "leg_ext_lbs", "grip_lbs"
)
tp <- meta$timepoint

# Primary: baseline (T1) samples, one per subject, vs pre_ outcomes.
base_ids <- meta$sample_id[tp == "T1"]
pre_cols <- c("age", paste0("pre_", outcome_bases))
trait_b <- as.matrix(cvh[tp == "T1", pre_cols])
rownames(trait_b) <- base_ids
cor_primary <- module_trait_cor(eig[base_ids, , drop = FALSE], trait_b)
write_csv(cor_primary, file.path(DAT, "module_trait_cor.csv"))

# Comparison arms: all samples, T1 -> pre_ and T2 -> post_ matched outcomes.
matched <- vapply(outcome_bases, function(b) {
  ifelse(tp == "T1", cvh[[paste0("pre_", b)]], cvh[[paste0("post_", b)]])
}, numeric(nrow(cvh)))
matched <- cbind(age = cvh$age, matched)
rownames(matched) <- meta$sample_id

cor_matched <- module_trait_cor(eig, matched)
write_csv(cor_matched, file.path(DAT, "module_trait_cor_matched.csv"))

eig_long <- as.data.frame(eig) |>
  rownames_to_column("sample_id") |>
  pivot_longer(-sample_id, names_to = "module", values_to = "eigengene")
trait_long <- as.data.frame(matched) |>
  rownames_to_column("sample_id") |>
  pivot_longer(-sample_id, names_to = "trait", values_to = "value")
block_vec <- setNames(meta$subject, meta$sample_id)
lmm_pheno <- module_trait_lmm(eig_long, trait_long, block_vec)
write_csv(lmm_pheno, file.path(DAT, "module_trait_lmm.csv"))

message(sprintf(
  "01 module stats: %d modules | fry rho %.3f | NES rows %d | pheno baseline n=%d",
  length(idx), set_rho, nrow(nes), length(base_ids)
))
