# Per-module member response on the pooled CR-vs-Ctl design. fry (self-contained
# rotation) asks whether a module's member proteins move together under each
# contrast; the eigengene omnibus F flags modules carrying any group effect (the
# main-figure gate). Same design and subject block as the module-trait lmm.

setwd(here::here())
pacman::p_load(limma, dplyr, tibble, readr)

DAT <- "04_Figures/F04_WGCNA/c_data"
expr <- readRDS(file.path(DAT, "imp_mat.rds"))
mc <- readRDS(file.path(DAT, "module_colors.rds"))
names(mc) <- rownames(expr)
MEs <- readRDS(file.path(DAT, "MEs.rds"))
meta <- read_csv(file.path(DAT, "meta.csv"), show_col_types = FALSE)

gt_levels <- intersect(
  c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1"), unique(meta$group_time)
)
samp <- colnames(expr)
grp <- factor(meta$group_time[match(samp, meta$sample_id)], levels = gt_levels)
block <- meta$subject[match(samp, meta$sample_id)]
stopifnot(!anyNA(grp))

design <- model.matrix(~ 0 + grp)
colnames(design) <- levels(grp)
contrast_vecs <- list(
  CRvH_Baseline = c(CRE_T1 = 0.25, CRE_T2 = 0.25, PLA_T1 = 0.25, PLA_T2 = 0.25, H_T1 = -1),
  CR_Training = c(CRE_T1 = -0.5, CRE_T2 = 0.5, PLA_T1 = -0.5, PLA_T2 = 0.5, H_T1 = 0)
)
cm <- vapply(contrast_vecs, function(v) v[colnames(design)], numeric(ncol(design)))
rownames(cm) <- colnames(design)

non_grey <- setdiff(unique(mc), "grey")
idx <- split(seq_len(nrow(expr)), mc[rownames(expr)])[non_grey]
idx <- idx[lengths(idx) >= 3]
set_rho <- max(duplicateCorrelation(expr, design, block = block)$consensus, 0)
settests <- bind_rows(lapply(colnames(cm), function(cn) {
  fr <- fry(expr, idx, design, contrast = cm[, cn], block = block, correlation = set_rho)
  tibble(
    module_color = rownames(fr), contrast = cn,
    direction = fr$Direction, fry_p = fr$PValue, fry_fdr = fr$FDR
  )
}))
write_csv(settests, file.path(DAT, "module_set_tests.csv"))

# Eigengene omnibus F: any group deviation from healthy, per module.
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

message(sprintf(
  "member response: %d modules x %d contrasts | set rho %.3f",
  length(idx), ncol(cm), set_rho
))
print(arrange(omnibus, p))
