# Module-level statistics for F05_WGCNA: fGSEA NES with modules as gene sets,
# and module-eigengene associations with clinical outcomes.

pacman::p_load(fgsea, dplyr, tibble, tidyr, purrr, WGCNA, lme4)

run_module_fgsea <- function(rank_wide, module_genes, contrasts,
                             min_size = 10, max_size = 5000, eps = 1e-10) {
  one <- function(contrast) {
    stats <- rank_wide[[paste0("t_", contrast)]]
    names(stats) <- rank_wide$gene
    stats <- stats[!is.na(stats) & !is.na(names(stats)) & names(stats) != ""]
    if (anyDuplicated(names(stats))) stats <- tapply(stats, names(stats), mean)
    stats <- sort(stats)
    set.seed(42)
    res <- fgsea::fgseaMultilevel(
      pathways = module_genes, stats = stats,
      minSize = min_size, maxSize = max_size, eps = eps
    )
    tibble(
      module = res$pathway, contrast = contrast,
      NES = res$NES, pval = res$pval, size = res$size
    )
  }
  bind_rows(lapply(contrasts, one)) |>
    group_by(contrast) |>
    mutate(padj = p.adjust(pval, "BH")) |>
    ungroup()
}

module_trait_cor <- function(eig, traits, n_samples = NULL) {
  r <- stats::cor(eig, traits, use = "pairwise.complete.obs")
  n_mat <- matrix(rep(colSums(!is.na(traits)), each = nrow(r)),
    nrow = nrow(r), dimnames = dimnames(r)
  )
  if (!is.null(n_samples)) n_mat[] <- n_samples
  p <- WGCNA::corPvalueStudent(r, n_mat)
  as_tibble(as.data.frame(as.table(r)), .name_repair = "minimal") |>
    setNames(c("module", "trait", "r")) |>
    mutate(p = as.vector(as.table(p)), n = as.vector(as.table(n_mat)), padj = p.adjust(p, "BH"))
}

module_trait_lmm <- function(eig_long, trait_long, block) {
  grid <- expand_grid(
    module = unique(eig_long$module), trait = unique(trait_long$trait)
  )
  fit_one <- function(module, trait) {
    d <- eig_long |>
      filter(module == !!module) |>
      inner_join(filter(trait_long, trait == !!trait), by = "sample_id") |>
      mutate(subject = block[sample_id], z = as.numeric(scale(value)))
    ok <- sum(!is.na(d$z) & !is.na(d$eigengene))
    if (ok < 4) {
      return(tibble(beta = NA_real_, se = NA_real_, t = NA_real_, p = NA_real_))
    }
    fit <- suppressWarnings(
      lme4::lmer(eigengene ~ z + (1 | subject), data = d, REML = TRUE)
    )
    cf <- summary(fit)$coefficients["z", ]
    df <- nrow(d) - 2
    tibble(
      beta = unname(cf["Estimate"]), se = unname(cf["Std. Error"]),
      t = unname(cf["t value"]), p = 2 * stats::pt(-abs(cf["t value"]), df)
    )
  }
  pmap_dfr(grid, function(module, trait) {
    bind_cols(tibble(module = module, trait = trait), fit_one(module, trait))
  }) |>
    mutate(padj = p.adjust(p, "BH"))
}
