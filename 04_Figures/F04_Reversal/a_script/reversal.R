# reversal.R  --  rejuvenation / signature-reversal engine for CvH.
#
# Pure analytical functions (no plotting, no file IO) consumed by the 04_Figures
# Reversal panels (via reversal_inputs.R + 00_build_fgsea_cache.R). Quantifies
# whether the training response (T = CR_post - CR_pre) reverses the disease
# deviation (D = CR_pre - H_pre), with a residual axis (R = CR_post - H_pre = D + T)
# tracking what remains.
#
# Design note (shared-baseline circularity): D and T share the CR_pre samples,
# so a structural negative cor(D, T) is mathematically guaranteed (Smyth & Altman
# 2013, PMID 23705896). Every directional claim is therefore tested against a
# protein-label permutation null that asks whether the disease-DEP set reverses
# MORE than random proteins under the SAME shared-baseline structure -- not
# whether reversal merely exceeds zero. fry (rotation, within the limma model)
# and RRHO2 (rank-based) corroborate from frameworks that are model-aware /
# non-parametric.
#
# Contrast names expected in the long DEP table (combined_results_pi.csv):
#   D = "CRvH_Baseline"   disease axis   (CR_pre - H_pre)
#   T = "CR_Training"     training axis  (CR_post - CR_pre)
#   R = "Resid"          residual axis  (CR_post - H_pre)
#
# Method lineage: Melov 2007 (PMID 17520024, proportion + permutation), Robinson
# 2017 (PMID 28273480, logFC correlation), Wu & Smyth 2010/2012 (PMID 20610611 /
# 22638577, ROAST/CAMERA), Cahill 2018 (RRHO2), Smyth & Altman 2013 (PMID
# 23705896, shared-baseline). Full citation catalog: see CvH_pipeline.qmd §5.

pacman::p_load(dplyr, tidyr, tibble, purrr, readr)

# Default rejuvenation-fraction band: |phi| < 0.25 = unchanged, >= 0.25 = moved.
REVERSAL_PHI_BAND <- 0.25
REVERSAL_CONTRASTS <- c(D = "CRvH_Baseline", T = "CR_Training", R = "Resid")

# 1. Reshape long DEP -> one row per protein with D/T/R columns
#' Pivot the long combined_results_pi.csv into a per-protein wide frame.
#' Keeps logFC, t, P.Value, pi_score and the Pi-significance flag for each axis.
load_reversal_table <- function(combined_pi_path,
                                contrasts = REVERSAL_CONTRASTS) {
  long <- readr::read_csv(combined_pi_path, show_col_types = FALSE)
  stopifnot(all(c("uniprot_id", "logFC", "t", "P.Value", "pi_score",
                  "sig_pi", "contrast", "gene") %in% names(long)))
  long <- long |> filter(contrast %in% contrasts)

  # dplyr:: qualified -- AnnotationDbi (pulled in by RRHO2/GO.db) masks select()
  ann <- long |>
    dplyr::select(any_of(c("uniprot_id", "gene", "protein", "description"))) |>
    distinct(uniprot_id, .keep_all = TRUE)

  vals <- long |>
    mutate(axis = names(contrasts)[match(contrast, contrasts)]) |>
    dplyr::select(uniprot_id, axis, logFC, t, P.Value, pi_score, sig_pi) |>
    pivot_wider(names_from = axis,
                values_from = c(logFC, t, P.Value, pi_score, sig_pi),
                names_sep = "_")

  ann |> inner_join(vals, by = "uniprot_id")
}

# 2. Rejuvenation fraction phi = -T / D + 3-way classification
#' phi > 0  : training opposes the disease change (reversal toward healthy)
#' phi = 1  : exact reversal (Resid = 0); phi > 1 = overshoot past healthy
#' phi < 0  : training amplifies the disease change (exacerbation)
#' Classification on the disease-signature set (Pi < 0.05 on the D axis):
#'   Normalized  phi >=  band   (training reverses >= band of disease deviation)
#'   Persistent  |phi| < band   (disease deviation largely unchanged by training)
#'   Exacerbated phi <= -band   (training pushes further from healthy)
compute_phi <- function(wide, band = REVERSAL_PHI_BAND, signature_axis = "D") {
  sig_col <- paste0("sig_pi_", signature_axis)
  stopifnot(all(c("logFC_D", "logFC_T", "logFC_R", sig_col) %in% names(wide)))

  wide |>
    mutate(
      # sig_pi encodes direction (+1 up / -1 down / 0 ns) -- signature = non-zero
      disease_sig = .data[[sig_col]] != 0,
      disease_dir = if_else(logFC_D > 0, "Disease Up", "Disease Down"),
      phi         = -logFC_T / logFC_D,
      reversed    = sign(logFC_T) != sign(logFC_D),   # phi > 0
      # residual fraction remaining relative to disease deviation: R = D * (1 - phi)
      resid_frac  = logFC_R / logFC_D,
      reversal_class = factor(case_when(
        phi >=  band ~ "Normalized",
        phi <= -band ~ "Exacerbated",
        TRUE         ~ "Persistent"
      ), levels = c("Normalized", "Persistent", "Exacerbated"))
    )
}

# 3. Directional asymmetry (disease-down vs disease-up reversal rate)
#' Two-proportion test on the disease-signature set: do disease-DOWN proteins
#' reverse at a different rate than disease-UP proteins? (Novel proteome-level
#' framing; mechanistic basis Murgia 2023 PMID 36517414.)
directional_asymmetry <- function(phi_tbl) {
  sig <- phi_tbl |> filter(disease_sig, is.finite(phi))
  by_dir <- sig |>
    group_by(disease_dir) |>
    summarise(n = n(), n_reversed = sum(reversed),
              pct_reversed = 100 * mean(reversed), .groups = "drop")

  up <- by_dir |> filter(disease_dir == "Disease Up")
  dn <- by_dir |> filter(disease_dir == "Disease Down")
  test <- if (nrow(up) == 1 && nrow(dn) == 1) {
    pt <- suppressWarnings(prop.test(c(dn$n_reversed, up$n_reversed),
                                     c(dn$n, up$n)))
    tibble(comparison = "Disease Down vs Disease Up",
           pct_down = dn$pct_reversed, pct_up = up$pct_reversed,
           chisq = unname(pt$statistic), p_value = pt$p.value)
  } else tibble(comparison = "Disease Down vs Disease Up",
                pct_down = NA_real_, pct_up = NA_real_,
                chisq = NA_real_, p_value = NA_real_)
  list(by_direction = by_dir, test = test)
}

# 4. Protein-label permutation null (shared-baseline circularity control)
#' Observed reversal fraction of the disease-DEP set vs a null where DEP
#' membership is reassigned to random proteins of the same size, holding ALL
#' fold changes (and the shared-baseline covariance) fixed. p = P(null >= obs).
#' Tests excess reversal in disease-DEPs over background, not reversal vs zero.
reversal_permutation_null <- function(phi_tbl, n_perm = 2000, seed = 42) {
  usable <- phi_tbl |> filter(is.finite(phi))
  n_dep  <- sum(usable$disease_sig)
  obs    <- mean(usable$reversed[usable$disease_sig])
  rev_all <- usable$reversed
  set.seed(seed)
  null <- replicate(n_perm, mean(rev_all[sample.int(length(rev_all), n_dep)]))
  tibble(
    n_dep = n_dep, n_background = nrow(usable),
    obs_reversal_frac = obs,
    null_mean = mean(null), null_sd = sd(null),
    z = (obs - mean(null)) / sd(null),
    p_value = (sum(null >= obs) + 1) / (n_perm + 1),
    n_perm = n_perm
  )
}

# 5. fry + camera rotation test of the disease signature on the T axis
#' Tests, within the limma model (block = Subject_ID, duplicateCorrelation),
#' whether the disease-UP / disease-DOWN sets move in the REVERSING direction
#' under the training contrast. fry = self-contained rotation (ROAST family,
#' robust to inter-gene correlation & repeated measures); camera = competitive
#' counterpart. Reversal => disease-Up set Direction "Down", disease-Down "Up".
reversal_rotation_test <- function(dal, phi_tbl,
                                   design_formula = ~ 0 + group_time,
                                   training_contrast = c(CR_post = 1, CR_pre = -1),
                                   block_col = "Subject_ID") {
  requireNamespace("limma", quietly = TRUE)
  y    <- as.matrix(dal$data)
  if (anyNA(y)) {                       # rotation tests need a complete matrix
    keep <- rowSums(is.na(y)) == 0
    warning(sprintf("rotation test: dropping %d/%d proteins with NAs (pass the imputed DAList to avoid this)",
                    sum(!keep), length(keep)))
    y <- y[keep, , drop = FALSE]
  }
  meta <- dal$metadata
  gt   <- factor(meta$group_time)
  design <- model.matrix(design_formula, data = meta)
  colnames(design) <- sub("^group_time", "", colnames(design))

  # training contrast vector aligned to design columns
  con <- numeric(ncol(design)); names(con) <- colnames(design)
  con[names(training_contrast)] <- training_contrast

  block <- meta[[block_col]]
  dc <- limma::duplicateCorrelation(y, design, block = block)
  fit_args <- list(y = y, design = design, contrast = con,
                   block = block, correlation = dc$consensus.correlation)

  sig <- phi_tbl |> filter(disease_sig)
  idx <- list(
    Disease_Up   = which(rownames(y) %in% sig$uniprot_id[sig$logFC_D > 0]),
    Disease_Down = which(rownames(y) %in% sig$uniprot_id[sig$logFC_D < 0])
  )
  idx <- idx[vapply(idx, length, integer(1)) >= 3]   # fry needs >= ~3 members

  fry_res <- do.call(limma::fry, c(list(index = idx), fit_args)) |>
    as.data.frame() |> rownames_to_column("set")
  # camera (competitive) does not support duplicateCorrelation -> no block here;
  # fry above is the repeated-measures-aware primary, camera is a robustness check
  cam_res <- limma::camera(y = y, index = idx, design = design, contrast = con) |>
    as.data.frame() |> rownames_to_column("set")

  # reversal expectation: Up set should fall (Down), Down set should rise (Up)
  expect <- c(Disease_Up = "Down", Disease_Down = "Up")
  fry_res$reversing <- fry_res$Direction == expect[fry_res$set]
  list(fry = fry_res, camera = cam_res,
       consensus_correlation = dc$consensus.correlation)
}

# 6. RRHO2 (threshold-free 4-quadrant concordance/discordance)
#' Ranks proteins by signed significance on D and T; the off-diagonal quadrants
#' (disease-Up x training-Down, disease-Down x training-Up) capture reversal.
#' Returns the RRHO2 object plus tidy hotspot genes per quadrant. Wrapped in
#' tryCatch because RRHO2 is sensitive to ties / tiny inputs.
reversal_rrho2 <- function(phi_tbl, labels = c("Disease", "Training"),
                           stepsize = NULL) {
  if (!requireNamespace("RRHO2", quietly = TRUE)) return(NULL)
  d <- phi_tbl |>
    filter(is.finite(logFC_D), is.finite(logFC_T),
           is.finite(P.Value_D), is.finite(P.Value_T))
  mk <- function(lfc, p) -log10(p) * sign(lfc)
  l1 <- data.frame(gene = d$uniprot_id, score = mk(d$logFC_D, d$P.Value_D))
  l2 <- data.frame(gene = d$uniprot_id, score = mk(d$logFC_T, d$P.Value_T))
  if (is.null(stepsize)) stepsize <- floor(sqrt(nrow(d)))

  tryCatch({
    obj <- RRHO2::RRHO2_initialize(l1, l2, labels = labels,
                                   log10.ind = TRUE, boundary = 0.02)
    obj
  }, error = function(e) {
    message("RRHO2 failed: ", conditionMessage(e)); NULL
  })
}

# 7. One-call summary
#' Runs phi classification, asymmetry, permutation null and (optionally) the
#' rotation test; returns a named list of tidy tables for the caller to persist.
run_reversal_analysis <- function(combined_pi_path, dal = NULL,
                                   band = REVERSAL_PHI_BAND, n_perm = 2000,
                                   seed = 42) {
  wide  <- load_reversal_table(combined_pi_path)
  phi   <- compute_phi(wide, band = band)
  asym  <- directional_asymmetry(phi)
  null  <- reversal_permutation_null(phi, n_perm = n_perm, seed = seed)

  class_counts <- phi |> filter(disease_sig) |>
    count(disease_dir, reversal_class, name = "n") |>
    group_by(disease_dir) |> mutate(pct = 100 * n / sum(n)) |> ungroup()

  rho <- phi |> filter(is.finite(logFC_D), is.finite(logFC_T))
  global_cor <- tibble(
    set = c("All proteins", "Disease signature"),
    n = c(nrow(rho), sum(rho$disease_sig)),
    spearman = c(cor(rho$logFC_D, rho$logFC_T, method = "spearman"),
                 cor(rho$logFC_D[rho$disease_sig], rho$logFC_T[rho$disease_sig],
                     method = "spearman"))
  )

  out <- list(phi_table = phi, class_counts = class_counts,
              asymmetry_by_dir = asym$by_direction, asymmetry_test = asym$test,
              permutation_null = null, global_correlation = global_cor)
  if (!is.null(dal)) {
    rot <- reversal_rotation_test(dal, phi)
    out$fry <- rot$fry; out$camera <- rot$camera
    out$rotation_meta <- tibble(consensus_correlation = rot$consensus_correlation)
  }
  out
}
