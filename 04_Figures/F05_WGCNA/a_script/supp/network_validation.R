# F05 supplement: is the co-expression network real, and are its parameters defensible?
#
# Five checks, in the order a sceptical reader would ask them:
#   1. Are any samples outliers that could drive the correlation structure?
#   2. Is the soft power justified, and how much does the choice matter?
#   3. Do the modules survive different clustering parameters, or are they an artifact
#      of minModuleSize and mergeCutHeight?
#   4. Does each module survive resampling of the subjects?
#   5. Does a matched noise dataset produce modules too?
# Nothing here gates a result. Each check reports a number the caption can carry.

setwd(here::here())
pacman::p_load(WGCNA, dplyr, tidyr, tibble, readr, purrr, ggplot2, patchwork)
source("04_Figures/shared/style.R")
source("04_Figures/shared/prediction_utils.R")

allowWGCNAThreads()
set.seed(42)

DAT <- "04_Figures/F05_WGCNA/c_data"
SUPP_PNG <- "04_Figures/F05_WGCNA/b_reports/supp/png"
SUPP_PDF <- "04_Figures/F05_WGCNA/b_reports/supp/pdf"
pdf_device <- get_pdf_device()

REF_POWER <- readRDS(file.path(DAT, "wgcna_network.rds"))$chosen_power
REF_MIN_SIZE <- 30L
REF_MERGE <- 0.25
N_SUBSAMPLE <- 100L
N_NULL <- 25L
SUBSAMPLE_FRAC <- 0.8

datExpr <- readRDS(file.path(DAT, "datExpr.rds"))
ref_colors <- readRDS(file.path(DAT, "module_colors.rds"))
names(ref_colors) <- colnames(datExpr)
meta <- read_csv(file.path(DAT, "meta.csv"), show_col_types = FALSE)
meta <- meta[match(rownames(datExpr), meta$sample_id), ]
ref_modules <- setdiff(unique(ref_colors), "grey")

build_net <- function(expr, power = REF_POWER, min_size = REF_MIN_SIZE,
                      merge_height = REF_MERGE, cor_type = "bicor") {
  net <- blockwiseModules(
    expr,
    power = power, networkType = "signed", TOMType = "signed",
    corType = cor_type, maxPOutliers = 0.05,
    minModuleSize = min_size, mergeCutHeight = merge_height,
    numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = FALSE, verbose = 0
  )
  out <- labels2colors(net$colors)
  names(out) <- colnames(expr)
  out
}

cor <- WGCNA::cor

# 1. Sample outliers. A sample whose standardised connectivity sits far below the others
# can manufacture correlation structure on its own.
adj <- adjacency(t(datExpr), type = "distance")
k <- colSums(adj) - 1
sample_qc <- tibble(
  sample_id = rownames(datExpr),
  group = meta$cancer_time,
  connectivity = k,
  z_k = as.numeric(scale(k))
) |>
  mutate(outlier = z_k < -2.5)
write_csv(sample_qc, file.path(DAT, "validation_sample_qc.csv"))

# 2. Soft power. Report the scale-free fit for both correlation functions so the choice
# of bicor is visible, not assumed.
powers <- 1:20
sft_curve <- map_dfr(c("bicor", "pearson"), function(ct) {
  sft <- pickSoftThreshold(datExpr,
    powerVector = powers, networkType = "signed",
    corFnc = if (ct == "bicor") bicor else stats::cor,
    corOptions = if (ct == "bicor") list(maxPOutliers = 0.05) else list(),
    verbose = 0
  )
  tibble(
    cor_type = ct, power = sft$fitIndices$Power,
    signed_r2 = -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq,
    mean_k = sft$fitIndices$mean.k., median_k = sft$fitIndices$median.k.
  )
})
write_csv(sft_curve, file.path(DAT, "validation_softpower.csv"))

# 3. Parameter grid. For each setting, how many modules appear, how much lands in grey,
# and how well the reference modules are recovered.
grid <- expand_grid(
  power = c(10L, 12L, 14L, 16L),
  min_size = c(15L, 20L, 30L, 40L),
  merge_height = c(0.15, 0.25, 0.35)
)
grid_res <- pmap_dfr(grid, function(power, min_size, merge_height) {
  mc <- build_net(datExpr, power, min_size, merge_height)
  matched <- match_modules(mc, ref_colors)
  tibble(
    power, min_size, merge_height,
    n_modules = length(setdiff(unique(mc), "grey")),
    pct_grey = 100 * mean(mc == "grey"),
    mean_recovery = mean(matched$jaccard),
    n_recovered = sum(matched$jaccard >= 0.5)
  )
})
write_csv(grid_res, file.path(DAT, "validation_grid.csv"))

grid_by_module <- pmap_dfr(grid, function(power, min_size, merge_height) {
  mc <- build_net(datExpr, power, min_size, merge_height)
  match_modules(mc, ref_colors) |>
    mutate(power = power, min_size = min_size, merge_height = merge_height)
})
write_csv(grid_by_module, file.path(DAT, "validation_grid_by_module.csv"))

# 4. Subject-level subsampling. Both timepoints of a subject leave together, so a module
# cannot be propped up by one subject appearing twice.
subjects <- unique(meta$subject)
stability <- map_dfr(seq_len(N_SUBSAMPLE), function(i) {
  keep <- sample(subjects, floor(SUBSAMPLE_FRAC * length(subjects)))
  rows <- which(meta$subject %in% keep)
  if (length(rows) < 15) {
    return(tibble())
  }
  mc <- build_net(datExpr[rows, , drop = FALSE])
  match_modules(mc, ref_colors) |> mutate(rep = i)
})
stability_summary <- stability |>
  group_by(full) |>
  summarise(
    mean_jaccard = mean(jaccard), median_jaccard = median(jaccard),
    min_jaccard = min(jaccard), pct_above_50 = 100 * mean(jaccard >= 0.5),
    n_reps = n(), .groups = "drop"
  ) |>
  mutate(verdict = case_when(
    mean_jaccard >= 0.85 ~ "stable",
    mean_jaccard >= 0.60 ~ "pattern",
    TRUE ~ "dissolved"
  ))
write_csv(stability_summary, file.path(DAT, "validation_stability.csv"))

# 5. Noise null. Permuting each protein independently across samples destroys every
# inter-protein correlation while preserving each protein's own distribution. Any module
# found here is an artifact of the clustering, not of biology.
null_res <- map_dfr(seq_len(N_NULL), function(i) {
  permuted <- apply(datExpr, 2, sample)
  dimnames(permuted) <- dimnames(datExpr)
  mc <- build_net(permuted)
  tibble(
    rep = i, n_modules = length(setdiff(unique(mc), "grey")),
    pct_grey = 100 * mean(mc == "grey")
  )
})
write_csv(null_res, file.path(DAT, "validation_null.csv"))

# 6. Borderline samples. Stage 01 drops a sample on >= 3/4 QC heuristics; these two
# scored 2/4 and were retained, then flagged again here on network connectivity. Rather
# than lower a pre-specified threshold after seeing this, rebuild without them and report
# whether the modules depend on them.
borderline <- sample_qc$sample_id[sample_qc$outlier]
keep_rows <- which(!rownames(datExpr) %in% borderline)
sens_colors <- build_net(datExpr[keep_rows, , drop = FALSE])
outlier_sens <- match_modules(sens_colors, ref_colors) |>
  mutate(dropped = paste(borderline, collapse = ", "), n_kept = length(keep_rows))
write_csv(outlier_sens, file.path(DAT, "validation_outlier_sensitivity.csv"))

cor <- stats::cor

mod_sizes <- table(ref_colors)[ref_modules]
stability_summary <- stability_summary |>
  mutate(n_proteins = as.integer(mod_sizes[full]))

p_power <- ggplot(sft_curve, aes(power, signed_r2, colour = cor_type)) +
  geom_hline(yintercept = 0.90, linetype = "dashed", colour = "grey55") +
  geom_vline(xintercept = REF_POWER, colour = "#D6604D", linewidth = 0.4) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = c(bicor = "#1B7837", pearson = "#8073AC"), name = NULL) +
  labs(
    title = "Soft-power selection", x = "Soft threshold (power)",
    y = expression(signed ~ R^2),
    subtitle = sprintf("dashed = 0.90 criterion; red = chosen power %d", REF_POWER)
  ) +
  FIG_THEME +
  theme(legend.position = "bottom")

p_grid <- ggplot(grid_res, aes(factor(min_size), factor(power), fill = n_modules)) +
  geom_tile(colour = "grey85", linewidth = 0.3) +
  geom_text(aes(label = n_modules), size = 2.4, colour = "grey15") +
  facet_wrap(~merge_height, labeller = label_both) +
  scale_fill_gradient(low = "#F7F7F7", high = "#2166AC", name = "modules") +
  labs(
    title = "Modules found across parameters", x = "minModuleSize", y = "soft power",
    subtitle = "the 5-module result is a parameter choice, not a fixed property"
  ) +
  FIG_THEME +
  theme(strip.text = element_text(size = 6.5))

p_stab <- stability_summary |>
  mutate(full = factor(full, levels = full[order(mean_jaccard)])) |>
  ggplot(aes(mean_jaccard, full, fill = verdict)) +
  geom_col(width = 0.7, colour = "grey30", linewidth = 0.3) +
  geom_vline(xintercept = c(0.5, 0.6, 0.85), linetype = "dashed", colour = "grey55") +
  geom_text(aes(label = sprintf("n=%d  %.0f%% of reps >0.5", n_proteins, pct_above_50)),
    hjust = -0.05, size = 1.9, colour = "grey25"
  ) +
  scale_fill_manual(
    values = c(stable = "#1B7837", pattern = "#F4A582", dissolved = "#D6604D"),
    name = NULL
  ) +
  scale_x_continuous(limits = c(0, 1.6), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(
    title = "Module survival under subject resampling",
    subtitle = sprintf(
      "%d subsamples at %.0f%% of subjects; Hennig thresholds 0.5 / 0.6 / 0.85",
      N_SUBSAMPLE, 100 * SUBSAMPLE_FRAC
    ),
    x = "mean Jaccard vs full-sample module", y = NULL
  ) +
  FIG_THEME

p_null <- ggplot(null_res, aes(n_modules)) +
  geom_histogram(binwidth = 1, fill = "grey70", colour = "grey30", linewidth = 0.3) +
  geom_vline(xintercept = length(ref_modules), colour = "#1B7837", linewidth = 1) +
  labs(
    title = "Modules recovered from noise",
    subtitle = sprintf(
      "%d permuted datasets; green line = observed (%d)",
      N_NULL, length(ref_modules)
    ),
    x = "modules found", y = "permutations"
  ) +
  FIG_THEME

validation_fig <- (p_power | p_grid) / (p_stab | p_null) +
  plot_annotation(
    title = "WGCNA network validation",
    subtitle = paste(
      "Parameter sensitivity, module stability under subject resampling, and a matched noise null.",
      sprintf("%d samples, %d proteins.", nrow(datExpr), ncol(datExpr))
    ),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, colour = "grey10"),
      plot.subtitle = element_text(face = "italic", size = 8, colour = "grey40")
    )
  )

ggsave(file.path(SUPP_PNG, "SUPP_F05_network_validation.png"), validation_fig,
  width = 280, height = 190, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(SUPP_PDF, "SUPP_F05_network_validation.pdf"), validation_fig,
  width = 280, height = 190, units = "mm", device = pdf_device
)

message(sprintf(
  "validation @ power %d: %d borderline samples | grid %d-%d modules | null median %d modules",
  REF_POWER, sum(sample_qc$outlier), min(grid_res$n_modules), max(grid_res$n_modules),
  stats::median(null_res$n_modules)
))
print(as.data.frame(stability_summary))
message("outlier sensitivity (rebuild without borderline samples):")
print(as.data.frame(outlier_sens[, c("full", "train", "jaccard")]))
