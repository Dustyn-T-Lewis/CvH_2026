# F06 stage 00: build the three feature spaces and the outcome vectors they are scored
# against. The spaces differ in how far their definition sits from this dataset:
# raw proteins are fully data-driven, module eigengenes are data-derived structure, and
# Hallmark singscores come from an external collection that never saw these samples.
# That gradient is what the figure compares.

setwd(here::here())
pacman::p_load(singscore, msigdbr, dplyr, tidyr, tibble, readr, purrr)

DAT <- "04_Figures/F06_Prediction/c_data"
F05 <- "04_Figures/F05_WGCNA/c_data"
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)
set.seed(42)

imp_mat <- readRDS(file.path(F05, "imp_mat.rds"))
ann <- read_csv(file.path(F05, "imp_annotations.csv"), show_col_types = FALSE)
mes <- readRDS(file.path(F05, "MEs.rds"))
meta <- read_csv(file.path(F05, "meta.csv"), show_col_types = FALSE)
stopifnot(identical(rownames(imp_mat), ann$uniprot_id))

# Protein space: z-score per protein so no single high-variance protein dominates the
# in-fold correlation ranking, then orient samples as rows.
x_protein <- t(scale(t(imp_mat)))
x_protein <- x_protein[stats::complete.cases(x_protein), , drop = FALSE]
x_protein <- t(x_protein)

# Module space: grey collects proteins with no assignment, so it is not a module.
x_module <- as.matrix(mes[, colnames(mes) != "MEgrey", drop = FALSE])
colnames(x_module) <- sub("^ME", "", colnames(x_module))
x_module <- x_module[rownames(x_protein), , drop = FALSE]

# Hallmark space. singscore ranks within each sample, so it needs one row per gene
# symbol; proteins mapping to the same symbol are averaged first.
gene_mat <- imp_mat[!is.na(ann$gene) & ann$gene != "", , drop = FALSE]
gene_key <- ann$gene[!is.na(ann$gene) & ann$gene != ""]
gene_sums <- rowsum(gene_mat, gene_key)
gene_mat <- gene_sums / as.vector(table(gene_key)[rownames(gene_sums)])

hallmark <- msigdbr(species = "Homo sapiens", collection = "H") |>
  distinct(gs_name, gene_symbol)
sets <- split(hallmark$gene_symbol, hallmark$gs_name)
sets <- keep(sets, ~ length(intersect(.x, rownames(gene_mat))) >= 10)

ranked <- rankGenes(gene_mat)
x_hallmark <- vapply(
  sets,
  function(s) simpleScore(ranked, upSet = intersect(s, rownames(gene_mat)))$TotalScore,
  numeric(ncol(gene_mat))
)
rownames(x_hallmark) <- colnames(gene_mat)
x_hallmark <- x_hallmark[rownames(x_protein), , drop = FALSE]

# Outcomes. Baseline is cross-sectional and unpaired; training is within-subject and
# only defined for the survivors sampled at both timepoints.
meta <- meta |> filter(sample_id %in% rownames(x_protein))
baseline_ids <- meta$sample_id[meta$timepoint == "T1"]
paired_subjects <- meta |>
  filter(cancer == "SURV") |>
  count(subject) |>
  filter(n == 2) |>
  pull(subject)
training_ids <- meta$sample_id[meta$subject %in% paired_subjects & meta$cancer == "SURV"]

outcomes <- list(
  baseline = list(
    ids = baseline_ids,
    y = as.integer(meta$cancer[match(baseline_ids, meta$sample_id)] == "SURV"),
    subject = meta$subject[match(baseline_ids, meta$sample_id)],
    paired = FALSE,
    label = "CR vs Ctl (baseline)"
  ),
  training = list(
    ids = training_ids,
    y = as.integer(meta$timepoint[match(training_ids, meta$sample_id)] == "T2"),
    subject = meta$subject[match(training_ids, meta$sample_id)],
    paired = TRUE,
    label = "CR pre vs post"
  )
)

features <- list(
  protein = x_protein,
  module = x_module,
  hallmark = x_hallmark
)

saveRDS(
  list(features = features, outcomes = outcomes, meta = meta),
  file.path(DAT, "feature_bundle.rds")
)

write_csv(
  tibble(
    space = names(features),
    n_features = vapply(features, ncol, integer(1)),
    n_samples = vapply(features, nrow, integer(1))
  ),
  file.path(DAT, "feature_summary.csv")
)

message(sprintf(
  "F06 features: %d proteins | %d modules | %d Hallmark sets (of 50) | baseline n=%d (%d CR) | training n=%d (%d subjects)",
  ncol(x_protein), ncol(x_module), ncol(x_hallmark),
  length(baseline_ids), sum(outcomes$baseline$y),
  length(training_ids), length(paired_subjects)
))
