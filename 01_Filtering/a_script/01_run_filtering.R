#!/usr/bin/env Rscript

pacman::p_load(proteoDA, here, readxl, readr, dplyr, tidyr, stringr, openxlsx, ggplot2, forcats, patchwork)
set.seed(42)

data_dir <- here("01_Filtering", "c_data") # outputs; cleared here, written at the end
report_dir <- here("01_Filtering", "b_reports")

# stage owns its outputs, so wipe stale runs (keep .gitkeep)
clear_dir <- function(d) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  unlink(setdiff(list.files(d, full.names = TRUE), file.path(d, ".gitkeep")), recursive = TRUE)
}
clear_dir(data_dir)
clear_dir(report_dir)

#### Load matrix + metadata ####
# Read the protein matrix, derive the 3-group pre/post scheme, align by sample, collapse duplicates.

raw <- read_excel(here("00_input", "CvH_raw.xlsx"))
annot_cols <- c("uniprot_id", "protein", "gene", "description", "n_seq")
annotation <- raw[, annot_cols]
intensity <- raw[, setdiff(names(raw), annot_cols)]

# CvH_meta.csv is the sample sheet, taken as given. Pool supplement arms (CRE/PLA)
# into the 3-group pre/post scheme and align the matrix to it by Col_ID.
metadata <- as.data.frame(read_csv(here("00_input", "CvH_meta.csv"), show_col_types = FALSE))
metadata$group_time <- dplyr::case_when(
  metadata$Group_Time %in% c("CRE_T1", "PLA_T1") ~ "CR_pre",
  metadata$Group_Time %in% c("CRE_T2", "PLA_T2") ~ "CR_post",
  metadata$Group_Time == "H_T1" ~ "H_pre"
)
stopifnot(!anyNA(metadata$group_time))
rownames(metadata) <- metadata$Col_ID
intensity <- intensity[, metadata$Col_ID] # errors if matrix/metadata IDs disagree
n_raw <- nrow(annotation)

if (any(duplicated(annotation$uniprot_id))) { # guard only; keep highest-mean row per accession
  rm_mean <- rowMeans(data.matrix(intensity), na.rm = TRUE)
  keep_idx <- tibble(i = seq_along(rm_mean), id = annotation$uniprot_id, m = rm_mean) |>
    group_by(id) |>
    slice_max(m, n = 1, with_ties = FALSE) |>
    pull(i)
  annotation <- annotation[keep_idx, ]
  intensity <- intensity[keep_idx, ]
}

#### Contaminant removal ####

# Remove a protein iff blood-derived AND NOT myofiber-expressed:
#   blood-derived = HPA "Secreted to blood" | immunoglobulin | erythrocyte-high
#   keep-override = high single-cell Myonuclei signal (genuine muscle protein)
# Myonuclei single-cell RNA is blood-free, so only the dominant hemoglobin band
# (ery_cut) is trusted; low-level platelet RNA is noisy in enucleate cells.
# Source: Human Protein Atlas (Uhlen 2015 Science; 2019 blood atlas).
bl <- read_tsv(here("00_input", "HPA_annotations.tsv"), show_col_types = FALSE) |>
  transmute(
    gene          = Gene,
    uniprot       = Uniprot,
    protein_class = `Protein class`,
    secretome     = `Secretome location`,
    blood_conc    = suppressWarnings(as.numeric(`Blood concentration - Conc. blood MS [pg/L]`)),
    ery           = suppressWarnings(as.numeric(`Single Cell Type RNA - Erythrocytes [nCPM]`)),
    myo           = suppressWarnings(as.numeric(`Single Cell Type RNA - Myonuclei [nCPM]`))
  ) |>
  mutate(
    secreted_blood = secretome == "Secreted to blood" & !is.na(secretome),
    is_ig = str_detect(coalesce(protein_class, ""), "Immunoglobulin genes"),
    is_erythrocyte = !is.na(ery) & ery >= 5000, # nCPM above this = hemoglobin-class red-cell protein
    muscle_keep = !is.na(myo) & myo >= 50 & (is.na(blood_conc) | blood_conc < 1e9), # rescue only if myonuclei-expressed AND not a true high-abundance plasma protein
    reason = case_when(
      muscle_keep & (secreted_blood | is_ig | is_erythrocyte) ~ "keep: muscle-expressed (rescued)",
      is_erythrocyte ~ "remove: erythrocyte (hemoglobin)",
      is_ig ~ "remove: immunoglobulin",
      secreted_blood ~ "remove: secreted-to-blood (plasma)",
      TRUE ~ "keep: not blood-associated"
    ),
    verdict = if_else(str_starts(reason, "remove"), "remove", "keep")
  )
hpa_genes <- unique(bl$gene)
remove_genes <- bl$gene[bl$verdict == "remove"]

# proteoDA::filter_proteins_by_annotation() errors under R >= 4.4 (its internal
# length-2 class() check), so annotation-based removal is done by direct subsetting;
# the group/sample filters below stay native.
flog <- tibble(step = "Raw input", n_after = nrow(annotation), n_removed = NA_integer_)
keep <- annotation$gene %in% hpa_genes
flog <- bind_rows(flog, tibble(step = "HPA presence", n_after = sum(keep), n_removed = sum(!keep)))
annotation <- annotation[keep, ]
intensity <- intensity[keep, ]
keep <- !(annotation$gene %in% remove_genes)
flog <- bind_rows(flog, tibble(step = "Blood contaminant removal", n_after = sum(keep), n_removed = sum(!keep)))
annotation <- annotation[keep, ]
intensity <- intensity[keep, ]

#### Build DAList ####
# Package the cleaned matrix/annotation/metadata into the proteoDA object.

int_mat <- as.data.frame(data.matrix(intensity))
rownames(int_mat) <- annotation$uniprot_id
annot_df <- as.data.frame(annotation)
rownames(annot_df) <- annotation$uniprot_id
meta_df <- as.data.frame(metadata)
rownames(meta_df) <- metadata$Col_ID
dal <- zero_to_missing(DAList(data = int_mat, annotation = annot_df, metadata = meta_df)) # DIA-NN 0 = non-detection

#### Missingness filter ####

n0 <- nrow(dal$data)
# keep a protein quantified in >= 5 samples in at least one group_time level
dal <- filter_proteins_by_group(dal, min_reps = 5, min_groups = 1, grouping_column = "group_time")
flog <- bind_rows(flog, tibble(
  step = "Missingness (>=5 in >=1 group)",
  n_after = nrow(dal$data), n_removed = n0 - nrow(dal$data)
))
flog <- flog |> mutate(pct_of_raw = round(n_after / n_raw * 100, 1))

#### Outlier consensus ####

# Four independent QC heuristics; a sample is dropped only on >= 3/4 agreement.
# Missingness adds a paired arm: a within-subject pre->post jump in NA rate is
# itself suspicious, which a single absolute fence would miss.
m <- dal$metadata
pct_missing <- colMeans(is.na(dal$data)) * 100
paired <- m |>
  filter(Group != "PPS") |>
  count(Subject_ID) |>
  filter(n == 2) |>
  pull(Subject_ID)
delta <- setNames(rep(NA_real_, nrow(m)), m$Col_ID)
for (s in paired) {
  r <- m[m$Subject_ID == s, ]
  t1 <- r$Col_ID[r$Timepoint == "T1"]
  t2 <- r$Col_ID[r$Timepoint == "T2"]
  if (length(t1) == 1 && length(t2) == 1) delta[c(t1, t2)] <- abs(pct_missing[t2] - pct_missing[t1])
}
miss_thr <- quantile(pct_missing, 0.75) + 1.5 * IQR(pct_missing) # Tukey upper fence
dv <- delta[!is.na(delta)]
delta_thr <- if (length(dv) > 2) quantile(dv, 0.75) + 1.5 * IQR(dv) else Inf
miss_flag <- pct_missing > miss_thr | (!is.na(delta) & delta > delta_thr)
sm <- apply(log2(dal$data), 2, median, na.rm = TRUE)
mad_flag <- abs(sm - median(sm)) > 3 * mad(sm) # Hampel 3xMAD on median intensity
pc3 <- prcomp(t(log2(dal$data[rowSums(is.na(dal$data)) == 0, ])), center = TRUE, scale. = TRUE)$x[, 1:3]
pca_flag <- mahalanobis(pc3, colMeans(pc3), cov(pc3)) > qchisq(0.99, df = 3) # chi-sq tail on PC1-3 position
cm <- cor(log2(dal$data), use = "pairwise.complete.obs")
mc <- apply(cm, 2, function(x) median(x[x < 1], na.rm = TRUE)) # x<1 drops the self-correlation
cor_flag <- mc < median(mc) - 3 * mad(mc)
outlier_diag <- tibble(
  Col_ID = colnames(dal$data), miss_flag = miss_flag[Col_ID], mad_flag = mad_flag[Col_ID],
  pca_flag = pca_flag[Col_ID], cor_flag = cor_flag[Col_ID]
) |>
  mutate(n_flags = miss_flag + mad_flag + pca_flag + cor_flag, consensus_outlier = n_flags >= 3)
outlier_ids <- outlier_diag$Col_ID[outlier_diag$consensus_outlier]
cat(sprintf("Outliers (>=3/4): %s\n", if (length(outlier_ids)) paste(outlier_ids, collapse = ", ") else "none"))
if (length(outlier_ids)) dal <- filter_samples(dal, !(Col_ID %in% outlier_ids))

#### Export ####
# Write the filtered DAList for Stage 02, plus a review workbook and the blood-contaminant plot.

saveRDS(dal, file.path(data_dir, "DAList_filtered.rds")) # un-normalized handoff to Stage 02

# gene -> MS protein name + description, so each drop is human-reviewable
ann_lookup <- raw |>
  as.data.frame() |>
  distinct(gene, .keep_all = TRUE) |>
  select(gene, protein, description)

removed <- bl |>
  filter(verdict == "remove", gene %in% raw$gene) |>
  left_join(ann_lookup, by = "gene") |>
  transmute(gene, uniprot, protein, description, reason,
    secretome, blood_conc,
    erythrocyte_nCPM = ery, myonuclei_nCPM = myo
  ) |>
  arrange(reason, desc(blood_conc))

# blood-derived proteins KEPT by the myonuclei rescue, for review
rescued <- bl |>
  filter(str_detect(reason, "rescued"), gene %in% raw$gene) |>
  left_join(ann_lookup, by = "gene") |>
  transmute(gene, uniprot, protein, description, reason,
    secretome, blood_conc,
    erythrocyte_nCPM = ery, myonuclei_nCPM = myo
  ) |>
  arrange(desc(myonuclei_nCPM))

write.xlsx(
  list(
    filter_log = flog,
    contaminants_removed = removed,
    contaminants_rescued = rescued,
    outlier_diagnostics = outlier_diag,
    blood_classification = bl |> select(
      gene, uniprot, protein_class, secretome,
      blood_conc, ery, myo, reason, verdict
    )
  ),
  file.path(data_dir, "filtering_report.xlsx"),
  overwrite = TRUE
)

pal <- c(
  "remove: secreted-to-blood (plasma)" = "#D6604D", "remove: immunoglobulin" = "#F4A582",
  "remove: erythrocyte (hemoglobin)" = "#9970AB"
)
bars <-
  (count(removed, reason) |> ggplot(aes(reorder(reason, n), n, fill = reason)) +
    geom_col(show.legend = FALSE) +
    geom_text(aes(label = n), hjust = -0.2, size = 3.5) +
    scale_fill_manual(values = pal) +
    coord_flip(clip = "off") +
    labs(title = "A  Filtered out, by reason", x = NULL, y = "proteins")) /
    (removed |> filter(!is.na(blood_conc)) |> slice_max(blood_conc, n = 25) |>
      ggplot(aes(fct_reorder(gene, blood_conc), log10(blood_conc + 1), fill = reason)) +
      geom_col(show.legend = FALSE) +
      scale_fill_manual(values = pal) +
      coord_flip() +
      labs(
        title = "B  Top filtered-out plasma proteins by blood concentration",
        x = NULL, y = "log10 blood conc [pg/L]"
      )) +
    plot_layout(heights = c(1, 3)) & theme_minimal(base_size = 11)
ggsave(file.path(report_dir, "blood_contaminants.pdf"), bars, width = 8, height = 9)

if (file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
cat(sprintf("Done: %d proteins x %d samples -> %s/\n", nrow(dal$data), ncol(dal$data), data_dir))
print(as.data.frame(flog))
