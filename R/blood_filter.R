# =============================================================================
# blood_filter.R  --  HPA-derived blood-contaminant classification (skeletal muscle)
#
# A protein is removed iff it is blood-derived AND not a genuine muscle protein:
#   keep    if antibody-detected in myocytes (IH low/medium/high)
#   remove  if in the canonical erythrocyte (RBC) set
#   remove  if HPA Secretome location == "Secreted to blood" (plasma)
#   keep    otherwise (incl. "Secreted to extracellular matrix" muscle ECM)
#
# Secretome location is the clean discriminator; the blood/muscle MS ratio is NOT
# used (HPA "blood vessel" MS is vessel-WALL tissue and over-removes muscle ECM).
# Source: HPA tissue MS + IH + secretome (Uhlen 2015 Science; 2019 blood atlas).
# =============================================================================

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr)
})

# Canonical erythrocyte (RBC) proteins -- small, unambiguous. CA3 excluded (muscle CA).
ERYTHROCYTE_GENES <- c("HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBQ1", "CA1", "CA2",
                       "PRDX2", "CAT", "BLVRB", "BPGM", "AHSP", "SLC4A1", "SPTA1", "SPTB",
                       "ANK1", "EPB41", "EPB42")

# Classify every protein in the HPA blood-annotation TSV.
# Returns: gene, uniprot, secretome, ih_myocyte, ms_blood, reason, verdict.
classify_blood_reference <- function(hpa_blood_path) {
  read_tsv(hpa_blood_path, show_col_types = FALSE) |>
    transmute(
      gene       = Gene,
      uniprot    = Uniprot,
      secretome  = `Secretome location`,
      ih_myocyte = str_remove_all(`Tissue Cell type Annotation (IH) - skeletal muscle - myocytes`, '"'),
      ms_blood   = suppressWarnings(as.numeric(`Tissue protein MS - blood vessel [Intensity]`))
    ) |>
    mutate(
      reason = case_when(
        ih_myocyte %in% c("low", "medium", "high") ~ "keep: myocyte-detected",
        gene %in% ERYTHROCYTE_GENES                ~ "remove: erythrocyte (RBC)",
        secretome == "Secreted to blood"           ~ "remove: secreted-to-blood (plasma)",
        TRUE                                       ~ "keep: not blood-associated"
      ),
      verdict = if_else(str_starts(reason, "remove"), "remove", "keep")
    )
}

# Convenience: the unique vector of gene symbols to remove.
blood_contaminant_genes <- function(hpa_blood_path) {
  classify_blood_reference(hpa_blood_path) |>
    filter(verdict == "remove") |>
    pull(gene) |>
    unique()
}
