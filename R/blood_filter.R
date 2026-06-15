# =============================================================================
# blood_filter.R  --  HPA-derived blood-contaminant classification (skeletal muscle)
#
# Input: a single HPA custom export (00_input/HPA_annotations.tsv) with columns
#   Gene, Uniprot, "Protein class", "Secretome location",
#   "Single Cell Type RNA - Erythrocytes [nCPM]", "... - Myonuclei [nCPM]",
#   (Blood concentration + Skeletal myofibers carried for reporting).
#
# Rule:  remove iff  (blood-derived)  AND NOT (myofiber-expressed)
#   blood-derived  = Secretome "Secreted to blood"            (plasma)
#                  | Protein class "Immunoglobulin genes"     (Ig)
#                  | Erythrocytes nCPM >= ERY_CUT             (hemoglobins)
#   keep-override  = Myonuclei nCPM >= MYO_KEEP               (myofiber / myokine)
#
# Why these signals: Secretome cleanly separates plasma from muscle ECM; the
# single-cell *Myonuclei* signal is blood-free (bulk tissue MS/RNA is biopsy-
# blood-contaminated). Platelet and low-level erythrocyte single-cell RNA are
# noisy (enucleate cells -> housekeeping genes inflate), so only the dominant
# hemoglobin band (ERY_CUT) is used; secreted platelet factors (PPBP/PF4) are
# caught by the plasma arm. Source: HPA (Uhlen 2015 Science; 2019 blood atlas).
# =============================================================================

suppressPackageStartupMessages({ library(readr); library(dplyr); library(stringr) })

ERY_CUT  <- 5000   # Erythrocyte nCPM above this = hemoglobin-class red-cell protein
MYO_KEEP <- 50     # Myonuclei nCPM at/above this = genuine muscle protein (rescued)

classify_blood_reference <- function(hpa_path) {
  read_tsv(hpa_path, show_col_types = FALSE) |>
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
      is_ig          = str_detect(coalesce(protein_class, ""), "Immunoglobulin genes"),
      is_erythrocyte = !is.na(ery) & ery >= ERY_CUT,
      muscle_keep    = !is.na(myo) & myo >= MYO_KEEP,
      reason = case_when(
        muscle_keep & (secreted_blood | is_ig | is_erythrocyte) ~ "keep: muscle-expressed (rescued)",
        is_erythrocyte                                          ~ "remove: erythrocyte (hemoglobin)",
        is_ig                                                   ~ "remove: immunoglobulin",
        secreted_blood                                          ~ "remove: secreted-to-blood (plasma)",
        TRUE                                                    ~ "keep: not blood-associated"
      ),
      verdict = if_else(str_starts(reason, "remove"), "remove", "keep")
    )
}

blood_contaminant_genes <- function(hpa_path) {
  classify_blood_reference(hpa_path) |>
    filter(verdict == "remove") |>
    pull(gene) |>
    unique()
}
