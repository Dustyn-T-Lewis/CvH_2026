# =============================================================================
# 00_build_blood_reference.R  --  HPA-derived blood-contaminant filter for muscle
#
# A protein is removed iff it is blood-derived AND not a genuine muscle protein:
#   keep    if antibody-detected in myocytes (IH low/medium/high)
#   remove  if in the canonical erythrocyte (RBC) set
#   remove  if HPA Secretome location == "Secreted to blood" (plasma)
#   keep    otherwise  (incl. "Secreted to extracellular matrix" muscle ECM)
#
# Secretome location is the clean signal (separates plasma from muscle ECM). The
# blood/muscle MS ratio is intentionally NOT used: HPA "blood vessel" MS reflects
# vessel-wall tissue, so it over-removes muscle ECM/cytoskeleton.
# Source: HPA tissue MS + IH + secretome (Uhlen 2015 Science; 2019 blood atlas).
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(readxl); library(dplyr)
  library(ggplot2); library(forcats); library(patchwork)
})
source(here("R", "blood_filter.R"))   # classify_blood_reference(), ERYTHROCYTE_GENES

# 1. Classify the HPA blood reference -----------------------------------------
ref <- classify_blood_reference(here("00_input", "HPA_blood_annotations.tsv"))

# validate: muscle markers kept, blood markers removed
stopifnot(
  ref$verdict[match(c("MB", "CKM", "CA3"), ref$gene)]      == "keep",
  ref$verdict[match(c("HBB", "HBA1", "TF", "C3"), ref$gene)] == "remove"
)

# 2. Filter the observed proteins ---------------------------------------------
observed <- read_excel(here("00_input", "CvH_raw.xlsx")) |> distinct(gene)
removed   <- inner_join(observed, ref, by = "gene") |> filter(verdict == "remove")

write_tsv(ref, here("01_Filtering", "c_data", "blood_reference.tsv"))
write_csv(select(removed, gene, uniprot, reason, ms_blood),
          here("01_Filtering", "c_data", "blood_contaminants.csv"))
message(sprintf("Filtered %d of %d observed proteins (%d plasma, %d RBC).",
                nrow(removed), nrow(observed),
                sum(grepl("plasma", removed$reason)), sum(grepl("RBC", removed$reason))))

# 3. Simple bars of what was filtered out -------------------------------------
pal <- c("remove: secreted-to-blood (plasma)" = "#D6604D", "remove: erythrocyte (RBC)" = "#9970AB")
p <-
  (count(removed, reason) |>
     ggplot(aes(reorder(reason, n), n, fill = reason)) +
     geom_col(show.legend = FALSE) + geom_text(aes(label = n), hjust = -0.2, size = 3.5) +
     scale_fill_manual(values = pal) + coord_flip(clip = "off") +
     labs(title = "A  Filtered out, by reason", x = NULL, y = "proteins")) /
  (slice_max(removed, ms_blood, n = 25) |>
     ggplot(aes(fct_reorder(gene, ms_blood), log10(ms_blood + 1), fill = reason)) +
     geom_col() + scale_fill_manual(values = pal, name = NULL) + coord_flip() +
     labs(title = "B  Top filtered out by blood-vessel MS abundance",
          x = NULL, y = "log10 blood MS")) +
  plot_layout(heights = c(1, 3)) & theme_minimal(base_size = 11)

ggsave(here("01_Filtering", "b_reports", "blood_filter_bars.png"), p, width = 8, height = 9, dpi = 150)
ggsave(here("01_Filtering", "b_reports", "blood_filter_bars.pdf"), p, width = 8, height = 9)
