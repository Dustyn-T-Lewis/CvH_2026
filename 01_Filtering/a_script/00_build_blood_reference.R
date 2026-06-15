# =============================================================================
# 00_build_blood_reference.R  --  classify the HPA blood reference + bar viz
# Thin wrapper over R/blood_filter.R; intersects with observed CvH proteins and
# shows what gets filtered out. See R/blood_filter.R for the rule + rationale.
# =============================================================================

suppressPackageStartupMessages({
  library(here); library(readxl); library(dplyr)
  library(readr); library(ggplot2); library(forcats); library(patchwork)
})
source(here("R", "blood_filter.R"))

ref <- classify_blood_reference(here("00_input", "HPA_annotations.tsv"))
write_tsv(ref, here("01_Filtering", "c_data", "blood_reference.tsv"))

observed <- read_excel(here("00_input", "CvH_raw.xlsx")) |> distinct(gene)
removed   <- inner_join(observed, ref, by = "gene") |> filter(verdict == "remove")
write_csv(select(removed, gene, uniprot, reason, blood_conc, ery, myo),
          here("01_Filtering", "c_data", "blood_contaminants.csv"))
message(sprintf("Filtered %d of %d observed proteins.", nrow(removed), nrow(observed)))

pal <- c("remove: secreted-to-blood (plasma)" = "#D6604D",
         "remove: immunoglobulin"             = "#F4A582",
         "remove: erythrocyte (hemoglobin)"   = "#9970AB")
p <-
  (count(removed, reason) |>
     ggplot(aes(reorder(reason, n), n, fill = reason)) +
     geom_col(show.legend = FALSE) + geom_text(aes(label = n), hjust = -0.2, size = 3.5) +
     scale_fill_manual(values = pal) + coord_flip(clip = "off") +
     labs(title = "A  Filtered out, by reason", x = NULL, y = "proteins")) /
  (removed |> filter(!is.na(blood_conc)) |> slice_max(blood_conc, n = 25) |>
     ggplot(aes(fct_reorder(gene, blood_conc), log10(blood_conc + 1), fill = reason)) +
     geom_col() + scale_fill_manual(values = pal, name = NULL) + coord_flip() +
     labs(title = "B  Top filtered-out plasma proteins by measured blood concentration",
          x = NULL, y = "log10 blood conc [pg/L]")) +
  plot_layout(heights = c(1, 3)) & theme_minimal(base_size = 11)

ggsave(here("01_Filtering", "b_reports", "blood_filter_bars.png"), p, width = 8, height = 9, dpi = 150)
ggsave(here("01_Filtering", "b_reports", "blood_filter_bars.pdf"), p, width = 8, height = 9)
