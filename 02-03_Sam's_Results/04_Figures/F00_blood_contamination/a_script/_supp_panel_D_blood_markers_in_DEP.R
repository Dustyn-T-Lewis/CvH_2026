# Supp D — 5 blood markers in our (ours-on-Sam's-data) Cancer_vs_Healthy DEPs.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(scales)
})

setwd(rprojroot::find_rstudio_root_file())
source("02-03_Sam's_Results/04_Figures/build_data_index.R")
source("02-03_Sam's_Results/04_Figures/shared/style.R")

dep <- read_csv(file.path(sam_idx$ours_on_his$per_contrast, "Cancer_vs_Healthy.csv"),
                show_col_types = FALSE)

blood_markers <- c("HBB", "HBA1", "MB", "ALB", "CKM")

gene_col <- intersect(c("gene", "Gene", "gene_symbol", "Symbol"), names(dep))[1]
if (is.na(gene_col)) stop("No gene-symbol column found in DEP table; saw: ", paste(names(dep), collapse = ", "))

supp_D_data <- dep |>
  filter(.data[[gene_col]] %in% blood_markers) |>
  transmute(
    gene = .data[[gene_col]],
    uniprot_id = if ("uniprot_id" %in% names(dep)) uniprot_id
                 else if ("Protein" %in% names(dep)) Protein else NA_character_,
    logFC, P.Value, FDR = adj.P.Val,
    pi_score = P.Value ^ abs(logFC),
    sig_FDR_10 = FDR < 0.10,
    sig_pi_05  = pi_score < 0.05
  ) |>
  arrange(match(gene, blood_markers))

missing <- setdiff(blood_markers, supp_D_data$gene)
if (length(missing) > 0) {
  supp_D_data <- bind_rows(supp_D_data, tibble(gene = missing,
    uniprot_id = NA_character_, logFC = NA, P.Value = NA, FDR = NA,
    pi_score = NA, sig_FDR_10 = NA, sig_pi_05 = NA))
}
supp_D_data$gene <- factor(supp_D_data$gene, levels = blood_markers)

write_csv(supp_D_data,
          "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/supp_D_blood_markers_DEP.csv")

p_supp_D <- ggplot(supp_D_data, aes(logFC, gene)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(aes(size = -log10(P.Value), shape = sig_FDR_10,
                 color = sig_pi_05), na.rm = TRUE) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1, `NA` = 4),
                     name = "FDR < 0.10",
                     labels = c(`FALSE` = "no", `TRUE` = "yes")) +
  scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey40"),
                     name = "π-score < 0.05") +
  scale_size_continuous(range = c(2, 6), name = "-log10 P") +
  labs(x = "logFC (Cancer vs Healthy, ours on Sam's data)",
       y = NULL,
       title = "Blood markers in Cancer_vs_Healthy DEPs",
       caption = "Open circles: not significant at FDR<0.10. X: marker missing from DEP table.") +
  theme_minimal(base_size = 10)

ggsave("02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/supp/png/panels/supp_D.png",
       p_supp_D, width = 5.5, height = 4, dpi = 300, bg = "white")

supp_D <- p_supp_D
