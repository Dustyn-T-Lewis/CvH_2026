# Supp A — per-protein log-intensity scatter for 5 blood markers, Sam vs ours.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

setwd(rprojroot::find_rstudio_root_file())
source("02-03_Sam's_Results/04_Figures/build_data_index.R")
source("02-03_Sam's_Results/04_Figures/shared/style.R")
source("02-03_Sam's_Results/04_Figures/F02_blood_contamination/a_script/_shared_keys.R")

sam     <- readRDS(sam_idx$sam$dalist_rds)
our_dal <- readRDS("01_normalization/c_data/03_DAList_normalized.rds")

blood_markers <- c("HBB", "HBA1", "MB", "ALB", "CKM")

# Sam's annotation has uniprot_id as both rownames and a column. Use the column.
sam_ann <- sam$annotation |> as.data.frame() |>
  filter(gene %in% blood_markers) |>
  select(uniprot_id, gene)

extract_long <- function(dal, ann, source_lab) {
  m <- dal$data
  hits <- intersect(ann$uniprot_id, rownames(m))
  if (length(hits) == 0) return(tibble(uniprot_id = character(), sample_id = character(),
                                       intensity = numeric(), source = character(), gene = character()))
  sub <- m[hits, , drop = FALSE]
  df <- as.data.frame(sub) |>
    tibble::rownames_to_column("uniprot_id") |>
    pivot_longer(-uniprot_id, names_to = "sample_id", values_to = "intensity") |>
    mutate(sample_id = normalize_sample_id(sample_id),
           source = source_lab)
  left_join(df, ann, by = "uniprot_id")
}

sam_long <- extract_long(sam,     sam_ann, "sam")
our_long <- extract_long(our_dal, sam_ann, "ours")

paired <- inner_join(
  sam_long |> select(uniprot_id, gene, sample_id, intensity_sam = intensity),
  our_long |> select(uniprot_id, sample_id, intensity_ours = intensity),
  by = c("uniprot_id", "sample_id")
) |>
  filter(!is.na(intensity_sam), !is.na(intensity_ours))

write_csv(paired,
          "02-03_Sam's_Results/04_Figures/F02_blood_contamination/c_data/supp_A_intensity_pairs.csv")

cors <- paired |>
  group_by(gene) |>
  summarise(r = cor(intensity_ours, intensity_sam, use = "complete.obs"),
            n = n(),
            .groups = "drop") |>
  mutate(label = sprintf("r = %.2f (n=%d)", r, n))

p_supp_A <- ggplot(paired, aes(intensity_ours, intensity_sam)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_text(data = cors, aes(label = label),
            x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = 3, inherit.aes = FALSE) +
  facet_wrap(~ gene, ncol = 3, scales = "free") +
  labs(x = "Our normalized log-intensity",
       y = "Sam's normalized log-intensity",
       title = "Per-protein blood-marker intensity, Sam vs ours") +
  theme_minimal(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey95", color = NA))

ggsave("02-03_Sam's_Results/04_Figures/F02_blood_contamination/b_reports/supp/png/panels/supp_A.png",
       p_supp_A, width = 7, height = 5, dpi = 300, bg = "white")

supp_A <- p_supp_A
