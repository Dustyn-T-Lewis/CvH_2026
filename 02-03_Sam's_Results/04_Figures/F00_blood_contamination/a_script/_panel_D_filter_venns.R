# Panel D — two sub-panels side by side:
#   D(i)  Euler/Venn of kept-protein UniProt sets (Sam vs ours)
#   D(ii) Fate-tile of the 5 named blood markers across each pipeline
#
# Blood markers (HBB, HBA1, ALB) that were filtered before DAList construction
# are detected from the raw input, not from the annotation, so all 5 are always
# represented in D(ii) regardless of which pipeline kept them.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(eulerr)
  library(readxl)
})

setwd(rprojroot::find_rstudio_root_file())
source("02-03_Sam's_Results/04_Figures/build_data_index.R")
source("02-03_Sam's_Results/04_Figures/shared/style.R")

# ---------------------------------------------------------------------------
# 1. Load both filtered protein sets
# ---------------------------------------------------------------------------
sam     <- readRDS(sam_idx$sam$dalist_rds)
our_dal <- readRDS("01_normalization/c_data/03_DAList_normalized.rds")

sam_kept <- rownames(sam$annotation)
our_kept <- rownames(our_dal$annotation)

# ---------------------------------------------------------------------------
# 2. Build protein-set universe and membership table
# ---------------------------------------------------------------------------
universe <- union(sam_kept, our_kept)

# Resolve gene names: Sam first, fall back to ours.
# Both annotations already have uniprot_id as a column (matching rownames).
sam_gene <- setNames(as.data.frame(sam$annotation)$gene,     sam_kept)
our_gene <- setNames(as.data.frame(our_dal$annotation)$gene, our_kept)

panel_D_proteins <- tibble(
  uniprot_id = universe,
  in_sam     = universe %in% sam_kept,
  in_ours    = universe %in% our_kept
) |>
  mutate(
    set_membership = case_when(
       in_sam &  in_ours ~ "both",
       in_sam & !in_ours ~ "sam_only",
      !in_sam &  in_ours ~ "ours_only",
      TRUE               ~ "neither"
    ),
    gene = coalesce(sam_gene[uniprot_id], our_gene[uniprot_id])
  ) |>
  select(uniprot_id, gene, in_sam, in_ours, set_membership)

write_csv(panel_D_proteins,
  "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/panel_D_venn_protein_sets.csv")

# ---------------------------------------------------------------------------
# 3. Build blood-marker fate table
#    Use the raw input as the universe so HBB/HBA1/ALB (filtered before
#    DAList) still appear as rows with in_sam=FALSE / in_ours=FALSE.
# ---------------------------------------------------------------------------
blood_markers <- c("HBB", "HBA1", "MB", "ALB", "CKM")

raw <- read_excel("00_input/CvH_raw.xlsx") |>
  select(uniprot_id, gene) |>
  filter(gene %in% blood_markers)

# Add Sam's blacklist flag where the protein reached Sam's annotation.
# uniprot_id is already a column, so no rownames_to_column needed.
sam_ann <- as.data.frame(sam$annotation) |>
  select(uniprot_id, in_blood_blacklist)

panel_D_markers <- raw |>
  left_join(sam_ann, by = "uniprot_id") |>
  mutate(
    in_sam  = uniprot_id %in% sam_kept,
    in_ours = uniprot_id %in% our_kept
  ) |>
  select(uniprot_id, gene, in_sam, in_ours, in_blood_blacklist) |>
  arrange(match(gene, blood_markers))

# Guarantee exactly 5 rows (one per marker) even if raw lookup missed any
missing_markers <- setdiff(blood_markers, panel_D_markers$gene)
if (length(missing_markers) > 0) {
  panel_D_markers <- bind_rows(
    panel_D_markers,
    tibble(uniprot_id     = NA_character_,
           gene           = missing_markers,
           in_sam         = FALSE,
           in_ours        = FALSE,
           in_blood_blacklist = NA)
  )
}
stopifnot(nrow(panel_D_markers) == 5)

write_csv(panel_D_markers,
  "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/panel_D_venn_blood_markers.csv")

# ---------------------------------------------------------------------------
# 4. D(i) — Euler diagram of kept-protein sets
# ---------------------------------------------------------------------------
n_both     <- sum(panel_D_proteins$set_membership == "both")
n_sam_only <- sum(panel_D_proteins$set_membership == "sam_only")
n_our_only <- sum(panel_D_proteins$set_membership == "ours_only")

euler_fit <- euler(
  c(Sam = n_sam_only, Ours = n_our_only, "Sam&Ours" = n_both)
)

p_venn_proteins <- plot(
  euler_fit,
  fills     = list(fill = c("#1f78b4", "#33a02c"), alpha = 0.5),
  quantities = list(cex = 0.9),
  labels     = list(cex = 0.9),
  main       = list(label = "Kept proteins after filtering", cex = 0.85)
)

# ---------------------------------------------------------------------------
# 5. D(ii) — blood-marker fate tile chart
# ---------------------------------------------------------------------------
marker_long <- panel_D_markers |>
  pivot_longer(c(in_sam, in_ours),
               names_to  = "pipeline",
               values_to = "kept") |>
  mutate(
    pipeline = factor(pipeline,
                      levels = c("in_sam", "in_ours"),
                      labels = c("Sam", "Ours")),
    gene = factor(gene, levels = blood_markers)
  )

p_markers <- ggplot(marker_long, aes(pipeline, gene, fill = kept)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(kept, "✓", "✗")), size = 5) +
  scale_fill_manual(
    values = c(`TRUE` = "#33a02c", `FALSE` = "#e31a1c"),
    guide  = "none"
  ) +
  labs(x = NULL, y = NULL,
       title = "5 blood markers — fate per pipeline") +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid   = element_blank(),
    axis.text    = element_text(size = 9),
    plot.title   = element_text(size = 9)
  )

# ---------------------------------------------------------------------------
# 6. Compose D(i) + D(ii) side by side and save
# ---------------------------------------------------------------------------
p_venn_grob <- cowplot::as_grob(p_venn_proteins)
p_D <- cowplot::plot_grid(
  cowplot::ggdraw(p_venn_grob),
  p_markers,
  ncol       = 2,
  rel_widths = c(2, 1),
  labels     = c("D(i)", "D(ii)"),
  label_size = 10
)

out_png <- "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/png/panels/panel_D.png"
ggsave(out_png, p_D, width = 7, height = 4, dpi = 300, bg = "white")
message("Saved: ", out_png)

# Export for use by a composite driver
panel_D <- p_D
