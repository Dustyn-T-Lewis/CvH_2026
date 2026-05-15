# F02 Blood Contamination — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the F02 blood-contamination QC figure (4 main + 4 supp panels) from Sam Cooper's per-sample blood-marker metadata, with cross-validation against our outlier consensus and a filter-intersection panel.

**Architecture:** R + ggplot2 + cowplot + proteoDA. Each panel is its own `_panel_*.R` script that produces (i) a `c_data/*.csv` data artifact and (ii) a `b_reports/{main,supp}/png/panels/*.png` rendering. Two driver scripts (`01_main_panels.R`, `02_supp_panels.R`) source the individual panel scripts. A stitcher (`90_stitch_F02.R`) composites panels into `MAIN_F02_*.pdf` + `SUPP_F02_*.pdf`. ID normalization is centralized in `_shared_keys.R`.

**Tech Stack:** R 4.x, proteoDA, ggplot2, cowplot, ggrepel, eulerr (for Venns), dplyr/tidyr, readr.

**Working directory for all `Rscript` commands:** `/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026/` (the `.Rproj` root that `rprojroot::find_rstudio_root_file()` resolves to).

**Companion docs:** `INTENT.md` (high-level intent), `SPEC.md` (design spec).

---

## File Structure

```
02-03_Sam's_Results/04_Figures/F00_blood_contamination/
├── a_script/
│   ├── INTENT.md                          (exists)
│   ├── SPEC.md                            (exists)
│   ├── PLAN.md                            (this file)
│   ├── _shared_keys.R                     ← Task 1
│   ├── _panel_A_blood_stacked.R           ← Task 2
│   ├── _panel_B_bm_ratio_dotplot.R        ← Task 3
│   ├── _panel_C_bm_vs_mahalanobis.R       ← Task 4
│   ├── _panel_D_filter_venns.R            ← Task 5
│   ├── _supp_panel_A_marker_intensity.R   ← Task 6
│   ├── _supp_panel_B_cutoff_sensitivity.R ← Task 7
│   ├── _supp_panel_C_corr_heatmap.R       ← Task 8
│   ├── _supp_panel_D_blood_markers_in_DEP.R ← Task 9
│   ├── 01_main_panels.R                   ← Task 10
│   ├── 02_supp_panels.R                   ← Task 11
│   └── 90_stitch_F02.R                    ← Task 12
├── b_reports/
│   ├── main/{pdf,png/panels}/
│   └── supp/{pdf,png/panels}/
└── c_data/
    ├── DATA_DICTIONARY.md                 ← Task 13
    └── *.csv (9 files, one per panel)
```

Each `_panel_*.R` is a sourceable script that:
1. Sources `04_Figures/build_data_index.R` for `sam_idx`.
2. Sources `04_Figures/shared/style.R` for palettes/themes.
3. Sources `_shared_keys.R` if it needs ID normalization.
4. Builds the panel's data frame and writes to `c_data/`.
5. Builds a `ggplot` object and saves a PNG to `b_reports/.../png/panels/`.
6. Exposes the plot as an object via `assign("panel_X", p, envir = .GlobalEnv)` so the stitcher can compose.

---

## Task 1: ID-normalization helper (`_shared_keys.R`)

**Files:**
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_shared_keys.R`

- [ ] **Step 1: Write the helper**

Create file with:

```r
# Sample-ID normalization for F02.
# Sam uses zero-padded 3-digit subject IDs (CR006_T1); our pipeline uses
# CR6_T1. Always canonicalize to the short form before joining.

normalize_sample_id <- function(x) {
  sub("^CR0*([0-9]+)_T([12])$", "CR\\1_T\\2", x)
}

# Quick self-check (only runs if sourced as main).
if (sys.nframe() == 0) {
  stopifnot(
    normalize_sample_id("CR006_T1") == "CR6_T1",
    normalize_sample_id("CR6_T1")   == "CR6_T1",
    normalize_sample_id("CR10_T2")  == "CR10_T2",
    normalize_sample_id("CR017_T2") == "CR17_T2"
  )
  message("_shared_keys.R: self-check passed")
}
```

- [ ] **Step 2: Run self-check**

```bash
Rscript "A_CvH_2026/02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_shared_keys.R"
```

Expected output: `_shared_keys.R: self-check passed`

- [ ] **Step 3: Commit**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
git add "02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_shared_keys.R"
git commit -m "feat(F02): add sample-ID normalization helper"
```

---

## Task 2: Panel A — per-sample blood-marker stacked bar

**Files:**
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_panel_A_blood_stacked.R`
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/panel_A_stacked_data.csv`
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/png/panels/panel_A.png`

- [ ] **Step 1: Write the panel script**

```r
# Panel A — per-sample stacked bar of blood-marker percentages.
# Bars sorted ascending by B_M_ratio (left = cleanest, right = most contaminated).
# Stack 5 markers only (HBB/HBA1/MB/ALB/CKM); muscle/other left implicit.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
  library(RColorBrewer)
})

setwd(rprojroot::find_rstudio_root_file())
source("02-03_Sam's_Results/04_Figures/build_data_index.R")
source("02-03_Sam's_Results/04_Figures/shared/style.R")

sam <- readRDS(sam_idx$sam$dalist_rds)
meta <- as.data.frame(sam$metadata)

# Build Group_Time from supp_time, falling back to cancer_time for control
# samples (Sam's CTL rows have supp = NA, so supp_time is NA_T1; map to
# the cancer_time value instead, e.g. CTL_T1).
meta$Group_Time <- ifelse(
  is.na(meta$supp) | meta$supp == "" | grepl("^NA_", meta$supp_time),
  as.character(meta$cancer_time),
  as.character(meta$supp_time)
)

panel_A_data <- meta |>
  select(sample_id, Group_Time, HBB_pct, HBA1_pct, MB_pct,
         ALB_pct, CKM_pct, B_M_ratio) |>
  pivot_longer(c(HBB_pct, HBA1_pct, MB_pct, ALB_pct, CKM_pct),
               names_to = "marker", values_to = "pct") |>
  mutate(marker = sub("_pct$", "", marker))

bm_order <- meta |> arrange(B_M_ratio) |> pull(sample_id)
panel_A_data$sample_id <- factor(panel_A_data$sample_id, levels = bm_order)

gt_levels <- intersect(c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2",
                         "CTL_T1", "CTL_T2"), unique(panel_A_data$Group_Time))
panel_A_data$Group_Time <- factor(panel_A_data$Group_Time, levels = gt_levels)

write_csv(panel_A_data, "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/panel_A_stacked_data.csv")

marker_palette <- setNames(brewer.pal(5, "Set2"),
                           c("HBB", "HBA1", "MB", "ALB", "CKM"))

p_A <- ggplot(panel_A_data, aes(sample_id, pct, fill = marker)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = marker_palette, name = "Marker") +
  labs(x = "Sample (sorted by B_M_ratio, ascending)",
       y = "% of total signal",
       title = "Per-sample blood-marker share") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        legend.position = "right")

ggsave("02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/png/panels/panel_A.png",
       p_A, width = 7, height = 4, dpi = 300, bg = "white")

panel_A <- p_A
```

- [ ] **Step 2: Run the panel script**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis
Rscript "A_CvH_2026/02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_panel_A_blood_stacked.R"
```

Expected: no error; PNG written.

- [ ] **Step 3: Verify outputs**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
test -f "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/panel_A_stacked_data.csv" && echo CSV_OK
test -f "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/png/panels/panel_A.png" && echo PNG_OK
Rscript -e 'd <- read.csv("02-03_Sam'"'"'s_Results/04_Figures/F00_blood_contamination/c_data/panel_A_stacked_data.csv"); stopifnot(nrow(d) == 35 * 5, all(c("sample_id","Group_Time","marker","pct","B_M_ratio") %in% names(d))); cat("rows:", nrow(d), "cols:", ncol(d), "\n")'
```

Expected: `CSV_OK`, `PNG_OK`, `rows: 175 cols: 5`.

- [ ] **Step 4: Commit**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
git add "02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_panel_A_blood_stacked.R" \
        "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/panel_A_stacked_data.csv" \
        "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/png/panels/panel_A.png"
git commit -m "feat(F02): panel A — per-sample blood-marker stacked bar"
```

---

## Task 3: Panel B — B_M_ratio dotplot by Group_Time

**Files:**
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_panel_B_bm_ratio_dotplot.R`
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/panel_B_dotplot_data.csv`
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/png/panels/panel_B.png`

- [ ] **Step 1: Write the panel script**

```r
# Panel B — B_M_ratio dotplot by Group_Time, cohort top-decile (top 4 of 35)
# colored red and labeled with sample_id.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(ggrepel)
})

setwd(rprojroot::find_rstudio_root_file())
source("02-03_Sam's_Results/04_Figures/build_data_index.R")
source("02-03_Sam's_Results/04_Figures/shared/style.R")

sam <- readRDS(sam_idx$sam$dalist_rds)
meta <- as.data.frame(sam$metadata)

# Group_Time derivation mirrors panel A: supp_time for SURV samples,
# cancer_time for CTL fallback.
meta$Group_Time <- ifelse(
  is.na(meta$supp) | meta$supp == "" | grepl("^NA_", meta$supp_time),
  as.character(meta$cancer_time),
  as.character(meta$supp_time)
)

bm_cut <- quantile(meta$B_M_ratio, 0.90, na.rm = TRUE)
gt_levels <- intersect(c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2",
                         "CTL_T1", "CTL_T2"), unique(meta$Group_Time))

panel_B_data <- meta |>
  transmute(sample_id,
            Group_Time = factor(Group_Time, levels = gt_levels),
            B_M_ratio,
            top_decile_flag = B_M_ratio >= bm_cut)

write_csv(panel_B_data,
          "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/panel_B_dotplot_data.csv")

kw_p <- kruskal.test(B_M_ratio ~ Group_Time, data = panel_B_data)$p.value
kw_label <- sprintf("Kruskal–Wallis p = %.3g", kw_p)

# Extend GROUP_COLORS with CTL_T1/CTL_T2 = the same green as H_T1.
group_palette <- c(GROUP_COLORS,
                   CTL_T1 = unname(GROUP_COLORS["H_T1"]),
                   CTL_T2 = unname(GROUP_COLORS["H_T1"]))

p_B <- ggplot(panel_B_data, aes(Group_Time, B_M_ratio)) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1),
               geom = "errorbar", width = 0.2, color = "grey50") +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, color = "grey30") +
  geom_jitter(aes(color = Group_Time, alpha = top_decile_flag,
                  size = top_decile_flag), width = 0.15, height = 0) +
  geom_text_repel(data = filter(panel_B_data, top_decile_flag),
                  aes(label = sample_id), size = 2.6, max.overlaps = Inf,
                  box.padding = 0.5, color = "firebrick") +
  scale_color_manual(values = group_palette, guide = "none") +
  scale_alpha_manual(values = c(`FALSE` = 0.7, `TRUE` = 1), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 1.8, `TRUE` = 2.6), guide = "none") +
  annotate("text", x = 0.7, y = max(panel_B_data$B_M_ratio, na.rm = TRUE),
           label = kw_label, hjust = 0, size = 3) +
  labs(x = "Group × Timepoint",
       y = "B_M_ratio (blood-to-muscle signal ratio)",
       title = "B_M_ratio by Group_Time (top-decile labeled)") +
  theme_minimal(base_size = 10)

ggsave("02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/png/panels/panel_B.png",
       p_B, width = 5.5, height = 4, dpi = 300, bg = "white")

panel_B <- p_B
```

- [ ] **Step 2: Run the panel script**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis
Rscript "A_CvH_2026/02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_panel_B_bm_ratio_dotplot.R"
```

Expected: no error.

- [ ] **Step 3: Verify outputs**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
Rscript -e 'd <- read.csv("02-03_Sam'"'"'s_Results/04_Figures/F00_blood_contamination/c_data/panel_B_dotplot_data.csv"); stopifnot(nrow(d) == 35, sum(d$top_decile_flag) %in% 3:5); cat("rows:", nrow(d), "top_decile:", sum(d$top_decile_flag), "\n")'
test -f "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/png/panels/panel_B.png" && echo PNG_OK
```

Expected: `rows: 35 top_decile: 4` (or 3/5 if ties), `PNG_OK`.

- [ ] **Step 4: Commit**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
git add "02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_panel_B_bm_ratio_dotplot.R" \
        "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/panel_B_dotplot_data.csv" \
        "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/png/panels/panel_B.png"
git commit -m "feat(F02): panel B — B_M_ratio dotplot with top-decile highlight"
```

---

## Task 4: Panel C — B_M_ratio vs our Mahalanobis distance

**Files:**
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_panel_C_bm_vs_mahalanobis.R`
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/panel_C_join_data.csv`
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/png/panels/panel_C.png`

- [ ] **Step 1: Write the panel script**

```r
# Panel C — cross-pipeline scatter: Sam's B_M_ratio vs our Mahalanobis distance.
# Color = consensus_outlier from our 4-method 01_normalization outlier_diag.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(ggrepel)
})

setwd(rprojroot::find_rstudio_root_file())
source("02-03_Sam's_Results/04_Figures/build_data_index.R")
source("02-03_Sam's_Results/04_Figures/shared/style.R")
source("02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_shared_keys.R")

sam <- readRDS(sam_idx$sam$dalist_rds)
sam_meta <- as.data.frame(sam$metadata) |>
  mutate(sample_id = normalize_sample_id(sample_id))

our_intermediates <- readRDS("01_normalization/c_data/00_report_intermediates.rds")
our_outlier <- our_intermediates$outlier_diag |>
  mutate(sample_id = normalize_sample_id(Col_ID))

joined <- inner_join(
  sam_meta |> select(sample_id, Group_Time = supp_time, B_M_ratio),
  our_outlier |> select(sample_id, mahal_dist, n_flags, consensus_outlier),
  by = "sample_id"
)

bm_cut <- quantile(joined$B_M_ratio, 0.90, na.rm = TRUE)
joined$label_this <- joined$consensus_outlier | joined$B_M_ratio >= bm_cut

write_csv(joined,
          "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/panel_C_join_data.csv")

sp <- suppressWarnings(cor.test(joined$B_M_ratio, joined$mahal_dist,
                                method = "spearman", exact = FALSE))
pe <- cor.test(joined$B_M_ratio, joined$mahal_dist, method = "pearson")
stats_label <- sprintf("Spearman ρ = %.2f (p = %.2g)\nPearson r = %.2f (p = %.2g)\nn = %d",
                      sp$estimate, sp$p.value, pe$estimate, pe$p.value, nrow(joined))

p_C <- ggplot(joined, aes(B_M_ratio, mahal_dist,
                          color = consensus_outlier,
                          shape = consensus_outlier)) +
  geom_point(size = 2.5) +
  geom_text_repel(data = filter(joined, label_this),
                  aes(label = sample_id), size = 2.6, max.overlaps = Inf,
                  box.padding = 0.5) +
  scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey40"),
                     name = "Our 4-method\nconsensus_outlier") +
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16),
                     name = "Our 4-method\nconsensus_outlier") +
  annotate("text", x = -Inf, y = Inf, label = stats_label,
           hjust = -0.1, vjust = 1.2, size = 3) +
  labs(x = "Sam's B_M_ratio",
       y = "Our Mahalanobis distance (PC1–2)",
       title = "Cross-pipeline contamination agreement") +
  theme_minimal(base_size = 10)

ggsave("02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/png/panels/panel_C.png",
       p_C, width = 5.5, height = 4, dpi = 300, bg = "white")

panel_C <- p_C
```

- [ ] **Step 2: Run the panel script**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis
Rscript "A_CvH_2026/02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_panel_C_bm_vs_mahalanobis.R"
```

Expected: no error.

- [ ] **Step 3: Verify outputs**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
Rscript -e 'd <- read.csv("02-03_Sam'"'"'s_Results/04_Figures/F00_blood_contamination/c_data/panel_C_join_data.csv"); stopifnot(nrow(d) >= 30, all(c("sample_id","B_M_ratio","mahal_dist","n_flags","consensus_outlier","label_this") %in% names(d))); cat("joined n:", nrow(d), "consensus_outliers in join:", sum(d$consensus_outlier), "\n")'
test -f "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/png/panels/panel_C.png" && echo PNG_OK
```

Expected: `joined n: 32` to `35`, `consensus_outliers in join: 0` to `3`, `PNG_OK`.

If joined n < 30, halt — investigate the ID-normalization or cohort mismatch before continuing.

- [ ] **Step 4: Commit**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
git add "02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_panel_C_bm_vs_mahalanobis.R" \
        "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/panel_C_join_data.csv" \
        "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/png/panels/panel_C.png"
git commit -m "feat(F02): panel C — B_M_ratio vs our Mahalanobis distance"
```

---

## Task 5: Panel D — filter intersection (two Venns)

**Files:**
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_panel_D_filter_venns.R`
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/panel_D_venn_protein_sets.csv`
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/panel_D_venn_blood_markers.csv`
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/png/panels/panel_D.png`

- [ ] **Step 1: Write the panel script**

```r
# Panel D — two Venns side by side:
#   D(i)  kept-protein UniProt sets (Sam vs ours)
#   D(ii) fate of the 5 named blood markers per pipeline

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(cowplot)
  library(eulerr)
})

setwd(rprojroot::find_rstudio_root_file())
source("02-03_Sam's_Results/04_Figures/build_data_index.R")
source("02-03_Sam's_Results/04_Figures/shared/style.R")

sam <- readRDS(sam_idx$sam$dalist_rds)
sam_kept <- rownames(sam$annotation)

our_dal <- readRDS("01_normalization/c_data/03_DAList_normalized.rds")
our_kept <- rownames(our_dal$annotation)

universe <- union(sam_kept, our_kept)
panel_D_proteins <- tibble(
  uniprot_id = universe,
  in_sam     = universe %in% sam_kept,
  in_ours    = universe %in% our_kept
) |>
  mutate(set_membership = case_when(
    in_sam &  in_ours ~ "both",
    in_sam & !in_ours ~ "sam_only",
   !in_sam &  in_ours ~ "ours_only",
    TRUE              ~ "neither"
  )) |>
  left_join(
    sam$annotation |> as.data.frame() |>
      tibble::rownames_to_column("uniprot_id") |>
      select(uniprot_id, gene_sam = gene),
    by = "uniprot_id"
  ) |>
  left_join(
    our_dal$annotation |> as.data.frame() |>
      tibble::rownames_to_column("uniprot_id") |>
      select(uniprot_id, any_of(c("gene_symbol", "Gene", "gene"))),
    by = "uniprot_id"
  )

# Resolve gene name (Sam first, then ours if missing).
panel_D_proteins$gene <- coalesce(
  panel_D_proteins$gene_sam,
  panel_D_proteins[[intersect(c("gene_symbol", "Gene", "gene"),
                              names(panel_D_proteins))[1]]]
)

write_csv(panel_D_proteins |> select(uniprot_id, gene, in_sam, in_ours, set_membership),
          "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/panel_D_venn_protein_sets.csv")

n_both     <- sum(panel_D_proteins$set_membership == "both")
n_sam_only <- sum(panel_D_proteins$set_membership == "sam_only")
n_our_only <- sum(panel_D_proteins$set_membership == "ours_only")

euler_fit <- euler(c(Sam = n_sam_only, Ours = n_our_only, "Sam&Ours" = n_both))
p_venn_proteins <- plot(euler_fit, fills = c("#1f78b4", "#33a02c"),
                        quantities = TRUE, alpha = 0.5,
                        main = "Kept proteins after filtering")

# D(ii) — 5 blood markers
blood_markers <- c("HBB", "HBA1", "MB", "ALB", "CKM")
sam_ann <- sam$annotation |> as.data.frame() |>
  tibble::rownames_to_column("uniprot_id")

marker_rows <- sam_ann |>
  filter(gene %in% blood_markers) |>
  select(uniprot_id, gene, in_blood_blacklist)

panel_D_markers <- marker_rows |>
  mutate(in_sam  = uniprot_id %in% sam_kept,
         in_ours = uniprot_id %in% our_kept) |>
  select(uniprot_id, gene, in_sam, in_ours, in_blood_blacklist)

# If a marker doesn't appear in Sam's annotation universe, append a stub row
# so the table is always 5 rows.
missing_markers <- setdiff(blood_markers, panel_D_markers$gene)
if (length(missing_markers) > 0) {
  panel_D_markers <- bind_rows(
    panel_D_markers,
    tibble(uniprot_id = NA_character_, gene = missing_markers,
           in_sam = FALSE, in_ours = FALSE, in_blood_blacklist = NA)
  )
}

write_csv(panel_D_markers,
          "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/panel_D_venn_blood_markers.csv")

# Render markers fate as a small indicator chart (5 points are too few for a Venn).
marker_long <- panel_D_markers |>
  tidyr::pivot_longer(c(in_sam, in_ours), names_to = "pipeline", values_to = "kept") |>
  mutate(pipeline = factor(pipeline, levels = c("in_sam", "in_ours"),
                           labels = c("Sam", "Ours")),
         gene = factor(gene, levels = blood_markers))

p_markers <- ggplot(marker_long, aes(pipeline, gene, fill = kept)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(kept, "✓", "✗")), size = 5) +
  scale_fill_manual(values = c(`TRUE` = "#33a02c", `FALSE` = "#e31a1c"),
                    guide = "none") +
  labs(x = NULL, y = NULL, title = "5 blood markers — fate") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 9))

# Convert the eulerr base plot into a ggplot-compatible grob.
p_venn_grob <- cowplot::ggdraw(p_venn_proteins)
p_D <- cowplot::plot_grid(p_venn_grob, p_markers, ncol = 2,
                          rel_widths = c(2, 1))

ggsave("02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/png/panels/panel_D.png",
       p_D, width = 7, height = 4, dpi = 300, bg = "white")

panel_D <- p_D
```

- [ ] **Step 2: Run the panel script**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis
Rscript "A_CvH_2026/02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_panel_D_filter_venns.R"
```

Expected: no error. If `eulerr` is not installed: `R -e 'install.packages("eulerr", repos="https://cloud.r-project.org")'`.

- [ ] **Step 3: Verify outputs**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
Rscript -e 'd <- read.csv("02-03_Sam'"'"'s_Results/04_Figures/F00_blood_contamination/c_data/panel_D_venn_protein_sets.csv"); cat("rows:", nrow(d), "  both:", sum(d$set_membership=="both"), "  sam_only:", sum(d$set_membership=="sam_only"), "  ours_only:", sum(d$set_membership=="ours_only"), "\n")'
Rscript -e 'm <- read.csv("02-03_Sam'"'"'s_Results/04_Figures/F00_blood_contamination/c_data/panel_D_venn_blood_markers.csv"); stopifnot(nrow(m) == 5); print(m)'
test -f "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/png/panels/panel_D.png" && echo PNG_OK
```

Expected: protein-set rows ≥ 2000 with non-zero counts in all three set categories; 5 marker rows printed.

- [ ] **Step 4: Commit**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
git add "02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_panel_D_filter_venns.R" \
        "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/panel_D_venn_protein_sets.csv" \
        "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/panel_D_venn_blood_markers.csv" \
        "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/png/panels/panel_D.png"
git commit -m "feat(F02): panel D — filter intersection (two Venns)"
```

---

## Task 6: Supp A — 5-marker intensity scatter (Sam vs ours)

**Files:**
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_supp_panel_A_marker_intensity.R`
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/supp_A_intensity_pairs.csv`
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/supp/png/panels/supp_A.png`

- [ ] **Step 1: Write the panel script**

```r
# Supp A — per-protein log-intensity scatter for 5 blood markers, Sam vs ours.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

setwd(rprojroot::find_rstudio_root_file())
source("02-03_Sam's_Results/04_Figures/build_data_index.R")
source("02-03_Sam's_Results/04_Figures/shared/style.R")
source("02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_shared_keys.R")

sam     <- readRDS(sam_idx$sam$dalist_rds)
our_dal <- readRDS("01_normalization/c_data/03_DAList_normalized.rds")

blood_markers <- c("HBB", "HBA1", "MB", "ALB", "CKM")

sam_ann <- sam$annotation |> as.data.frame() |>
  tibble::rownames_to_column("uniprot_id") |>
  filter(gene %in% blood_markers) |>
  select(uniprot_id, gene)

extract_long <- function(dal, ann, source_lab) {
  m <- dal$data
  hits <- intersect(ann$uniprot_id, rownames(m))
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
          "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/supp_A_intensity_pairs.csv")

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

ggsave("02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/supp/png/panels/supp_A.png",
       p_supp_A, width = 7, height = 5, dpi = 300, bg = "white")

supp_A <- p_supp_A
```

- [ ] **Step 2: Run the panel script**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis
Rscript "A_CvH_2026/02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_supp_panel_A_marker_intensity.R"
```

Expected: no error.

- [ ] **Step 3: Verify outputs**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
Rscript -e 'd <- read.csv("02-03_Sam'"'"'s_Results/04_Figures/F00_blood_contamination/c_data/supp_A_intensity_pairs.csv"); cat("rows:", nrow(d), "  unique genes:", length(unique(d$gene)), "\n"); print(table(d$gene))'
test -f "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/supp/png/panels/supp_A.png" && echo PNG_OK
```

Expected: rows ≥ 100, unique genes between 3 and 5 (some markers may be absent from our DAList).

- [ ] **Step 4: Commit**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
git add "02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_supp_panel_A_marker_intensity.R" \
        "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/supp_A_intensity_pairs.csv" \
        "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/supp/png/panels/supp_A.png"
git commit -m "feat(F02): supp A — 5-marker intensity scatter Sam vs ours"
```

---

## Task 7: Supp B — B_M_ratio cutoff sensitivity on DEP count

**Files:**
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_supp_panel_B_cutoff_sensitivity.R`
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/supp_B_cutoff_sensitivity.csv`
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/supp/png/panels/supp_B.png`

- [ ] **Step 1: Inspect Sam's DAList contrast configuration**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
Rscript -e '
sam <- readRDS("02-03_Sam'"'"'s_Results/00_input/01_normalized_DAList_SURV_stringent_muscle.RDS")
cat("Has $design:", !is.null(sam$design), "\n")
if (!is.null(sam$design)) { cat("Design cols:\n"); print(colnames(sam$design)) }
cat("Has $eBayes_fit:", !is.null(sam$eBayes_fit), "\n")
cat("Has $results:", !is.null(sam$results), "\n")
if (!is.null(sam$results)) print(names(sam$results))
'
```

Expected: confirms whether refit needs design rebuild or can reuse Sam's existing design.

- [ ] **Step 2: Write the panel script**

```r
# Supp B — B_M_ratio cutoff sensitivity for Cancer_vs_Healthy DEP count.
# Sweep top-X% B_M_ratio sample drops and refit limma.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(proteoDA)
})

setwd(rprojroot::find_rstudio_root_file())
source("02-03_Sam's_Results/04_Figures/build_data_index.R")
source("02-03_Sam's_Results/04_Figures/shared/style.R")

sam <- readRDS(sam_idx$sam$dalist_rds)
meta <- as.data.frame(sam$metadata)
fdr_thresh <- 0.10

# Inspect cancer_time levels to build a correct contrast for Sam's cohort.
ct_levels <- sort(unique(as.character(meta$cancer_time)))
surv_levels <- grep("^SURV", ct_levels, value = TRUE)
ctl_levels  <- grep("^CTL",  ct_levels, value = TRUE)
stopifnot(length(surv_levels) >= 1, length(ctl_levels) >= 1)
surv_expr <- if (length(surv_levels) == 1) surv_levels
             else sprintf("(%s)/%d", paste(surv_levels, collapse = " + "),
                                     length(surv_levels))
ctl_expr  <- if (length(ctl_levels) == 1) ctl_levels
             else sprintf("(%s)/%d", paste(ctl_levels, collapse = " + "),
                                     length(ctl_levels))
contrast_expr <- sprintf("%s - %s", surv_expr, ctl_expr)
message("Supp B contrast: Cancer_vs_Healthy = ", contrast_expr)

fit_for_subset <- function(dal_sub) {
  dal_sub <- add_design(dal_sub, design_formula = ~0 + cancer_time + sex)
  dal_sub <- add_contrasts(dal_sub,
                           contrasts_vector = c(Cancer_vs_Healthy = contrast_expr))
  dal_sub <- fit_limma_model(dal_sub, robust = TRUE, trend = TRUE)
  res <- dal_sub$results$statistical_results$Cancer_vs_Healthy
  res$FDR <- p.adjust(res$P.Value, method = "BH")
  sum(res$FDR < fdr_thresh, na.rm = TRUE)
}

cutoffs <- c(0, 5, 10, 15, 20)
sweep <- lapply(cutoffs, function(p) {
  n_drop <- ceiling(nrow(meta) * p / 100)
  drop_ids <- if (n_drop == 0) character(0) else
    meta |> arrange(desc(B_M_ratio)) |> slice_head(n = n_drop) |> pull(sample_id)
  keep_ids <- setdiff(meta$sample_id, drop_ids)
  dal_sub <- sam
  dal_sub$data <- dal_sub$data[, keep_ids, drop = FALSE]
  dal_sub$metadata <- dal_sub$metadata[keep_ids, , drop = FALSE]
  dal_sub$design <- NULL
  dal_sub$eBayes_fit <- NULL
  dal_sub$results <- NULL
  n_DEP <- tryCatch(fit_for_subset(dal_sub),
                    error = function(e) { warning(conditionMessage(e)); NA_integer_ })
  data.frame(
    cutoff_pct = p,
    n_samples_dropped = n_drop,
    n_proteins_used = nrow(dal_sub$data),
    n_DEP_fdr10 = n_DEP,
    dropped_sample_ids = paste(drop_ids, collapse = ";")
  )
})
supp_B_data <- do.call(rbind, sweep)

write_csv(supp_B_data,
          "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/supp_B_cutoff_sensitivity.csv")

baseline <- supp_B_data$n_DEP_fdr10[supp_B_data$cutoff_pct == 0]

p_supp_B <- ggplot(supp_B_data, aes(cutoff_pct, n_DEP_fdr10)) +
  geom_hline(yintercept = baseline, linetype = "dashed", color = "grey50") +
  geom_line(color = "steelblue") +
  geom_point(color = "steelblue", size = 2.5) +
  geom_text(aes(label = n_DEP_fdr10), vjust = -1, size = 3) +
  labs(x = "% of top-B_M_ratio samples dropped",
       y = "Cancer_vs_Healthy DEPs at FDR < 0.10",
       title = "DEP count sensitivity to contamination cutoff",
       caption = paste0("Dashed = baseline (cutoff = 0). Refit via proteoDA::fit_limma_model with robust+trend.")) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20)) +
  theme_minimal(base_size = 10)

ggsave("02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/supp/png/panels/supp_B.png",
       p_supp_B, width = 5.5, height = 4, dpi = 300, bg = "white")

supp_B <- p_supp_B
```

- [ ] **Step 3: Run the panel script**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis
Rscript "A_CvH_2026/02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_supp_panel_B_cutoff_sensitivity.R"
```

Expected: no error. Runtime ~30 s for 5 refits.

If `cancer_time` factor levels differ from `SURV_T1/T2/CTL_T1/T2` (e.g., levels are `SURV` and `CTL` with no time), inspect `levels(meta$cancer_time)` and adjust the `add_contrasts()` vector accordingly before re-running.

- [ ] **Step 4: Verify outputs**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
Rscript -e 'd <- read.csv("02-03_Sam'"'"'s_Results/04_Figures/F00_blood_contamination/c_data/supp_B_cutoff_sensitivity.csv"); print(d); stopifnot(nrow(d) == 5, all(d$n_DEP_fdr10 >= 0, na.rm = TRUE))'
test -f "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/supp/png/panels/supp_B.png" && echo PNG_OK
```

Expected: 5 rows, baseline `n_DEP_fdr10` between 100 and 800 (per Sam's known result range).

- [ ] **Step 5: Commit**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
git add "02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_supp_panel_B_cutoff_sensitivity.R" \
        "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/supp_B_cutoff_sensitivity.csv" \
        "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/supp/png/panels/supp_B.png"
git commit -m "feat(F02): supp B — B_M_ratio cutoff sensitivity on DEP count"
```

---

## Task 8: Supp C — correlation heatmap of blood markers

**Files:**
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_supp_panel_C_corr_heatmap.R`
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/supp_C_corr_matrix.csv`
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/supp/png/panels/supp_C.png`

- [ ] **Step 1: Write the panel script**

```r
# Supp C — 6x6 Spearman correlation heatmap of blood markers + B_M_ratio.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

setwd(rprojroot::find_rstudio_root_file())
source("02-03_Sam's_Results/04_Figures/build_data_index.R")
source("02-03_Sam's_Results/04_Figures/shared/style.R")

sam <- readRDS(sam_idx$sam$dalist_rds)
meta <- as.data.frame(sam$metadata)

vars <- c("HBB_pct", "HBA1_pct", "MB_pct", "ALB_pct", "CKM_pct", "B_M_ratio")
m <- meta[, vars]

cor_mat <- cor(m, method = "spearman", use = "pairwise.complete.obs")
p_mat <- outer(seq_along(vars), seq_along(vars), Vectorize(function(i, j) {
  if (i == j) return(NA_real_)
  suppressWarnings(cor.test(m[[i]], m[[j]], method = "spearman", exact = FALSE)$p.value)
}))
rownames(p_mat) <- colnames(p_mat) <- vars

cor_long <- as.data.frame(as.table(cor_mat)) |>
  rename(var1 = Var1, var2 = Var2, rho = Freq) |>
  mutate(p = as.vector(p_mat))

write_csv(cor_long,
          "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/supp_C_corr_matrix.csv")

cor_long$var1 <- factor(cor_long$var1, levels = vars)
cor_long$var2 <- factor(cor_long$var2, levels = rev(vars))

p_supp_C <- ggplot(cor_long, aes(var1, var2, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman ρ") +
  labs(x = NULL, y = NULL,
       title = "Inter-marker + B_M_ratio correlation (Spearman, N=35)") +
  coord_equal() +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

ggsave("02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/supp/png/panels/supp_C.png",
       p_supp_C, width = 5, height = 4.5, dpi = 300, bg = "white")

supp_C <- p_supp_C
```

- [ ] **Step 2: Run the panel script**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis
Rscript "A_CvH_2026/02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_supp_panel_C_corr_heatmap.R"
```

Expected: no error.

- [ ] **Step 3: Verify outputs**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
Rscript -e 'd <- read.csv("02-03_Sam'"'"'s_Results/04_Figures/F00_blood_contamination/c_data/supp_C_corr_matrix.csv"); stopifnot(nrow(d) == 36); cat("rows:", nrow(d), "\n")'
test -f "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/supp/png/panels/supp_C.png" && echo PNG_OK
```

Expected: 36 rows (6×6), `PNG_OK`.

- [ ] **Step 4: Commit**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
git add "02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_supp_panel_C_corr_heatmap.R" \
        "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/supp_C_corr_matrix.csv" \
        "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/supp/png/panels/supp_C.png"
git commit -m "feat(F02): supp C — blood-marker correlation heatmap"
```

---

## Task 9: Supp D — blood markers in our DEP results

**Files:**
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_supp_panel_D_blood_markers_in_DEP.R`
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/supp_D_blood_markers_DEP.csv`
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/supp/png/panels/supp_D.png`

- [ ] **Step 1: Write the panel script**

```r
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

gene_col <- intersect(c("Gene", "gene", "gene_symbol", "Symbol"), names(dep))[1]
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
```

- [ ] **Step 2: Run the panel script**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis
Rscript "A_CvH_2026/02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_supp_panel_D_blood_markers_in_DEP.R"
```

Expected: no error. If gene column isn't named `Gene`, the script will halt with a clear "No gene-symbol column found" message — inspect the DEP CSV header and update `gene_col` resolution.

- [ ] **Step 3: Verify outputs**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
Rscript -e 'd <- read.csv("02-03_Sam'"'"'s_Results/04_Figures/F00_blood_contamination/c_data/supp_D_blood_markers_DEP.csv"); stopifnot(nrow(d) == 5); print(d)'
test -f "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/supp/png/panels/supp_D.png" && echo PNG_OK
```

Expected: 5 rows printed (some may have NA values if missing from DEP).

- [ ] **Step 4: Commit**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
git add "02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/_supp_panel_D_blood_markers_in_DEP.R" \
        "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/supp_D_blood_markers_DEP.csv" \
        "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/supp/png/panels/supp_D.png"
git commit -m "feat(F02): supp D — blood markers in Cancer_vs_Healthy DEPs"
```

---

## Task 10: Main driver (`01_main_panels.R`)

**Files:**
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/01_main_panels.R`

- [ ] **Step 1: Write the driver**

```r
# F02 main-panels driver. Sources each _panel_*.R in order; each script
# leaves its ggplot object in .GlobalEnv as panel_A/B/C/D.

setwd(rprojroot::find_rstudio_root_file())
base_dir <- "02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script"

source(file.path(base_dir, "_panel_A_blood_stacked.R"))
source(file.path(base_dir, "_panel_B_bm_ratio_dotplot.R"))
source(file.path(base_dir, "_panel_C_bm_vs_mahalanobis.R"))
source(file.path(base_dir, "_panel_D_filter_venns.R"))

stopifnot(exists("panel_A"), exists("panel_B"),
          exists("panel_C"), exists("panel_D"))
message("F02 main panels: all 4 ggplot objects built.")
```

- [ ] **Step 2: Run the driver end-to-end**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis
Rscript "A_CvH_2026/02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/01_main_panels.R"
```

Expected output: `F02 main panels: all 4 ggplot objects built.`

- [ ] **Step 3: Commit**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
git add "02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/01_main_panels.R"
git commit -m "feat(F02): main-panels driver"
```

---

## Task 11: Supp driver (`02_supp_panels.R`)

**Files:**
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/02_supp_panels.R`

- [ ] **Step 1: Write the driver**

```r
# F02 supp-panels driver. Sources each _supp_panel_*.R; objects:
# supp_A / supp_B / supp_C / supp_D.

setwd(rprojroot::find_rstudio_root_file())
base_dir <- "02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script"

source(file.path(base_dir, "_supp_panel_A_marker_intensity.R"))
source(file.path(base_dir, "_supp_panel_B_cutoff_sensitivity.R"))
source(file.path(base_dir, "_supp_panel_C_corr_heatmap.R"))
source(file.path(base_dir, "_supp_panel_D_blood_markers_in_DEP.R"))

stopifnot(exists("supp_A"), exists("supp_B"),
          exists("supp_C"), exists("supp_D"))
message("F02 supp panels: all 4 ggplot objects built.")
```

- [ ] **Step 2: Run the driver**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis
Rscript "A_CvH_2026/02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/02_supp_panels.R"
```

Expected: `F02 supp panels: all 4 ggplot objects built.`

- [ ] **Step 3: Commit**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
git add "02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/02_supp_panels.R"
git commit -m "feat(F02): supp-panels driver"
```

---

## Task 12: Composite stitcher (`90_stitch_F02.R`)

**Files:**
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/90_stitch_F02.R`
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/pdf/MAIN_F00_blood_contamination.pdf`
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/supp/pdf/SUPP_F00_blood_contamination.pdf`

- [ ] **Step 1: Write the stitcher**

```r
# F02 composite stitcher.
# Composes MAIN (2x2: A,B / C,D) and SUPP (2x2: A,B / C,D) into PDFs + PNGs.

suppressPackageStartupMessages({
  library(cowplot); library(ggplot2)
})

setwd(rprojroot::find_rstudio_root_file())
base_dir <- "02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script"

source(file.path(base_dir, "01_main_panels.R"))
source(file.path(base_dir, "02_supp_panels.R"))

main_composite <- plot_grid(
  panel_A, panel_B,
  panel_C, panel_D,
  ncol = 2, labels = c("A", "B", "C", "D"),
  label_size = 14, align = "hv"
)

supp_composite <- plot_grid(
  supp_A, supp_B,
  supp_C, supp_D,
  ncol = 2, labels = c("A", "B", "C", "D"),
  label_size = 14, align = "hv"
)

out_main_pdf <- "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/pdf/MAIN_F00_blood_contamination.pdf"
out_main_png <- "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/png/MAIN_F00_blood_contamination.png"
out_supp_pdf <- "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/supp/pdf/SUPP_F00_blood_contamination.pdf"
out_supp_png <- "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/supp/png/SUPP_F00_blood_contamination.png"

ggsave(out_main_pdf, main_composite, width = 12, height = 9, units = "in")
ggsave(out_main_png, main_composite, width = 12, height = 9, units = "in", dpi = 300, bg = "white")
ggsave(out_supp_pdf, supp_composite, width = 12, height = 9, units = "in")
ggsave(out_supp_png, supp_composite, width = 12, height = 9, units = "in", dpi = 300, bg = "white")

message("F02 composites written:\n  ", out_main_pdf, "\n  ", out_supp_pdf)
```

- [ ] **Step 2: Run the stitcher**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis
Rscript "A_CvH_2026/02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/90_stitch_F02.R"
```

Expected: composite PDFs and PNGs written.

- [ ] **Step 3: Verify composites**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
for f in \
  "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/pdf/MAIN_F00_blood_contamination.pdf" \
  "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/png/MAIN_F00_blood_contamination.png" \
  "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/supp/pdf/SUPP_F00_blood_contamination.pdf" \
  "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/supp/png/SUPP_F00_blood_contamination.png"; do
  test -f "$f" && echo "OK $f" || echo "MISSING $f"
done
```

Expected: 4 `OK` lines.

- [ ] **Step 4: Commit**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
git add "02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/90_stitch_F02.R" \
        "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/main/" \
        "02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/supp/"
git commit -m "feat(F02): composite stitcher (MAIN + SUPP)"
```

---

## Task 13: Data dictionary (`DATA_DICTIONARY.md`)

**Files:**
- Create: `02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/DATA_DICTIONARY.md`

- [ ] **Step 1: Write the data dictionary**

```markdown
# F02 c_data dictionary

One line per CSV. Producer script in parens. Columns + types in inline list.

## Main

- **panel_A_stacked_data.csv** (`_panel_A_blood_stacked.R`)
  `sample_id <chr>` | `Group_Time <fct>` | `marker <chr>` | `pct <dbl>` | `B_M_ratio <dbl>` — 175 rows (35 samples × 5 markers).
- **panel_B_dotplot_data.csv** (`_panel_B_bm_ratio_dotplot.R`)
  `sample_id <chr>` | `Group_Time <fct>` | `B_M_ratio <dbl>` | `top_decile_flag <lgl>` — 35 rows.
- **panel_C_join_data.csv** (`_panel_C_bm_vs_mahalanobis.R`)
  `sample_id <chr>` | `Group_Time <chr>` | `B_M_ratio <dbl>` | `mahal_dist <dbl>` | `n_flags <int>` | `consensus_outlier <lgl>` | `label_this <lgl>` — ≥30 rows (joined Sam ∩ ours).
- **panel_D_venn_protein_sets.csv** (`_panel_D_filter_venns.R`)
  `uniprot_id <chr>` | `gene <chr>` | `in_sam <lgl>` | `in_ours <lgl>` | `set_membership <chr>` — union of Sam-kept (1944) and our-kept (~2582) UniProt IDs.
- **panel_D_venn_blood_markers.csv** (`_panel_D_filter_venns.R`)
  `uniprot_id <chr>` | `gene <chr>` | `in_sam <lgl>` | `in_ours <lgl>` | `in_blood_blacklist <lgl>` — 5 rows (HBB, HBA1, MB, ALB, CKM).

## Supp

- **supp_A_intensity_pairs.csv** (`_supp_panel_A_marker_intensity.R`)
  `uniprot_id <chr>` | `gene <chr>` | `sample_id <chr>` | `intensity_sam <dbl>` | `intensity_ours <dbl>` — per-protein per-sample paired intensities for the 5 markers.
- **supp_B_cutoff_sensitivity.csv** (`_supp_panel_B_cutoff_sensitivity.R`)
  `cutoff_pct <int>` | `n_samples_dropped <int>` | `n_proteins_used <int>` | `n_DEP_fdr10 <int>` | `dropped_sample_ids <chr>` — 5 rows (cutoffs 0/5/10/15/20%).
- **supp_C_corr_matrix.csv** (`_supp_panel_C_corr_heatmap.R`)
  `var1 <chr>` | `var2 <chr>` | `rho <dbl>` | `p <dbl>` — 36 rows (long form of 6×6 Spearman matrix).
- **supp_D_blood_markers_DEP.csv** (`_supp_panel_D_blood_markers_in_DEP.R`)
  `gene <chr>` | `uniprot_id <chr>` | `logFC <dbl>` | `P.Value <dbl>` | `FDR <dbl>` | `pi_score <dbl>` | `sig_FDR_10 <lgl>` | `sig_pi_05 <lgl>` — 5 rows.

## Conventions

- All `sample_id` values are short-form (`CR6_T1`, not `CR006_T1`), normalized via `_shared_keys.R::normalize_sample_id`.
- FDR threshold: 0.10 (BH).
- π-score threshold: 0.05 (Xiao 2014, `P.Value ^ |logFC|`).
- All correlations are Spearman unless stated otherwise.
```

- [ ] **Step 2: Commit**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
git add "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/DATA_DICTIONARY.md"
git commit -m "docs(F02): data dictionary for c_data artifacts"
```

---

## Task 14: End-to-end acceptance verification

**Files:** read-only.

- [ ] **Step 1: Run all panel scripts from scratch**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis
Rscript "A_CvH_2026/02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script/90_stitch_F02.R"
```

Expected: no errors; composites overwrite cleanly.

- [ ] **Step 2: Verify all acceptance criteria from SPEC.md**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
ROOT="02-03_Sam's_Results/04_Figures/F00_blood_contamination"

echo "--- Panel PNGs ---"
for p in panel_A panel_B panel_C panel_D; do
  test -f "$ROOT/b_reports/main/png/panels/$p.png" && echo "OK $p" || echo "MISS $p"
done
for s in supp_A supp_B supp_C supp_D; do
  test -f "$ROOT/b_reports/supp/png/panels/$s.png" && echo "OK $s" || echo "MISS $s"
done

echo "--- c_data CSVs (expect 9) ---"
ls "$ROOT/c_data/"*.csv | wc -l

echo "--- DATA_DICTIONARY ---"
test -f "$ROOT/c_data/DATA_DICTIONARY.md" && echo OK || echo MISS

echo "--- Composites ---"
test -f "$ROOT/b_reports/main/pdf/MAIN_F00_blood_contamination.pdf" && echo "OK main pdf"
test -f "$ROOT/b_reports/supp/pdf/SUPP_F00_blood_contamination.pdf" && echo "OK supp pdf"

echo "--- Panel C join n ---"
Rscript -e 'd <- read.csv("'"$ROOT"'/c_data/panel_C_join_data.csv"); cat("joined n:", nrow(d), "\n"); stopifnot(nrow(d) >= 30)'

echo "--- Supp B sweep n ---"
Rscript -e 'd <- read.csv("'"$ROOT"'/c_data/supp_B_cutoff_sensitivity.csv"); print(d); stopifnot(nrow(d) == 5)'

echo "--- Normalize-ID round-trip ---"
Rscript -e 'source("'"$ROOT"'/a_script/_shared_keys.R")'
```

Expected: 8 `OK` panel lines, `9` c_data CSVs, `OK` for DATA_DICTIONARY and both composites, `joined n` ≥ 30, supp B with 5 rows, `_shared_keys.R: self-check passed`.

- [ ] **Step 3: Tag completion**

If all acceptance checks pass, no further commit needed — the work is already in the history. Otherwise iterate on the failing task.

---

## Cross-task notes

### Risk mitigations from SPEC.md

| Risk | Plan-time handling |
|---|---|
| Supp B refit needs design rebuild | Task 7 Step 1 explicitly inspects `sam$design` before writing. If contrast names differ, adjust `add_contrasts()` per the diagnostic output. |
| Panel D(ii) Venn with 5 points | Plan uses an indicator heatmap (`geom_tile + ✓/✗`) instead of a Venn — visually cleaner for n=5. |
| Panel C join sparsity | Task 4 Step 3 halts if join < 30 samples. |
| `consensus_outlier == TRUE` samples missing from Sam's cohort | Panel C is a join; missing samples simply don't appear. Caption notes joined `n` to make this explicit. |

### Package dependencies

The plan assumes installed: `dplyr`, `tidyr`, `readr`, `ggplot2`, `cowplot`, `ggrepel`, `RColorBrewer`, `eulerr`, `proteoDA`, `rprojroot`, `tibble`, `scales`. If any are missing the panel scripts will fail at `library()` — install one-off via `install.packages("<pkg>")`.

### Git commit cadence

One commit per task (panel + data + PNG + commit). Tasks 10–13 commit individually. Task 14 is verification-only with no commit. This produces ~13 commits in sequence — keep the cadence; do not batch.
