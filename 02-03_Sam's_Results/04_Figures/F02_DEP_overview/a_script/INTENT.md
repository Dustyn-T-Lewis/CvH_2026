# F02_DEP_overview — Sam-parallel DEP overview (6 main + 3 supp panels)

Mirrors YvO 2026 F02 applied to Sam's CRvH cohort (N=35, 1944 proteins).
All panels use the canonical pipeline rerun of DEP on Sam's normalized data
(`03_DEP/c_data/`).

## Main panels (3×2 composite, 178 × 115 mm)

| Panel | Content | Key inputs | Plot type |
|-------|---------|------------|-----------|
| A | Sample PCA + PERMANOVA | `00_input/01_normalized_DAList_SURV_stringent_muscle.RDS` | PCA biplot, 5-group palette |
| B | logFC density histograms | combined_results_CRvH + combined_results_CR | Histogram + density, 4 contrasts |
| C | DEPs per contrast (stacked bar, 3 thresholds) | per-contrast CSVs (pi_score) | Horizontal stacked bar |
| D | UpSet contrast overlap | per-contrast CSVs (sig_pi) | Custom ggplot2 UpSet (dual bar + dot matrix) |
| E | fGSEA pathway enrichment | `04_Figures/shared/fgsea_cache/*.rds` | Stacked dodged bar (Up/Down × DB) |
| F | DEP rank barcode | combined results + per-contrast pi_score | Barcode + density traces per contrast |

## Supplementary panels (composite, 178 × 115 mm)

| Panel | Content | Inputs | Plot type |
|-------|---------|--------|-----------|
| SA | Per-protein CV% scatter (group pairs) | normalized DAList | Scatter faceted by group pair |
| SB | Per-group CV% violin | normalized DAList | Violin + boxplot |
| SC | Intra-individual variability (CRvH only) | normalized DAList + per-contrast logFC | Boxplot per subject |

## Contrast selection rationale (4 for UpSet/Barcode)

6 contrasts total. Two biological axes:
- **Cancer axis**: `Cancer_vs_Healthy`, `Training_CR` (primary disease recovery story)
- **Supplement axis**: `Training_CRE`, `Training_PLA`, `Baseline_Supplement`, `Supplement_Interaction`

For UpSet (panel D), we use 4: `Cancer_vs_Healthy`, `Training_CR`, `Training_CRE`, `Training_PLA`.
This captures both the disease-recovery axis and the supplement-stratified training effects.
`Baseline_Supplement` and `Supplement_Interaction` are secondary; they can appear in barcode (F)
and fGSEA (E) using all 6 contrasts.

## Run order

```
Rscript 02-03_Sam's_Results/04_Figures/F02_DEP_overview/a_script/90_stitch_F02.R
```

or individually:
```
Rscript 02-03_Sam's_Results/04_Figures/F02_DEP_overview/a_script/02_supp_panels.R
Rscript 02-03_Sam's_Results/04_Figures/F02_DEP_overview/a_script/01_main_panels.R
```
