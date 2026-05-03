# CvH F1 Figure Design — Proteomics Overview

## Context

Adapt YvO F1 (8 panels + 1 supplementary) for CvH dataset.

**CvH design:** CR_CRE (n=7, T1/T2), CR_PLA (n=8, T1/T2), PPS (n=10, T1 only).
2582 proteins, 39 samples. Two limma models: CRvH (all samples) and CR (cancer only).

**5 contrasts:** Cancer_vs_Healthy, Training_CR, Training_CRE, Training_PLA, Supplement_Interaction.

## Directory Structure

```
A_CvH_2026/04_Figures/
├── shared/style.R, pathway_utils.R
└── F1/a_script/ (style.R + panels A-H), b_reports/, c_data/
```

## Color Palette

- CRE_T1=#2166AC, CRE_T2=#67A9CF, PLA_T1=#D6604D, PLA_T2=#F4A582, H_T1=#4DAF4A
- Contrasts: Cancer_vs_Healthy=#4CAF50, Training_CR=#9C27B0, Training_CRE=#2166AC, Training_PLA=#D6604D, Supplement_Interaction=#FF8F00

## Panel Specifications

| Panel | Content | Groups/Contrasts | Size (mm) | Notes |
|-------|---------|-----------------|-----------|-------|
| A | CV% violins | CRE, PLA, Healthy facets | 160×120 | PPS=T1 only baseline |
| B | CV scatter triptych | CRE, PLA (paired) + ΔCV | 300×120 | PPS excluded |
| C | Intra-individual variability | CRE vs PLA subjects | 160×90 | Wilcoxon CRE vs PLA |
| D | logFC density | 5 contrasts | 140×160 | Median |logFC| + CIs |
| E | PCA biplot | 5 groups + PERMANOVA | 145×100 | Group+Time+Supp terms |
| F | DEP counts (pseudo-log bars) | 5 contrasts | 200×70 | p/FDR(0.10)/Pi tiers |
| G | UpSet overlap | 5 contrast sig sets | 200×120 | Pi<0.05 threshold |
| H | fGSEA grouped bars | 5 contrasts × 4 DBs | 160×variable | Hallmark/KEGG/Reactome/GO:BP |

## Key Adaptations from YvO

1. 5 contrasts (was 3)
2. PPS as T1-only baseline reference in panels A, E
3. CRE/PLA replace Young/Old for paired analyses (panels A-C)
4. FDR 0.10 exploratory threshold (was 0.05)
5. Wider panels for 5 contrasts (F, G)

## Data Sources

- Normalized: `01_normalization/c_data/02_normalized.csv`
- Imputed: `02_Imputation/c_data/01_imputed.csv`
- DEP CRvH: `03_DEP/c_data/03_combined_results_CRvH.csv`
- DEP CR: `03_DEP/c_data/03_combined_results_CR.csv`
- Per-contrast: `03_DEP/c_data/04_per_contrast_results/*.csv`
- Metadata: `00_input/CvH_meta.csv`
