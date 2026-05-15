# F03 — Volcano Rings (Sam parallel)

YvO F03 transposed to Sam's CvH cohort. Four polar-coordinate volcano-ring
panels, one per headline contrast, arranged in a 2×2 composite.

## Main panels (A–D)

| Tag | Contrast             | Biological question                          |
|-----|----------------------|----------------------------------------------|
| A   | Cancer_vs_Healthy    | Disease proteome shift (136 sig pathways)    |
| B   | Training_CR          | Combined training response (36 sig pathways) |
| C   | Training_CRE         | Creatine arm training (44 sig pathways)      |
| D   | Training_PLA         | Placebo arm training (1 sig pathway)         |

## Supplementary panels (A–E)

| Tag | Panel                      | Content                                      |
|-----|----------------------------|----------------------------------------------|
| A   | Raw p-value histograms     | P.Value per contrast (4-facet)               |
| B   | Pi-score distributions     | pi_score per contrast with Π < 0.05 threshold|
| C   | FDR distributions          | adj.P.Val per contrast                       |
| D   | MA plots                   | Mean intensity vs logFC per contrast         |
| E   | Outlier sensitivity        | DEP retention with/without outliers          |

## Notes

- pi-score threshold: < 0.05 (primary significance criterion)
- Ring terms: top 12 from Hallmark + GO Slim, FDR < 0.05
- Training_PLA (panel D) has only 1 sig pathway — ring will show 1 arc
- Supp E reads `03_DEP/c_data/11_outlier_sensitivity.csv` (main CvH pipeline)
- fGSEA inputs: `04_Figures/shared/fgsea_cache/*.rds` (6 per-contrast RDS files)

## Outputs

- `b_reports/main/pdf/MAIN_F03_composite.pdf`  (178 × 180 mm)
- `b_reports/main/png/MAIN_F03_composite.png`
- `b_reports/supp/pdf/SUPP_F03_composite.pdf`  (178 × 225 mm)
- `b_reports/supp/png/SUPP_F03_composite.png`
- `c_data/F03_supplementary.xlsx`
