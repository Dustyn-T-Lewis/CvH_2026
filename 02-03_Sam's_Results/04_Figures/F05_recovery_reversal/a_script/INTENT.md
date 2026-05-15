# F05 Recovery Reversal — INTENT

## Biological Question

Does resistance training reverse the proteomic signature of cancer (vs healthy)?

We cross two contrasts from Sam Cooper's CvH muscle proteomics:
- **Cancer_vs_Healthy** (SURV − CTL): proteins that differ between cancer survivors and healthy controls
- **Training_CR** (Pre → Post combined across CRE + PLA arms): proteins that change with training in cancer survivors

Concordance interpretation:
- **Reversal**: opposite direction (Cancer-Up + Training-Down, or Cancer-Down + Training-Up)
- **Non-reversal** (Exacerbated): same direction (both up or both down)

## Panel Plan

### Main (5 panels, inline in 01_main_panels.R)

| Panel | Name                       | Content                                               |
|-------|----------------------------|-------------------------------------------------------|
| A     | Quadrant ORA               | Scatter (Cancer logFC vs Training logFC) + flanking ORA bars per quadrant |
| B     | Pattern Heatmap            | Per-protein classification heatmap with GO Slim categories |
| C     | fry Barcode                | Rotation test: do Cancer-Up DEPs train Down and vice versa? |
| D     | NES Scatter                | Pathway-level reversal: fGSEA NES (Cancer) vs NES (Training_CR) |
| E     | RRHO2                      | Threshold-free rank-rank hypergeometric overlap       |

### Supplementary (6 panels, inline in 02_supp_panels.R)

| Panel | Name                            | Defends       |
|-------|---------------------------------|---------------|
| A     | ORA Dedup Sensitivity           | Main A        |
| B     | Pearson r Bootstrap             | Main A + D    |
| C     | Circularity Diagnostic          | Main C (fry)  |
| D     | Reversal Threshold Sensitivity  | Main B        |
| E     | GO Slim Category Distribution   | Main B        |
| F     | fry Leading-Edge Proteins       | Main C        |

## Contrast Mapping (YvO → Sam)

| YvO contrast    | Sam contrast      | Axis |
|-----------------|-------------------|------|
| Aging           | Cancer_vs_Healthy | X    |
| Training_Old    | Training_CR       | Y    |

## Key Files

- DEP CSVs: `02-03_Sam's_Results/03_DEP/c_data/04_per_contrast_results/{Cancer_vs_Healthy,Training_CR}.csv`
- fGSEA cache: `02-03_Sam's_Results/04_Figures/shared/fgsea_cache/{Cancer_vs_Healthy,Training_CR}.rds`
- Protein columns: `logFC`, `t`, `pi_score`, `adj.P.Val`
