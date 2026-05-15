# F04 Training Concordance — Intent

## Question
Does resistance training produce the same proteomic signature in creatine-supplemented (CRE)
vs placebo (PLA) participants? I.e., is the training response supplement-independent?

## Contrast Mapping (YvO analogue)
| YvO F04 | Sam CvH F04 |
|---------|-------------|
| Training_Young | Training_CRE |
| Training_Old   | Training_PLA |

## Concordance Interpretation
- **Same direction (Concordant Up / Down)** = supplement-independent training effect
- **Opposite direction (Discordant)** = training effect is supplement-modified

## Pre-Warning: Sparse DEPs
Both contrasts have 0 DEPs at Pi < 0.05 (Training_CRE: 0 sig, Training_PLA: 0 sig).
- Panel A scatter will show all proteins NS; quadrant ORA will be gene-level (unsorted by pi)
- Panel C fry will use empty gene sets — barcode shows full ranked distribution
- Panel B heatmap will include proteins at relaxed pi < 0.10 or padj < 0.10 thresholds
- **Pathway-level (Panel D: fGSEA NES) is the primary signal carrier**
  - Training_CRE: 44 sig pathways; Training_PLA: 1 sig pathway

## Panels
### Main (5 panels)
- **A**: Quadrant ORA scatter (logFC CRE vs logFC PLA) + flanking ORA bars
- **B**: Pattern heatmap (per-protein concordance classification; relaxed threshold)
- **C**: fry rotation test (CRE-direction proteins → PLA t-stats)
- **D**: fGSEA NES scatter (pathway-level concordance) — **primary signal**
- **E**: RRHO2 threshold-free rank-rank overlap

### Supplementary (6 panels)
- **A**: ORA dedup sensitivity (Jaccard cutoffs; concordant quadrant genes)
- **B**: Pearson r bootstrap (logFC concordance CI)
- **C**: Circularity diagnostic (permuted null for r)
- **D**: Threshold sensitivity (concordant/discordant % by logFC cutoff)
- **E**: GO Slim category distribution by concordance quadrant
- **F**: Leading-edge proteins (top by |t_PLA|)

## Data Sources
- DEP: `02-03_Sam's_Results/03_DEP/c_data/04_per_contrast_results/{Training_CRE,Training_PLA}.csv`
- fGSEA: `02-03_Sam's_Results/04_Figures/shared/fgsea_cache/{Training_CRE,Training_PLA}_fgsea.rds`
