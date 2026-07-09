# CvH Proteomics Analysis

Label-free DIA-MS skeletal-muscle proteomics comparing cancer-recovery (`CR`)
against healthy (`H`) participants. The CR arm carries pre- and post-training
timepoints; the CRE and PLA training sub-arms are pooled into a three-cell design:

- `H_pre`: healthy baseline
- `CR_pre`: cancer-recovery baseline
- `CR_post`: cancer-recovery after training

## Design and Contrasts

DEP fits the means model `~ 0 + group_time + (1 | Subject_ID)` over the three
cells, with the within-subject pre-to-post pairing carried by
`duplicateCorrelation` on `Subject_ID`. Three contrasts, each a linear
combination of the three cell means:

- `CRvH_Baseline = CR_pre - H_pre` — disease deviation from healthy (D)
- `CR_Training = CR_post - CR_pre` — effect of training (T)
- `Resid = CR_post - H_pre` — what remains after training (R = D + T)

The three form a closed triangle: `Resid = CRvH_Baseline + CR_Training`.
Significance is the Pi-score (`Pi = P.Value^|log2FC|`, Xiao et al. 2014) at
`Pi < 0.05`, with Benjamini-Hochberg FDR as a secondary criterion.

## Pipeline Overview

| Stage | Directory | Canonical logic |
| --- | --- | --- |
| `00` | `00_input/` | Raw DIA-NN intensity matrix, sample sheet, HPA annotations |
| `01` | `01_Filtering/` | HPA blood-contaminant removal with a myonuclei rescue, missingness filter, consensus outlier detection -> `DAList_filtered.rds` |
| `02` | `02_Normalization/` | `cycloess` normalization; `imputation/` holds three arms (`imp4p`, MsCoreUtils hybrid, `missForest`), each writing a method-tagged `DAList_imputed_<method>.rds`. `imp4p` is the canonical imputed arm |
| `03` | `03_DEP/` | `a_non_imputed/`: primary `limma + duplicateCorrelation`, three contrasts, Pi-score. `b_imputed/`: concordance DEP on the imputed matrices (`imp4p` canonical; `missForest` and MsCoreUtils as comparison) |
| `04` | `04_Figures/` | F01-F04 and F06 figure trees, plus a `Reversal/` signature-reversal analysis |

## Canonical Run Order

```sh
Rscript 01_Filtering/a_script/01_run_filtering.R
Rscript 02_Normalization/a_script/01_run_normalization.R
Rscript 02_Normalization/imputation/a_script/a_imp4p.R          # canonical imputed arm
Rscript 02_Normalization/imputation/a_script/b_mscoreutils.R    # comparison
Rscript 02_Normalization/imputation/a_script/c_missforest.R     # comparison

Rscript 03_DEP/a_non_imputed/a_script/01_run_dep.R             # primary DEP
Rscript 03_DEP/b_imputed/a_script/01_run_dep_imputed.R         # imputed concordance
```

Figures under `04_Figures/` read the DEP outputs; rerun them after stages 01-03
complete. `CvH_pipeline.qmd` is the frozen end-to-end walkthrough.

## Repository Conventions

- `a_script/`: scripts and narrative notebooks
- `b_reports/`: generated figure renders and QC reports
- `c_data/`: stage outputs read by downstream steps
- primary DEP uses the non-imputed normalized matrix; imputation feeds only the
  concordance check and figures

## Reproducibility Rules

- Path resolution uses `here::here()` from the project root (figure scripts
  anchor the working directory with `setwd(here::here())`)
- stochastic steps use `set.seed(42)`
- primary DEP uses the non-imputed matrix; `imp4p` is the canonical imputed arm
- repeated-measures blocking uses `Subject_ID`
