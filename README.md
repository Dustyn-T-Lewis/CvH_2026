# CvH Proteomics Analysis

Label-free DIA-MS skeletal-muscle proteomics comparing cancer-recovery (`CR`)
against healthy (`H`) participants. The CR arm carries pre- and post-training
timepoints; healthy controls give a single baseline biopsy.

39 samples survive QC: 10 healthy at T1, 15 CR at T1, 14 CR at T2.

## Design and Contrasts

DEP fits the means model `~ 0 + model_cell + (1 | Subject_ID)` over five cells
(`H_pre` plus `CRE`/`PLA` × `pre`/`post`), with the within-subject pre-to-post
pairing carried by `duplicateCorrelation` on `Subject_ID` (consensus rho = 0.248).
The CRE and PLA arms enter the model so estimates are supplement-adjusted.

Three reported contrasts average the two supplement arms 50:50:

- `CRvH_Baseline = ½(CRE_pre + PLA_pre) - H_pre` — disease deviation from healthy (D)
- `CR_Training = ½(CRE_post + PLA_post) - ½(CRE_pre + PLA_pre)` — effect of training (T)
- `Resid = ½(CRE_post + PLA_post) - H_pre` — what remains after training (R = D + T)

These form a closed triangle: `Resid = CRvH_Baseline + CR_Training`.

Four further contrasts test the supplement arms directly and are written to the
same results table, so `combined_results_pi.csv` carries seven contrasts, not three:
`Baseline_Supplement`, `Training_CRE`, `Training_PLA`, `Supplement_Interaction`.

Significance is the Pi-score (`Pi = P.Value^|log2FC|`, Xiao et al. 2014) at
`Pi < 0.05`, with Benjamini-Hochberg FDR as a secondary criterion. Proteins whose
contrast is not estimable carry `sig_pi = NA`, not `0`.

## Pipeline Overview

| Stage | Directory | Canonical logic |
| --- | --- | --- |
| `00` | `00_input/` | DIA-NN intensity matrix, sample sheet, HPA annotations, RBC proteome reference |
| `01` | `01_Filtering/` | HPA presence, blood-contaminant removal with myonuclei rescue, red-cell tracking removal, missingness filter, consensus outlier detection -> `DAList_filtered.rds` |
| `02` | `02_Normalization/` | `cycloess` normalization; `imputation/` holds three arms (`imp4p`, MsCoreUtils hybrid, `missForest`), each writing `DAList_imputed_<method>.rds`. **missForest is the arm the figures read** |
| `03` | `03_DEP/` | `a_non_imputed/`: primary `limma + duplicateCorrelation`, seven contrasts, Pi-score. `b_imputed/`: concordance DEP on all three imputed matrices |
| `04` | `04_Figures/` | `F01`–`F06` figure trees plus `shared/` |

Filter cascade, from `01_Filtering/c_data/filtering_report.xlsx`:

| Step | Proteins remaining | Removed |
| --- | --- | --- |
| Raw input | 3172 | — |
| HPA presence | 2962 | 210 |
| Blood contaminant removal | 2811 | 151 |
| Red-cell tracking removal | 2437 | 374 |
| Missingness (>=5 in >=1 group) | 2176 | 261 |

## Figures

| Figure | Contents |
| --- | --- |
| `F01_Phenotype` | Phenotype and strength outcomes; supplement covers LBM, chest press, leg extension, grip |
| `F02_Proteome_Overview` | PCA, DEP counts, effect sizes, overlap, direction, pathways; QC supplement (CV, ICC, dbRDA) |
| `F03_Enrich_Volcanoes` | `enrichVolcano` ring grid over the Model-1 contrast trio |
| `F04_Reversal` | Reversal landscape, trajectory clustering, fry rotation; diagnostics and methods supplements |
| `F05_WGCNA` | Module card: counts, member response (fGSEA NES over fry), eigengene trajectory, ORA |
| `F06_Prediction` | Three feature spaces x two outcomes, with a circularity ladder |

WGCNA network: bicor, signed, soft power 14, five modules — turquoise 675,
blue 447, brown 408, yellow 182, green 59, plus 405 unassigned (grey).

## Canonical Run Order

```sh
Rscript 01_Filtering/a_script/01_run_filtering.R
Rscript 02_Normalization/a_script/01_run_normalization.R
Rscript 02_Normalization/imputation/a_script/a_imp4p.R
Rscript 02_Normalization/imputation/a_script/b_mscoreutils.R
Rscript 02_Normalization/imputation/a_script/c_missforest.R      # figures read this arm

Rscript 03_DEP/a_non_imputed/a_script/01_run_dep.R               # primary DEP
Rscript 03_DEP/b_imputed/a_script/01_run_dep_imputed.R           # imputed concordance

Rscript 04_Figures/shared/build_fgsea_cache.R                    # pathway cache
Rscript 04_Figures/F01_Phenotype/a_script/90_stitch_F01.R
Rscript 04_Figures/F01_Phenotype/a_script/supp/90_stitch_F01_supp.R
Rscript 04_Figures/F02_Proteome_Overview/a_script/90_stitch_F02.R
Rscript 04_Figures/F02_Proteome_Overview/a_script/supp/90_stitch_F02_supp.R
Rscript 04_Figures/F03_Enrich_Volcanoes/a_script/01_enrich_volcanoes.R
Rscript 04_Figures/F04_Reversal/a_script/90_stitch_F04.R
Rscript 04_Figures/F05_WGCNA/a_script/00_run_F05.R
Rscript 04_Figures/F06_Prediction/a_script/00_run_F06.R          # needs F04 and F05 first
```

`04_Figures/F05_WGCNA/a_script/supp/network_validation.R` is a long-running
parameter sweep, run on demand rather than from the driver.

## Repository Conventions

- `a_script/`: scripts
- `b_reports/`: generated figure renders and QC reports
- `c_data/`: stage outputs read by downstream steps
- stage outputs are committed as a full mirror, so a fresh clone can read any
  stage without re-running the one above it

## Reproducibility Rules

- paths resolve from the project root; figure scripts anchor with `setwd(here::here())`
- stochastic steps use `set.seed(42)`
- primary DEP uses the non-imputed matrix; imputation feeds the concordance check,
  PCA, and WGCNA only
- repeated-measures blocking uses `Subject_ID`
- `shared/wgcna_stats.R` and `shared/prediction_utils.R` call WGCNA namespace-qualified
  and never attach it: `WGCNA::cor()` returns a 1x1 matrix and silently breaks any bare
  `cor()` in a script sourced afterwards

## Known Environment Issue

PDF renders drop `Π`, `ρ`, and `✱` on machines without XQuartz: `cairo_pdf` cannot
load, `get_pdf_device()` falls back to base `pdf()`, and that device cannot encode
those glyphs. PNG output is unaffected. Installing XQuartz restores them.
