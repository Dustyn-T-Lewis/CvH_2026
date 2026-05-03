# CvH 2026 -- Cancer Recovery Skeletal Muscle Proteomics

This repository contains the CvH proteomics analysis adapted from the validated
`A_YvO_2026` workflow, but rederived for the actual CvH experimental design
rather than copied from YvO assumptions.

The current pipeline is organized as:

`00_input -> 01_normalization -> 02_Imputation -> 03_DEP -> 05_WGCNA -> 04_Figures`

`04_Figures` consumes outputs from the upstream analysis stages; `05_WGCNA`
feeds the module-level figure streams.

## Study design

### Factors present in CvH

| Factor | Levels |
| --- | --- |
| Cancer recovery group | `CR_CRE`, `CR_PLA` |
| Healthy reference | `PPS` / `H_T1` |
| Timepoint | `T1`, `T2` for cancer-recovery subjects only |
| Supplement | `CRE`, `PLA` for cancer-recovery subjects; healthy controls baseline-only |

### Design consequences

- Cancer-recovery subjects support repeated-measures analyses.
- Healthy controls are baseline-only and support cross-sectional baseline comparisons.
- There is no healthy `T2`, so YvO's symmetric repeated-measures 2 x 2 design does not transfer directly.
- CvH therefore uses two valid limma/LMM model families:
  - `CRvH`: all subjects, for baseline cancer-vs-healthy and pooled survivor training response
  - `CR`: cancer-recovery subjects only, for creatine-vs-placebo baseline and training contrasts

### Planned contrasts

#### CRvH model

- `Cancer_vs_Healthy = (CRE_T1 + PLA_T1)/2 - H_T1`
- `Training_CR = (CRE_T2 + PLA_T2)/2 - (CRE_T1 + PLA_T1)/2`

#### CR model

- `Baseline_Supplement = CRE_T1 - PLA_T1`
- `Training_CRE = CRE_T2 - CRE_T1`
- `Training_PLA = PLA_T2 - PLA_T1`
- `Supplement_Interaction = (CRE_T2 - CRE_T1) - (PLA_T2 - PLA_T1)`

Primary significance framing remains YvO-style exploratory proteomics:

- nominal `p <= 0.10`
- Benjamini-Hochberg adjusted `p`
- Pi-score threshold `Pi < 0.05`

## Inputs

| File | Role |
| --- | --- |
| `00_input/CvH_raw.xlsx` | DIA-MS protein matrix with annotation |
| `00_input/CvH_meta.csv` | Analysis metadata aligned to matrix column names |
| `00_input/CRm_meta.csv` | Source phenotype/clinical metadata used for cross-checking and phenotype figures |
| `00_input/HPA_skeletal_muscle_annotations.tsv` | Skeletal muscle tissue reference |

## Shared validation

`R/cvh_design.R` is the shared entrypoint for:

- sample-ID harmonization (`CR006_T1 -> CR6_T1`)
- consistency checks between `CvH_meta.csv` and `CRm_meta.csv`
- repeated-measures design validation for cancer-recovery subjects
- enforcement of healthy baseline-only structure

Core analysis stages should source this helper instead of carrying local,
inconsistent metadata assumptions.

## Stage 01 -- Normalization

Script:

- `01_normalization/a_script/01_run_normalization.R`

Main logic:

- HPA skeletal-muscle filter
- blood/immunoglobulin contaminant removal
- UniProt deduplication
- missingness filtering by `Group_Time`
- 4-method outlier consensus
- cycloess normalization

Key outputs:

- `01_normalization/c_data/02_normalized.csv`
- `01_normalization/c_data/03_DAList_normalized.rds`
- `01_normalization/c_data/05_normalization_supp.xlsx`
- `01_normalization/b_reports/01_norm_comparison.pdf`
- `01_normalization/b_reports/02_qc_pre.pdf`
- `01_normalization/b_reports/03_qc_post.pdf`
- `01_normalization/b_reports/04_diagnostics.pdf`

## Stage 02 -- Imputation

Scripts:

- `02_Imputation/a_script/apply_missforest.R`
- `02_Imputation/a_script/02_imputation_reports.R`

Main logic:

- YvO-style 3-method MAR/MNAR classification
- missForest imputation
- low-confidence imputation flagging
- workbook/report generation from the active script outputs

Key outputs:

- `02_Imputation/c_data/01_imputed.csv`
- `02_Imputation/c_data/01_DAList_imputed.rds`
- `02_Imputation/c_data/02_imputation.xlsx`
- `02_Imputation/c_data/02_mar_mnar_classification.csv`
- `02_Imputation/c_data/07_imputation_mask.csv`
- `02_Imputation/c_data/08_mnar_imputation_audit.csv`
- `02_Imputation/c_data/09_imputation_summary.txt`
- `02_Imputation/b_reports/01_missingness_report.pdf`
- `02_Imputation/b_reports/02_imputation_report.pdf`

## Stage 03 -- Differential abundance

Scripts:

- `03_DEP/a_script/01_run_dep.R`
- `03_DEP/a_script/02_dep_reports.R`
- `03_DEP/a_script/03_dep_robustness.R`
- `03_DEP/a_script/04_dep_overview.R`

Main logic:

- limma + duplicateCorrelation blocking
- one `CRvH` model and one `CR` model
- per-contrast result tables, Pi-scores, reports, robustness summaries

Key outputs:

- `03_DEP/c_data/03_combined_results_CRvH.csv`
- `03_DEP/c_data/03_combined_results_CR.csv`
- `03_DEP/c_data/04_per_contrast_results/*.csv`
- `03_DEP/c_data/05_results_CRvH.xlsx`
- `03_DEP/c_data/05_results_CR.xlsx`
- `03_DEP/c_data/10_DEP_supplementary.xlsx`
- `03_DEP/b_reports/02_dep_overview.pdf`

## Stage 05 -- WGCNA

Script:

- `05_WGCNA/a_script/01_run_wgcna.R`

Main logic:

- signed WGCNA network on the imputed matrix
- module assignments and hub proteins
- module-level enrichment
- LMM contrast testing aligned to the same `CRvH` and `CR` model logic

Key outputs:

- `05_WGCNA/c_data/wgcna/*`
- `05_WGCNA/b_reports/soft_threshold_SUPP.pdf`
- `04_Figures/F06/c_data/*`

## Figures

Figure streams already present in this repository are CvH-specific analogues, not
blind YvO copies.

### Direct or near-direct analogues

- `04_Figures/F01`: phenotype-level summaries
- `04_Figures/F03/CRvH` and `04_Figures/F03/CR`: DEP overview, overlap, rank, enrichment views
- `04_Figures/Reversal`: cancer-vs-healthy baseline signal versus pooled CR training response
- `04_Figures/F08`: WGCNA/module summaries

### Important non-transfers

- YvO age-group language does not transfer.
- YvO healthy post-training comparisons do not transfer.
- YvO classifier-style phenotype prediction does not have a default one-to-one CvH equivalent.

See `docs/yvo_to_cvh_method_mapping.md` for the full transfer matrix and rerun order.

## Reproducibility notes

- All stochastic steps set `set.seed(42)`.
- `01_normalization/c_data/02_normalized.csv` is the canonical handoff to stages 02 and 03 for deterministic float serialization.
- Metadata validation is part of the active pipeline, not a manual pre-step.
- Raw inputs are not overwritten by stage scripts.

## Recommended rerun order

1. `01_normalization/a_script/01_run_normalization.R`
2. `01_normalization/a_script/02_norm_reports.R`
3. `02_Imputation/a_script/apply_missforest.R`
4. `02_Imputation/a_script/02_imputation_reports.R`
5. `03_DEP/a_script/01_run_dep.R`
6. `03_DEP/a_script/02_dep_reports.R`
7. `03_DEP/a_script/03_dep_robustness.R`
8. `03_DEP/a_script/04_dep_overview.R`
9. `05_WGCNA/a_script/01_run_wgcna.R`
10. Figure scripts or stitchers that consume refreshed outputs
