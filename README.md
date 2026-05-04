# CvH 2026 -- Cancer Recovery Skeletal Muscle Proteomics

Adapted from the validated `A_YvO_2026` workflow, rederived for the actual CvH
experimental design rather than copied from YvO assumptions.

## Pipeline layout

```
00_input        raw matrix + metadata + HPA annotations
01_normalization HPA filter, blood removal, dedup, missingness, outliers, cycloess
02_Imputation    MAR/MNAR consensus + missForest (live); 17-method benchmark (manual)
03_DEP           two-model proteoDA/limma + sensitivity arms in xlsx supplement
04_Figures       per-figure F01..F08 panel scripts (live); F00 supp QC TBD
02-03_Sam's_Results  collaborator outputs + our DEP rerun on his data
```

Each stage uses the `a_script/` + `b_reports/` + `c_data/` triple to mirror YvO.

## Study design

### Factors

| Factor | Levels |
| --- | --- |
| Cancer recovery group | `CR_CRE`, `CR_PLA` |
| Healthy reference | `PPS` / `H_T1` |
| Timepoint | `T1`, `T2` for cancer-recovery subjects only |
| Supplement | `CRE`, `PLA` for cancer-recovery subjects; healthy controls baseline-only |

### Two-model rationale

There is no healthy `T2`, so YvO's symmetric repeated-measures 2x2 design does
not transfer. CvH uses two limma/LMM model families:

- **CRvH** (all subjects): baseline cancer-vs-healthy and pooled survivor training response
- **CR** (cancer-recovery only): creatine-vs-placebo baseline, per-arm training, supplement-by-time interaction

Diverges from YvO's single full-cohort model -- intentional and documented.

#### CRvH model contrasts

- `Cancer_vs_Healthy = (CRE_T1 + PLA_T1)/2 - H_T1`
- `Training_CR = (CRE_T2 + PLA_T2)/2 - (CRE_T1 + PLA_T1)/2`

#### CR model contrasts

- `Baseline_Supplement = CRE_T1 - PLA_T1`
- `Training_CRE = CRE_T2 - CRE_T1`
- `Training_PLA = PLA_T2 - PLA_T1`
- `Supplement_Interaction = (CRE_T2 - CRE_T1) - (PLA_T2 - PLA_T1)`

Significance framing (YvO-style exploratory proteomics):

- nominal `p < 0.10`
- Benjamini-Hochberg adjusted `p`
- Pi-score threshold `Pi < 0.05` (Xiao et al. 2014, PMID 24478644)

## Inputs

| File | Role |
| --- | --- |
| `00_input/CvH_raw.xlsx` | DIA-MS protein matrix with annotation |
| `00_input/CvH_meta.csv` | Analysis metadata aligned to matrix column names |
| `00_input/CRm_meta.csv` | Source phenotype/clinical metadata used for cross-checking and phenotype figures |
| `00_input/HPA_skeletal_muscle_annotations.tsv` | Skeletal muscle tissue reference |

## Shared utilities

`R/cvh_design.R` is the shared entrypoint for:

- sample-ID harmonization (`CR006_T1 -> CR6_T1`)
- consistency checks between `CvH_meta.csv` and `CRm_meta.csv`
- repeated-measures design validation
- enforcement of healthy baseline-only structure
- `BLOOD_CONTAMINANTS` plasma-protein gene list (Geyer 2016, PMID 27135364)

Core analysis stages source this helper instead of carrying local copies.

## Stage 01 -- Normalization

Scripts: `01_normalization/a_script/01_run_normalization.R` and `02_norm_reports.R`

Pipeline: HPA skeletal-muscle filter -> blood/Ig contaminant removal ->
UniProt deduplication -> missingness filter by `Group_Time` -> 4-method
outlier consensus -> cycloess (Bolstad 2003).

Key outputs:

- `c_data/02_normalized.csv` (canonical CSV handoff to stages 02 and 03)
- `c_data/03_DAList_normalized.rds` (DAList object for downstream)
- `c_data/05_normalization_supp.xlsx`
- `b_reports/01-04_*.pdf` (4 PDFs: norm comparison, QC pre, QC post, diagnostics)

## Stage 02 -- Imputation

Live scripts: `02_Imputation/a_script/apply_missforest.R` and `02_imputation_reports.R`

Logic: 3-method MAR/MNAR consensus (kmeans + global-logistic + left-tail) ->
missForest on full matrix (Stekhoven 2012, PMID 22039212) -> low-confidence
flagging (>50 percent missing).

Optional 17-method benchmark in `a_script/benchmark/`, aligned with YvO's
registry. Run manually first time:

```
Rscript 02_Imputation/a_script/benchmark/_run_all.R
```

Outputs go to `c_data/benchmark/`. Once present, `apply_missforest.R` and
`02_imputation_reports.R` automatically pick up the composite ranking. Six
CvH-only hybrid methods (BPCA_QRILC, KNN_QRILC, imp4p_mixed, msImpute_v2_mnar,
MAI, RF_MsCoreUtils) remain in `methods/` as a library; re-register them in
`benchmark/_common.R::BASE_METHODS` to include them in the run.

Key outputs:

- `c_data/01_imputed.csv`, `01_DAList_imputed.rds`
- `c_data/02_imputation.xlsx`, `02_mar_mnar_classification.csv`
- `c_data/07_imputation_mask.csv`, `08_mnar_imputation_audit.csv`
- `b_reports/01_missingness_report.pdf`, `02_imputation_report.pdf`

## Stage 03 -- Differential abundance

Scripts: `03_DEP/a_script/01_run_dep.R`, `02_dep_reports.R`, `03_dep_robustness.R`, `04_dep_overview.R`

Logic: limma + duplicateCorrelation blocking on `subject` (Smyth 2005),
`robust=TRUE, trend=TRUE` eBayes, BH multiple testing, Pi-score
(Xiao 2014). Sensitivity arms (response differential, bootstrap BCa CIs,
power, imputation Spearman) appended as sheets in
`c_data/10_DEP_supplementary.xlsx` (mirrors YvO's pattern).

Multi-threshold overview (FDR < 0.05/0.10, Pi < 0.05, P < 0.05/0.01,
with/without outlier removal) in `c_data/13_DEP_overview.csv` /
`14_DEP_overview.xlsx`.

Key outputs:

- `c_data/01_limma_DAList_{CRvH,CR}.rds`
- `c_data/03_combined_results_{CRvH,CR}.csv`
- `c_data/04_per_contrast_results/<contrast>.csv`
- `c_data/05_results_{CRvH,CR}.xlsx`
- `c_data/10_DEP_supplementary.xlsx`
- `b_reports/01_proteoDA_CRvH/`, `02_proteoDA_CR/`
- `b_reports/02_dep_overview.pdf`, `03_contrast_summaries/`

## Stage 04 -- Figures

`04_Figures/F01..F08` per-figure directories, each with `a_script/`,
`b_reports/`, `c_data/`. Per-figure orchestrator named
`90_stitch_F0x[_stream].R` (YvO convention). Shared infrastructure under
`shared/`: `style.R`, `volcano_ring.R`, `pathway_utils.R`,
`go_slim_categories.R`, `figure_supplement_helpers.R`, frozen `fgsea_CRvH.csv`.

`Reversal/` is a self-contained mini-stage exploring signature-reversal
methods between cancer-vs-healthy baseline and pooled CR training response.

`04_Figures/keys/` holds per-figure gene/protein key lists used by panel scripts.

`04_Figures/archive/` is a gitignored on-disk safety net containing
pre-reorg WIP (supp/ subdirs not yet tracked); triage when convenient.

Deferred for follow-up:

- `F00/` pipeline-QC supplementary figure (per YvO)
- `shared/comparison_panels/` and `shared/print_scale_apply_380mm.R`
- WGCNA module-trait analysis (YvO has it as F06 supp; not yet ported to CvH)

## Collaborator results: `02-03_Sam's_Results/`

Sam's pre-DEP DAList, his limma xlsx (8 contrasts), and our two-model DEP
rerun on his data live here. Three-way comparison artifact and narrative:

- `our_rerun/c_data/comparison_3way.xlsx` -- per-contrast 3-way table
- `our_rerun/c_data/comparison_3way_summary.csv` -- summary with rho, sig counts, overlap
- `SAM_VS_CVH_DIFF.md` -- methodology diff and interpretation

Headline: Spearman rho >= 0.795 between Sam's logFC and ours-on-his-data
across all 5 math-equivalent contrasts (peaks at 0.955 for Cancer_vs_Healthy).
Significance counts diverge driven by our `robust+trend` eBayes vs Sam's
default; underlying biology agrees.

## Reproducibility

- All stochastic steps set `set.seed(42)`.
- `01_normalization/c_data/02_normalized.csv` is the canonical handoff for
  deterministic float serialization across stages.
- Metadata validation is part of the active pipeline, not a manual pre-step.
- Raw inputs are not overwritten by stage scripts.
- Stage 01 cleans up the stray `Rplots.pdf` produced by proteoDA QC routines.

## Rerun order

1. `01_normalization/a_script/01_run_normalization.R`
2. `01_normalization/a_script/02_norm_reports.R`
3. `02_Imputation/a_script/apply_missforest.R`
4. `02_Imputation/a_script/02_imputation_reports.R`
5. `03_DEP/a_script/01_run_dep.R`
6. `03_DEP/a_script/02_dep_reports.R`
7. `03_DEP/a_script/03_dep_robustness.R`
8. `03_DEP/a_script/04_dep_overview.R`
9. Figure scripts and `90_stitch_F0x.R` orchestrators that consume refreshed outputs

Optional: `02_Imputation/a_script/benchmark/_run_all.R` to refresh the
17-method benchmark.

Optional: `02-03_Sam's_Results/our_rerun/a_script/run_dep_on_sam.R` followed
by `compare_3way.R` to refresh the collaborator comparison.

## Cross-pipeline crosswalk vs YvO_2026

| Element | YvO_2026 | CvH_2026 |
| --- | --- | --- |
| Stage triple | a/b/c convention | matches |
| Stage 01 normalization | cycloess + 4-method outlier | matches |
| Stage 02 imputation | 17-method benchmark + missForest | matches (subset to YvO's 17 in audit 2026-05-04; CvH-extra hybrids in `methods/` library) |
| Stage 03 DEP model | single full-cohort | two-model (CRvH + CR-only) |
| Sensitivity arms | embedded in 03 (xlsx sheets) | embedded in 03 (xlsx sheets in `10_DEP_supplementary.xlsx`) |
| Stage 04 figure layout | F00..F07 + 90_stitch_F0x.R | F01..F08 + 90_stitch_F0x.R (F00 deferred) |
| Shared figure infra | `shared/style.R`, volcano_ring, pathway_utils, frozen fgsea cache | matches |
| WGCNA (YvO F06) | dedicated supp + permutation cache | not yet ported |
| Collaborator inputs | n/a | `02-03_Sam's_Results/` |

See `docs/yvo_to_cvh_method_mapping.md` for the full transfer matrix.
