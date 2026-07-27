# CvH Proteomics Analysis

Label-free DIA-MS on skeletal muscle, comparing cancer-recovery participants (CR)
against healthy controls (H). The CR arm has pre- and post-training biopsies; the
healthy controls give one baseline biopsy each.

39 samples survive QC: 10 healthy at T1, 15 CR at T1, 14 CR at T2.

## The model

DEP fits `~ 0 + model_cell + (1 | Subject_ID)` over five cells: `H_pre`, plus
`CRE`/`PLA` crossed with `pre`/`post`. The pre-to-post pairing is carried by
`duplicateCorrelation` on `Subject_ID` (consensus rho = 0.248). Both supplement
arms enter the model, so every estimate is supplement-adjusted.

Three contrasts average the two supplement arms 50:50:

- `CRvH_Baseline` — how far disease sits from healthy (D)
- `CR_Training` — what training does (T)
- `Resid` — what is left afterwards (R)

They close a triangle: `Resid = CRvH_Baseline + CR_Training`.

Four more contrasts test the supplement arms directly (`Baseline_Supplement`,
`Training_CRE`, `Training_PLA`, `Supplement_Interaction`) and land in the same
table. So `combined_results_pi.csv` carries seven contrasts, not three.

Significance is the Pi-score, `Pi = P.Value^|log2FC|` (Xiao 2014), at `Pi < 0.05`,
with BH-FDR as a secondary criterion. A protein whose contrast cannot be estimated
gets `sig_pi = NA`, never `0`.

## Stages

Each stage keeps scripts in `a_script/`, renders in `b_reports/`, and outputs that
later stages read in `c_data/`. Stage outputs are committed, so a fresh clone can
start anywhere without re-running the stage above.

`00_input/` holds the DIA-NN intensity matrix, the sample sheet, HPA annotations,
and the red-blood-cell proteome reference.

`01_Filtering/` applies HPA presence, blood-contaminant removal with myonuclei
rescue, red-cell tracking removal, a missingness filter, and consensus outlier
detection, ending at `DAList_filtered.rds`. The cascade runs 3172 proteins down to
2176:

| Step | Remaining | Removed |
| --- | --- | --- |
| Raw input | 3172 | — |
| HPA presence | 2962 | 210 |
| Blood contaminants | 2811 | 151 |
| Red-cell tracking | 2437 | 374 |
| Missingness (>=5 in >=1 group) | 2176 | 261 |

`02_Normalization/` runs cycloess. Its `imputation/` subfolder holds three arms:
`imp4p`, an MsCoreUtils hybrid, and `missForest`, each writing
`DAList_imputed_<method>.rds`. The figures read the missForest arm.

`03_DEP/` splits in two. `a_non_imputed/` is the primary analysis: limma with
`duplicateCorrelation`, seven contrasts, Pi-score. `b_imputed/` re-runs the same
DEP on all three imputed matrices as a concordance check.

`04_Figures/` holds `F01`–`F06` plus `shared/`.

## Figures

| Figure | What it shows |
| --- | --- |
| `F01_Phenotype` | Phenotype and strength; supplement adds LBM, chest press, leg extension, grip |
| `F02_Proteome_Overview` | PCA, DEP counts, effect sizes, overlap, direction, pathways; QC supplement covers CV, ICC, dbRDA |
| `F03_Enrich_Volcanoes` | `enrichVolcano` ring grid over the Model-1 contrast trio |
| `F04_Reversal` | Reversal landscape, trajectory clustering, fry rotation, plus diagnostics |
| `F05_WGCNA` | Module card: counts, member response (fGSEA NES over fry), eigengene trajectory, ORA |
| `F06_Prediction` | Three feature spaces by two outcomes, with a circularity ladder |

The WGCNA network uses bicor, signed, soft power 14, and yields five modules:
turquoise 675, blue 447, brown 408, yellow 182, green 59, plus 405 unassigned.

## Running it

```sh
Rscript 01_Filtering/a_script/01_run_filtering.R
Rscript 02_Normalization/a_script/01_run_normalization.R
Rscript 02_Normalization/imputation/a_script/a_imp4p.R
Rscript 02_Normalization/imputation/a_script/b_mscoreutils.R
Rscript 02_Normalization/imputation/a_script/c_missforest.R      # figures read this arm

Rscript 03_DEP/a_non_imputed/a_script/01_run_dep.R               # primary DEP
Rscript 03_DEP/b_imputed/a_script/01_run_dep_imputed.R           # concordance check

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

`04_Figures/F05_WGCNA/a_script/supp/network_validation.R` is a long parameter
sweep. Run it on demand, not from the driver.

Three test suites cover the shared statistics helpers (59 tests):

```sh
Rscript -e "testthat::test_file('04_Figures/F05_WGCNA/tests/test-wgcna_stats.R')"
Rscript -e "testthat::test_file('04_Figures/F06_Prediction/tests/test-prediction_utils.R')"
Rscript -e "testthat::test_file('04_Figures/shared/tests/test-stats.R')"
```

`docs/decisions.md` records the choices behind the pipeline and why they were
made. Read it before changing anything statistical.

## Rules that matter

Paths resolve from the project root; figure scripts anchor with
`setwd(here::here())`. Anything stochastic uses `set.seed(42)`, and that includes
plot jitter — F01 pins it with `position_jitter(seed = 42)` so the panels render
byte-identically twice in a row.

Every figure writes to `b_reports/{main,supp}/{pdf,png}/`, per-panel output one
level down in `panels/`. Style lives in `shared/style.R`; the statistics the
panels share live in `shared/stats.R`, which is what the third test suite covers.

The primary DEP uses the non-imputed matrix. Imputation feeds only the
concordance check, PCA, and WGCNA.

Repeated measures block on `Subject_ID`.

`shared/wgcna_stats.R` and `shared/prediction_utils.R` call WGCNA
namespace-qualified and never attach it. `WGCNA::cor()` returns a 1x1 matrix and
will silently break any bare `cor()` in a script sourced afterwards.

## Known environment issue

Without XQuartz, `cairo_pdf` cannot load, `get_pdf_device()` falls back to base
`pdf()`, and that device cannot encode `Π`, `ρ`, or `✱`, so PDF renders drop
them. PNG output is fine. Installing XQuartz restores the glyphs.
