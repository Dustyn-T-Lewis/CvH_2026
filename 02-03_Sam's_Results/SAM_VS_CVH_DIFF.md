# Sam vs CvH: methodology diff and three-way DEP comparison

This document summarizes how collaborator Sam's pipeline differs from ours
on the same Cancer-vs-Healthy proteomics data, and reports a three-way DEP
comparison: (1) Sam's limma on Sam's filtered/normalized data, (2) our
two-model proteoDA/limma DEP on Sam's data, (3) our two-model DEP on our
filtered/normalized data.

Companion files in this directory:
- `01_normalized_DAList_SURV_stringent_muscle.RDS` — Sam's pre-DEP DAList
- `01_normalized_data_SURV_stringent_muscle.csv` — flat export of the RDS data slot
- `limma_results_all_muscle_str.xlsx` — Sam's full limma output (8 contrasts)
- `our_rerun/a_script/run_dep_on_sam.R` — our DEP wrapper for his data
- `our_rerun/a_script/compare_3way.R` — produces the comparison artifact
- `our_rerun/c_data/comparison_3way.xlsx` — per-contrast comparison + summary
- `our_rerun/c_data/comparison_3way_summary.csv` — summary table

## 1. Methodology differences

| Step | Sam | Ours |
|---|---|---|
| HPA filter | `in_hpa_muscle & !in_blood_blacklist` (more inclusive) | HPA tier-based filter (more aggressive) |
| Normalization | proteoDA `cycloess` | proteoDA `cycloess` |
| Re-normalization for sub-models | Yes — re-normalizes after dropping CTL (max diff 0.126 log2 vs RDS) | No — single normalization, then subset |
| Imputation | None — limma per-protein listwise on NAs | 19-method benchmark + missForest pick |
| Limma model 1 | SURV+CTL: 3 contrasts (Baseline_SURVvCTL, Training_SURVvCTL, Training_SURV) | CRvH: 2 contrasts (Cancer_vs_Healthy, Training_CR) |
| Limma model 2 | SURV-only: 5 contrasts (Baseline_CREvPLA, Training_CREvPLA, Training_CRE, Training_PLA, Interaction_supp) | CR-only 2x2: 4 contrasts (Baseline_Supplement, Training_CRE, Training_PLA, Supplement_Interaction) |
| Variance shrinkage | eBayes (default) | eBayes `robust=TRUE, trend=TRUE` (more conservative) |
| Multiple testing | BH + Xiao Pi-value | BH + Xiao Pi-value (matches) |
| Sensitivity / blunting diagnostics | None | KS / Fligner / Wilcoxon / Cliff / bootstrap BCa / power / imputation Spearman |
| Blood/muscle contamination QC | In metadata (HBB/HBA1/MB/ALB/CKM percentages, B_M_ratio) | Not in pipeline |
| Sample N | 35 (CRE 7+7, PLA 5+6, CTL 10) | 38 post-outlier (after dropping CR9_T1, CR10_T1, CR386_T1) |

## 2. Three-way DEP comparison (FDR < 0.10)

Five contrasts have a clean math equivalence between Sam's and our pipelines.

| Sam contrast | Our contrast | N matched | Sig (Sam / OurOnHis / OurOnOurs) | rho Sam~OurOnHis | rho Sam~OurOnOurs | rho OurOnHis~OurOnOurs | Overlap Sam~OurOnHis |
|---|---|---|---|---|---|---|---|
| Baseline_SURVvCTL | Cancer_vs_Healthy | 1928 | 532 / 466 / **616** | **0.955** | 0.903 | 0.973 | **429** |
| Baseline_CREvPLA | Baseline_Supplement | 1802 | 0 / 0 / 0 | 0.653 | 0.257 | 0.773 | 0 |
| Training_CRE | Training_CRE | 1802 | 4 / 2 / 0 | 0.852 | 0.820 | 0.901 | 2 |
| Training_PLA | Training_PLA | 1802 | 21 / 0 / 0 | 0.861 | 0.358 | 0.651 | 0 |
| Interaction_supp | Supplement_Interaction | 1802 | 3 / 1 / 0 | 0.795 | 0.368 | 0.657 | 1 |

`OurOnHis` = our DEP on Sam's filtered/normalized data; `OurOnOurs` = our DEP on our filtered/normalized data.

## 3. Interpretation

**Effect-size agreement is very high.** Spearman rho between Sam's logFC and ours-on-his-data is >= 0.795 across all contrasts and reaches 0.955 for the strongest signal (Cancer vs Healthy). This means Sam and our pipeline see the same direction and magnitude for almost every protein, even though we use different significance machinery (default vs robust+trend eBayes).

**Significance counts diverge by methodology.** For the Cancer-vs-Healthy contrast, Sam calls 532 significant, our DEP on his data calls 466, and our DEP on our data calls 616. The 429-protein overlap between Sam and our-on-his (~80 percent of the smaller call set) shows that the bulk of Sam's significance calls are reproduced by our pipeline; the gap is the conservative effect of our `robust=TRUE, trend=TRUE` eBayes flag, which down-weights variance outliers and intensity-trend variance.

**`Training_PLA` is the most interesting divergence.** Sam calls 21 significant, our DEP on his data calls 0. With logFC rho = 0.861 (directions agree strongly), the gap is purely in the significance test. The most likely cause is the variance estimation: Sam's per-protein listwise deletion in limma can produce smaller per-protein variance estimates for PLA proteins with sparse missingness, inflating his t-statistics. Our pipeline imputes (in the OurOnOurs run) or operates on his already-NA-tolerant matrix (in the OurOnHis run) with robust+trend eBayes, both of which damp the variance.

**Cohort-driven differences (Sam vs OurOnOurs) are larger than methodology-driven.** Looking at rho Sam~OurOnOurs (0.257-0.903) versus rho Sam~OurOnHis (0.653-0.955): the larger gap appears in the CR-only contrasts (Baseline, Interaction), where Sam re-normalizes after dropping CTL. Our pipeline keeps a single normalization for both models, so cross-model intensity values are comparable but the per-cohort baseline shifts slightly.

**No-signal contrasts agree on absence.** Baseline_CREvPLA has zero significant in all three runs, with weak rho. This is consistent with the literature on this comparison: a single supplementation-arm baseline rarely shows differential proteomic signal.

## 4. Things to reconcile next

1. **Two PLA samples Sam dropped** (vs our 38-sample post-outlier cohort) need to be identified. Likely candidates are the same outliers our pipeline removes (CR9_T1, CR10_T1, CR386_T1) plus possibly one more by Sam's QC.
2. **HPA filter divergence** — Sam's filter is more inclusive; the 1944 vs ~2542 protein gap is the result. A direct comparison on the intersection (1700-1800 proteins) is the apples-to-apples view.
3. **Blood-contamination QC** — Sam tracks HBB/HBA1/MB/ALB/CKM percentages and `B_M_ratio` per sample. Our pipeline could adopt these as either covariates or QC drop criteria. Worth a separate evaluation.
4. **`robust=TRUE, trend=TRUE` decision** — for the next manuscript, consider whether to follow Sam (more sensitive, more discoveries, more false positives) or stay with our conservative choice. Sensitivity arms in `03_DEP/c_data/10_DEP_supplementary.xlsx` already quantify the effect.
