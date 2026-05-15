# F06_WGCNA — Sam-parallel WGCNA Figure

Sam's N=35 normalized DAList (`01_normalized_DAList_SURV_stringent_muscle.RDS`,
1944 proteins) run through the same WGCNA pipeline used for YvO F06 and main
CvH F06. Signed Pearson, blockwiseModules, LMM contrasts per module eigengene.

## Script execution order

```
Rscript 02-03_Sam's_Results/04_Figures/F06_WGCNA/a_script/00_run_wgcna.R
Rscript 02-03_Sam's_Results/04_Figures/F06_WGCNA/a_script/01_main_panels.R
Rscript 02-03_Sam's_Results/04_Figures/F06_WGCNA/a_script/02_supp_panels.R
Rscript 02-03_Sam's_Results/04_Figures/F06_WGCNA/a_script/90_stitch_F06.R
```

Run all from the project root (`A_CvH_2026/`).

## WGCNA parameters (mirrors YvO + main CvH)

- networkType = "signed", TOMType = "signed"
- corType = "pearson"
- minModuleSize = 30
- mergeCutHeight = 0.25
- deepSplit = 2
- pamRespectsDendro = FALSE
- Soft power: auto-selected (R^2 > 0.87 or elbow)

## Main panels (2)

- **Panel A**: Module-trait LMM heatmap
  Rows = modules (ME), cols = design contrasts (Cancer_vs_Healthy, Training_CR,
  Training_CRE, Training_PLA, Supplement_Interaction + Baseline_Supplement).
  Overlay: r_equiv, BH-corrected stars, nominal p dashed border.
  Right flange: module protein count bar.
- **Panel B**: Eigengene scatter for top 2 modules vs strongest trait.
  Pre/post boxplot per Group_Time (CRE_T1, CRE_T2, PLA_T1, PLA_T2, H_T1)
  with paired lines (subject_key), faceted by module.

## Supplementary panels (~5)

- **Supp A**: Soft-power selection plot (R^2 and mean connectivity vs power)
- **Supp B**: WGCNA dendrogram + module color bar
- **Supp C**: Module size distribution (bar chart, labelled)
- **Supp D**: Module-trait simple Pearson heatmap (full, all 10 traits)
- **Supp E**: Top-module hub network (igraph, top-15 kME hub proteins,
              connectivity weighted edges, coloured by module)

## Inputs

- Sam's DAList: `02-03_Sam's_Results/00_input/01_normalized_DAList_SURV_stringent_muscle.RDS`
- Sam's DEP results for LMM contrasts: `02-03_Sam's_Results/03_DEP/c_data/01_limma_DAList_CRvH.rds`
                                       `02-03_Sam's_Results/03_DEP/c_data/01_limma_DAList_CR.rds`

## Outputs

- `c_data/wgcna/`: network object, SFT summary, module assignments, hub
  proteins, enrichment, module-trait CSVs, LMM CSV
- `c_data/`: panel-level CSVs, RDS helpers (MEs, datExpr, meta)
- `b_reports/main/`: MAIN_F06_composite.{pdf,png}
- `b_reports/supp/`: SUPP_F06_composite.{pdf,png}
