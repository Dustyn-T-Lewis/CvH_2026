# F01 — DEP landscape on Sam's data

Main figure summarising DEP across all 8 Sam contrasts + our 6 mapped contrasts.

## Intended panels (main)
- A. volcano grid: Sam's 5 math-equivalent contrasts (Baseline_SURVvCTL, Baseline_CREvPLA, Training_CRE, Training_PLA, Interaction_supp)
- B. UpSet of significant sets (Sam FDR<0.10) across the 5 contrasts
- C. DEP count bar (Sam vs OurOnHis vs OurOnOurs)
- D. logFC heatmap of top-50 Cancer_vs_Healthy proteins

## Intended panels (supp)
- A. volcanos for Sam's other 3 contrasts (Training_SURV etc)
- B. pi-score histograms
- C. per-contrast direction tables

## Inputs
`sam_idx$sam$limma_xlsx` (Sam's 8 contrasts), `sam_idx$ours_on_his$combined_CRvH/CR`,
`sam_idx$comparison$xlsx` (counts/overlaps).

## Composite
`90_stitch_F01.R` -> `b_reports/main/pdf/MAIN_F01_DEP_landscape.pdf` + supp.
