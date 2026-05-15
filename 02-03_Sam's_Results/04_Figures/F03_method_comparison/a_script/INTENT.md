# F03 — 3-way method comparison (Sam vs OurOnHis vs OurOnOurs)

Converts the existing SAM_VS_CVH_DIFF.md table + comparison_3way.xlsx into
publishable figures.

## Intended panels (main)
- A. Spearman rho matrix (3x3) per contrast (5 panels)
- B. logFC concordance scatters: Sam vs OurOnHis (5 contrasts as facets)
- C. DEP count bar: Sam / OurOnHis / OurOnOurs side-by-side per contrast
- D. UpSet of significant-protein sets (FDR<0.10) across 3 methods, Cancer_vs_Healthy contrast

## Intended panels (supp)
- A. logFC scatters Sam vs OurOnOurs (the cohort-driven gap)
- B. methodology diff table (rendered as gt::gt() from SAM_VS_CVH_DIFF.md section 1)
- C. Training_PLA divergence deep-dive: per-protein variance estimate Sam vs OurOnHis (explains the 21 vs 0 sig gap)

## Inputs
`sam_idx$comparison$xlsx`, `sam_idx$comparison$summary`,
`sam_idx$sam$limma_xlsx`, `sam_idx$ours_on_his$combined_CRvH/CR`.

## Composite
`90_stitch_F03.R` -> `b_reports/main/pdf/MAIN_F03_method_comparison.pdf`.
