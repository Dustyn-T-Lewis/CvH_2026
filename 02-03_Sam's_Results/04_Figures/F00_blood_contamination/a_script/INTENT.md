# F00 — Blood/muscle contamination QC (Sam-unique)

The most novel figure in this subdir — Sam tracks per-sample blood-marker
percentages (HBB, HBA1, MB, ALB, CKM) and a B_M_ratio that our main pipeline
doesn't have. This figure puts those QC metrics into a publishable form.

## Intended panels (main)
- A. per-sample stacked bar of HBB/HBA1/MB/ALB/CKM percent of total signal
- B. B_M_ratio dotplot by Group_Time with sample IDs (highlight outliers)
- C. correlation of B_M_ratio with PCA outlier consensus from our pipeline
- D. how Sam's HPA filter (`in_hpa_muscle & !in_blood_blacklist`) vs our HPA tier-based filter intersect: Venn + retention rates

## Intended panels (supp)
- A. per-protein blood-marker concentration scatter (Sam vs ours)
- B. impact of B_M_ratio cutoff on Cancer_vs_Healthy DEP count

## Inputs
Sam's DAList metadata slot (B_M_ratio etc.) -> `sam_idx$sam$dalist_rds` then
`as.data.frame(sam$metadata)`.

## Composite
`90_stitch_F02.R` -> `b_reports/main/pdf/MAIN_F00_blood_contamination.pdf`.

## Why this figure matters
Contaminant filtering is a high-stakes upstream decision. Showing Sam's
metric and how it correlates with our outlier-consensus framework gives the
manuscript a methodological transparency story — and answers the "did blood
contamination drive the Cancer_vs_Healthy signal?" reviewer question.
