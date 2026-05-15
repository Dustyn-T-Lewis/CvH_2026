# F01 — Phenotype Figure (Sam-parallel YvO F01 port)

YvO F01 transposed to Sam's CvH muscle proteomics cohort (N=25 SURV).
Displays pre/post phenotype change by supplement arm (CRE vs PLA) with
paired delta bars, mirroring the Young vs Old structure of YvO F01.

## Data source
`02-03_Sam's_Results/00_input/01_normalized_DAList_SURV_stringent_muscle.RDS`
— `$metadata` slot, wide-format pre_/post_ phenotype columns.

## Group axis
CRE vs PLA within SURV. CTL excluded (no training intervention pre/post).

## Main panels (3)
- A. Leg Extension 1RM pre/post (strength recovery metric)
- B. DXA Lean Body Mass pre/post (direct YvO B match)
- C. Appendicular Lean Mass pre/post (closest to YvO VL thickness)

## Supp panels (3)
- A. Chest Press 1RM pre/post
- B. Grip Strength pre/post
- C. Sit-to-Stand Max Power pre/post

## Composite outputs
`90_stitch_F01.R` -> `b_reports/main/{pdf,png}/MAIN_F01_composite.{pdf,png}`
                  -> `b_reports/supp/{pdf,png}/SUPP_F01_composite.{pdf,png}`
