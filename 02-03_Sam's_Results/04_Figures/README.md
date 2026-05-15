# 02-03_Sam's_Results/04_Figures/

Self-contained figure pipeline driven from Sam Cooper's pre-processed muscle
proteomics outputs + our DEP rerun on his data + the 3-way comparison.

Layout mirrors `A_YvO_2026/04_Figures/` (canonical YvO pattern), one level
deeper because everything here is scoped to Sam's analysis stream.

## Figures
- `F00_dataset_overview/` — Sam's filter cascade + cohort QC
- `F01_DEP_landscape/` — DEP volcanos + count summaries across Sam's 8 + our mapped contrasts
- `F02_blood_contamination/` — Sam-unique blood-marker QC (HBB/HBA1/MB/ALB/CKM, B_M_ratio)
- `F03_method_comparison/` — Sam vs OurOnHis vs OurOnOurs visualization

Each figure has `a_script/INTENT.md` outlining intended panels and inputs.

## Conventions
- `a_script/` — panel scripts + `90_stitch_FXX.R` composite assembler
- `b_reports/main/{pdf,png/panels}/` — main figure outputs (MAIN_FXX_*.pdf)
- `b_reports/supp/{pdf,png/panels}/` — supplementary outputs (SUPP_FXX_*.pdf)
- `c_data/` — per-figure data dictionaries + intermediate artifacts

## Shared utilities
- `build_data_index.R` — provides `sam_idx` (named list of file paths)
- `shared/style.R` — sources main pipeline style (palettes, themes, sizing)

## Reading order for someone new
1. `SAM_VS_CVH_DIFF.md` (the analytical narrative)
2. `build_data_index.R` (the data map)
3. Any `FXX_*/a_script/INTENT.md` (per-figure intent)
