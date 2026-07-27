# Decisions

## 2026-07-26 — Section-D review of `docs/cvh-pipeline-redesign`

| Question | Decision |
|---|---|
| ORA multiple-testing family | Keep per-database BH; rank top-N within database instead of across a mixed pile |
| Non-estimable contrasts in `sig_pi` | Code as `NA`, not `0`; make `overview_panels.R:99` NA-safe |
| Re-run cascade | Run in place, git as the safety net |
| Stale `CvH_pipeline.qmd` numbers | Banner as stale now, rewrite once after the cascade |
| Orphaned figure code | Delete `shared/volcano_ring.R` and the `style.R` orphans, each verified first |
| Convention sweep | `styler` only. Pipes stay mixed; `setwd(here::here())` stays as the house convention |
| Reproducibility pinning | No `.lintr` (conflicts with no-config-files preference, catches only cosmetics). Revisit `renv` after the numbers settle |
| Test gaps | Value tests for `fisher_z_ci`, `classify_proteins_f4`, `perm_p_unpaired`, `deduplicate_enrichment` |

### Applied in this review

Verified not to change results: `cleanup_after_workbook` preserve-guard (stale dir names,
unguarded recursive `unlink`, bare-dirname miss); `03_dep.qmd` read path; double `log2()` in
`02_normalization.qmd`; `SUPP_FILL` clobber in `panel_icc.R`.

Result-changing: `direction = "<"` pinned at the three *classifier* ROC sites; QRILC
`MARGIN = c(1L, 2L)`; F05 `CRvH_Baseline` weights aligned to `03_DEP`; `fry_camera` sheet label.

### F06 result changed

Unpinned ROC direction was corrupting the permutation null, not just the point estimates.
The null refits on permuted labels; k=1 selection inverts out-of-fold, and `auto` folded every
anti-predictive draw back above 0.5 (measured null mean 0.708, P(null ≥ 0.82) = 0.263). Pinned:
null mean 0.364, P = 0.003.

Module baseline is now the one positive result: AUC 0.82, perm p = 0.0030, q_BH = 0.018.
Training arm is below chance (0.389 / 0.458 / 0.354). Reversal arm perm p 0.114 → 0.0149.

Leave `uni_auc` / `perm_p_unpaired` / `perm_p_paired` on `direction = "auto"` — eigengene sign
is arbitrary and those deliberately fold observed and null together via `fold_auc()`.

### Pending

Re-run order: mscoreutils imputation → `03_DEP` both arms → F02, F03, F04, F05 → F06.
`03_DEP/b_imputed/` was already stale before this review (3 contrasts committed, 7 in code).
Stages 01–03 call `clear_dir()` at script top, ~200 lines before the first write.

## 2026-07-27 — dead-weight sweep

Scope settled up front: roots are the README run order plus the two `tests/`
directories and `network_validation.R`. Cut tiers 0-3 (OS cruft, dead source,
orphan outputs, `docs/` scratch) and left tier 4, the three-arm imputed DEP, since
the README defines it as the concordance check.

Write-only artifacts keep their write calls. Only files with no writer at all were
removed, so the 280 proteoDA `static_plots/` and the DA report HTML stay.

### What went

- 9 unreferenced symbols, 252 lines: `matrix_to_df`, `read_matrix_sheet`,
  `export_slim_mapping`, `classify_pathway_func`, `run_enrichment_pipeline`,
  `perm_p_unpaired`, `GROUP_COLORS`, `reversal_rrho2`, `run_reversal_analysis`.
  Each appeared exactly once, at its own definition.
- Two forwarding `style.R` shims in F05 and F06; their three callers now source
  `shared/style.R` directly.
- Five F05 `c_data` CSVs stranded by the three-stage driver rewrite (76ef03a):
  `01_panel_A_heatmap_data`, `03_panel_B_eigengene_data`,
  `03_panel_B_heatmap_zscores`, `03_panel_B_triptych_enrichment`,
  `04_panel_C_hub_proteins`.
- 42 files of `docs/` process scratch, `Rplots.pdf`, 16 `.DS_Store`.

### What stayed, and why

All 69 R files were reachable, so no script was dead. The two `construction.R`
files are distinct and both sourced by full path from `02_clustering.R`, so there
was no duplicate to merge. The proteoDA `static_plots/` looked orphaned only
because proteoDA builds its filenames internally.

### Verification

Full cold re-run, all 16 commands, zero failures. Of 110 pre-sweep data files, 97
byte-identical and 9 differing, all `.xlsx`. Comparing those 9 sheet by sheet gave
78 sheets with zero value differences, so the byte drift is zip timestamp
metadata. F05 re-ran again after its five orphans were deleted: none regenerated,
every surviving output byte-identical.

`docs/` was gitignored, so that deletion has no git backup and cannot be undone.
`decisions.md` is now tracked to stop the same thing happening to this log.
