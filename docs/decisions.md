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

## 2026-07-27 — trim to the minimum that runs

Scope: cut anything not feeding a figure or a documented result, collapse
duplicated logic, strip AI tells. Roots unchanged from the dead-weight sweep.
All 67 R files were reachable again, so nothing was cut wholesale; the dead
weight was inside live files.

### Result-changing, each verified against a before/after run

| Change | Effect |
|---|---|
| `circularity_ladder.R` read `chosen_power` off the committed network instead of hardcoding `12` | Tier 2 refits now run at 14, the power the network was actually built at. Mean Jaccard 0.451 -> 0.455, failed folds 71 -> 66. Tier 1 byte-identical. `baseline/turquoise` unmoved (0.960 in-sample, 0.953 out-of-fold), so the one positive F06 result stands. `training/turquoise` crossed 1 -> 0 failed folds and so gained a stability box in F06 panel C. |
| F01 jitter pinned with `position_jitter(seed = 42)` | F01 was nondeterministic: two runs of identical code produced different PNGs, because every panel jittered unseeded against the README's own `set.seed(42)` rule. Point positions are now fixed. All seven F01 audit CSVs stayed byte-identical, so no statistic moved. |
| F01 panels D-G collapsed onto `pre_post_panel()` | Adopted panel D's behaviour per the brief: `max(abs(delta))` for the bracket position (robust when every delta is negative) and jitter alpha 0.35. E/F/G previously used `max(delta)` and alpha 0.5. 504 lines -> 196. Audit CSVs byte-identical. |
| F02 contrast palettes derived from `CONTRAST_COLORS` | F02 was painting `CRvH_Baseline` in `#D6604D`, the same hex every other figure uses for `Training_PLA`, and `Resid` in `Baseline_Supplement`'s teal. `SUPP_PAL` derived identically, so only the main figure recoloured. No count changed. |

### Cut

- `reversal.R` 209 -> 58 lines. `compute_phi`, `directional_asymmetry`,
  `reversal_permutation_null` and `reversal_rotation_test` had zero callers;
  the analyses they implement were re-implemented inline in
  `SUPP_directional_asymmetry.R`, `SUPP_melov_proportion.R`,
  `panel_F_trajectory.R` and `panel_D_fry.R`. Only `load_reversal_table` is
  live, via `f04_data.R:6`. The four inline copies are left alone: folding them
  back onto one engine would change results and is a separate decision.
- 35 top-level objects computed and never read, found by iterating "name occurs
  once repo-wide" to a fixpoint over five rounds. Every one was a pure
  expression, so no write or other side effect was lost.
- `style.R:add_tag` (shadowed by a different local definition in
  `90_stitch_F02.R`) and `figure_supplement_helpers.R:read_sheet_df`.
- Four F05 artifacts whose contents already live inside the consolidated
  `wgcna_network.rds`: the raw `net`, `sft_summary`, a byte-identical second
  copy of `wgcna_module_assignments.csv`, plus `kME_all.rds` and
  `key_modules.txt`. `c_data/wgcna/wgcna_lmm_contrast_audit.csv` had no writer
  at all. `c_data/wgcna/sft_fitIndices.rds` stays: `supp/construction.R:7`
  reads it.

### Fixed

- `SUPP_fry_leading.R` read `panel_D_fry/driving_proteins.csv` with an xlsx
  fallback, but `90_stitch_F04.R:183` deletes that directory at `:180` and only
  sources the script at `:211`. The CSV branch could never run. It now reads the
  workbook directly, which is what always happened.
- `reversal.R:4` named `reversal_inputs.R` and `00_build_fgsea_cache.R`, neither
  of which exists, and claimed "no file IO" while defining `load_reversal_table`.
- `overview_panels.R` and `01_enrich_volcanoes.R` read the fgsea cache with no
  guard. `panel_B_nes_scatter.R:22` already had one; they now match it.

### Conventions settled

`# === SECTION ===` was recorded as the house header form but appeared **zero**
times in the repo. What existed was 160 banner separators in four other forms
(`# --- x ---` x43, `# ── x ──` x81, `# -- x ----` x15, `# ═══` x12). All 160 are
gone, rewritten to plain one-line comments or deleted where they carried no
words. Treat a plain `# label` as the house form; there is no rule-decorated
variant any more.

23 ceremonial `message("...done")` / `cat("...done\n")` calls removed. One
terminal message per top-level driver stays, as does every message carrying a
count, p-value or dimension.

### Gene dedup divergence — documented, deliberately not unified

The same operation, collapsing multiple rows per gene, is done five different
ways:

| Rule | Site |
|---|---|
| `slice_max(abs(t))` | `build_fgsea_cache.R:30` |
| `pivot_wider(values_fn = mean)` | `01_module_stats.R:115`, `:160` |
| `tapply(mean)` | `wgcna_stats.R:15` |
| `slice_min(pi_score)` | `01_enrich_volcanoes.R:44` |
| bare `distinct(gene, .keep_all = TRUE)` | `panel_E_rrho2.R:23`, `panel_G_resid_volcano.R:50` |

`slice_max(abs(t))` is the most defensible and `distinct()` the least, but
switching any caller changes that figure's numbers. Left as-is on purpose;
settle it as a statistical question, not a refactor.

### Reversal design note, moved here from `reversal.R`

Shared-baseline circularity: D and T share the CR_pre samples, so a structural
negative `cor(D, T)` is mathematically guaranteed (Smyth & Altman 2013, PMID
23705896). Every directional claim is therefore tested against a protein-label
permutation null asking whether the disease-DEP set reverses *more* than random
proteins under the same shared-baseline structure, not merely whether reversal
exceeds zero. fry (rotation, within the limma model) and RRHO2 (rank-based)
corroborate from model-aware and non-parametric frameworks.

Method lineage: Melov 2007 (PMID 17520024, proportion + permutation), Robinson
2017 (PMID 28273480, logFC correlation), Wu & Smyth 2010/2012 (PMID 20610611 /
22638577, ROAST/CAMERA), Cahill 2018 (RRHO2), Smyth & Altman 2013 (PMID
23705896, shared baseline).

### Still duplicated, not yet collapsed

Measured, left for a follow-up because each changes many call sites at once:

- The PNG+PDF `ggsave` pair, 88 calls across 44 sites. Flags are inconsistent:
  16 of 44 PNG calls set `bg = "white"`, 28 do not.
- The F04 panel prologue (`setwd` + `source` + path constants + `dir.create` +
  `get_pdf_device`), 16 files, ~165 lines. 78 `dir.create` calls repo-wide.
- `clear_dir` defined four times identically across stages 01-03, which have no
  shared file to hold it.
- `CONSOLIDATED_PATHWAY_ORDER` / `CONSOLIDATED_COLORS` duplicated between
  `pathway_utils.R:316` and `go_slim_categories.R:42`, the latter behind an
  `if (!exists(...))` load-order guard.
- `compute_cv` (twice), `boot_median_ci` (twice, byte-identical).
