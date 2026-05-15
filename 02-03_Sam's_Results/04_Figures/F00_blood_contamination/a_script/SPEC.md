# F02 — Blood/muscle contamination QC (design spec)

Date: 2026-05-15
Status: approved (brainstorming) → ready for implementation plan
Author: claude + dtl0018

Companion: `INTENT.md` (high-level intent). This file is the implementation spec.

## Purpose

Use Sam Cooper's per-sample blood-marker tracking (`HBB_pct`, `HBA1_pct`,
`MB_pct`, `ALB_pct`, `CKM_pct`, `blood_pct`, `muscle_pct`, `B_M_ratio`) to
build a publication-grade QC figure that:

1. Surfaces blood contamination as a per-sample, per-group story.
2. Cross-validates Sam's `B_M_ratio` against our 4-method outlier consensus.
3. Documents how Sam's stringent HPA + blood-blacklist filter compares with
   our HPA tier-based filter at the protein-set level.
4. Reassures the reviewer that `Cancer_vs_Healthy` is robust to contamination
   via a B_M_ratio cutoff sensitivity analysis.

## Inputs (canonical paths via `sam_idx`)

| Slot | Path | Used in |
|---|---|---|
| `sam_idx$sam$dalist_rds` | `02-03_Sam's_Results/00_input/01_normalized_DAList_SURV_stringent_muscle.RDS` | panels A–D, supp A, supp C |
| Our normalization intermediates | `A_CvH_2026/01_normalization/c_data/00_report_intermediates.rds` | panel C (`outlier_diag$mahal_dist`, `consensus_outlier`) |
| Our normalized DAList | `A_CvH_2026/01_normalization/c_data/03_DAList_normalized.rds` | panel D Venn (our kept-protein universe), supp A |
| Our DEP results (Cancer_vs_Healthy) | `02-03_Sam's_Results/03_DEP/c_data/04_per_contrast_results/Cancer_vs_Healthy.csv` | supp D |
| Our imputed matrix | `A_CvH_2026/02_Imputation/c_data/01_imputed.csv` | supp B (refit base) |

## Sample-ID normalization

Sam: `CR006_T1` (zero-padded 3-digit); Ours: `CR6_T1`.

Helper in `_shared_keys.R`:

```r
normalize_sample_id <- function(x) {
  sub("^CR0*([0-9]+)_T([12])$", "CR\\1_T\\2", x)
}
```

Apply to *both* sides before any join. Returns the canonical short form
(`CR6_T1`). All `c_data/*.csv` outputs use the short form.

## Main figure — 2×2 grid

### Panel A — per-sample blood-marker stacked bar

- **Data:** `as.data.frame(sam$metadata)`, columns
  `c("sample_id", "Group_Time", "HBB_pct", "HBA1_pct", "MB_pct",
    "ALB_pct", "CKM_pct", "B_M_ratio")`.
- **Pivot:** long on the 5 `*_pct` columns → `(sample_id, marker, pct)`.
- **Order:** `sample_id` as a factor with levels = sample_id sorted by
  ascending `B_M_ratio`.
- **Geom:** `geom_col(position = "stack")`, fill = marker.
  Palette = 5 distinguishable hues (e.g. `RColorBrewer::brewer.pal(5, "Set2")`).
- **Decoration:**
  - Thin colored band along x-axis encoding `Group_Time`
    (use `GROUP_COLORS` from `shared/style.R`).
  - Secondary annotation row above showing `B_M_ratio` as small numeric labels
    (rotated 90°).
- **Axes:** y = "% of total signal", x = "Sample (sorted by B_M_ratio)".
- **Title:** "Per-sample blood-marker share (sorted by B_M_ratio, ascending)".
- **Output data:** `c_data/panel_A_stacked_data.csv`
  (`sample_id, Group_Time, marker, pct, B_M_ratio`).

### Panel B — B_M_ratio dotplot by Group_Time

- **Data:** same as panel A, one row per sample.
- **Groups:** 5 `Group_Time` levels — `CRE_T1, CRE_T2, PLA_T1, PLA_T2, CTL_T1`
  (CTL has no T2 in Sam's data; confirm at scaffold).
- **Geom:** `geom_jitter(width = 0.15)` + `stat_summary(fun = mean, geom = "crossbar")`
  + `stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar")`.
- **Outlier highlight:** cohort-wide top decile of `B_M_ratio` across N=35
  (top 4 samples) colored red and labeled with `sample_id` via `ggrepel::geom_text_repel`.
- **Annotation:** `kruskal.test()` p-value in upper-right corner.
- **Axes:** y = "B_M_ratio (blood-to-muscle signal ratio)", x = "Group_Time".
- **Color:** points fill = `GROUP_COLORS[Group_Time]`.
- **Output data:** `c_data/panel_B_dotplot_data.csv`
  (`sample_id, Group_Time, B_M_ratio, top_decile_flag`).

### Panel C — cross-pipeline scatter (B_M_ratio vs Mahalanobis)

- **Data join:**
  1. `sam_meta <- as.data.frame(sam$metadata)` (N=35)
  2. `our_outlier <- readRDS(".../00_report_intermediates.rds")$outlier_diag` (N=41)
  3. `our_outlier <- our_outlier |> mutate(sample_id = normalize_sample_id(Col_ID))`
  4. `sam_meta   <- sam_meta   |> mutate(sample_id = normalize_sample_id(sample_id))`
  5. `joined <- inner_join(sam_meta, our_outlier, by = "sample_id")`
  6. Expected join size: ~32–35 paired samples.
- **Geom:** `geom_point(aes(B_M_ratio, mahal_dist, color = consensus_outlier, shape = consensus_outlier))`.
  - FALSE → grey circle; TRUE → red triangle (larger).
- **Labels:** `geom_text_repel` on rows where `consensus_outlier == TRUE` OR
  `B_M_ratio >= quantile(B_M_ratio, 0.9)`.
- **Annotation:**
  - Spearman ρ + p (use `cor.test(method = "spearman", exact = FALSE)`).
  - Pearson r + p.
  - n joined samples.
  - Place all in upper-left text block.
- **Axes:** x = "Sam's B_M_ratio", y = "Our Mahalanobis distance (PC1–2)".
- **Output data:** `c_data/panel_C_join_data.csv`
  (`sample_id, Group_Time, B_M_ratio, mahal_dist, n_flags, consensus_outlier, label_this`).

### Panel D — filter intersection (two Venns)

- **Data sources:**
  - Sam-kept: `rownames(sam$annotation)` (UniProt IDs, length 1944).
  - Our-kept: `rownames(our_dal$annotation)` (length ~2582) from
    `03_DAList_normalized.rds`.
  - 5 blood-marker UniProts: pull from Sam's annotation by gene symbol
    `c("HBB", "HBA1", "MB", "ALB", "CKM")`.
- **D(i): kept-protein Venn (2-set).**
  - Use `ggVennDiagram::ggVennDiagram` or `eulerr::euler` (decide at impl;
    eulerr handles 2 sets cleanly).
  - Labels: counts only (no protein names).
  - Fill: light tint of two contrasting hues.
  - Title: "Kept proteins after filtering".
- **D(ii): blood-marker fate.**
  - If 2-set Venn is clean for n=5: render as a small Venn with point labels.
  - Else (more likely): 5-row indicator chart, rows = `c("HBB", "HBA1", "MB", "ALB", "CKM")`, two columns = "Sam-kept", "Ours-kept", cell = ✓/✗.
  - Note: `in_blood_blacklist == TRUE` typically excludes plasma proteins
    (HBB/HBA1/ALB likely dropped; MB/CKM kept).
  - Title: "5 named blood markers — fate per pipeline".
- **Composition:** `cowplot::plot_grid(p_venn_proteins, p_venn_markers, ncol = 2, rel_widths = c(2, 1))`
  inside the panel D slot.
- **Output data:**
  - `c_data/panel_D_venn_protein_sets.csv` columns:
    `uniprot_id, gene, in_sam, in_ours, set_membership` where `set_membership ∈ {both, sam_only, ours_only}`.
  - `c_data/panel_D_venn_blood_markers.csv` columns:
    `uniprot_id, gene, in_sam, in_ours, in_blood_blacklist`.

## Supplementary figure — 2×2 grid

### Supp A — 5-marker intensity scatter (Sam vs ours)

- **Data:**
  - Sam normalized matrix: `sam$data` (proteins × samples).
  - Our normalized matrix: `our_dal$data`.
  - For each of `c("HBB", "HBA1", "MB", "ALB", "CKM")`:
    - find UniProt in both annotations,
    - subset to shared sample IDs (after `normalize_sample_id`),
    - extract log-intensity vectors.
- **Geom:** 5-facet (`facet_wrap(~ gene, ncol = 3, scales = "free")`) scatter:
  x = our intensity, y = Sam intensity, point per shared sample.
- **Decoration:** identity line `geom_abline(slope = 1, intercept = 0, linetype = "dashed")`,
  Pearson r annotation per facet via `ggpubr::stat_cor` or manual.
- **Output data:** `c_data/supp_A_intensity_pairs.csv`
  (`gene, uniprot_id, sample_id, intensity_ours, intensity_sam`).

### Supp B — B_M_ratio cutoff sensitivity

- **Inputs:** Sam's DAList (`sam$data`, `sam$metadata`, `sam$design`,
  `sam$annotation`).
- **Procedure:**
  1. Define cutoffs = `c(0, 5, 10, 15, 20)` % top-B_M_ratio samples dropped.
  2. For each cutoff:
     - Identify samples to drop (highest B_M_ratio first).
     - Subset DAList in place (data, metadata, design rows).
     - Re-run `proteoDA::fit_limma_model(dal, contrasts = "Cancer_vs_Healthy", robust = TRUE, trend = TRUE)`.
     - Apply BH FDR; count proteins with `FDR < 0.10`.
- **No re-imputation needed.** Sam's data is already normalized + cycloess,
  and limma handles residual missingness via `lmFit`. If MAR concerns arise,
  document but do not re-impute (justification: this is a sensitivity
  check on cohort composition, not on imputation method).
- **Geom:** `geom_line() + geom_point()`, x = % dropped, y = n DEPs at FDR<0.10.
  Horizontal dashed reference line at "main result" (cutoff = 0).
- **Output data:** `c_data/supp_B_cutoff_sensitivity.csv`
  (`cutoff_pct, n_samples_dropped, n_proteins_used, n_DEP_fdr10, dropped_sample_ids`).

### Supp C — correlation heatmap of blood markers

- **Data:** N=35 × 6 matrix (5 `*_pct` + `B_M_ratio`).
- **Compute:** Spearman correlation matrix.
- **Geom:** `ggcorrplot::ggcorrplot` or manual `geom_tile`, diverging palette
  (`scale_fill_gradient2(low, mid = "white", high, limits = c(-1, 1))`),
  cell-label ρ to 2 decimals.
- **Output data:** `c_data/supp_C_corr_matrix.csv` (long form: var1, var2, rho, p).

### Supp D — blood markers in our DEP results

- **Data:** `read_csv("...Cancer_vs_Healthy.csv")` from our 03_DEP outputs.
  Filter to 5 blood-marker UniProts.
- **Geom:** dot/bar showing `logFC` per marker, color/size by `−log10(P.Value)`,
  shape encoding `FDR < 0.10` (filled vs open). Add π-score column
  (`P.Value^abs(logFC)`) — flag if `pi_score < 0.05`.
- **Layout:** small table-style figure, 5 rows.
- **Output data:** `c_data/supp_D_blood_markers_DEP.csv`
  (`gene, uniprot_id, logFC, P.Value, FDR, pi_score, sig_FDR_10, sig_pi_05`).

## Composite (`90_stitch_F02.R`)

- Main: `cowplot::plot_grid(pA, pB, pC, pD, ncol = 2, labels = "AUTO",
  label_size = 14, align = "hv")` → ggsave PDF (180×180 mm) + PNG (300 dpi).
- Supp: same grid, labels A–D → `SUPP_F00_blood_contamination.{pdf,png}`.
- Both write to `b_reports/main/` and `b_reports/supp/` respectively.

## Data dictionary

`c_data/DATA_DICTIONARY.md` — one entry per CSV listing:
- producing script,
- columns + types,
- N rows expected.

## Style + conventions

- Source `04_Figures/build_data_index.R` at top of every script for `sam_idx`.
- Source `04_Figures/shared/style.R` for `GROUP_COLORS`, `DIR_COLORS`,
  sizing constants.
- `_shared_keys.R` centralizes ID normalization.
- snake_case for variables; concise annotations (no AI artifacts).
- Each panel script is sourceable and idempotent — re-running overwrites
  `c_data/` + `b_reports/png/panels/` outputs deterministically.

## Acceptance criteria

- [ ] All 8 panels render without warning (4 main + 4 supp).
- [ ] All 9 `c_data/*.csv` artifacts produced + listed in `DATA_DICTIONARY.md`
      (panel A/B/C, panel D × 2, supp A/B/C/D).
- [ ] `MAIN_F00_blood_contamination.pdf` and `SUPP_F00_blood_contamination.pdf`
      exist in correct `b_reports/{main,supp}/pdf/`.
- [ ] Panel C join hits ≥ 30 samples; report `n` in caption.
- [ ] Supp B cutoff sweep shows monotonic or near-monotonic DEP trend
      (document if non-monotonic).
- [ ] Sample-ID normalization tested against both formats (CR6_T1 ↔ CR006_T1).

## Risks + mitigations

| Risk | Mitigation |
|---|---|
| Supp B refit fails if Sam's DAList lacks `$design` slot or contrast spec | Inspect at scaffold; if missing, rebuild from `metadata` columns (`supp`, `cancer`, `timepoint`) using the same formula as `01_run_dep.R` |
| Panel D(ii) Venn looks empty with 5 points | Fall back to 2-column indicator chart; decide at implementation |
| Panel C n drops below 30 if outliers + Sam's exclusions overlap heavily | Caption explicitly notes joined n; if < 25, demote panel to supp and promote one supp panel up |
| `consensus_outlier == TRUE` samples don't appear in Sam's cohort | Possible — Sam excluded different samples. Panel C will mark them as "ours-only" via shape or annotate as missing from join |

## Out of scope

- Cross-tissue (serum) integration with the Box-hosted serum DIA-NN data.
- Pipeline-level changes to the contamination filter in our main pipeline.
- Touching `A_CvH_2026/04_Figures/F02/` (CR/CRvH split deviation, deferred Phase A).
