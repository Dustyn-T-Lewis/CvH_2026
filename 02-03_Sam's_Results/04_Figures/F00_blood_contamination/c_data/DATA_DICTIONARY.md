# F00 c_data dictionary

One line per CSV. Producer script in parens. Columns + types in inline list.

## Main

- **panel_A_stacked_data.csv** (`_panel_A_blood_stacked.R`)
  `sample_id <chr>` | `Group_Time <fct>` | `B_M_ratio <dbl>` | `marker <chr>` | `pct <dbl>` — 175 rows (35 samples × 5 markers).
- **panel_B_dotplot_data.csv** (`_panel_B_bm_ratio_dotplot.R`)
  `sample_id <chr>` | `Group_Time <fct>` | `B_M_ratio <dbl>` | `top_decile_flag <lgl>` — 35 rows.
- **panel_C_join_data.csv** (`_panel_C_bm_vs_mahalanobis.R`)
  `sample_id <chr>` | `Group_Time <chr>` | `B_M_ratio <dbl>` | `mahal_dist <dbl>` | `n_flags <int>` | `consensus_outlier <lgl>` | `label_this <lgl>` — 35 rows (joined Sam ∩ ours after PPS ID normalization).
- **panel_D_venn_protein_sets.csv** (`_panel_D_filter_venns.R`)
  `uniprot_id <chr>` | `gene <chr>` | `in_sam <lgl>` | `in_ours <lgl>` | `set_membership <chr>` — 2558 rows (union of Sam-kept 1944 and our-kept 2542 UniProt IDs; 1928 both, 16 sam_only, 614 ours_only).
- **panel_D_venn_blood_markers.csv** (`_panel_D_filter_venns.R`)
  `uniprot_id <chr>` | `gene <chr>` | `in_sam <lgl>` | `in_ours <lgl>` | `in_blood_blacklist <lgl>` — 5 rows (HBB, HBA1, MB, ALB, CKM). HBB/HBA1/ALB filtered by both pipelines (in_sam = in_ours = FALSE, in_blood_blacklist = NA); MB and CKM kept by both.

## Supp

- **supp_A_intensity_pairs.csv** (`_supp_panel_A_marker_intensity.R`)
  `uniprot_id <chr>` | `gene <chr>` | `sample_id <chr>` | `intensity_sam <dbl>` | `intensity_ours <dbl>` — 70 rows (2 surviving markers MB+CKM × 35 paired samples). HBB/HBA1/ALB absent from both DALists (filtered upstream).
- **supp_B_cutoff_sensitivity.csv** (`_supp_panel_B_cutoff_sensitivity.R`)
  `cutoff_pct <int>` | `n_samples_dropped <int>` | `n_samples_kept <int>` | `n_proteins_used <int>` | `n_DEP_fdr10 <int>` | `dropped_sample_ids <chr>` — 5 rows (cutoffs 0/5/10/15/20%). Baseline 653 DEPs at cutoff 0; monotonically declines to 475 at cutoff 20.
- **supp_C_corr_matrix.csv** (`_supp_panel_C_corr_heatmap.R`)
  `var1 <chr>` | `var2 <chr>` | `rho <dbl>` | `p <dbl>` — 36 rows (long form of 6×6 Spearman matrix over HBB_pct/HBA1_pct/MB_pct/ALB_pct/CKM_pct/B_M_ratio).
- **supp_D_blood_markers_DEP.csv** (`_supp_panel_D_blood_markers_in_DEP.R`)
  `gene <chr>` | `uniprot_id <chr>` | `logFC <dbl>` | `P.Value <dbl>` | `FDR <dbl>` | `pi_score <dbl>` | `sig_FDR_10 <lgl>` | `sig_pi_05 <lgl>` — 5 rows (MB and CKM with stats, HBB/HBA1/ALB rows are NA — filtered from DEP table upstream).

## Conventions

- All `sample_id` values are short-form (`CR6_T1`, not `CR006_T1`; `PPS2`, not `PPS02`), normalized via `_shared_keys.R::normalize_sample_id`. The helper strips leading zeros from both `CR##_T#` and `PPS##` IDs.
- FDR threshold: 0.10 (BH).
- π-score threshold: 0.05 (Xiao 2014, `P.Value ^ |logFC|`).
- All correlations are Spearman unless stated otherwise.
- Cohort: 35 samples = 25 CR (cancer-recovery survivors, with timepoints T1/T2) + 10 PPS (healthy controls, no timepoint).
