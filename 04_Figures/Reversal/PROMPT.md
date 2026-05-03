# CvH Cancer Recovery Reversal Figure — Build Prompt

## What to build

A 5-panel "Cancer Recovery Reversal" figure asking: **does exercise training in cancer recovery patients reverse the cancer-vs-healthy proteome signature?** Plus 6 supplementary diagnostic panels.

This recapitulates the YvO F05 "Aging Reversal" figure from `/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_YvO_2025/04_Figures/F05/` but for CvH data with different contrasts. Use the YvO F05 scripts as **reference architecture** — read them to understand the logic, then write new self-contained CvH scripts. Do NOT copy YvO code verbatim or import YvO utilities.

## Project root

`/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026`

All scripts use `setwd(rprojroot::find_rstudio_root_file())` and paths relative to this root.

## Contrast mapping (YvO → CvH)

| Role | YvO contrast | CvH contrast | CvH column suffix |
|------|-------------|--------------|-------------------|
| Disease/baseline signature (X-axis) | `Aging` (Old_Pre − Young_Pre) | `Cancer_vs_Healthy` ((CRE_T1+PLA_T1)/2 − H_T1) | `_Cancer_vs_Healthy` |
| Training recovery (Y-axis) | `Training_Old` (Old_Post − Old_Pre) | `Training_CR` ((CRE_T2+PLA_T2)/2 − (CRE_T1+PLA_T1)/2) | `_Training_CR` |

Reversal = opposite-sign logFC between Cancer_vs_Healthy and Training_CR.

**Circularity**: Cancer_vs_Healthy and Training_CR share CRE_T1/PLA_T1 baselines with opposite signs → structural negative correlation (same issue as YvO F05). This must be acknowledged and tested in supplementary.

## Data sources (all verified to exist)

| Data | Path | Notes |
|------|------|-------|
| DEP results | `03_DEP/c_data/03_combined_results_CRvH.csv` | Wide format. Key columns: `logFC_Cancer_vs_Healthy`, `logFC_Training_CR`, `t_Cancer_vs_Healthy`, `t_Training_CR`, `pi_score_Cancer_vs_Healthy`, `pi_score_Training_CR`, `P.Value_*`, `adj.P.Val_*`, `gene` |
| Imputation flags | `02_Imputation/c_data/02_mar_mnar_classification.csv` | Column `classification`; "Complete" = not imputed |
| fGSEA cache | `04_Figures/shared/fgsea_CRvH.csv` | **LONG format** (4,700 rows). Columns: `pathway`, `pval`, `padj`, `NES`, `size`, `leadingEdge`, `database`, `contrast`. Contrast values: `"Cancer_vs_Healthy"`, `"Training_CR"`. Must pivot wider for Panel B scatter. |
| Imputed DAList | `02_Imputation/c_data/01_DAList_imputed.rds` | For fry design matrix + expression matrix |
| Normalized DAList | `01_normalization/c_data/03_DAList_normalized.rds` | Metadata: `dal$metadata` has `Col_ID`, `Group_Time`, `Subject_ID`, `Timepoint`, `Supplement` |

## Existing CvH shared utilities (USE THESE — do not rewrite)

- `04_Figures/shared/style.R` — palettes, theme, helpers. **Already has**:
  - `classify_proteins_f4(pi_CvH, pi_TR)` → "Sig Both" / "Sig Cancer only" / "Sig Training only" / "NS"
  - `SIG_COLORS_F4`, `SIG_LABEL_FILL_F4`, `SIG_LABEL_TEXT_F4`
  - `ORA_QUAD_COLORS_F4` = "Reversed (Cancer Up)", "Reversed (Cancer Down)", "Exacerbated Up", "Exacerbated Down"
  - `CONTRAST_COLORS` (Cancer_vs_Healthy = green, Training_CR = purple)
  - `FIG_THEME`, `get_pdf_device()`, `clean_pathway_name()`, `fisher_z_ci()`, `sig_stars()`, `make_sigmoid_ribbon()`, `scale_text()`
  - Text sizes: `FIG_TITLE_SIZE=12`, `FIG_SUBTITLE_SIZE=9`, `FIG_AXIS_TEXT=8.5`, `PANEL_MD=180`
  - **No PRINT_SCALE** — CvH uses direct sizing, not print-scale compression
- `04_Figures/shared/pathway_utils.R` — **Full dedup pipeline**:
  - `build_pathway_collection(min_size, max_size, include_goslim, exclude_variants)` — MSigDB Hallmark + KEGG Medicus + Reactome + GO:BP, with disease/cancer term filtering
  - `run_ora_deduplicated(genes, universe, pathways, jaccard_cutoff=0.5)` — per-database `fgsea::fora()` → database-stratified Jaccard dedup (Reimand et al. 2019). Returns only significant deduplicated terms.
  - `run_fgsea_deduplicated(ranks, pathways, jaccard_cutoff=0.5)` — `fgseaMultilevel` → stratified Jaccard dedup on significant results
  - `run_enrichment_pipeline(stats_list, pw_list, jaccard_cutoff=0.35)` — multi-contrast fGSEA with `collapsePathways` + Jaccard dedup
  - `deduplicate_enrichment(results, pathways, jaccard_cutoff)` — two-pass: within-database dedup first, then cross-database interleaving
  - `deduplicate_enrichment_flat(results, pathways, jaccard_cutoff)` — greedy single-pass Jaccard dedup
  - `classify_database()`, `classify_pathway_func()`, `CONSOLIDATED_PATHWAY_ORDER`, `CONSOLIDATED_COLORS`
- `04_Figures/shared/go_slim_categories.R` — `assign_go_slim_consolidated()` maps genes to 15 consolidated GO Slim BP categories via GO.db GOBPOFFSPRING hierarchy

## Target directory (already created)

```
04_Figures/Reversal/
├── a_script/
│   ├── 90_stitch_figure.R              # master orchestrator
│   ├── main/
│   │   ├── 90_stitch_main.R            # 5-panel composite + xlsx workbook
│   │   └── panels/
│   │       ├── panel_A_ORA.R           # quadrant ORA scatter + flanking bars (~400 lines)
│   │       ├── panel_B_nes_scatter.R   # pathway NES scatter (~280 lines)
│   │       ├── panel_C_pattern_heatmap.R  # pattern heatmap + sankey (~400 lines)
│   │       ├── panel_D_fry.R           # fry rotation test (~500 lines)
│   │       └── panel_E_rrho2.R         # RRHO2 heatmap (~380 lines)
│   └── supp/
│       ├── 90_stitch_supp.R            # 6-panel supp composite
│       └── panels/
│           ├── enrichment_heatmap.R     # pathway enrichment heatmap by reversal pattern
│           ├── SUPP_fry_circularity.R  # permutation test for circularity bias
│           ├── SUPP_fry_leading.R      # top 25 fry driving proteins dotplot
│           ├── SUPP_goslim_bars.R      # GO Slim distribution by quadrant
│           ├── SUPP_ora_dedup.R        # ORA sensitivity across Jaccard cutoffs
│           └── SUPP_r_bootstrap.R      # Pearson r bootstrap CI
├── b_reports/main/{pdf,png}/panels/    # individual panel outputs
├── b_reports/supp/{pdf,png}/panels/    # supp panel outputs
└── c_data/                             # intermediate CSVs → xlsx workbook
```

## Panel-by-panel specification

### Panel A — Quadrant ORA Scatter (self-contained, ~400 lines)

**Data loading**:
```r
dep_df <- read_csv("03_DEP/c_data/03_combined_results_CRvH.csv", show_col_types = FALSE)
imputation_df <- read_csv("02_Imputation/c_data/02_mar_mnar_classification.csv", show_col_types = FALSE) %>%
  transmute(gene, imputed = classification != "Complete")
```

**Scatter data frame**:
```r
scatter_df <- dep_df %>%
  transmute(gene,
            logFC_CvH = logFC_Cancer_vs_Healthy,
            logFC_TR  = logFC_Training_CR,
            pi_CvH    = pi_score_Cancer_vs_Healthy,
            pi_TR     = pi_score_Training_CR) %>%
  filter(!is.na(logFC_CvH), !is.na(logFC_TR)) %>%
  left_join(imputation_df, by = "gene") %>%
  mutate(
    imputed   = replace_na(imputed, FALSE),
    sig_class = classify_proteins_f4(pi_CvH, pi_TR),
    is_sig    = sig_class != "NS",
    quadrant  = case_when(
      logFC_CvH > 0 & logFC_TR < 0 ~ "Reversed (Cancer Up / Training Down)",
      logFC_CvH < 0 & logFC_TR > 0 ~ "Reversed (Cancer Down / Training Up)",
      logFC_CvH > 0 & logFC_TR > 0 ~ "Exacerbated Up",
      TRUE                          ~ "Exacerbated Down"))
```

**ORA per quadrant**: `build_pathway_collection(min_size=15, max_size=500, include_goslim=FALSE, exclude_variants=TRUE)` then `run_ora_deduplicated()` per quadrant. Top 5 per quadrant.

**Scatter plot**:
- Quadrant background: blue tint for reversed (off-diagonal), red tint for exacerbated (diagonal)
- Reference line: slope = −1 (reversal diagonal), dashed
- NS points: grey80, small, low alpha
- Sig points: colored by sig_class, shape 21 (filled circle), imputed = black border
- Gene labels: top 5 per sig_class ranked by `|logFC_CvH| + |logFC_TR|`, `ggrepel::geom_label_repel()`
- Corner labels: quadrant name + sig_count/total_count
- Axis titles inside plot: `log₂FC (Cancer vs Healthy)` and `log₂FC (Training CR)`
- Custom significance key below scatter

**Flanking half-bars**: 4 panels (UL, LL, UR, LR) showing top 5 ORA terms per quadrant as horizontal bars, `−log₁₀(p_adj)` axis, asterisks for significance, pathway labels inside or outside bars depending on bar width.

**Composite**: patchwork `area()` layout — bars flanking scatter, key below. Save PNG + PDF.

**Exports**: `c_data/panel_A/ora_quadrant.csv`

### Panel B — Pathway NES Scatter (~280 lines)

**fGSEA cache loading** — critical format conversion:
```r
fgsea_long <- read_csv("04_Figures/shared/fgsea_CRvH.csv", show_col_types = FALSE)
# Pivot to wide format for scatter
fgsea_wide <- fgsea_long %>%
  select(pathway, NES, padj, database, contrast) %>%
  pivot_wider(names_from = contrast,
              values_from = c(NES, padj),
              names_glue = "{.value}_{contrast}") %>%
  filter(!is.na(NES_Cancer_vs_Healthy), !is.na(NES_Training_CR))
```

Then filter to Hallmark + GO Slim pathways (or all databases — match what makes biological sense for a scatter). Compute Spearman ρ with Fisher z CI. Plot NES_Cancer_vs_Healthy vs NES_Training_CR. Reference slope = −1. Quadrant backgrounds: off-diagonal blue (reversed), diagonal red (exacerbated). Top 12 pathways labeled. `clean_pathway_name()` for display.

**Exports**: `c_data/panel_B/nes_scatter.csv`

### Panel C — Pattern Heatmap (~400 lines)

**Filter**: proteins with `pi_score_Cancer_vs_Healthy < 0.05 | pi_score_Training_CR < 0.05`

**Classification**:
```r
quadrant = case_when(
  logFC_Cancer_vs_Healthy > 0 & logFC_Training_CR < 0 ~ "Reversed Up",
  logFC_Cancer_vs_Healthy < 0 & logFC_Training_CR > 0 ~ "Reversed Down",
  TRUE ~ "Non-reversed")
sig_cat = case_when(
  pi_score_Cancer_vs_Healthy < 0.05 & pi_score_Training_CR < 0.05 ~ "Both",
  pi_score_Cancer_vs_Healthy < 0.05  ~ "Cancer",
  pi_score_Training_CR < 0.05        ~ "Tr.(CR)",
  TRUE ~ "NS")
```

**Layout** (Y-axis = proteins, sorted by logFC_Cancer_vs_Healthy within quadrant):
- Left: quadrant color strip (Reversed Up = red, Reversed Down = blue, Non-reversed = green)
- Center: 2-column heatmap (logFC_CvH, logFC_TR) with diverging blue-white-red scale
- Right: significance category strip
- Far right: GO Slim consolidated category bars (stacked) + Sankey ribbons from quadrant → category

Column headers: "Cancer" and "Tr.(CR)" colored by `CONTRAST_COLORS`

**Exports**: `c_data/panel_C_heatmap/pattern_classification.csv`, `sankey_links.csv`, `bar_data.csv`

### Panel D — fry Rotation Test (~500 lines)

**Purpose**: test whether cancer-significant protein sets collectively reverse with training.

**Design matrix construction** — must subset to CR subjects only:
```r
dal <- readRDS("02_Imputation/c_data/01_DAList_imputed.rds")
mat <- dal$data
meta <- as.data.frame(dal$metadata)
# Subset to CR subjects (exclude H_T1 healthy controls)
cr_idx <- which(meta$Group_Time != "H_T1")
mat_cr <- mat[, cr_idx]
meta_cr <- meta[cr_idx, ]
meta_cr$Group_Time <- factor(meta_cr$Group_Time,
                              levels = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2"))
design <- model.matrix(~ 0 + Group_Time, data = meta_cr)
colnames(design) <- gsub("^Group_Time", "", colnames(design))
```

**Within-subject correlation**: `limma::duplicateCorrelation(mat_cr, design, block = meta_cr$Subject_ID)`

**Gene sets**: cancer-significant DEPs from `03_combined_results_CRvH.csv` at Pi < 0.05:
- `cancer_up`: genes with `pi_score_Cancer_vs_Healthy < 0.05 & logFC_Cancer_vs_Healthy > 0`
- `cancer_down`: genes with `pi_score_Cancer_vs_Healthy < 0.05 & logFC_Cancer_vs_Healthy < 0`

**Contrast for testing**: `Training_CR = (CRE_T2 + PLA_T2)/2 - (CRE_T1 + PLA_T1)/2`

**fry call**: `limma::fry(mat_cr, index = list(cancer_up_idx, cancer_down_idx), design = design, contrast = Training_CR_contrast, block = meta_cr$Subject_ID, correlation = dupCor$consensus)`

**Expected reversal**: cancer-up set → negative direction in Training_CR; cancer-down set → positive direction.

**Circularity**: `has_circularity = TRUE`. Report `r = cor(t_Cancer_vs_Healthy, t_Training_CR)` in subtitle.

**Visualization**: Running enrichment score curves (ES, like GSEA barcode plots) for both gene sets against Training_CR ranked t-statistics. Barcode tick marks at gene set positions. Flanking ORA bars for driving proteins (those in set whose t-stat sign matches expected reversal direction).

**Exports**: `c_data/panel_D_fry/fry_results_all.csv`, `driving_proteins.csv`

### Panel E — RRHO2 Heatmap (~380 lines)

**Rank lists**: proteins ranked by `t_Cancer_vs_Healthy` and `t_Training_CR` (descending).

**RRHO2 implementation**: pure-R `phyper()` stratified hypergeometric test (Cahill et al. 2018). NOT the RedRibbon package (it segfaults with ~2000+ genes). Build step-size based on `floor(n/100)` or similar. Compute −log₁₀(p) in each quadrant of the rank-rank grid.

**4-quadrant heatmap**:
- UU = "Exacerbated Up" (both up)
- DD = "Exacerbated Down" (both down)
- UD = "Reversed (Cancer↑ Tr↓)"
- DU = "Reversed (Cancer↓ Tr↑)"
- JET colormap, −log₁₀(p) intensity

**Hotspot gene extraction**: genes in the max-signal region of each quadrant.

**Per-quadrant ORA**: `run_ora_deduplicated()` on hotspot genes per quadrant.

**Supplementary ORA barplot**: horizontal bars per quadrant, top enriched terms.

**Exports**: `c_data/panel_E/rrho2_summary.csv`, `rrho2_hotspot_genes.csv`, `rrho2_ora_concordant.csv`, `rrho2_ora_discordant.csv`

### 90_stitch_main.R — Composite + Workbook

**Source order**: A → B → D → E → C (C last — loads AnnotationDbi which masks dplyr::select; also preserves stat snapshots from earlier panels).

**Layout**: 3-column patchwork grid matching YvO F05:
- Top row: A (ORA, 8 cols) + B (heatmap, 6 cols) — B is the pattern heatmap (relabeled from C)
- Bottom row: C (fry, 6 cols) + D (NES scatter, 4 cols) + E (RRHO2, 4 cols)
- Note: visual panel labels (A-E) don't match script names — the stitcher maps them

**Panel tags/titles/subtitles**: cowplot `draw_label()` overlay with stat summaries captured between panel sourcing.

**Excel workbook**: `build_workbook()` from `04_Figures/shared/figure_supplement_helpers.R` (you'll need to write a minimal version or source from shared if it exists in CvH — check first). Sheets: `panel_A_ora_quadrant`, `panel_B_pattern_class`, `panel_B_sankey`, `panel_B_bar`, `panel_C_fry_results`, `panel_C_fry_driving`, `panel_D_nes_scatter`, `panel_E_rrho2_summary`, `panel_E_rrho2_hotspot`, `panel_E_rrho2_ora_*`, plus all SUPP sheets.

**Cleanup**: remove consumed CSVs from `c_data/` subdirectories after workbook build.

**Output**: `b_reports/main/pdf/MAIN_Reversal_composite.pdf`, `b_reports/main/png/MAIN_Reversal_composite.png`

### Supplementary Panels (6)

All write to `c_data/panel_supp/` and expose plot objects for the supp stitcher.

1. **SUPP_ora_dedup.R** — Run ORA at Jaccard cutoffs 0.3, 0.5, 0.7, 1.0 for reversed quadrants. Grouped bar chart showing pathway count stability across cutoffs.

2. **SUPP_r_bootstrap.R** — 1000 bootstrap replicates of Pearson r between logFC_CvH and logFC_TR. Histogram + observed r + 95% percentile CI.

3. **SUPP_fry_circularity.R** — 1000 protein-label permutations. Shuffles gene labels, recomputes correlation. Demonstrates observed r is more extreme than circularity bias alone.

4. **SUPP_reversal_threshold.R** — Line plot: % Reversed / Exacerbated / Negligible across |logFC| thresholds 0.05–0.30. Shows classification stability.

5. **SUPP_goslim_bars.R** — Stacked horizontal bars showing GO Slim category distribution by reversal quadrant (Reversed Up/Down, Non-reversed). Filter on Pi < 0.05.

6. **SUPP_fry_leading.R** — Dotplot of top 25 fry driving proteins ranked by |t_Training_CR|. Colors by cancer-set direction (cancer-up reversed = blue, cancer-down reversed = red).

### 90_stitch_supp.R — Supplementary Composite

Sources all 6 supp panel scripts, stitches into 3×2 cowplot grid with panel tags (A–F). Output: `b_reports/supp/pdf/SUPP_Reversal_diagnostics.pdf`, `b_reports/supp/png/SUPP_Reversal_diagnostics.png`

### 90_stitch_figure.R — Master Orchestrator

Sources enrichment_heatmap.R → 90_stitch_supp.R → 90_stitch_main.R. Final cleanup.

## Critical implementation notes

1. **ORA deduplication**: Use `run_ora_deduplicated()` from `pathway_utils.R` for ALL ORA calls. It runs per-database `fgsea::fora()` then database-stratified Jaccard dedup (Reimand et al. 2019). Default cutoff 0.5.

2. **fGSEA redundancy correction**: The fGSEA cache (`fgsea_CRvH.csv`) was already generated with `run_enrichment_pipeline()` which applies `collapsePathways` + Jaccard dedup. For Panel B, you're reading pre-computed NES values — no re-running needed. For any new fGSEA runs (e.g., enrichment_heatmap.R supp), use `run_fgsea_deduplicated()` which applies `fgseaMultilevel` + stratified Jaccard dedup.

3. **fGSEA cache format**: CvH's cache is LONG (one row per pathway×contrast), not WIDE. You MUST `pivot_wider()` before plotting Panel B scatter.

4. **No PRINT_SCALE**: CvH uses direct font sizing (FIG_TITLE=12pt, FIG_AXIS=8.5pt). Do not import YvO's `print_scale_380.R` or use `PRINT_SCALE` multipliers.

5. **Gene column**: Use `gene` (not `gene_symbol`) from combined results — this matches pathway gene sets.

6. **fry design**: Subset to CR subjects only (exclude H_T1) because Training_CR is defined only within CR subjects. Block on `Subject_ID`.

7. **RRHO2**: Use pure-R `phyper()`, NOT RedRibbon package. RedRibbon segfaults at ~2000+ genes.

8. **PDF device**: Always use `get_pdf_device()` from style.R (cairo_pdf → pdf fallback).

9. **figure_supplement_helpers.R**: Check if CvH has this at `04_Figures/shared/figure_supplement_helpers.R`. If not, write a minimal version with `add_sheet()`, `safe_read()`, `build_workbook()`, `cleanup_after_workbook()` — or build the xlsx inline with `openxlsx`.

10. **All outputs use `MAIN_` or `SUPP_` prefix** in filenames.

## Reference implementation

Read the YvO F05 scripts for architecture reference:
- `A_YvO_2025/04_Figures/F05/a_script/main/panels/panel_A_ORA.R` (417 lines — self-contained quadrant ORA)
- `A_YvO_2025/04_Figures/shared/comparison_panels/panel_B_nes_scatter.R` (278 lines — NES scatter template)
- `A_YvO_2025/04_Figures/shared/comparison_panels/panel_C_pattern_heatmap.R` (408 lines — heatmap template)
- `A_YvO_2025/04_Figures/shared/comparison_panels/panel_D_fry.R` (576 lines — fry template)
- `A_YvO_2025/04_Figures/shared/comparison_panels/panel_E_rrho2.R` (381 lines — RRHO2 template)
- `A_YvO_2025/04_Figures/F05/a_script/main/90_stitch_main.R` (255 lines — stitcher)

Adapt the logic — do not copy verbatim. Remap all contrast names, column names, color palettes, and text labels to CvH equivalents.

## Execution plan

Build in this order:
1. Panel A (self-contained, most straightforward)
2. Panel B (depends on fGSEA cache pivot)
3. Panel D (depends on DAList + limma fry)
4. Panel E (depends on RRHO2 implementation)
5. Panel C (depends on GO Slim — source last due to AnnotationDbi masking)
6. 90_stitch_main.R (sources all 5, builds composite + xlsx)
7. Supp panels (6 scripts)
8. 90_stitch_supp.R
9. 90_stitch_figure.R (orchestrator)

Test each panel individually before stitching. Run with `Rscript --no-init-file` (CvH has `.Rprofile` that may source nonexistent renv).
