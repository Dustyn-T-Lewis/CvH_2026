# CvH F1 Figure Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Adapt YvO F1 proteomics overview figure (8 panels + 1 supplementary) for the CvH cancer recovery dataset.

**Architecture:** Mirror YvO's `04_Figures/` structure with `shared/` utilities and per-figure directories. Each panel is a standalone R script sourcing a figure-specific `style.R`. Panels F→G share state (sig_sets). Panel H uses `pathway_utils.R` from YvO (copied).

**Tech Stack:** R (ggplot2, dplyr, tidyr, patchwork, vegan, fgsea, msigdbr, ggrepel, cowplot)

**Source reference:** All panels adapt from `/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_YvO_2025/04_Figures/F1/a_script/`

---

### Task 1: Create directory structure and shared utilities

**Files:**
- Create: `04_Figures/shared/style.R`
- Create: `04_Figures/F1/a_script/style.R`
- Copy: `04_Figures/shared/pathway_utils.R` (from YvO)
- Create: `04_Figures/F1/b_reports/` (empty)
- Create: `04_Figures/F1/c_data/` (empty)

**Step 1: Create directories**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
mkdir -p 04_Figures/shared 04_Figures/F1/a_script 04_Figures/F1/b_reports 04_Figures/F1/c_data
```

**Step 2: Create `04_Figures/shared/style.R`**

Adapt from YvO `shared/style.R`. Key changes:
- Replace `AGE_COLORS` with CvH group colors from `00_shared/helpers.R` (`PAL_GT`)
- Replace `GROUP_FILL` with CvH group fills (5 groups instead of 4)
- Keep `DB_COLORS`, `FIG_THEME`, all utility functions unchanged
- Add `DIR_COLORS` (Up/Down/NS)

```r
# CvH Shared Figure Style
# Adapted from YvO shared/style.R

library(ggplot2)

# ── Palettes ──
GROUP_COLORS <- c(
  CRE_T1 = "#2166AC", CRE_T2 = "#67A9CF",
  PLA_T1 = "#D6604D", PLA_T2 = "#F4A582",
  H_T1   = "#4DAF4A")

DIR_COLORS <- c(Up = "#D6604D", Down = "#4393C3", NS = "grey70")

GROUP_FILL <- c(
  CRE_T1 = scales::alpha("#2166AC", 0.7),
  CRE_T2 = scales::alpha("#67A9CF", 0.7),
  PLA_T1 = scales::alpha("#D6604D", 0.7),
  PLA_T2 = scales::alpha("#F4A582", 0.7),
  H_T1   = scales::alpha("#4DAF4A", 0.7))

DB_COLORS <- c(
  Hallmark    = "#E41A1C", KEGG      = "#377EB8",
  Reactome    = "#4DAF4A", WikiPathways = "#984EA3",
  `GO:BP`     = "#FF7F00", BioCarta  = "#A65628",
  PID         = "#F781BF")

# ── Sizing ──
PANEL_MD   <- 180
BASE_GENE  <- 3.2
BASE_STAT  <- 3.5

scale_text <- function(base, panel_w, ref = PANEL_MD)
  base * sqrt(ref / panel_w)

# ── Theme ──
FIG_THEME <- theme_bw(base_size = 10) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(face = "bold.italic", size = 9, colour = "grey30"),
    strip.text       = element_text(face = "bold", size = 10),
    axis.title       = element_text(face = "bold", size = 10),
    axis.text        = element_text(size = 8.5),
    legend.text      = element_text(size = 8.5),
    legend.key.size  = unit(3, "mm"),
    panel.grid.minor = element_blank())

# ── Utility functions ──
get_pdf_device <- function() {
  if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
}

fmt_p <- function(p) {
  ifelse(p < 0.001, "< 0.001",
    ifelse(p < 0.01, sprintf("= %.3f", p),
      sprintf("= %.2f", p)))
}

sig_stars <- function(padj) {
  ifelse(padj < 0.001, "***",
    ifelse(padj < 0.01, "**",
      ifelse(padj < 0.05, "*", "ns")))
}

reorder_within <- function(x, by, within, fun = mean, sep = "___") {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}

darken_color <- function(col, factor = 0.7) {
  r <- grDevices::col2rgb(col)
  grDevices::rgb(r[1]*factor, r[2]*factor, r[3]*factor, maxColorValue = 255)
}

fisher_z_ci <- function(r, n, level = 0.95) {
  z  <- atanh(r)
  se <- 1 / sqrt(n - 3)
  q  <- qnorm((1 + level) / 2)
  lo <- tanh(z - q * se)
  hi <- tanh(z + q * se)
  c(lo = lo, hi = hi)
}
```

**Step 3: Create `04_Figures/F1/a_script/style.R`**

Adapt from YvO `F1/a_script/style.R`. Key changes:
- Replace `CONTRAST_COLORS` with 5 CvH contrasts
- Replace `PCA_COLORS` and `PCA_SHAPES` for 5 groups
- Update contrast label maps

```r
# CvH F1 Style — extends shared/style.R

source(file.path(rprojroot::find_rstudio_root_file(),
                 "04_Figures/shared/style.R"))

# ── Contrast palette (5 contrasts) ──
CONTRAST_COLORS <- c(
  Cancer_vs_Healthy      = "#4CAF50",
  Training_CR            = "#9C27B0",
  Training_CRE           = "#2166AC",
  Training_PLA           = "#D6604D",
  Supplement_Interaction = "#FF8F00")

# ── PCA palette (5 groups) ──
PCA_COLORS <- c(
  CRE_T1 = "#2166AC", CRE_T2 = "#67A9CF",
  PLA_T1 = "#D6604D", PLA_T2 = "#F4A582",
  H_T1   = "#4DAF4A")

PCA_SHAPES <- c(
  CRE_T1 = 16, CRE_T2 = 17,
  PLA_T1 = 16, PLA_T2 = 17,
  H_T1   = 15)

# ── Contrast labels ──
CTR_SHORT <- c(
  Cancer_vs_Healthy      = "CR vs H",
  Training_CR            = "Tr.(CR)",
  Training_CRE           = "Tr.(CRE)",
  Training_PLA           = "Tr.(PLA)",
  Supplement_Interaction = "CRE-PLA")

CTR_FACET <- CTR_SHORT
CTR_AXIS  <- CTR_SHORT

# ── Supplement group labels ──
SUPP_LABELS <- c(CRE = "Creatine", PLA = "Placebo", H = "Healthy")

# ── F1 sizing ──
BASE_COUNT <- 4.0
BASE_GENE  <- 3.8
BASE_STAT  <- 4.0

FIG_TITLE_SIZE    <- 12
FIG_SUBTITLE_SIZE <- 9
FIG_STRIP_SIZE    <- 10
FIG_AXIS_TEXT     <- 9.5
FIG_LEGEND_TITLE  <- 10.5
FIG_LEGEND_TEXT   <- 9.5

FIG_THEME <- theme_bw(base_size = 10) +
  theme(
    plot.title       = element_text(face = "bold", size = FIG_TITLE_SIZE),
    plot.subtitle    = element_text(face = "bold.italic", size = FIG_SUBTITLE_SIZE,
                                    colour = "grey30"),
    strip.text       = element_text(face = "bold", size = FIG_STRIP_SIZE),
    axis.title       = element_text(face = "bold", size = 10),
    axis.text        = element_text(size = FIG_AXIS_TEXT),
    legend.title     = element_text(face = "bold", size = FIG_LEGEND_TITLE),
    legend.text      = element_text(size = FIG_LEGEND_TEXT),
    legend.key.size  = unit(3, "mm"),
    panel.grid.minor = element_blank())
```

**Step 4: Copy pathway_utils.R from YvO**

```bash
cp /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_YvO_2025/04_Figures/shared/pathway_utils.R \
   /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026/04_Figures/shared/pathway_utils.R
```

**Step 5: Verify by sourcing style.R in R**

```bash
Rscript -e 'setwd("/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026"); source("04_Figures/F1/a_script/style.R"); cat("OK\n")'
```

Expected: `OK` with no errors.

---

### Task 2: Panel A — CV% Violins (3 facets: CRE, PLA, Healthy)

**Files:**
- Create: `04_Figures/F1/a_script/panel_A.R`
- Reference: YvO `04_Figures/F1/a_script/panel_A.R`

**Adaptation from YvO:**
- YvO: 2 facets (Young, Old), each with Pre/Post violins
- CvH: 3 facets (Creatine, Placebo, Healthy). CRE/PLA have T1/T2 violins with paired Wilcoxon. Healthy has T1-only single violin as baseline reference.
- Parse sample IDs: `CR*_T1/T2` → CR group, `PPS*` → Healthy
- Compute CV% on linear scale (2^normalized) per protein within each sample group
- Bootstrap 95% CI on median CV, Wilcoxon + Cliff's delta for CRE/PLA T1 vs T2

**Step 1: Write panel_A.R**

Read the YvO version first (`/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_YvO_2025/04_Figures/F1/a_script/panel_A.R`), then adapt:

- Replace age/time parsing with supplement/timepoint parsing from CvH metadata
- 3 facets: `facet_wrap(~ supplement, scales = "free_x")`
- CRE and PLA facets show T1 vs T2 (paired Wilcoxon)
- Healthy facet shows single T1 violin
- Keep: bootstrap CI, Cliff's delta, audit CSV exports
- Output: `04_Figures/F1/b_reports/panel_A_cv.pdf` (160×120 mm)

**Step 2: Run panel_A.R**

```bash
Rscript /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026/04_Figures/F1/a_script/panel_A.R
```

Expected: PDF + PNG in `b_reports/`, audit CSVs in `c_data/`.

---

### Task 3: Panel B — CV Scatter Triptych (CRE, PLA, ΔCV)

**Files:**
- Create: `04_Figures/F1/a_script/panel_B.R`
- Reference: YvO `04_Figures/F1/a_script/panel_B.R`

**Adaptation from YvO:**
- YvO: Young Pre vs Post, Old Pre vs Post, ΔCV Young vs Old
- CvH: CRE Pre vs Post, PLA Pre vs Post, ΔCV CRE vs PLA
- PPS excluded (no paired data)
- Keep: Pearson correlation, top 15 labels, 98th percentile cap, audit CSVs
- Output: `panel_B_cv_scatter.pdf` (300×120 mm)

**Step 1: Write panel_B.R** — Read YvO version, adapt group parsing.

**Step 2: Run and verify output.**

---

### Task 4: Panel C — Intra-Individual Variability

**Files:**
- Create: `04_Figures/F1/a_script/panel_C.R`
- Reference: YvO `04_Figures/F1/a_script/panel_C.R`

**Adaptation from YvO:**
- YvO: log2FC per subject, faceted by Young/Old, Wilcoxon Young vs Old
- CvH: log2FC per subject (T2-T1), faceted by CRE/PLA, Wilcoxon CRE vs PLA
- PPS excluded (no T2)
- Keep: ordered by median logFC within facets, audit CSVs
- Output: `panel_C_intra_variability.pdf` (160×90 mm)

**Step 1: Write panel_C.R** — Read YvO version, adapt group parsing.

**Step 2: Run and verify output.**

---

### Task 5: Panel D — logFC Density Histograms (5 contrasts)

**Files:**
- Create: `04_Figures/F1/a_script/panel_D.R`
- Reference: YvO `04_Figures/F1/a_script/panel_D.R`

**Adaptation from YvO:**
- YvO: 3 contrasts (Aging, Training_Young, Training_Old) from single combined results file
- CvH: 5 contrasts from TWO results files (CRvH: Cancer_vs_Healthy + Training_CR; CR: Training_CRE + Training_PLA + Supplement_Interaction)
- Read both `03_combined_results_CRvH.csv` and `03_combined_results_CR.csv`
- Extract `logFC_<contrast>` columns from each
- Stack into long format with contrast as factor
- 5 vertically stacked density histograms
- Keep: median |logFC| with bootstrap CI, n(>0.5), KS/Fligner tests
- Output: `panel_D_logfc_density.pdf` (140×200 mm — taller for 5 panels)

**Step 1: Write panel_D.R** — Read YvO version, adapt to read two result files and handle 5 contrasts.

**Step 2: Run and verify output.**

---

### Task 6: Panel E — PCA Biplot + PERMANOVA

**Files:**
- Create: `04_Figures/F1/a_script/panel_E.R`
- Reference: YvO `04_Figures/F1/a_script/panel_E.R`

**Adaptation from YvO:**
- YvO: 4 groups (Young_Pre/Post, Old_Pre/Post), PERMANOVA: age × time × interaction
- CvH: 5 groups (CRE_T1, CRE_T2, PLA_T1, PLA_T2, H_T1), use `PCA_COLORS` and `PCA_SHAPES` from style.R
- PERMANOVA terms: Group (CR vs H) + Timepoint + Supplement (CRE vs PLA)
- Note: PPS has no T2, so PERMANOVA design needs careful handling — use Group_Time as single factor
- Alternative: use `adonis2(dist ~ Group_Time, ...)` as one-way, then report pairwise
- Keep: 80% confidence ellipses, bootstrap PC variance CIs, betadisper, audit CSVs
- Output: `panel_E_pca.pdf` (145×100 mm)

**Step 1: Write panel_E.R** — Read YvO version, adapt for 5 groups and unbalanced design.

**Step 2: Run and verify output.**

---

### Task 7: Panel F — DEP Counts (pseudo-log stacked bars, 5 contrasts)

**Files:**
- Create: `04_Figures/F1/a_script/panel_F.R`
- Reference: YvO `04_Figures/F1/a_script/panel_F.R`

**Adaptation from YvO:**
- YvO: 3 contrasts, reads single combined results, FDR 0.05 tier
- CvH: 5 contrasts from TWO result files, FDR 0.10 tier (exploratory)
- Read both combined results, extract per-contrast: logFC, P.Value, adj.P.Val, pi_score, sig_pi
- Build sig_sets, dir_map, all_genes for Panel G handoff
- Tiers: p<0.05 (alpha 0.25), FDR<0.10 (alpha 0.55), Pi<0.05 (alpha 1.0)
- Output: `panel_F_dep_counts.pdf` (200×70 mm — wider for 5 contrasts)
- IMPORTANT: Must populate `sig_sets`, `dir_map`, `all_genes`, `SET_LABELS`, `SET_DISPLAY_COLORS`, `pi_total`, `fdr_total` in the calling environment for Panel G

**Step 1: Write panel_F.R** — Read YvO version, adapt for 5 contrasts and FDR 0.10.

**Step 2: Run and verify output.**

---

### Task 8: Panel G — UpSet Overlap Plot

**Files:**
- Create: `04_Figures/F1/a_script/panel_G.R`
- Reference: YvO `04_Figures/F1/a_script/panel_G.R`

**Adaptation from YvO:**
- YvO: 3 contrast sets
- CvH: 5 contrast sets (more complex intersection matrix)
- Depends on Panel F for `sig_sets`, `dir_map`, `all_genes`, etc.
- Keep: dual bar + dot matrix, up/down stratification, Fisher's exact enrichment
- Output: `panel_G_upset.pdf` (200×140 mm — wider/taller for 5 sets)

**Step 1: Write panel_G.R** — Read YvO version, adapt dimensions. Logic is mostly generic (operates on sig_sets).

**Step 2: Source panel_F.R then panel_G.R together to verify.**

---

### Task 9: Panel H — fGSEA Grouped Bar Chart

**Files:**
- Create: `04_Figures/F1/a_script/panel_H.R`
- Reference: YvO `04_Figures/F1/a_script/panel_H.R`

**Adaptation from YvO:**
- YvO: 5 contrasts (Aging, Training_Young/Old, Interaction, Reversal), reads single combined results
- CvH: 5 contrasts from two models, needs to read both combined results and merge t-statistics
- Ranking metric: moderated t-statistic (`t_<contrast>` columns)
- Databases: Hallmark, KEGG, Reactome, GO:BP (same as YvO)
- Uses `run_fgsea_perdb()` from `pathway_utils.R`
- Cache to `04_Figures/F1/c_data/fgsea_tstat_all.csv` (no F2/F3 yet)
- Output: `panel_H_fgsea.pdf` (160 mm width, variable height)

**Step 1: Write panel_H.R** — Read YvO version, adapt to merge t-stats from two model files.

**Step 2: Run and verify output.**

---

### Task 10: Supplementary S1 — Pi-score Distributions

**Files:**
- Create: `04_Figures/F1/a_script/supp_S1.R`
- Reference: YvO `04_Figures/F1/a_script/supp_S1_7.R`

**Adaptation from YvO:**
- YvO: 4 contrasts in 2×2 layout
- CvH: 5 contrasts — use 3×2 layout (6 panels, one empty or use for legend)
- Read both combined results, extract `pi_score_<contrast>` columns
- Histogram + ranked scatter per contrast
- Output: `supp_S1_pi_scores.pdf`

**Step 1: Write supp_S1.R** — Read YvO version, adapt layout for 5 contrasts.

**Step 2: Run and verify output.**

---

### Task 11: Final verification and cleanup

**Step 1: Run all panels sequentially**

```bash
cd /Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026
for script in 04_Figures/F1/a_script/panel_{A,B,C,D,E,F}.R; do
  echo "Running $script..."
  Rscript "$script" || echo "FAILED: $script"
done
# F and G must run together (F populates sig_sets for G)
Rscript -e 'setwd("/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026"); source("04_Figures/F1/a_script/panel_F.R"); source("04_Figures/F1/a_script/panel_G.R")'
Rscript 04_Figures/F1/a_script/panel_H.R
Rscript 04_Figures/F1/a_script/supp_S1.R
```

**Step 2: Verify all outputs exist**

```bash
ls -la 04_Figures/F1/b_reports/*.pdf
ls -la 04_Figures/F1/c_data/*.csv
```

Expected: 8 panel PDFs + supplementary PDF, multiple audit CSVs.
