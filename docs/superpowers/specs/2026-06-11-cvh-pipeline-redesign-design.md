# CvH Pipeline Redesign + Documentation — Design Spec

**Date:** 2026-06-11
**Status:** Draft for user review
**Scope:** A_CvH_2026 (primary, build first) + A_Mito_2026 (phase 2 outline)
**Author:** brainstorming session (Claude + Dustyn)

---

## 1. Purpose & goals

Rationalize the CvH (Cancer-Recovery vs Healthy skeletal-muscle proteomics) pipeline into a clean,
logical, ground-up reproducible workflow, and document it in a Quarto methods+narrative document.
Same treatment for Mito (H9c2 mito-transplant) in phase 2, plus a cell-culture white paper.

Driving requirements from the user:

- **Ground-up per-pipeline Quarto** explaining design rationale and every stage (code shown, frozen
  outputs embedded).
- **Trim/validate scripts** so they are logical and human-readable; use proteoDA cleanly.
- **Blood-contaminant removal** in CvH must be **systematic, unbiased, literature-backed, HPA-derived**
  — not a hand-typed gene list.
- **Clean, clear statistical models** with explicit, interpretable contrasts.
- **Clean, clear figures** consolidated to a core set.
- Mito gets a **cell-culture white paper** (goals + placeholder academic budget).

### Working method (locked)
- Brainstorm both now; **build CvH end-to-end first**, then Mito.
- Quartos: narrate + show code (`eval: false`) around **frozen** figures/tables.
- Script changes: **proposed in this spec/Quarto first**, applied only after approval, on a feature branch.
- Refactor: extract duplicated/buried logic into documented `R/` helpers **with testthat tests**
  (repo rule: 80%+ on statistical functions).

---

## 2. New directory structure (CvH)

Filtering is split from normalization; imputation is a dedicated numbered stage.

```
00_input/         raw matrix, metadata, HPA downloads (muscle + blood atlas)
01_Filtering/     HPA muscle filter -> blood C=B\M removal -> UniProt dedup -> missingness -> outlier consensus
02_Normalization/ cycloess + QC reports
03_Imputation/    MAR/MNAR consensus + missForest (imputed matrix for figures/WGCNA only)
04_DEP/           single 3-group limma fit; primary + reversal contrasts; sensitivity arms
05_Figures/       F00_QC, F01_Overview, F02_Baseline_Training, F03_Reversal, F04_WGCNA
```

- **DEP uses the non-imputed cycloess matrix** (limma handles per-protein NAs). 03_Imputation output
  feeds only figures/WGCNA. This split is preserved and documented.
- Outlier **sample** detection lives in 01_Filtering (it is sample filtering).
- Each stage keeps the `a_script/ b_reports/ c_data/` triple.

> Migration note: current repo has `01_normalization` (filter+norm combined) and `02_Imputation`.
> The redesign separates these and renumbers `03_DEP -> 04_DEP`, figures -> `05_Figures`. Old stage
> outputs are regenerated; raw inputs untouched.

---

## 3. Blood-contaminant removal (the centerpiece)

Replaces the hand-typed `BLOOD_CONTAMINANTS` vector (Geyer 2016 list) with a reproducible,
HPA-derived, set-theoretic rule. All sets regenerate from the HPA download; nothing hand-typed.

### 3.1 Three HPA-derived sets
- **M = muscle-expressed set** — genes HPA reports as expressed in skeletal muscle (nTPM above
  detection in HPA "skeletal muscle" tissue, and/or HPA tissue-specificity ∈ {tissue-enriched,
  group-enriched, tissue-enhanced} that includes skeletal muscle).
- **B = blood/plasma candidate set** — genes HPA classifies as blood-borne: protein classes
  {Plasma proteins, Immunoglobulin genes, Blood group antigen proteins, Secreted to blood}
  ∪ hemoglobin/erythrocyte genes ∪ blood-cell-enriched genes from the HPA Blood Atlas.
- **C = B \ M** — blood proteins with **no** legitimate skeletal-muscle expression. This is the
  contaminant filter list.

### 3.2 Application
`observed_proteins ∩ C` is removed. The `B \ M` construction self-corrects the hard cases:
- **Kept (B∩M):** myoglobin (MB), creatine kinase (CKM), carbonic anhydrase (CA3) — muscle proteins
  that also circulate; real signal.
- **Removed (B\M):** hemoglobins (HBA1/HBA2/HBB), albumin (ALB), immunoglobulins, complement,
  apolipoproteins — blood with no muscle expression.

### 3.3 Quantitative QC + sensitivity (so contamination is *measured*, not assumed)
- **Per-sample blood-contamination index** from canonical erythrocyte/plasma markers
  (hemoglobins ± broader plasma panel), reported in F00_QC. High-contamination samples flagged into
  the outlier-consensus step.
- **Removal-sensitivity arm** in DEP: rerun key contrasts with vs without `C` removal; confirm
  conclusions hold (mirrors existing sensitivity pattern).
- **Geyer 2016** retained only as a cross-check that `C` recovers the known plasma panel.

### 3.4 Literature anchors (verify exact PMIDs/DOIs at build)
- HPA tissue atlas: Uhlén et al. 2015, *Science*.
- HPA blood atlas / secretome: Uhlén et al. 2019.
- Plasma contamination panel cross-check: Geyer et al. 2016, PMID 27135364.

> Build-time verification: confirm HPA's actual class names and tissue-expression columns against the
> downloaded annotation file before coding the set logic.

---

## 4. Statistical model (CvH)

### 4.1 Simplified 3-group design
Collapse the supplement (CRE/PLA) distinction; pool CR. Three `group_time` levels:
`H_pre`, `CR_pre`, `CR_post`.

**Single limma fit:**
```r
dal <- add_design(dal, "~ 0 + group_time + (1|subject)")
dal <- add_contrasts(dal, c(
  "CRvH_Baseline = CR_pre  - H_pre",   # cross-sectional: cancer-recovered vs healthy at baseline
  "CR_Training   = CR_post - CR_pre",  # within-subject (paired): recovery/training effect
  "Resid         = CR_post - H_pre"    # residual deviation after training (reversal readout)
))
dal <- fit_limma_model(dal)
```

This replaces the old two-model (CRvH + CR) structure and its 6 contrasts with **one fit and two
primary contrasts** (+ a derived reversal contrast). 6-across-2-models -> 2-in-1-model.

### 4.2 How proteoDA handles the random effect (verified from source)
proteoDA's `add_design` parses `(1|subject)`, **strips it from the design matrix**, and stores
`subject` as `random_factor`. `fit_limma_model` runs `limma::duplicateCorrelation(block = subject)`,
gets one **consensus intra-block correlation** ρ̂, and fits
`lmFit(..., block = subject, correlation = ρ̂)` then `eBayes(robust = TRUE)`.

Implications (documented in the Quarto):
- It is **not** a full mixed model — it is one blocking factor + one global ρ̂ (limma duplicateCorrelation;
  Smyth, Michaud & Scott 2005).
- **Exactly one** blocking factor is supported. `subject` is the one we need (the pre/post pairing),
  so **supplement cannot also be a random effect** — hence supplement is a CR-only sensitivity covariate.
- For a **2-timepoint** design, random-intercept LMM ≡ duplicateCorrelation compound symmetry, so
  this is the appropriate (not compromised) tool for a design with **both** within-subject
  (CR_pre→CR_post) and between-subject (CR vs Healthy) contrasts. A fixed-`subject` model could not
  estimate the cross-sectional contrast.

### 4.3 Interpretation (for the Quarto)
- `~ 0 + group_time`: each coefficient is a group mean log2 abundance; contrasts are differences of means.
- ρ̂: average within-subject correlation after removing group effects; a clearly positive ρ̂ is the
  pairing showing up and is what makes `CR_Training` a proper paired test. Report ρ̂.
- eBayes shrinks per-protein variance toward a global trend (power at small n).
- Pi-score `= P.Value^|logFC|` ranks by significance and effect together (Xiao 2014).

### 4.4 eBayes trend note
proteoDA's `fit_limma_model` sets `robust = TRUE` but **not `trend = TRUE`** (verified in source),
despite the README claiming robust+trend. Reconcile during build: for label-free intensities,
`trend = TRUE` is often the better default — document the choice and, if changed, show it as a
sensitivity comparison.

### 4.5 Sensitivity arms
- **Paired-model check for Training:** explicit `~ subject + timepoint` on CR subjects only; confirm
  `CR_Training` logFCs match the duplicateCorrelation version.
- **Supplement covariate:** `~ 0 + group_time + supplement + (1|subject)` on CR subjects; confirm
  `CRvH_Baseline` unmoved by creatine.
- **Blood-removal arm:** key contrasts with vs without `C` removal.

---

## 5. Reversal / rejuvenation analysis (CvH)

Goal: does training (`CR_Training`) reverse the cancer signature (`CRvH_Baseline`)? Use **all three
tiers**, unified by a rotation null that neutralizes shared-baseline circularity.

### 5.1 Validity threat (must address): shared-baseline circularity
`D = CR_pre − H_pre` and `T = CR_post − CR_pre` share `CR_pre` with opposite signs, inducing a
**structural negative covariance** even under the null. Any naive D–T correlation / RRHO2 / NES
anti-correlation is biased toward "reversal". (Circular-analysis principle: Kriegeskorte et al. 2009.)

**Fix (the backbone):** a **sample-label rotation/permutation null** that recomputes both contrasts
(and gene-set membership) under permuted labels, so the artifact is in the null. Plus the
**`Resid = CR_post − H_pre`** framing, which shares `H_pre` at the *same* sign (no manufactured reversal)
and is the interpretable "distance to healthy" readout.

### 5.2 Three tiers
- **Tier 1 — display:** logFC–logFC quadrant scatter, NES–NES pathway scatter, RRHO2
  (Plaisier 2010; Cahill 2018). Shown in the figure; claims require the rotation null.
- **Tier 2 — primary test (sensitive):** gene-set rotation tests. Define disease-UP/DOWN sets from `D`;
  test whether `T` moves them oppositely with **`fry`/`roast`** (self-contained, directional;
  Wu et al. 2010) and **`camera`** (competitive, inter-gene-correlation corrected; Wu & Smyth 2012).
- **Tier 3 — explicit per-protein model:** decomposition `D`, `T`, `Resid = D + T`.
  - Restrict to disease-dysregulated proteins (D FDR-significant).
  - **Rejuvenation fraction φ = −T/D** (1 = normalized to healthy, 0 = untouched, <0 = exacerbated,
    >1 = overshoot). (Signature-reversal lineage: Lamb 2006 CMap; Sirota 2011.)
  - Classify: *Normalized* (D sig, Resid n.s./|Resid|≪|D|, T opposes D), *Persistent* (D sig, Resid sig,
    same sign), *Exacerbated* (T same sign as D).
  - Global statistic: proportion normalized, or projection `T·D/‖D‖`, tested via the rotation null.

### 5.3 Citations to verify at build
ROAST (Wu et al. 2010); CAMERA (Wu & Smyth 2012); RRHO/RRHO2 (Plaisier 2010 / Cahill 2018);
signature reversal (Lamb 2006 / Sirota 2011); fgsea (Korotkevich 2021); GSEA (Subramanian 2005);
circular analysis (Kriegeskorte 2009).

---

## 6. Figures (CvH, 05_Figures)

Consolidated core set; everything else -> `supp/`.

- **F00_QC** — proteoDA QC/norm reports + filtering waterfall + missingness + 4-method outlier
  consensus + **per-sample blood-contamination index**.
- **F01_Proteome_Overview** — PCA, DEP counts, rank plots (port the YvO-style overview).
- **F02_Baseline_and_Training** — **volcano-rings** for `CRvH_Baseline` and `CR_Training` independently
  (DEPs in the core + enrichment NES arcs) + enrichment summary panels. Built on the existing-but-unused
  `shared/volcano_ring.R` / `pathway_utils.R`.
- **F03_Reversal** — single composite: logFC scatter + NES scatter + RRHO2 (display), `fry`/`camera`
  bar (primary test), per-protein φ classification, rotation null throughout; `Resid` readout.
- **F04_WGCNA** — modules like YvO (merge the overlapping current F06 + F08 into one figure).

Stubs (current F02/F03/F04/F05/F07) deleted or repurposed into the above. LIMPA-vs-limma comparison
-> supp. Current `Reversal/` mini-stage is promoted/folded into F03.

---

## 7. Refactor + helper extraction (CvH)

Extract duplicated/buried logic into documented `R/` helpers with testthat tests:
- `R/pi_score.R` — `compute_pi_scores()` (currently duplicated across DEP scripts).
- `R/outlier_consensus.R` — 4-method consensus (missingness, PCA-Mahalanobis, MAD intensity,
  inter-sample correlation; ≥3/4).
- `R/mar_mnar.R` — 3-method MAR/MNAR vote (kmeans + global-logistic + left-tail).
- `R/blood_filter.R` — HPA set builders (M, B) and `C = B \ M`, plus the per-sample contamination index.
- `R/reversal.R` — φ fraction, classification, rotation-null helpers, fry/camera wrappers.

Long scripts to modularize: `01_run_normalization.R` (380) -> filtering helpers; `03_dep_robustness.R`
(436) -> per-arm functions; `04_dep_overview.R` (394). Keep `R/cvh_design.R` as the shared design
entrypoint (now with the 3-group `group_time` derivation + supplement retained for sensitivity).

Tests target: 80%+ line coverage on statistical functions (Pi-score, outlier consensus, MAR/MNAR,
blood-set logic, reversal stats), real-DB fixtures over mocks per repo preference.

---

## 8. Quarto document (CvH)

`docs/CvH_pipeline.qmd` — ground-up narrative:
1. Study design (3-group rationale, why supplement is pooled, paired structure).
2. Stage 01 Filtering — HPA muscle filter, blood `C=B\M` algorithm (with the kept/removed examples),
   dedup, missingness, outlier consensus.
3. Stage 02 Normalization — cycloess, QC.
4. Stage 03 Imputation — MAR/MNAR consensus, missForest, why DEP avoids it.
5. Stage 04 DEP — the single model, proteoDA/duplicateCorrelation mechanics + interpretation, contrasts,
   Pi-score, sensitivity arms.
6. Stage 05 Figures — what each figure shows and the reversal method (3 tiers + rotation null).
7. Reproducibility (seeds, frozen handoffs).

Code chunks `eval: false`; frozen figures/tables embedded from `c_data/`/`b_reports/`.

---

## 9. Mito (phase 2 — outline, to be detailed in its own spec)

Lighter scope (model already clean):
- **Document, don't restructure** the `~ 0 + group + (1|Replicate)` duplicateCorrelation model.
- Quarto `docs/Mito_pipeline.qmd` mirroring the CvH structure.
- Figure consolidation: keep F01–F06 (all real); extract the ~85%-duplicated `_panel_A_quadrant.R`
  (F03/F04) into one parameterized engine; document the `04_Figures` (engines/caches) vs `05_Figures`
  (composites) split.
- Optionally apply the same HPA-class-derived approach to the FBS/keratin contaminant step (rat
  orthologs), lighter-touch since Mito's problem is FBS/keratin not blood.
- **Cell-culture white paper** `docs/Mito_white_paper.md`: scientific goals (disease/intervention/
  rescue/interaction), H9c2 rat-cardiomyocyte + PHE + mito-transplant design, deliverables, and a
  **placeholder academic budget** (MS instrument time, mito isolation, cell-culture consumables,
  reagents, labor) marked as estimates for core-facility quotes.

Mito will get its own brainstorm (esp. white-paper cost detail) and spec before building.

---

## 10. Out of scope (YAGNI)
- No new imputation methods (keep missForest; benchmark stays as-is).
- No change to Mito's statistical model.
- No live/executable Quarto (frozen outputs only).
- No supplement as a contrast of interest (sensitivity only).

## 11. Open items to resolve at build
- Verify HPA class names / tissue-expression columns against the actual download.
- Verify all method PMIDs/DOIs (§3.4, §4.2, §5.3) via literature tools before drafting methods prose.
- Confirm eBayes `trend` decision (§4.4) empirically.
- Confirm ρ̂ is clearly positive (else reconsider blocking).

## 12. Build sequence
1. CvH: 01_Filtering (incl. blood algorithm + tests) → 02_Normalization → 03_Imputation → 04_DEP
   (model + reversal + sensitivity) → 05_Figures (F00–F04) → Quarto.
2. Mito: own spec → document + figure consolidation → white paper.
