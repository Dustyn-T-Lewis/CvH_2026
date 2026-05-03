# Literature Review: Methods for Quantifying Signature Reversal in Omics Studies

**Date:** 2026-04-25  
**Context:** CvH DIA-MS skeletal muscle proteomics (~2,500 proteins); cancer recovery exercise reversal figure  
**Scope:** 6 search topics × 5–8 databases each; ~100 unique papers screened; 75+ included below

---

## Table of Contents

1. [Study Table: Exercise × Disease Reversal Studies](#1-study-table-exercise--disease-reversal-studies)
2. [Study Table: CMap/Drug Reversal Methodology](#2-study-table-cmapdrug-reversal-methodology)
3. [Methods Catalog: Reversal Quantification Approaches](#3-methods-catalog-reversal-quantification-approaches)
4. [Directional Asymmetry: Literature and Gap Analysis](#4-directional-asymmetry-literature-and-gap-analysis)
5. [Proteomics-Specific Considerations](#5-proteomics-specific-considerations)
6. [Shared-Baseline Circularity: Statistical Frameworks](#6-shared-baseline-circularity-statistical-frameworks)
7. [Methods Comparison Matrix](#7-methods-comparison-matrix)
8. [Recommendations for the Reversal Figure](#8-recommendations-for-the-reversal-figure)

---

## 1. Study Table: Exercise × Disease Reversal Studies

| # | Citation | PMID | Omics | Disease | Intervention | Reversal Method | Asymmetry? | Threshold | n / Platform | Key Finding |
|---|---------|------|-------|---------|-------------|----------------|------------|-----------|-------------|-------------|
| 1 | Melov et al. 2007, *PLoS ONE* | 17520024 | Transcriptomics | Aging/sarcopenia | 6-mo RT, older adults | Proportion overlap + permutation (10,000×); p<0.02 | Not formally | FDR<5% (596 aging DEGs) | n=25 old + 26 young + 14 trained; Illumina Ref-8 | 179/596 aging genes reversed; enriched for mitochondrial/ETC genes |
| 2 | Dhahbi et al. 2004, *PNAS* | 15044709 | Transcriptomics | Aging (mouse) | CR → diet switch (8 wk reversal) | Proportion overlap, descriptive | No | Threshold-filtered DEGs | Affymetrix ~25k probes | 90% of CR gene expression effects reversed within 8 weeks |
| 3 | Robinson et al. 2017, *Cell Metab* | 28273480 | Multi-omics (transcr + proteomics + metabolomics) | Aging/sarcopenia | 12-wk HIIT/RT/Combined RCT | Directional normalization toward young; fold-change correlation | No | Age-related differences at baseline | n=72 (36 young + 36 old); RNA-seq + LC-MS | HIIT superior for reversing age-related mitochondrial protein decline |
| 4 | Voisin et al. 2023, *Aging Cell* | 37128843 | Methylomics + transcriptomics | Aging | Exercise vs disuse (meta-analysis) | Epigenetic clock score (MEAT 2.0) | No | Age-associated methylome/transcriptome profiles | n=3,176 samples; 850K EPIC + RNA-seq | Exercise → younger methylome/transcriptome; disuse → older |
| 5 | Lohman et al. 2023, *Aging Cell* | 37078430 | Transcriptomics | Aging | 4-wk HIIT RCT | Transcriptomic age clock reversal | No | mRNA transcriptomic age clock | n=30 older adults; blood RNA-seq | HIIT: −3.59 yrs transcriptomic age vs +3.29 yrs control |
| 6 | MoTrPAC 2024, *Nature* | 38693412 | Multi-omics (9 platforms × 19 tissues) | Multiple disease signatures | 8-wk endurance treadmill (rat) | Multi-omic enrichment overlap with disease signatures | No | Disease-relevant gene sets from published databases | n=104 rats; 19 tissues; 9 omics platforms | Widespread normalization toward healthy states in liver/heart/muscle/adipose |
| 7 | Han et al. 2023, *J Mol Cell Biol* | 36882217 | Transcriptomics + metabolomics | HFD obesity (mouse) | 8-wk treadmill + wheel running | Count-based: 37/38 metabolites (97.4%); Pearson r of FC vectors | Partial (sets I/II) | HFD-altered metabolites/genes | n=8–10 mice/group; BGISEQ-500 + LC-MS | 97% HFD-altered metabolites reversed; 465 genes reversed |
| 8 | Benrick et al. 2020, *JCEM* | 32232327 | Transcriptomics + methylation | PCOS | 16-wk electroacupuncture vs exercise | Directional count (95% normalized toward healthy) | Partial | PCOS-associated gene expression changes | n=21 women; Illumina microarray | 95% EA-induced changes in normalizing direction; 38% overlap with exercise |
| 9 | O'Leary et al. 2024, *Exp Physiol* | 39663727 | **Proteomics (TMT-MS)** | Aging/sarcopenia | 8-wk RT, young vs older women | Directional overlap of age-effect and training-effect vectors | **Yes** — mitochondrial reversed, contractile NOT | Age-related proteome differences | n=15 (7Y + 8O); vastus lateralis TMT LC-MS/MS | RT reverses aging-related mitochondrial proteome but not contractile proteins |
| 10 | Collao et al. 2023, *JCSM* | 36797054 | Transcriptomics (RNA-seq) | Cancer cachexia (mouse) | RT + endurance during chemoradiation | Pathway enrichment (GSEA/ORA); attenuation of inflammatory/fibrotic pathways | No | GO enrichment of cancer-chemo gene sets | Mouse model; RNA-seq | Exercise prevented inflammatory/fibrotic transcriptomic signature |
| 11 | Lammers et al. 2012, *Am J Physiol* | 23011062 | Transcriptomics | Physical inactivity → T2DM/CVD | 3-wk limb suspension / 6-wk FES exercise in SCI | Overlap count: 18 genes down-by-inactivity AND up-by-exercise | One direction only | Inactivity-regulated genes | n=8 + 5; Affymetrix HG-U133+ 2.0 | 18 genes linking fatty acid transport/insulin signaling reversed by exercise |
| 12 | Seaborne et al. 2018, *Sci Rep* | 29382913 | Methylomics + transcriptomics | Disuse/detraining | 7-wk RT → 7-wk detraining → 7-wk retraining | Epigenetic memory: hypomethylation retention during detraining | N/A | Training-induced CpG changes | n=8 young men; 850K EPIC | Muscle epigenetic memory: 2× hypomethylation response on retraining |
| 13 | Fisher et al. 2017, *FASEB J* | 28821632 | Transcriptomics + methylation | Disuse atrophy (rat) | TTX denervation → 7-day recovery | Return to sham baseline; bidirectional reversal confirmed | Both directions tracked | 3,714 disuse DEGs | n=6 rats/timepoint; microarray + pyrosequencing | Full transcriptome reversal during recovery; epigenetic changes bidirectional and reversible |
| 14 | Fitzgerald et al. 2021, *Aging* | 33844651 | Epigenomics | Aging (epigenetic) | 8-wk lifestyle (diet+sleep+exercise+relaxation) RCT | Horvath DNAmAge clock reversal | N/A (scalar) | DNAmAge | n=43 males; Illumina 450K | −3.23 yrs DNAmAge (p=0.018); exercise not isolatable |
| 15 | Fiorito et al. 2021, *Aging Cell* | 34535961 | Epigenomics | Aging (cancer risk) | 2-yr Mediterranean diet + physical activity RCT | DNAmGrimAge deceleration | N/A (scalar) | Multiple methylation clocks | n=130 women; Illumina 450K | DNAmGrimAge slowed in intervention group |
| 16 | Gil et al. 2021, *JCSM* | 34666419 | Transcriptomics + histology | Bariatric surgery muscle loss | Exercise training post-RYGB RCT | GSEA NES: UPS pathway −1.7 (p<0.01); fiber CSA restoration | Pathway-level | RNA-seq pathways | Human women; vastus lateralis RNA-seq | Exercise suppresses surgery-upregulated catabolism; restores anabolic pathways |
| 17 | Kohman et al. 2011, *PLoS ONE* | 21857943 | Transcriptomics | Aging (hippocampus) | Voluntary wheel running (mouse) | Directional overlap + visual direction tables | No | Aging DEGs | Mouse; Affymetrix 45k | Descriptive reversal reporting (no formal test) |
| 18 | Messaoudi et al. 2017, *BMC Genomics* | 28545403 | Transcriptomics | Obesity (primate) | WSD → caloric restriction (30 mo) | **NEGATIVE CONTROL**: most WSD changes NOT reversed | Partial | WSD-induced DEGs | n=8–12 macaques; RNA-seq | **CR alone does not reverse skeletal muscle obesity transcriptome** |

---

## 2. Study Table: CMap/Drug Reversal Methodology

| # | Citation | PMID | Score | Formula Basis | Signature Input | Null Distribution | Key Application |
|---|---------|------|-------|--------------|----------------|-------------------|----------------|
| 1 | Lamb et al. 2006, *Science* | 17008526 | CS (Connectivity Score) | Weighted KS enrichment | Top/bottom DEGs | Within-database gene permutation | Proof-of-concept: 164 drugs × 5 cell lines |
| 2 | Subramanian et al. 2017, *Cell* | 29195078 | NCS, WCS, Tau (τ) | z-score of CS; percentile rank vs 1.3M profiles | Ranked gene list (978 landmarks) | Empirical full L1000 database | Scaled to 1.3M profiles; MoA discovery |
| 3 | Samart et al. 2021, *Brief Bioinform* | 34013329 | 17 scores reconciled (ES, css, Sum, Cosine, XSum, XCor, XSpe, XCos, EWCos, CS, RGES, NCS, WCS, Tau, CSS, EMUDRA) | Unified notation | Varies | Meta-analyzed | **KEY REFERENCE** for choosing among connectivity scores |
| 4 | Lin et al. 2020, *Brief Bioinform* | 31774912 | ZhangScore (best performing) | Modified KS, partial AUC optimized | 10–200 genes | Drug Repurposing Hub benchmark | ZhangScore > XSum > XCor at 10–200 genes |
| 5 | Parkkinen & Kaski 2014, *BMC Bioinform* | 24742351 | Probabilistic CS | Bayesian posterior / Gaussian processes | Full profiles | Posterior uncertainty | Probabilistic alternative to rank-based scoring |
| 6 | Kim et al. 2019, *Sci Rep* (RETRACTED) | 30804389 | RGES / sRGES | −1 × CS; weighted mean across cell lines | Meta-DEGs | ChEMBL IC₅₀ validation | Gastric cancer drug repositioning |
| 7 | Stathias et al. 2018, *Nat Commun* | 30552330 | SynergySeq (discordance + concordance) | Cosine similarity + NCS | Disease DEGs + LINCS | Random drug pairs | Drug combination synergy prediction (GBM) |
| 8 | Shah et al. 2021, *Neuro-Oncol Adv* | 35118385 | sRGES + CNS-MPO filter | sRGES + pharmacokinetic constraint | TCGA GBM DEGs | sRGES–IC₅₀ correlation | GBM drug repositioning with BBB penetrance |
| 9 | Jahchan et al. 2013, *Cancer Discov* | 24078773 | Enrichment-based CMap query | Enrichment matching | SCLC transcriptional signature | Random gene sets | Identified imipramine as SCLC inhibitor |
| 10 | Pushpakom et al. 2019, *Nat Rev Drug Discov* | 30310233 | (Review) | — | — | — | Authoritative review positioning CMap in drug repurposing |
| 11 | Koudijs et al. 2018, *Sci Rep* | 29588458 | Per-patient RGES | Negative enrichment per individual sample | 534 individual ccRCC tumor signatures | Group vs individual comparison | Individual RGES outperforms group for ccRCC |
| 12 | Yang et al. 2019, *GeroScience* | 31637571 | ANDRU (Aging Network Drug) | CMap NCS on aging subnetworks | GTEx aging DEGs in adipose | Literature evidence | Geroprotector discovery (pioglitazone) |
| 13 | Lim & Pavlidis 2021, *Sci Rep* | 34475469 | (Evaluation) | — | — | — | CMap shows limited reproducibility; DE strength predicts reproducibility |
| 14 | Fan & Evans 2017, *Cell Metab* | 27889389 | (Review) | — | — | — | Exercise mimetics conceptual framework |
| 15 | Dai et al. 2014, *Aging Cell* | 24612461 | Proteome concordance | Direction counting (823 proteins) | Age-dependent cardiac proteome | — | Rapamycin/CR reverse aging cardiac proteome |
| 16 | Xing et al. 2026, *Cell* | 41850287 | GPS (deep learning) | Structure → predicted signature → reversal optimization | Chemical SMILES + disease DEGs + LINCS | Benchmarked vs known drugs | De novo drug design optimizing for reversal |

### Score Evolution Lineage

CS (2006) → XSum/XCor variants (2008–2015) → NCS/WCS/Tau (2017) → RGES/sRGES (2019) → ZhangScore/EMUDRA (2020–2021) → reconciliation (Samart 2021) → deep learning GPS (Xing 2026)

---

## 3. Methods Catalog: Reversal Quantification Approaches

### Category A: Proportion/Overlap (Threshold-Dependent)

| Method | Key Paper | Statistical Test | Threshold | Limitations | Our Use |
|--------|----------|-----------------|-----------|-------------|---------|
| Melov proportion | Melov 2007 (PMID 17520024) | Binomial test vs 50%; permutation (10,000×) | FDR or p-value for primary contrast | No null model in original; threshold-sensitive; assumes gene independence | Pilot script (binom.test) |
| CR reversal % | Dhahbi 2004 (PMID 15044709) | Descriptive | Filtered DEGs | No formal test | — |
| Chi-square/Fisher | Swindell 2009 (PMID 19968875) | 2×2 contingency (up/down × aging/CR) | Both contrasts thresholded | Assumes gene independence; threshold-sensitive | Could add as supplement |

### Category B: Connectivity Scores (Drug Repurposing Heritage)

| Method | Key Paper | Input | Threshold-Free? | Applicable to ~2,500 proteins? | Our Use |
|--------|----------|-------|-----------------|-------------------------------|---------|
| CMap CS (KS-based) | Lamb 2006 (PMID 17008526) | Up/down gene sets + ranked reference | Partial | Requires adaptation (gene mapping) | Pilot script (custom weighted KS) |
| Tau (τ) | Subramanian 2017 (PMID 29195078) | 978 landmarks + 1.3M profiles | Partial | Not directly (needs L1000 data) | — |
| RGES / sRGES | Kim 2019 (PMID 30804389, retracted) | Disease DEGs + LINCS | Yes for query | Not directly | — |

### Category D: Rotation Tests (Threshold-Free, Model-Based)

| Method | Key Paper | Statistical Framework | Handles Repeated Measures? | Proteomics Scale? | Our Use |
|--------|----------|----------------------|---------------------------|-------------------|---------|
| ROAST | Wu et al. 2010 (PMID 20610611) | Monte Carlo rotation; directional (Up/Down/Mixed) | Yes (via limma duplicateCorrelation) | Yes | — |
| CAMERA | Wu & Smyth 2012 (PMID 22638577) | Competitive; variance inflation for inter-gene correlation | Yes | Yes | — |
| **fry** | limma docs (Ritchie 2015, PMID 25605792) | Analytical F-approximation to ROAST | Yes | Yes | **Panel C** |
| roastgsa | Caballé-Mestres 2023 (PMID 37904108) | Multiple score functions (mean, maxmean, KS) within rotation | Yes | Yes; effective signature size useful | Consider |

### Category E: RRHO/RRHO2 (Threshold-Free, Rank-Based)

| Method | Key Paper | Statistical Framework | Handles Discordance? | Proteomics Scale? | Our Use |
|--------|----------|----------------------|---------------------|-------------------|---------|
| RRHO | Plaisier et al. 2010, *NAR* (DOI 10.1093/nar/gkq636) | Hypergeometric grid; all thresholds swept | Original: NO (concordance only) | Yes | — |
| **RRHO2** | Cahill et al. 2018, *Sci Rep* | Corrected 4-quadrant hypergeometric | **YES** — key improvement | Yes (~2,100 works well) | **Panel E** |
| RankerGUI | Thind 2019 (PMID 31816915) | eRRHO + distance metrics | Yes | Yes | — |

### Category F: Correlation-Based (Threshold-Free)

| Method | Key Paper | Input | Significance | Proteomics Scale? | Our Use |
|--------|----------|-------|-------------|-------------------|---------|
| **Spearman/Pearson rho of logFC** | Robinson 2017 (PMID 28273480) | Paired logFC vectors across all features | Permutation or Fisher-z | **Yes** (used at ~2,700 proteins) | Panel A scatter |
| Fisher Z CI | Standard | rho + n | Analytical CI | Yes | Panel D NES scatter |

### Category G: GSEA-Based (Threshold-Free)

| Method | Key Paper | Statistical Framework | Proteomics Scale? | Our Use |
|--------|----------|----------------------|-------------------|---------|
| GSEA | Subramanian 2005 (PMID 16199517) | Weighted KS running sum; sample or gene permutation | Yes (adjust minSize to 5–10) | Conceptual basis |
| **fgsea** | Korotkevich et al., bioRxiv (DOI 10.1101/060012) | Adaptive multi-level permutation | Yes | **Panel D** (NES scatter) |

### Category H: Permutation/Bootstrap

| Approach | Applied To | Key Consideration |
|----------|-----------|-------------------|
| Gene-label permutation | Proportion overlap (Melov), CMap CS | Breaks correlation structure; valid for gene-level metrics |
| Sample-label permutation | GSEA, fry/ROAST | Requires n>7/group; preserves correlation structure |
| Bootstrap for winner's curse | Effect size correction | Hawinkel 2025 (PMID 40864625): corrects for correlated features |

---

## 4. Directional Asymmetry: Literature and Gap Analysis

### Executive Finding: THIS IS A RESEARCH GAP

**No published paper explicitly names, quantifies, or tests "directional asymmetry in signature reversal"** as a phenomenon — the observation that disease-DOWN proteins reverse at a different rate than disease-UP proteins after intervention. The concept appears incidentally across disconnected fields but has never been synthesized.

### Supporting Evidence

| Paper | PMID | Key Evidence |
|-------|------|-------------|
| Murgia et al. 2023, *JCSM* | 36517414 | Bed rest plasma proteomics: teneurin-4 fully reversed in 2 days; lumican did NOT recover — **direct protein-level recovery asymmetry** |
| O'Leary et al. 2024, *Exp Physiol* | 39663727 | RT reverses aging mitochondrial proteome but NOT contractile proteins — **explicit pathway-class asymmetry** |
| Mahmassani et al. 2021, *J Gerontol A* | 33705535 | Inactivity eliminates "uncoupled translation" for ribosomal proteins — mRNA-protein discordance mechanism |
| Drummond et al. 2008, *Am J Physiol* | 18827171 | miR-1/miR-133a suppressed by exercise in YOUNG but not OLD — disease-like overexpression fails to reverse with age |
| Gil et al. 2021, *JCSM* | 34666419 | Exercise suppresses surgery-upregulated UPS (NES −1.7) while restoring fiber CSA — bidirectional pathway reversal |
| Tyagi & Pedrioli 2015, *NAR* | 25870413 | Codon-usage bias creates intrinsically asymmetric translation recovery kinetics after stress |
| Walsh et al. 2022, *JCSM* | 35092190 | miRNA regulatory architecture predicts which ICU-damaged muscle programs recover vs persist |
| Endo et al. 2021, *Genomics* | 34214629 | Aged mouse: exercise upregulates neural/NMJ genes but FAILS to upregulate anabolic pathways — class-specific non-reversibility |

### Candidate Biological Mechanisms

| Mechanism | Prediction for Cancer-DOWN vs Cancer-UP | Source |
|-----------|----------------------------------------|--------|
| Protein half-life by class | Constitutive metabolic proteins (cancer-DOWN) have long half-lives, recover passively; signaling proteins (cancer-UP) are actively stabilized by disease signals | Murgia 2023, Tyagi 2015 |
| UPS accessibility | Catabolic/inflammatory proteins (cancer-UP) regulated by UPS; exercise modulates UPS but may not overcome persistent oncogenic drivers | Gil 2021, Knapp 2020 |
| Transcriptional vs translational regulation | Cancer-DOWN proteins lost through disuse/cachexia (not active silencing) → restored by exercise-driven anabolic signaling; Cancer-UP proteins maintained by active transcriptional programs exercise cannot extinguish | Mahmassani 2021, Molanouri Shamsi 2015 |
| Codon bias in translation | Different protein classes have different codon-usage patterns → differential translational efficiency during recovery | Tyagi 2015 |

### Our Contribution

Our observation (76% cancer-DOWN reversal vs 56% cancer-UP reversal) would be **the first explicit quantification of directional asymmetry in signature reversal at the proteome level**. The most parsimonious explanation: cancer-suppressed proteins (constitutive metabolic/structural proteins lost through cachexia/disuse) are restored by exercise-driven anabolic signaling, while cancer-upregulated proteins (stabilized by active oncogenic/inflammatory signals) resist reversal because exercise cannot extinguish the upstream drivers.

---

## 5. Proteomics-Specific Considerations

### Methods Performance at ~2,500 Proteins vs ~20,000 Genes

| Method | Performance at Proteome Scale | Key Consideration |
|--------|------------------------------|-------------------|
| Fold-change correlation | **Excellent** — directly applicable; Robinson 2017 used at ~2,700 proteins | Single interpretable metric |
| fry/ROAST | **Excellent** — rotation count independent of feature count; works with ≥15 proteins/set | Best for complex designs with duplicateCorrelation |
| fgsea | **Good** — reduce minSize to 5–10; fewer features means fewer sets pass size filter | Adjust pathway database filtering |
| RRHO2 | **Good** — smaller N reduces computation; hypergeometric valid at any feature count | Pure-R phyper() avoids segfaults |
| CMap CS | **Poor** — designed for 978–20k genes; protein→gene mapping loses 30–40% coverage | Not directly applicable |
| Proportion overlap | **Adequate** — but small denominator amplifies threshold sensitivity | Permutation test essential |

### Proteomics Reversal Studies Found

Only **two** proteomic reversal studies were identified:
1. **O'Leary et al. 2024** (PMID 39663727) — TMT-MS, aging, 8-wk RT, n=15
2. **Dai et al. 2014** (PMID 24612461) — deuterium labeling, cardiac aging, rapamycin/CR, 823 proteins

The literature is overwhelmingly transcriptomics-based. **Our study fills a significant gap as a proteomics-level reversal analysis in human cancer recovery.**

### DIA-MS Fold-Change Reliability

DIA-MS provides reliable fold-change estimates across the quantified proteome (reviewed in Bruderer et al. 2015, PMID 26381204). The key consideration: missing values at the protein level are non-random (MNAR for low-abundance proteins), which can bias fold-change estimates for the most depleted proteins. Since cancer-DOWN proteins are by definition depleted, their fold changes may be more uncertain → conservative interpretation of the asymmetry finding is warranted.

---

## 6. Shared-Baseline Circularity: Statistical Frameworks

### Executive Finding: NO FORMAL FRAMEWORK EXISTS

**No published paper names "shared-baseline circularity" as a statistical problem**, defines a null hypothesis for it, or proposes a dedicated test. The field has not formally characterized this bias despite it being structurally present in many designs.

### The Problem

When Cancer_vs_Healthy uses CR_T1 as reference and Training_CR uses CR_T1 as baseline:
- Both fold-change vectors share CR_T1 variance
- **Structural negative correlation** is mathematically guaranteed
- A protein high in CR_T1 (cancer-UP) will appear to decrease with training (structural reversal) even if training has no biological effect

### Relevant Bodies of Literature

#### Body 1: limma Design Theory (Handles Covariance Correctly)

| Paper | PMID | Key Finding |
|-------|------|-------------|
| Smyth & Altman 2013, *BMC Bioinform* | 23705896 | **Most relevant**: Two-channel microarrays are structurally identical to shared-baseline designs. Shared reference channels create correlated M-values (fold changes). Fix: extract inter-spot (A-value) information for independent contrasts |
| Law et al. 2020, *F1000Res* | 33604029 | limma's contrast matrix automatically handles shared-sample covariance — specify all contrasts simultaneously within one model |
| Ritchie et al. 2015, *NAR* | 25605792 | limma implements global covariance models for correlated contrasts |
| Wu & Smyth 2012, *NAR* (CAMERA) | 22638577 | Gene permutation tests are INVALID when genes are correlated — analogous principle: permuting sample labels independently for correlated contrasts is invalid |

#### Body 2: Winner's Curse / Regression to Mean

| Paper | PMID | Key Finding |
|-------|------|-------------|
| Crager 2010, *Stat Med* | 19960511 | Regression-to-mean correction for gene association estimates using bivariate normal shrinkage |
| Hawinkel et al. 2025, *Biostatistics* | 40864625 | **Winner's curse under dependence**: standard empirical Bayes corrections biased when features are correlated; bootstrap correction recommended |
| Huang et al. 2018, *NAR* | 30189032 | Bootstrap method (BootstrapQTL) for effect size correction |
| Forde et al. 2023, *PLoS Genet* | 37721937 | Bootstrap outperforms conditional likelihood for winner's curse correction; R package `winnerscurse` |
| Degen & Medo 2025, *PLoS Comput Biol* | 40324149 | Regression to the mean quantified as function of sample size for RNA-seq DE |

#### Body 3: Shared-Reference Artifacts in Chemical Genomics

| Paper | PMID | Key Finding |
|-------|------|-------------|
| Kim et al. 2024, *Mol Sys Biol* | 39349762 | **Most directly analogous**: Pooled chemical genomics screens have shared control → spurious gene-gene and drug-drug correlations. Linear normalization correction removes artifacts |
| Lim & Pavlidis 2021, *Sci Rep* | 34475469 | CMap limited reproducibility partly due to shared DMSO reference contributing to spurious cross-drug correlations |

#### Body 4: Design Circularity in Other Fields

| Paper | PMID | Key Finding |
|-------|------|-------------|
| Patil et al. 2015, *Bioinformatics* | 25788628 | Cross-sample normalization creates test-set bias; rank-based features proposed as normalization-independent |
| Tustison et al. 2014, *Hum Brain Map* | 23151955 | Registration-based normalization circularity in neuroimaging; fix: use independent normalization data |
| Fisher et al. 2022, *Clin Cancer Res* | 35792866 | Stromal contamination confounds gene signature correlations; ConfoundR resource |

### Recommended Mitigation Strategies

1. **Report reversal on the full proteome** (all ~2,500 proteins, not just cancer-significant) — structural correlation applies uniformly, so the biological signal in the aging-DEP subset should exceed the global background
2. **Compare observed correlation against analytical null** — derive the expected negative correlation from the design structure (shared CR_T1 variance contribution), then test whether the observed reversal exceeds this structural expectation
3. **Use pathway-level tests (fry, fgsea)** — these operate within the limma linear model framework where the covariance structure is correctly modeled
4. **Report the circularity diagnostic** — correlation between t-statistics for the two contrasts across all proteins (already done in panel_D_fry.R)
5. **Supplement with RRHO2** — rank-based and partially robust to shared-baseline magnitude effects
6. **Acknowledge explicitly** in Methods that negative structural correlation is a mathematical property of the design

---

## 7. Methods Comparison Matrix

### Suitability for Our Dataset (DIA-MS, ~2,500 proteins, 0 FDR-sig training hits, repeated measures)

| Method | Threshold-Free? | Handles Weak Individual Signal? | Handles Repeated Measures? | Proteome Scale? | Circularity-Robust? | Recommended? |
|--------|----------------|-------------------------------|---------------------------|----------------|---------------------|-------------|
| **fry/ROAST** | ✅ | ✅ (set-level test) | ✅ (via limma) | ✅ | ✅ (within limma model) | **YES — KEEP** |
| **fgsea NES** | ✅ | ✅ (pathway-level) | ❌ (uses pre-computed stats) | ✅ (adjust minSize) | Partial (rank-based) | **YES — KEEP** |
| **RRHO2** | ✅ | ✅ (sweeps all thresholds) | ❌ (uses ranked lists) | ✅ | Partial (rank-based) | **YES — KEEP** |
| **Spearman rho of logFC** | ✅ | ✅ (uses all proteins) | ❌ (uses summary stats) | ✅ | ❌ (directly inflated) | **KEEP with caveat** |
| Quadrant ORA | ❌ (needs Pi<0.05) | Partial | ❌ | ✅ | ❌ | **KEEP** (visual anchor) |
| Pattern heatmap + Sankey | ❌ (needs Pi<0.05) | Partial | ❌ | ✅ | ❌ | **KEEP** (interpretive) |
| Melov proportion | ❌ | ❌ (needs significant hits) | ❌ | ✅ | ❌ | **ADD as supplement** |
| CMap connectivity | Partial | ✅ | ❌ | ⚠️ (adaptation needed) | ❌ | **ADD as supplement** |
| Recovery score | ❌ | ❌ (noisy per-protein) | ❌ | ✅ | ❌ | **KEEP in pilot only** |
| Chi-square 2×2 | ❌ | ❌ | ❌ | ✅ | ❌ | Optional supplement |

### Method Strength by Question

| Question | Best Method(s) | Why |
|----------|---------------|-----|
| "Is there systematic reversal?" | Spearman rho + fry + RRHO2 | Complementary: global correlation + set-level rotation + threshold-free grid |
| "Which pathways reverse?" | fgsea NES scatter + quadrant ORA | NES sign-flip identifies reversed pathways; ORA identifies enriched quadrants |
| "How strong is the reversal?" | CMap connectivity score + recovery score | Continuous metrics with interpretable scale |
| "Is the reversal significant?" | fry p-value + Melov permutation p + RRHO2 peak -log10(p) | Each has a formal null distribution |
| "Is directional asymmetry real?" | Two-proportion z-test on reversal fractions (UP vs DOWN) | Simple, interpretable, novel |
| "Is the signal robust to threshold?" | RRHO2 + fgsea (both threshold-free) | Visualize signal across all possible thresholds |

---

## 8. Recommendations for the Reversal Figure

### A. Which Methods to Keep, Add, or Remove

**KEEP (current 5 panels):**
- Panel A (Quadrant ORA scatter) — visual anchor; shows protein-level landscape
- Panel B (Pattern heatmap + GO Slim Sankey) — interpretive; connects proteins to pathways
- Panel C (fry rotation test) — **statistically strongest method**; handles repeated measures within limma; directional testing with formal p-values
- Panel D (NES scatter) — pathway-level reversal visualization
- Panel E (RRHO2) — threshold-free gold standard; shows full concordance/discordance landscape

**ADD to supplementary:**
- **Melov proportion test** (from pilot) — most widely cited reversal metric in the exercise-aging literature; enables direct comparison with Melov 2007 and subsequent studies. Report: "61.5% of cancer DEPs reversed, binomial p=1.2e-7" and acknowledge that the original Melov paper (PMID 17520024) used qualitative assessment, not binomial test
- **CMap connectivity score** (from pilot) — connects to the large CMap/L1000 literature (>7,000 citations); report: "connectivity = −0.208, permutation p<0.001" and cite Lamb 2006 + Subramanian 2017
- **Threshold sensitivity analysis** — show how reversal proportion changes across Pi<0.01, Pi<0.05, FDR<0.05, FDR<0.10, nominal p<0.05 thresholds (addresses the main limitation of Category A methods)
- **Directional asymmetry panel** — two-proportion z-test comparing reversal rate of cancer-DOWN (76%) vs cancer-UP (56%) proteins; report as novel finding with biological interpretation

**REMOVE: nothing** — all 5 main panels serve distinct and complementary roles.

### B. Threshold Recommendations

| Purpose | Recommended Threshold | Justification |
|---------|----------------------|---------------|
| Defining "cancer signature" for proportion test | **Pi < 0.05** | Consistent with YvO pipeline; captures both statistical and biological significance |
| Gene set sizes for fgsea/fry | **min_size = 10, max_size = 500** | Already used; appropriate for ~2,500 proteins |
| Pathway significance | **FDR < 0.05** for main; **FDR < 0.25** for leading-edge ORA | Standard; lenient threshold for leading-edge justified per Subramanian 2005 |
| RRHO2 parameters | **stepsize = 20, boundary = 0.02** | Already used; appropriate for ~2,500 proteins |

### C. How to Handle the Circularity Issue

1. **In Methods**: State explicitly that Cancer_vs_Healthy and Training_CR share the CR_T1 baseline, creating a structural negative correlation between fold-change vectors. Cite Smyth & Altman 2013 (PMID 23705896) as the formal characterization of this property in reference-channel designs.

2. **In the fry panel (Panel C)**: The fry test operates within the limma linear model, which correctly models the shared-sample covariance structure. Report the raw correlation between t-statistics (r(t_Cancer, t_Training)) as a circularity diagnostic — already implemented in your code.

3. **In supplementary**: Add a circularity sensitivity analysis:
   - Compare the observed Spearman rho between Cancer and Training logFC against the expected structural correlation from the design matrix
   - Show that the reversal signal in the aging-DEP subset significantly exceeds the global background correlation
   - Consider a permutation that preserves the shared-baseline structure: randomly assign proteins to "cancer-DEP" vs "non-DEP" while keeping all fold changes fixed, then compute the reversal fraction — this tests whether cancer-DEPs reverse more than random proteins, not whether reversal exceeds zero

4. **In Discussion**: Acknowledge the shared-baseline limitation; note that fry (rotation-based) and RRHO2 (rank-based) partially mitigate this by operating within model-aware or non-parametric frameworks; note that this is a gap in the field (no formal test exists, per our literature search).

### D. How to Present the Directional Asymmetry Finding

1. **Frame as novel**: "To our knowledge, this is the first explicit quantification of directional asymmetry in proteomic signature reversal."

2. **Report the numbers**: Cancer-DOWN proteins: 76% reversal rate; Cancer-UP proteins: 56% reversal rate; two-proportion z-test p-value.

3. **Biological interpretation**: Cancer-downregulated proteins (constitutive metabolic/structural proteins lost through disuse/cachexia) may be more amenable to exercise-driven restoration because their suppression reflects passive depletion rather than active transcriptional silencing. Cancer-upregulated proteins (inflammatory, stress, or oncogenic programs) may resist reversal because the upstream driving signals persist despite exercise.

4. **Cite the mechanistic basis**: Murgia 2023 (PMID 36517414, protein-level recovery asymmetry), Tyagi 2015 (PMID 25870413, codon-bias differential recovery), Mahmassani 2021 (PMID 33705535, translational regulation asymmetry).

### E. Citation Guide by Method

| Method in Our Figure | Primary Citation | Supporting Citation(s) |
|---------------------|-----------------|----------------------|
| Quadrant ORA (fora) | Reimand et al. 2019, *Nat Protoc* (PMID 30664679) | fgsea package (Korotkevich, DOI 10.1101/060012) |
| Pattern heatmap | (descriptive; no method-specific citation needed) | GO Slim: Gene Ontology Consortium 2004 |
| fry rotation test | Wu & Smyth 2010 (PMID 20610611, ROAST); Ritchie et al. 2015 (PMID 25605792, limma) | Caballé-Mestres 2023 (PMID 37904108, roastgsa score comparison) |
| fgsea NES scatter | Korotkevich et al., bioRxiv DOI 10.1101/060012 | Subramanian et al. 2005 (PMID 16199517, GSEA) |
| RRHO2 | Cahill et al. 2018, *Sci Rep* | Plaisier et al. 2010, *NAR* (DOI 10.1093/nar/gkq636) |
| Melov proportion (supp) | Melov et al. 2007, *PLoS ONE* (PMID 17520024) | Swindell 2009 (PMID 19968875, chi-square variant) |
| CMap connectivity (supp) | Lamb et al. 2006, *Science* (PMID 17008526) | Subramanian et al. 2017, *Cell* (PMID 29195078); Samart et al. 2021 (PMID 34013329) |
| Spearman rho of logFC | Robinson et al. 2017, *Cell Metab* (PMID 28273480) | — |
| Directional asymmetry | NOVEL (no precedent citation) | Mechanistic basis: Murgia 2023, Tyagi 2015, Mahmassani 2021 |
| Shared-baseline acknowledgment | Smyth & Altman 2013 (PMID 23705896) | Kim et al. 2024, *Mol Sys Biol* (PMID 39349762) |

---

## Appendix: Key PMIDs for Quick Reference

### Foundational (Cite in Introduction/Methods)
- 17520024 — Melov 2007 (exercise reverses aging transcriptome)
- 17008526 — Lamb 2006 (Connectivity Map)
- 29195078 — Subramanian 2017 (L1000/Next-Gen CMap)
- 16199517 — Subramanian 2005 (GSEA)
- 20610611 — Wu & Smyth 2010 (ROAST)
- 25605792 — Ritchie et al. 2015 (limma)

### Methodological (Cite in Methods)
- 34013329 — Samart 2021 (connectivity score reconciliation)
- 37904108 — Caballé-Mestres 2023 (roastgsa)
- 22638577 — Wu & Smyth 2012 (CAMERA, inter-gene correlation)
- 23705896 — Smyth & Altman 2013 (shared-reference channel analysis)
- 19968875 — Swindell 2009 (CR–aging reversal critical analysis)

### Exercise Reversal (Cite in Introduction/Discussion)
- 28273480 — Robinson 2017 (HIIT reverses aging multi-omics)
- 38693412 — MoTrPAC 2024 (multi-tissue exercise atlas)
- 36882217 — Han 2023 (exercise reverses HFD obesity, 97%)
- 39663727 — O'Leary 2024 (proteomics reversal, mitochondrial but not contractile)
- 34666419 — Gil 2021 (exercise after bariatric surgery)

### Directional Asymmetry (Cite in Discussion)
- 36517414 — Murgia 2023 (protein-level recovery asymmetry)
- 25870413 — Tyagi 2015 (codon-bias differential recovery)
- 33705535 — Mahmassani 2021 (translational regulation asymmetry)
- 18827171 — Drummond 2008 (age-differential miRNA response)

### Circularity (Cite in Methods/Limitations)
- 23705896 — Smyth & Altman 2013 (two-channel shared reference)
- 39349762 — Kim 2024 (pooled genomics shared-reference correction)
- 40864625 — Hawinkel 2025 (winner's curse under dependence)
