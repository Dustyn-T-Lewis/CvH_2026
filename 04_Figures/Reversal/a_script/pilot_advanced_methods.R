# Pilot: Advanced Reversal Methods
#
# 1. Analytical structural null (pooled-screen correction analog)
#    Kim et al. 2024, Mol Sys Biol (PMID 39349762)
#    Smyth & Altman 2013, BMC Bioinform (PMID 23705896)
#
# 2. roastgsa multi-score comparison
#    Caballé-Mestres et al. 2023, BMC Bioinform (PMID 37904108)
#
# 3. Corrected correlation (observed - structural)
#    Hawinkel et al. 2025, Biostatistics (PMID 40864625): winner's curse context
#
# Design:
#   Cancer_vs_Healthy = mean(CR_T1) - mean(H_T1)
#   Training_CR       = mean(CR_T2) - mean(CR_T1)
#   Shared term: CR_T1 appears with +1 in CvH and -1 in TR
#   → structural Cov(logFC_CvH, logFC_TR) = -Var(mean(CR_T1))
#   → expected negative r even under H0: no biological reversal
setwd(here::here())
source("04_Figures/shared/style.R")
library(tidyverse)
library(limma)
library(patchwork)

PILOT_DIR <- "04_Figures/Reversal/b_reports/pilot"
DAT       <- "04_Figures/Reversal/c_data/pilot_advanced"
dir.create(PILOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT,       recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

# ═══════════════════════════════════════════════════════════════════════════════
# DATA
# ═══════════════════════════════════════════════════════════════════════════════
dal    <- readRDS("02_Imputation/c_data/01_DAList_imputed.rds")
dep_df <- read_csv("03_DEP/c_data/03_combined_results_CRvH.csv",
                    show_col_types = FALSE)

mat  <- dal$data
meta <- as.data.frame(dal$metadata)

dep_fc <- dep_df %>%
  filter(!is.na(logFC_Cancer_vs_Healthy), !is.na(logFC_Training_CR)) %>%
  mutate(reversed = sign(logFC_Cancer_vs_Healthy) != sign(logFC_Training_CR),
         cancer_sig = pi_score_Cancer_vs_Healthy < 0.05)

observed_r <- cor(dep_fc$logFC_Cancer_vs_Healthy, dep_fc$logFC_Training_CR,
                  method = "pearson")
n_proteins <- nrow(dep_fc)

cat(sprintf("\n=== DATA ===\nProteins: %d | Observed r: %.4f\n",
            n_proteins, observed_r))

# ═══════════════════════════════════════════════════════════════════════════════
# 1. ANALYTICAL STRUCTURAL NULL
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n═══ 1. ANALYTICAL STRUCTURAL CORRELATION ═══\n")

# Count samples per group
group_counts <- table(meta$Group_Time)
cat("Group sizes:\n"); print(group_counts)

# For a given protein with residual variance sigma^2:
#   logFC_CvH = mean(CR_T1) - mean(H_T1)
#   logFC_TR  = mean(CR_T2) - mean(CR_T1)
#
# Var(logFC_CvH) = sigma^2 * (1/n_CR_T1 + 1/n_H_T1)
# Var(logFC_TR)  = sigma^2 * (1/n_CR_T2 + 1/n_CR_T1)
# Cov(logFC_CvH, logFC_TR) = -sigma^2 / n_CR_T1
#
# Expected structural r = Cov / sqrt(Var1 * Var2)

# Our design pools CRE and PLA for the CR group
# CR_T1 = CRE_T1 + PLA_T1; CR_T2 = CRE_T2 + PLA_T2
n_cr_t1 <- sum(meta$Group_Time %in% c("CRE_T1", "PLA_T1"))
n_cr_t2 <- sum(meta$Group_Time %in% c("CRE_T2", "PLA_T2"))
n_h_t1  <- sum(meta$Group_Time == "H_T1")

cat(sprintf("  n(CR_T1) = %d, n(CR_T2) = %d, n(H_T1) = %d\n",
            n_cr_t1, n_cr_t2, n_h_t1))

# Analytical expected correlation (sigma^2 cancels)
var_cvh <- 1/n_cr_t1 + 1/n_h_t1
var_tr  <- 1/n_cr_t2 + 1/n_cr_t1
cov_structural <- -1/n_cr_t1

r_structural <- cov_structural / sqrt(var_cvh * var_tr)
cat(sprintf("  Var(logFC_CvH) ∝ %.4f\n", var_cvh))
cat(sprintf("  Var(logFC_TR)  ∝ %.4f\n", var_tr))
cat(sprintf("  Cov(structural) ∝ %.4f\n", cov_structural))
cat(sprintf("  Expected structural r = %.4f\n", r_structural))
cat(sprintf("  Observed r            = %.4f\n", observed_r))
cat(sprintf("  Excess r (biology)    = %.4f\n", observed_r - r_structural))

# Per-protein structural correlation from the actual limma fit
# Use the contrast covariance from the design matrix
meta_all <- meta
meta_all$Group_Time <- factor(meta_all$Group_Time,
                               levels = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1"))
design_all <- model.matrix(~ 0 + Group_Time, data = meta_all)
colnames(design_all) <- gsub("^Group_Time", "", colnames(design_all))

cm_both <- makeContrasts(
  Cancer_vs_Healthy = (CRE_T1 + PLA_T1) / 2 - H_T1,
  Training_CR = (CRE_T2 + PLA_T2) / 2 - (CRE_T1 + PLA_T1) / 2,
  levels = design_all
)

# The correlation between contrasts from the design matrix
# cov(contrasts) = t(cm) %*% solve(t(X) %*% X) %*% cm
XtX_inv <- solve(t(design_all) %*% design_all)
contrast_cov <- t(cm_both) %*% XtX_inv %*% cm_both
contrast_cor <- cov2cor(contrast_cov)
r_design <- contrast_cor[1, 2]
cat(sprintf("\n  Design-matrix contrast correlation: %.4f\n", r_design))
cat(sprintf("  (This is the exact structural r from limma's perspective)\n"))

# Bootstrap test: is observed r significantly more negative than structural r?
set.seed(42)
B <- 10000
# Simulate null: for each protein, generate logFC_CvH and logFC_TR from
# bivariate normal with the structural correlation but NO biological reversal
null_r <- replicate(B, {
  # Sample n_proteins from bivariate normal with r = r_design
  z1 <- rnorm(n_proteins)
  z2 <- rnorm(n_proteins)
  x1 <- z1
  x2 <- r_design * z1 + sqrt(1 - r_design^2) * z2
  cor(x1, x2)
})

p_excess <- (sum(null_r <= observed_r) + 1) / (B + 1)
cat(sprintf("  Null simulation p (observed r more negative than structural): %s\n",
            format.pval(p_excess, digits = 3)))

# ═══════════════════════════════════════════════════════════════════════════════
# 2. ROASTGSA MULTI-SCORE COMPARISON
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n═══ 2. ROASTGSA MULTI-SCORE COMPARISON ═══\n")

if (requireNamespace("roastgsa", quietly = TRUE)) {
  library(roastgsa)

  # Test 1: Cancer sets on Training_CR contrast (CR subjects only)
  cr_idx  <- which(meta$Group_Time != "H_T1")
  mat_cr  <- mat[, cr_idx]
  meta_cr <- meta[cr_idx, ]

  covar_cr <- data.frame(
    Group_Time = factor(meta_cr$Group_Time,
                        levels = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2")),
    Subject_ID = factor(meta_cr$Subject_ID),
    row.names = rownames(meta_cr)
  )

  # Cancer-sig gene sets (indices into mat_cr rows)
  sig_cancer <- dep_df %>%
    filter(pi_score_Cancer_vs_Healthy < 0.05, uniprot_id %in% rownames(mat_cr))
  idx_up <- which(rownames(mat_cr) %in%
                    sig_cancer$uniprot_id[sig_cancer$logFC_Cancer_vs_Healthy > 0])
  idx_dn <- which(rownames(mat_cr) %in%
                    sig_cancer$uniprot_id[sig_cancer$logFC_Cancer_vs_Healthy < 0])

  index_list <- list(cancer_up = idx_up, cancer_down = idx_dn)

  # Design for roastgsa (same as fry panel)
  design_cr <- model.matrix(~ 0 + Group_Time, data = covar_cr)
  colnames(design_cr) <- gsub("^Group_Time", "", colnames(design_cr))

  cm_tr <- makeContrasts(
    Training_CR = (CRE_T2 + PLA_T2) / 2 - (CRE_T1 + PLA_T1) / 2,
    levels = design_cr
  )

  scores <- c("mean", "maxmean", "absmean")
  roast_results <- list()

  for (sc in scores) {
    cat(sprintf("  Running roastgsa (score = %s)...\n", sc))
    tryCatch({
      rg <- roastgsa(mat_cr, covar = covar_cr,
                     form = as.formula("~ Group_Time"),
                     design = design_cr,
                     contrast = cm_tr[, "Training_CR"],
                     index = index_list,
                     set.statistic = sc,
                     self.contained = TRUE,
                     nrot = 999,
                     mccores = 1,
                     executation.info = FALSE)
      res <- rg$res
      res$score <- sc
      roast_results[[sc]] <- res
      cat(sprintf("    %s: cancer_up p=%.4f, cancer_down p=%.4f\n",
                  sc, res$pval[1], res$pval[2]))
    }, error = function(e) {
      cat(sprintf("    %s: ERROR — %s\n", sc, e$message))
    })
  }

  if (length(roast_results)) {
    roast_df <- bind_rows(roast_results) %>%
      mutate(set = rep(names(index_list), length(roast_results)))
    write_csv(roast_df, file.path(DAT, "roastgsa_multiscore.csv"))
    cat("\n  roastgsa results:\n")
    print(roast_df)

    # Effective signature size via varrotrand
    cat("\n  Computing effective signature sizes (varrotrand)...\n")
    tryCatch({
      rg_mean <- roastgsa(mat_cr, covar = covar_cr,
                          form = as.formula("~ Group_Time"),
                          design = design_cr,
                          contrast = cm_tr[, "Training_CR"],
                          index = index_list,
                          set.statistic = "mean",
                          self.contained = TRUE,
                          nrot = 999,
                          mccores = 1,
                          executation.info = FALSE)
      vr <- varrotrand(rg_mean, mat_cr,
                       testedsizes = c(seq(2, 20, by = 2), seq(25, 100, by = 10),
                                       seq(120, 400, by = 20)),
                       nrep = 50)
      # Save the varrotrand object for the effective size plot
      saveRDS(list(roastgsa_obj = rg_mean, varrot = vr),
              file.path(DAT, "roastgsa_varrot.rds"))
      cat("  varrotrand saved\n")
    }, error = function(e) {
      cat(sprintf("  varrotrand error: %s\n", e$message))
    })
  }
} else {
  cat("  roastgsa not installed — skipping\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# 3. CORRECTED CORRELATION
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n═══ 3. CORRECTED CORRELATION ═══\n")

# Approach: subtract structural r from observed r
r_corrected_simple <- observed_r - r_design
cat(sprintf("  Simple correction: %.4f - %.4f = %.4f\n",
            observed_r, r_design, r_corrected_simple))

# Cancer-DEP subset
dep_subset <- dep_fc %>% filter(cancer_sig)
r_dep <- cor(dep_subset$logFC_Cancer_vs_Healthy, dep_subset$logFC_Training_CR,
             method = "pearson")
r_dep_corrected <- r_dep - r_design
cat(sprintf("  Cancer-DEP:  observed = %.4f, corrected = %.4f\n",
            r_dep, r_dep_corrected))
cat(sprintf("  Non-DEP:     observed = %.4f\n",
            cor(dep_fc$logFC_Cancer_vs_Healthy[!dep_fc$cancer_sig],
                dep_fc$logFC_Training_CR[!dep_fc$cancer_sig],
                method = "pearson")))

# Fisher Z test: is corrected r significantly different from 0?
n_dep_n <- nrow(dep_subset)
fisher_z <- atanh(r_dep_corrected) * sqrt(n_dep_n - 3)
fisher_p <- 2 * pnorm(-abs(fisher_z))
cat(sprintf("  Fisher Z test (corrected r vs 0): z = %.3f, p = %s\n",
            fisher_z, format.pval(fisher_p, digits = 3)))

# ═══════════════════════════════════════════════════════════════════════════════
# SUMMARY TABLE
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n═══ SUMMARY ═══\n")

summary_df <- tibble(
  metric = c("observed_r_all", "observed_r_dep", "structural_r_analytical",
             "structural_r_design_matrix", "excess_r_all", "excess_r_dep",
             "null_sim_p", "fisher_z_corrected_p"),
  value = c(observed_r, r_dep, r_structural, r_design,
            observed_r - r_design, r_dep - r_design,
            p_excess, fisher_p),
  interpretation = c(
    "Pearson r across all proteins",
    "Pearson r within cancer-DEPs (Pi<0.05)",
    "Expected r from shared CR_T1 (simple formula)",
    "Expected r from limma contrast covariance (exact)",
    "Biology-only r (all proteins)",
    "Biology-only r (cancer-DEPs)",
    sprintf("Observed r more negative than structural? (n_sim=%d)", B),
    "Is corrected cancer-DEP r ≠ 0?"
  )
)
write_csv(summary_df, file.path(DAT, "structural_correction_summary.csv"))
print(summary_df)

# ═══════════════════════════════════════════════════════════════════════════════
# VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════════════

# Panel 1: Structural null simulation
p_struct <- ggplot(tibble(x = null_r), aes(x = x)) +
  geom_histogram(bins = 50, fill = "grey70", color = "grey40", linewidth = 0.2) +
  geom_vline(xintercept = r_design, color = "grey30", linewidth = 0.7,
             linetype = "dashed") +
  geom_vline(xintercept = observed_r, color = DIR_COLORS["Up"], linewidth = 0.8) +
  geom_vline(xintercept = r_dep, color = "#2563EB", linewidth = 0.8) +
  annotate("text", x = observed_r, y = Inf, vjust = 1.5, hjust = -0.1,
           label = sprintf("All r = %.3f", observed_r),
           color = DIR_COLORS["Up"], fontface = "bold", size = 3) +
  annotate("text", x = r_dep, y = Inf, vjust = 3.0, hjust = -0.1,
           label = sprintf("DEP r = %.3f", r_dep),
           color = "#2563EB", fontface = "bold", size = 3) +
  annotate("text", x = r_design, y = Inf, vjust = 4.5, hjust = -0.1,
           label = sprintf("Structural r = %.3f", r_design),
           color = "grey30", fontface = "bold", size = 3) +
  annotate("text", x = max(null_r), y = Inf, vjust = 1.5, hjust = 1,
           label = sprintf("p %s", fmt_p(p_excess)),
           color = "grey30", fontface = "bold.italic", size = 3.5) +
  labs(title = "Analytical Structural Null",
       subtitle = sprintf("Expected r from shared CR_T1 baseline = %.3f", r_design),
       x = "Simulated r (structural correlation only)", y = "Count") +
  FIG_THEME

# Panel 2: Corrected correlation decomposition
decomp_df <- tibble(
  group = factor(c("All proteins", "Cancer-DEPs", "Non-DEPs"),
                 levels = c("All proteins", "Cancer-DEPs", "Non-DEPs")),
  observed = c(observed_r,
               r_dep,
               cor(dep_fc$logFC_Cancer_vs_Healthy[!dep_fc$cancer_sig],
                   dep_fc$logFC_Training_CR[!dep_fc$cancer_sig])),
  structural = r_design,
  biological = c(observed_r - r_design, r_dep - r_design,
                 cor(dep_fc$logFC_Cancer_vs_Healthy[!dep_fc$cancer_sig],
                     dep_fc$logFC_Training_CR[!dep_fc$cancer_sig]) - r_design)
)

decomp_long <- decomp_df %>%
  pivot_longer(cols = c(structural, biological),
               names_to = "component", values_to = "r") %>%
  mutate(component = factor(component, levels = c("structural", "biological")))

p_decomp <- ggplot(decomp_long, aes(x = group, y = r, fill = component)) +
  geom_col(width = 0.55, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_point(data = decomp_df, aes(x = group, y = observed),
             inherit.aes = FALSE, shape = 18, size = 4, color = "black") +
  scale_fill_manual(values = c("structural" = "#94A3B8", "biological" = "#2563EB"),
                    labels = c("Structural (shared baseline)", "Biological (excess)"),
                    name = NULL) +
  geom_text(data = decomp_df, aes(x = group, y = observed,
                                    label = sprintf("r = %.3f", observed)),
            inherit.aes = FALSE, vjust = -0.8, size = 3, fontface = "bold") +
  labs(title = "Correlation Decomposition",
       subtitle = sprintf("Structural r = %.3f (from design matrix)", r_design),
       x = NULL, y = "Pearson r") +
  FIG_THEME +
  theme(legend.position = "bottom",
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 8))

# Panel 3: roastgsa score comparison (if available)
if (exists("roast_df") && nrow(roast_df) > 0) {
  roast_plot_df <- roast_df %>%
    mutate(neg_log_p = -log10(pmax(pval, 1e-10)),
           set = factor(set, levels = c("cancer_up", "cancer_down")),
           score = factor(score, levels = scores))

  p_roast <- ggplot(roast_plot_df, aes(x = score, y = neg_log_p, fill = set)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6,
             color = "black", linewidth = 0.2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               color = "grey40", linewidth = 0.3) +
    annotate("text", x = 3.4, y = -log10(0.05),
             label = "p = 0.05", hjust = 0, vjust = -0.5,
             size = 2.5, fontface = "italic", color = "grey40") +
    scale_fill_manual(values = c("cancer_up" = "#E57373", "cancer_down" = "#64B5F6"),
                      labels = c("Cancer-Up", "Cancer-Down"),
                      name = NULL) +
    labs(title = "roastgsa Score Comparison",
         subtitle = "Self-contained test: Cancer sets on Training_CR contrast",
         x = "Score function", y = expression(-log[10](p))) +
    FIG_THEME +
    theme(legend.position = "bottom",
          legend.key.size = unit(3, "mm"),
          legend.text = element_text(size = 8))
} else {
  p_roast <- ggplot() + theme_void() +
    annotate("text", x = 0.5, y = 0.5, label = "roastgsa not available",
             size = 4, color = "grey60")
}

# Combine
pilot_composite <- p_struct + p_decomp + p_roast +
  plot_layout(widths = c(1, 1, 1)) +
  plot_annotation(
    title = "Advanced Reversal Methods: Structural Correction & Score Comparison",
    subtitle = sprintf("Structural r = %.3f | Observed r = %.3f (all), %.3f (DEPs) | Excess r = %.3f (DEPs)",
                        r_design, observed_r, r_dep, r_dep - r_design),
    tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "grey30"),
      plot.tag = element_text(face = "bold", size = 13)))

ggsave(file.path(PILOT_DIR, "pilot_advanced_methods.png"), pilot_composite,
       width = 420, height = 140, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(PILOT_DIR, "pilot_advanced_methods.pdf"), pilot_composite,
       width = 420, height = 140, units = "mm", device = pdf_device)

cat(sprintf("\nPilot saved to: %s\n", PILOT_DIR))
