# Pilot: Cross-project structural artifact comparison
# Runs the analytical structural null on YvO_2025 data (F05 reversal)
# to cross-validate the CvH finding.
#
# YvO design:
#   Aging         = Old_Pre - Young_Pre  (uses Old_Pre as "disease" baseline)
#   Training_Old  = Old_Post - Old_Pre   (uses Old_Pre as pre-training)
#   Shared term: Old_Pre appears with +1 in Aging and -1 in Training_Old
#   → same structural negative correlation as CvH
#
# CvH design:
#   Cancer_vs_Healthy = CR_T1 - H_T1
#   Training_CR       = CR_T2 - CR_T1
#   Shared: CR_T1
#
# Question: does the structural artifact dominate in YvO too?

library(tidyverse)
library(limma)

cat("\n════════════════════════════════════════════════════════════\n")
cat("  CROSS-PROJECT STRUCTURAL ARTIFACT COMPARISON\n")
cat("════════════════════════════════════════════════════════════\n")

# ═══════════════════════════════════════════════════════════════════════════════
# YvO DATA
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n──── YvO_2025 (Aging × Training_Old) ────\n")

yvo_dal <- readRDS("/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_YvO_2025/02_Imputation/c_data/01_DAList_imputed.rds")
yvo_dep <- read_csv("/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_YvO_2025/03_DEP/c_data/03_combined_results.csv",
                     show_col_types = FALSE)

yvo_mat  <- yvo_dal$data
yvo_meta <- as.data.frame(yvo_dal$metadata)

cat(sprintf("  Proteins: %d | Samples: %d\n", nrow(yvo_mat), ncol(yvo_mat)))
cat("  Group sizes:\n")
print(table(yvo_meta$Group, yvo_meta$Timepoint))

# Group sizes for structural r calculation
# Aging = Old_Pre - Young_Pre
# Training_Old = Old_Post - Old_Pre
n_old_pre   <- sum(yvo_meta$Group == "Old" & yvo_meta$Timepoint == "Pre")
n_young_pre <- sum(yvo_meta$Group == "Young" & yvo_meta$Timepoint == "Pre")
n_old_post  <- sum(yvo_meta$Group == "Old" & yvo_meta$Timepoint == "Post")

cat(sprintf("  n(Old_Pre) = %d, n(Young_Pre) = %d, n(Old_Post) = %d\n",
            n_old_pre, n_young_pre, n_old_post))

# Analytical structural r
var_aging <- 1/n_old_pre + 1/n_young_pre
var_tr_old <- 1/n_old_post + 1/n_old_pre
cov_yvo <- -1/n_old_pre
r_structural_yvo <- cov_yvo / sqrt(var_aging * var_tr_old)

# Design matrix approach (exact)
yvo_meta$Group_Time <- factor(paste0(yvo_meta$Group, "_", yvo_meta$Timepoint))
design_yvo <- model.matrix(~ 0 + Group_Time, data = yvo_meta)
colnames(design_yvo) <- gsub("^Group_Time", "", colnames(design_yvo))

cm_yvo <- makeContrasts(
  Aging = Old_Pre - Young_Pre,
  Training_Old = Old_Post - Old_Pre,
  levels = design_yvo
)

XtX_inv_yvo <- solve(t(design_yvo) %*% design_yvo)
contrast_cov_yvo <- t(cm_yvo) %*% XtX_inv_yvo %*% cm_yvo
contrast_cor_yvo <- cov2cor(contrast_cov_yvo)
r_design_yvo <- contrast_cor_yvo["Aging", "Training_Old"]

# Observed r
yvo_fc <- yvo_dep %>%
  filter(!is.na(logFC_Aging), !is.na(logFC_Training_Old)) %>%
  mutate(reversed = sign(logFC_Aging) != sign(logFC_Training_Old),
         aging_sig = pi_score_Aging < 0.05)

r_obs_yvo <- cor(yvo_fc$logFC_Aging, yvo_fc$logFC_Training_Old, method = "pearson")

yvo_dep_sub <- yvo_fc %>% filter(aging_sig)
r_dep_yvo <- cor(yvo_dep_sub$logFC_Aging, yvo_dep_sub$logFC_Training_Old, method = "pearson")

pct_rev_yvo <- 100 * mean(yvo_dep_sub$reversed)

cat(sprintf("\n  Structural r (simple):   %.4f\n", r_structural_yvo))
cat(sprintf("  Structural r (design):   %.4f\n", r_design_yvo))
cat(sprintf("  Observed r (all %d):     %.4f\n", nrow(yvo_fc), r_obs_yvo))
cat(sprintf("  Observed r (DEPs %d):    %.4f\n", nrow(yvo_dep_sub), r_dep_yvo))
cat(sprintf("  Excess r (all):          %.4f\n", r_obs_yvo - r_design_yvo))
cat(sprintf("  Excess r (DEPs):         %.4f\n", r_dep_yvo - r_design_yvo))
cat(sprintf("  %% aging-DEPs reversed:   %.1f%%\n", pct_rev_yvo))

# ═══════════════════════════════════════════════════════════════════════════════
# CvH DATA (recalculate for side-by-side)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n──── CvH_2026 (Cancer × Training_CR) ────\n")

cvh_dal <- readRDS("/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026/02_Imputation/c_data/01_DAList_imputed.rds")
cvh_dep <- read_csv("/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026/03_DEP/c_data/03_combined_results_CRvH.csv",
                     show_col_types = FALSE)

cvh_mat  <- cvh_dal$data
cvh_meta <- as.data.frame(cvh_dal$metadata)

n_cr_t1 <- sum(cvh_meta$Group_Time %in% c("CRE_T1", "PLA_T1"))
n_cr_t2 <- sum(cvh_meta$Group_Time %in% c("CRE_T2", "PLA_T2"))
n_h_t1  <- sum(cvh_meta$Group_Time == "H_T1")

cat(sprintf("  Proteins: %d | Samples: %d\n", nrow(cvh_mat), ncol(cvh_mat)))
cat(sprintf("  n(CR_T1) = %d, n(CR_T2) = %d, n(H_T1) = %d\n",
            n_cr_t1, n_cr_t2, n_h_t1))

cvh_meta$Group_Time <- factor(cvh_meta$Group_Time,
                               levels = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1"))
design_cvh <- model.matrix(~ 0 + Group_Time, data = cvh_meta)
colnames(design_cvh) <- gsub("^Group_Time", "", colnames(design_cvh))

cm_cvh <- makeContrasts(
  Cancer_vs_Healthy = (CRE_T1 + PLA_T1) / 2 - H_T1,
  Training_CR = (CRE_T2 + PLA_T2) / 2 - (CRE_T1 + PLA_T1) / 2,
  levels = design_cvh
)

XtX_inv_cvh <- solve(t(design_cvh) %*% design_cvh)
contrast_cov_cvh <- t(cm_cvh) %*% XtX_inv_cvh %*% cm_cvh
contrast_cor_cvh <- cov2cor(contrast_cov_cvh)
r_design_cvh <- contrast_cor_cvh[1, 2]

cvh_fc <- cvh_dep %>%
  filter(!is.na(logFC_Cancer_vs_Healthy), !is.na(logFC_Training_CR)) %>%
  mutate(reversed = sign(logFC_Cancer_vs_Healthy) != sign(logFC_Training_CR),
         cancer_sig = pi_score_Cancer_vs_Healthy < 0.05)

r_obs_cvh <- cor(cvh_fc$logFC_Cancer_vs_Healthy, cvh_fc$logFC_Training_CR, method = "pearson")
cvh_dep_sub <- cvh_fc %>% filter(cancer_sig)
r_dep_cvh <- cor(cvh_dep_sub$logFC_Cancer_vs_Healthy, cvh_dep_sub$logFC_Training_CR, method = "pearson")
pct_rev_cvh <- 100 * mean(cvh_dep_sub$reversed)

cat(sprintf("\n  Structural r (design):   %.4f\n", r_design_cvh))
cat(sprintf("  Observed r (all %d):     %.4f\n", nrow(cvh_fc), r_obs_cvh))
cat(sprintf("  Observed r (DEPs %d):    %.4f\n", nrow(cvh_dep_sub), r_dep_cvh))
cat(sprintf("  Excess r (all):          %.4f\n", r_obs_cvh - r_design_cvh))
cat(sprintf("  Excess r (DEPs):         %.4f\n", r_dep_cvh - r_design_cvh))
cat(sprintf("  %% disease-DEPs reversed: %.1f%%\n", pct_rev_cvh))

# ═══════════════════════════════════════════════════════════════════════════════
# COMPARISON TABLE
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n──── SIDE-BY-SIDE COMPARISON ────\n")

comp_df <- tibble(
  metric = c("Proteins", "Shared baseline n",
             "Structural r (design)", "Observed r (all)", "Observed r (DEPs)",
             "Excess r (all)", "Excess r (DEPs)", "% DEPs reversed",
             "n DEPs"),
  YvO_2025 = c(nrow(yvo_fc), n_old_pre,
               round(r_design_yvo, 4), round(r_obs_yvo, 4), round(r_dep_yvo, 4),
               round(r_obs_yvo - r_design_yvo, 4), round(r_dep_yvo - r_design_yvo, 4),
               round(pct_rev_yvo, 1), nrow(yvo_dep_sub)),
  CvH_2026 = c(nrow(cvh_fc), n_cr_t1,
               round(r_design_cvh, 4), round(r_obs_cvh, 4), round(r_dep_cvh, 4),
               round(r_obs_cvh - r_design_cvh, 4), round(r_dep_cvh - r_design_cvh, 4),
               round(pct_rev_cvh, 1), nrow(cvh_dep_sub))
)

print(comp_df, n = 20)

write_csv(comp_df,
          "/Users/dtl0018/Desktop/A_Proteomics_Analysis/A_CvH_2026/04_Figures/Reversal/c_data/pilot_advanced/structural_comparison_yvo_cvh.csv")

cat("\n══ INTERPRETATION ══\n")
cat(sprintf("  YvO: structural r = %.3f, observed r = %.3f → excess = %.3f (%s)\n",
            r_design_yvo, r_obs_yvo, r_obs_yvo - r_design_yvo,
            ifelse(r_obs_yvo > r_design_yvo,
                   "LESS negative than structural = biological concordance signal",
                   "MORE negative than structural = biological reversal signal")))
cat(sprintf("  CvH: structural r = %.3f, observed r = %.3f → excess = %.3f (%s)\n",
            r_design_cvh, r_obs_cvh, r_obs_cvh - r_design_cvh,
            ifelse(r_obs_cvh > r_design_cvh,
                   "LESS negative than structural = biological concordance signal",
                   "MORE negative than structural = biological reversal signal")))
