# F04 CRvH Supplementary: Directional Asymmetry Analysis
# Q1: Is concordant-up overlap stronger than concordant-down?
# Q2: Are specific biological processes asymmetrically conserved?
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F04/a_script/style.R")

library(tidyverse)
library(patchwork)

RPT <- "04_Figures/F04/CRvH/b_reports/supp"
DAT <- "04_Figures/F04/CRvH/c_data/supp"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

dep_df <- readr::read_csv("03_DEP/a_non_imputed/c_data/combined_results_pi.csv",
                          show_col_types = FALSE) |>
  dplyr::mutate(contrast = dplyr::recode(contrast,
                                         CRvH_Baseline = "Cancer_vs_Healthy",
                                         CR_Training   = "Training_CR")) |>
  tidyr::pivot_wider(id_cols = c(uniprot_id, gene, protein, description),
                     names_from = contrast,
                     values_from = c(logFC, t, P.Value, adj.P.Val, pi_score, sig_pi),
                     names_glue = "{.value}_{contrast}")

rr <- dep_df %>%
  transmute(gene,
            t_cvh = t_Cancer_vs_Healthy, t_tr = t_Training_CR,
            lfc_cvh = logFC_Cancer_vs_Healthy, lfc_tr = logFC_Training_CR,
            pi_cvh = pi_score_Cancer_vs_Healthy, pi_tr = pi_score_Training_CR,
            sig_cvh = sig_pi_Cancer_vs_Healthy, sig_tr = sig_pi_Training_CR) %>%
  filter(!is.na(t_cvh), !is.na(t_tr))

n_total <- nrow(rr)

# -- 1. Protein-level directional overlap --

# Define directional gene sets by sign of t-statistic
cvh_up   <- rr$gene[rr$t_cvh > 0]
cvh_down <- rr$gene[rr$t_cvh < 0]
tr_up    <- rr$gene[rr$t_tr > 0]
tr_down  <- rr$gene[rr$t_tr < 0]

# Concordant overlap
conc_up   <- intersect(cvh_up, tr_up)      # both increase
conc_down <- intersect(cvh_down, tr_down)   # both decrease

# Expected by chance
exp_conc_up   <- length(cvh_up) * length(tr_up) / n_total
exp_conc_down <- length(cvh_down) * length(tr_down) / n_total

# Hypergeometric enrichment (one-sided, greater than expected)
p_conc_up <- phyper(length(conc_up) - 1, length(tr_up),
                    n_total - length(tr_up), length(cvh_up), lower.tail = FALSE)
p_conc_down <- phyper(length(conc_down) - 1, length(tr_down),
                      n_total - length(tr_down), length(cvh_down), lower.tail = FALSE)

# Enrichment ratios
or_up   <- length(conc_up) / exp_conc_up
or_down <- length(conc_down) / exp_conc_down

# Fisher's exact test: is concordance asymmetric between up and down?
fisher_tbl <- matrix(c(
  length(conc_up),
  length(cvh_up) - length(conc_up),
  length(conc_down),
  length(cvh_down) - length(conc_down)
), nrow = 2, byrow = TRUE,
dimnames = list(c("Up_direction", "Down_direction"),
                c("Concordant", "Discordant")))

fisher_res <- fisher.test(fisher_tbl)

cat("=== Protein-level directional asymmetry (F04 CRvH) ===\n")
cat(sprintf("Concordant Up:   %d / %d expected (OR = %.3f, p = %.2e)\n",
            length(conc_up), round(exp_conc_up), or_up, p_conc_up))
cat(sprintf("Concordant Down: %d / %d expected (OR = %.3f, p = %.2e)\n",
            length(conc_down), round(exp_conc_down), or_down, p_conc_down))
cat(sprintf("Fisher asymmetry test: OR = %.3f, p = %.4f\n",
            fisher_res$estimate, fisher_res$p.value))
cat(sprintf("  Interpretation: %s\n",
            ifelse(fisher_res$p.value < 0.05,
                   ifelse(fisher_res$estimate > 1,
                          "Up-direction MORE concordant than down",
                          "Down-direction MORE concordant than up"),
                   "No significant directional asymmetry")))

# -- 2. Magnitude-weighted asymmetry --

# For significant proteins only (Pi < 0.05 in either contrast)
sig_either <- rr %>% filter(abs(sig_cvh) > 0 | abs(sig_tr) > 0)
sig_conc_up   <- sig_either %>% filter(t_cvh > 0, t_tr > 0)
sig_conc_down <- sig_either %>% filter(t_cvh < 0, t_tr < 0)
sig_disc      <- sig_either %>% filter(sign(t_cvh) != sign(t_tr))

cat(sprintf("\nSignificant proteins (Pi < 0.05 in either):\n"))
cat(sprintf("  Concordant Up:   %d (mean |t_cvh| = %.2f, mean |t_tr| = %.2f)\n",
            nrow(sig_conc_up), mean(abs(sig_conc_up$t_cvh)), mean(abs(sig_conc_up$t_tr))))
cat(sprintf("  Concordant Down: %d (mean |t_cvh| = %.2f, mean |t_tr| = %.2f)\n",
            nrow(sig_conc_down), mean(abs(sig_conc_down$t_cvh)), mean(abs(sig_conc_down$t_tr))))
cat(sprintf("  Discordant:      %d\n", nrow(sig_disc)))

# Wilcoxon: is |t_tr| larger for concordant-up vs concordant-down proteins?
if (nrow(sig_conc_up) > 5 && nrow(sig_conc_down) > 5) {
  w_test <- wilcox.test(abs(sig_conc_up$t_tr), abs(sig_conc_down$t_tr))
  cat(sprintf("  Wilcoxon |t_TR|: up vs down p = %.4f (median %.2f vs %.2f)\n",
              w_test$p.value, median(abs(sig_conc_up$t_tr)), median(abs(sig_conc_down$t_tr))))
}

# -- 3. Pathway-level directional asymmetry --

# Load fGSEA results
fgsea_all <- read_csv(file.path("04_Figures/F04/CRvH/c_data/shared", "fgsea_tstat_CRvH.csv"),
                      show_col_types = FALSE)

fgsea_wide <- fgsea_all %>%
  filter(contrast %in% c("Cancer_vs_Healthy", "Training_CR")) %>%
  select(pathway, contrast, NES, padj) %>%
  pivot_wider(names_from = contrast, values_from = c(NES, padj),
              names_glue = "{.value}_{contrast}") %>%
  filter(!is.na(NES_Cancer_vs_Healthy), !is.na(NES_Training_CR))

# Classify pathways by direction
pw_class <- fgsea_wide %>%
  mutate(
    sig_cvh = padj_Cancer_vs_Healthy < 0.05,
    sig_tr  = padj_Training_CR < 0.05,
    conc_up   = NES_Cancer_vs_Healthy > 0 & NES_Training_CR > 0,
    conc_down = NES_Cancer_vs_Healthy < 0 & NES_Training_CR < 0,
    concordant = conc_up | conc_down,
    sig_both  = sig_cvh & sig_tr
  )

cat("\n=== Pathway-level directional asymmetry ===\n")
cat(sprintf("Total pathways: %d\n", nrow(pw_class)))
cat(sprintf("Concordant Up:   %d (%d sig both)\n",
            sum(pw_class$conc_up), sum(pw_class$conc_up & pw_class$sig_both)))
cat(sprintf("Concordant Down: %d (%d sig both)\n",
            sum(pw_class$conc_down), sum(pw_class$conc_down & pw_class$sig_both)))

# Compare NES magnitude between concordant-up and concordant-down pathways
cu <- pw_class %>% filter(conc_up)
cd <- pw_class %>% filter(conc_down)

cat(sprintf("Mean |NES_CvH| -- Up: %.2f, Down: %.2f\n",
            mean(abs(cu$NES_Cancer_vs_Healthy)), mean(abs(cd$NES_Cancer_vs_Healthy))))
cat(sprintf("Mean |NES_TR|  -- Up: %.2f, Down: %.2f\n",
            mean(abs(cu$NES_Training_CR)), mean(abs(cd$NES_Training_CR))))

# Concordance ratio: |NES_TR / NES_CvH| (closer to 1 = more concordant)
cu_ratio <- abs(cu$NES_Training_CR / cu$NES_Cancer_vs_Healthy)
cd_ratio <- abs(cd$NES_Training_CR / cd$NES_Cancer_vs_Healthy)

w_ratio <- wilcox.test(cu_ratio, cd_ratio)
cat(sprintf("\nConcordance ratio |NES_TR/NES_CvH| -- Up: %.2f, Down: %.2f (Wilcoxon p = %.4f)\n",
            median(cu_ratio, na.rm = TRUE), median(cd_ratio, na.rm = TRUE), w_ratio$p.value))
cat(sprintf("  Interpretation: %s direction is more concordant\n",
            ifelse(median(cu_ratio, na.rm = TRUE) > median(cd_ratio, na.rm = TRUE),
                   "Up", "Down")))

# -- 4. Visualization: Asymmetry summary panel --

# Panel A: Protein-level enrichment comparison
protein_df <- tibble(
  direction = c("Concordant Up\n(both increase)", "Concordant Down\n(both decrease)"),
  observed  = c(length(conc_up), length(conc_down)),
  expected  = c(exp_conc_up, exp_conc_down),
  OR        = c(or_up, or_down),
  p_val     = c(p_conc_up, p_conc_down)
)

pA <- ggplot(protein_df, aes(x = direction)) +
  geom_col(aes(y = observed, fill = direction), width = 0.6, alpha = 0.85) +
  geom_point(aes(y = expected), shape = 4, size = 4, stroke = 1.2) +
  geom_text(aes(y = observed + 15,
                label = sprintf("OR = %.2f\np = %.1e", OR, p_val)),
            size = 2.8, lineheight = 0.85) +
  scale_fill_manual(values = c("Concordant Up\n(both increase)" = unname(DIR_COLORS["Up"]),
                                "Concordant Down\n(both decrease)" = unname(DIR_COLORS["Down"])),
                    guide = "none") +
  labs(title = "Protein-Level Overlap",
       subtitle = sprintf("Fisher asymmetry: OR = %.2f, p = %.3f",
                          fisher_res$estimate, fisher_res$p.value),
       y = "Proteins with concordant sign", x = NULL) +
  FIG_THEME

# Panel B: Pathway-level concordance ratio
ratio_df <- bind_rows(
  tibble(direction = "Up", ratio = cu_ratio),
  tibble(direction = "Down", ratio = cd_ratio)
) %>% filter(is.finite(ratio), ratio < 5)  # cap outliers

pB <- ggplot(ratio_df, aes(x = direction, y = ratio, fill = direction)) +
  geom_boxplot(width = 0.5, alpha = 0.7, outlier.size = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("Up" = unname(DIR_COLORS["Up"]), "Down" = unname(DIR_COLORS["Down"])),
                    guide = "none") +
  annotate("text", x = 1.5, y = max(ratio_df$ratio) * 0.95,
           label = sprintf("Wilcoxon p = %.3f", w_ratio$p.value),
           size = 3) +
  labs(title = "Pathway Concordance Ratio",
       subtitle = "|NES_TR / NES_CvH| (1 = perfectly concordant)",
       y = "Concordance ratio", x = NULL) +
  FIG_THEME

# Panel C: Scatter of NES CvH vs TR colored by direction
pC <- ggplot(pw_class, aes(x = NES_Cancer_vs_Healthy, y = NES_Training_CR)) +
  annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
           fill = DIR_COLORS["Up"], alpha = 0.08) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = DIR_COLORS["Down"], alpha = 0.08) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = case_when(
    conc_up & sig_both ~ "Conc. Up (sig both)",
    conc_down & sig_both ~ "Conc. Down (sig both)",
    concordant ~ "Concordant (NS)",
    TRUE ~ "Discordant"
  ), size = -log10(pmin(padj_Cancer_vs_Healthy, padj_Training_CR))),
  alpha = 0.7) +
  scale_color_manual(values = c(
    "Conc. Up (sig both)" = unname(DIR_COLORS["Up"]),
    "Conc. Down (sig both)" = unname(DIR_COLORS["Down"]),
    "Concordant (NS)" = "grey60",
    "Discordant" = "grey80"
  ), name = NULL) +
  scale_size_continuous(range = c(1, 4), guide = "none") +
  labs(title = "Pathway NES: Cancer vs Healthy vs Training CR",
       subtitle = "Identity line = perfectly concordant; deviation = divergence",
       x = "NES (Cancer vs Healthy)", y = "NES (Training CR)") +
  coord_fixed() +
  FIG_THEME +
  theme(legend.position = "bottom", legend.text = element_text(size = 7))

p_composite <- pA + pB + pC + plot_layout(widths = c(1, 1, 1.5)) +
  plot_annotation(
    title = "Directional Asymmetry of Cancer Recovery Concordance",
    subtitle = sprintf("%d proteins | limma + dupCor | %d pathways",
                        n_total, nrow(pw_class)),
    theme = theme(plot.title = element_text(size = 12, face = "bold"),
                  plot.subtitle = element_text(size = 9, color = "grey30"))
  )

ggsave(file.path(RPT, "supp_directional_asymmetry_SUPP.pdf"), p_composite,
       width = 320, height = 130, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "supp_directional_asymmetry_SUPP.png"), p_composite,
       width = 320, height = 130, units = "mm", dpi = 300)

# -- 5. Export audit data --

write_csv(protein_df, file.path(DAT, "asymmetry_protein_level.csv"))
write_csv(
  pw_class %>% select(pathway, NES_Cancer_vs_Healthy, NES_Training_CR,
                      padj_Cancer_vs_Healthy, padj_Training_CR,
                      conc_up, conc_down, concordant, sig_both),
  file.path(DAT, "asymmetry_pathway_level.csv")
)

message("F04 CRvH directional asymmetry analysis done")
