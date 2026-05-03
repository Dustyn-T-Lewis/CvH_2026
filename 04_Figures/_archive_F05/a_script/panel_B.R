# Figure 5 — Panel B: Supplement-Specific Reversal Scatter
# Among cancer-signature proteins (P.Value_CvH < 0.05):
# X = logFC_Training_CRE, Y = logFC_Training_PLA
# Color by Cancer_vs_Healthy direction. Identity line + correlation.
# Asks: "Does creatine enhance or hinder reversal compared to placebo?"
# Outputs: panel_B_supp_reversal.pdf/png, supplement_reversal.csv

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
})

PB_W <- 200; PB_H <- 200

RPT <- "04_Figures/F05/b_reports"
DAT <- "04_Figures/F05/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# --- Load DEP results from both models ---
dep_crvh <- read_csv("03_DEP/c_data/03_combined_results_CRvH.csv",
                      show_col_types = FALSE)
dep_cr   <- read_csv("03_DEP/c_data/03_combined_results_CR.csv",
                      show_col_types = FALSE)

# --- Identify cancer-signature proteins (P.Value < 0.05 in CRvH model) ---
cancer_sig <- dep_crvh |>
  filter(P.Value_Cancer_vs_Healthy < 0.05) |>
  select(gene, uniprot_id,
         logFC_CvH = logFC_Cancer_vs_Healthy,
         pval_CvH  = P.Value_Cancer_vs_Healthy)

cat(sprintf("Cancer-signature proteins (P.Value < 0.05): %d\n", nrow(cancer_sig)))

# --- Match with CR model by gene ---
cr_subset <- dep_cr |>
  select(gene, uniprot_id,
         logFC_Training_CRE, logFC_Training_PLA,
         P.Value_Training_CRE, P.Value_Training_PLA,
         logFC_Supplement_Interaction,
         P.Value_Supplement_Interaction)

merged <- cancer_sig |>
  inner_join(cr_subset, by = c("gene", "uniprot_id")) |>
  distinct(uniprot_id, .keep_all = TRUE)

cat(sprintf("Matched across models: %d proteins\n", nrow(merged)))

# --- Classify by cancer direction ---
merged <- merged |>
  mutate(
    cancer_dir = ifelse(logFC_CvH > 0, "Cancer Up", "Cancer Down"),
    cancer_dir = factor(cancer_dir, levels = c("Cancer Up", "Cancer Down"))
  )

# --- Correlation statistics ---
cor_test <- cor.test(merged$logFC_Training_CRE, merged$logFC_Training_PLA,
                     method = "spearman")
cor_label <- sprintf("rho = %.2f\np %s  (n = %d)",
                     cor_test$estimate,
                     fmt_p(cor_test$p.value),
                     nrow(merged))

# --- Quadrant labels ---
q_labels <- merged |>
  mutate(
    quad = case_when(
      logFC_Training_CRE > 0 & logFC_Training_PLA > 0 ~ "Both Up",
      logFC_Training_CRE < 0 & logFC_Training_PLA < 0 ~ "Both Down",
      logFC_Training_CRE > 0 & logFC_Training_PLA < 0 ~ "CRE Up\nPLA Down",
      logFC_Training_CRE < 0 & logFC_Training_PLA > 0 ~ "CRE Down\nPLA Up"
    )
  ) |>
  count(quad)

# --- Axis limits ---
ax_max <- max(abs(c(merged$logFC_Training_CRE, merged$logFC_Training_PLA)),
              na.rm = TRUE) * 1.15
ax_lim <- c(-ax_max, ax_max)

# --- Top labels: proteins far from identity line (potential supplement effect) ---
merged <- merged |>
  mutate(
    delta_from_identity = logFC_Training_CRE - logFC_Training_PLA,
    abs_delta = abs(delta_from_identity)
  )

top_diff <- merged |>
  slice_max(abs_delta, n = 12, with_ties = FALSE)

# --- Quadrant annotation positions ---
quad_annot <- data.frame(
  x = c(ax_max * 0.9, -ax_max * 0.9, -ax_max * 0.9, ax_max * 0.9),
  y = c(ax_max * 0.9, ax_max * 0.9, -ax_max * 0.9, -ax_max * 0.9),
  hjust = c(1, 0, 0, 1),
  vjust = c(1, 1, 0, 0),
  label = vapply(c("Both Up", "CRE Down\nPLA Up", "Both Down", "CRE Up\nPLA Down"),
                 function(q) {
                   n <- q_labels$n[q_labels$quad == q]
                   if (length(n) == 0) n <- 0L else n <- n[1]
                   sprintf("%s\nn=%d", q, n)
                 }, character(1))
)

# --- Plot ---
pB <- ggplot(merged, aes(x = logFC_Training_CRE, y = logFC_Training_PLA)) +
  # Guidelines
  geom_abline(slope = 1, intercept = 0, linetype = "solid",
              color = "grey60", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  # Points
  geom_point(aes(color = cancer_dir), alpha = 0.5, size = 1.5) +
  # Labels for proteins with largest CRE-PLA divergence
  geom_text_repel(data = top_diff,
                  aes(label = gene, color = cancer_dir),
                  size = 2.3, max.overlaps = 15, seed = 42,
                  segment.size = 0.3, segment.color = "grey50",
                  fontface = "italic",
                  box.padding = 0.4, point.padding = 0.3,
                  show.legend = FALSE) +
  # Quadrant annotations
  geom_text(data = quad_annot,
            aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
            size = 2.5, color = "grey40", fontface = "bold",
            inherit.aes = FALSE) +
  # Correlation annotation
  annotate("text", x = -ax_max * 0.95, y = ax_max * 0.95,
           label = cor_label, hjust = 0, vjust = 1,
           size = 2.8, color = "black", fontface = "bold") +
  # Identity line annotation
  annotate("text", x = ax_max * 0.7, y = ax_max * 0.55,
           label = "y = x\n(equal response)", hjust = 0, vjust = 0,
           size = 2.2, color = "grey50", fontface = "italic") +
  scale_color_manual(values = CANCER_DIR_COLORS,
                     name = "Cancer\nsignature") +
  coord_cartesian(xlim = ax_lim, ylim = ax_lim) +
  labs(title = "Supplement-Specific Training Reversal",
       subtitle = "Cancer-signature proteins: creatine vs placebo training response",
       x = paste0("logFC — ", CTR_SHORT["Training_CRE"]),
       y = paste0("logFC — ", CTR_SHORT["Training_PLA"]),
       tag = "B") +
  FIG_THEME +
  theme(legend.position = "right",
        aspect.ratio = 1)

# --- Save ---
ggsave(file.path(RPT, "panel_B_supp_reversal.pdf"), pB,
       width = PB_W, height = PB_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_B_supp_reversal.png"), pB,
       width = PB_W, height = PB_H, units = "mm", dpi = 300)

# --- Export ---
export_df <- merged |>
  select(gene, uniprot_id, logFC_CvH, pval_CvH, cancer_dir,
         logFC_Training_CRE, logFC_Training_PLA,
         P.Value_Training_CRE, P.Value_Training_PLA,
         logFC_Supplement_Interaction, P.Value_Supplement_Interaction,
         delta_from_identity) |>
  arrange(desc(abs(delta_from_identity)))

write.csv(as.data.frame(export_df),
          file.path(DAT, "supplement_reversal.csv"), row.names = FALSE)

cat(sprintf("Panel B (supplement reversal scatter) done: %d proteins.\n",
            nrow(merged)))
