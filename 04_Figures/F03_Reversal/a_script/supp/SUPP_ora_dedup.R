# SUPP: ORA Deduplication Sensitivity
# Pathway count stability across Jaccard cutoffs for reversed quadrants
setwd(here::here())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")
pacman::p_load(tidyverse, fgsea)

RPT_PNG <- "04_Figures/F03_Reversal/b_reports/supp/png/panels"
RPT_PDF <- "04_Figures/F03_Reversal/b_reports/supp/pdf/panels"
DAT     <- "04_Figures/F03_Reversal/c_data/panel_supp"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT,     recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

# ── Data ─────────────────────────────────────────────────────────────────────
source("04_Figures/F03_Reversal/a_script/f03_data.R")
dep_df <- dep_df %>%
  transmute(gene,
            logFC_CvH = logFC_CRvH_Baseline,
            logFC_TR  = logFC_CR_Training) %>%
  filter(!is.na(logFC_CvH), !is.na(logFC_TR))

# Reversed quadrants: Cancer Up + Training Down, Cancer Down + Training Up
cancer_up_train_dn   <- dep_df %>% filter(logFC_CvH > 0, logFC_TR < 0) %>% pull(gene)
cancer_dn_train_up   <- dep_df %>% filter(logFC_CvH < 0, logFC_TR > 0) %>% pull(gene)
universe <- unique(dep_df$gene)

pathways <- build_pathway_collection(min_size = 15, max_size = 500,
                                     include_goslim = FALSE,
                                     exclude_variants = TRUE)

# ── Sweep Jaccard cutoffs ────────────────────────────────────────────────────
cutoffs <- c(0.3, 0.5, 0.7, 1.0)
quad_list <- list(
  `Cancer Up / Training Down` = cancer_up_train_dn,
  `Cancer Down / Training Up` = cancer_dn_train_up
)

sweep_res <- map_dfr(cutoffs, function(jc) {
  map_dfr(names(quad_list), function(q) {
    ora <- run_ora_deduplicated(quad_list[[q]], universe, pathways,
                                jaccard_cutoff = jc, padj_cutoff = 0.05)
    tibble(jaccard_cutoff = jc, quadrant = q, n_sig = nrow(ora))
  })
})

write_csv(sweep_res, file.path(DAT, "SUPP_ora_dedup_sensitivity.csv"))

# ── Plot ─────────────────────────────────────────────────────────────────────
QUAD_COLORS <- c(
  `Cancer Up / Training Down` = unname(DIR_COLORS["Up"]),
  `Cancer Down / Training Up` = unname(DIR_COLORS["Down"])
)

pS_ora_dedup <- ggplot(sweep_res,
                       aes(x = factor(jaccard_cutoff), y = n_sig,
                           fill = quadrant)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = n_sig),
            position = position_dodge(width = 0.7), vjust = -0.4, size = 3) +
  scale_fill_manual(values = QUAD_COLORS) +
  labs(title    = "ORA Deduplication Sensitivity",
       subtitle = "Significant pathways across Jaccard cutoffs",
       x = "Jaccard cutoff", y = "Significant pathways (padj < 0.05)",
       fill = "Quadrant") +
  FIG_THEME +
  theme(legend.position = "bottom")

ggsave(file.path(RPT_PNG, "SUPP_ora_dedup.png"), pS_ora_dedup,
       width = 140, height = 100, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_ora_dedup.pdf"), pS_ora_dedup,
       width = 140, height = 100, units = "mm", device = pdf_device)

message("Done: SUPP_ora_dedup")
