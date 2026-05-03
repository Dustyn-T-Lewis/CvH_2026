# F04 CRvH Supp: NES Scatter with collapsePathways (H + KEGG + Reactome + GO:BP)
# Detail-tier databases with fgsea::collapsePathways() to select independent terms.
# Contrast with main panel_D.R which uses the a priori GO Slim + Hallmark collection.
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F04/a_script/style.R")

library(tidyverse)
library(ggrepel)
library(fgsea)

PG_W <- 200
RPT <- "04_Figures/F04/CRvH/b_reports/supp"
DAT <- "04_Figures/F04/CRvH/c_data/supp"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# --- 1. Load DEP results and build ranked lists
dep_df <- read_csv("03_DEP/c_data/03_combined_results_CRvH.csv", show_col_types = FALSE)

contrasts <- c("Cancer_vs_Healthy", "Training_CR")
stats_list <- setNames(lapply(contrasts, function(ctr) {
  col <- paste0("t_", ctr)
  s <- setNames(dep_df[[col]], dep_df$gene)
  s[!is.na(s)]
}), contrasts)

# --- 2. Build detail-tier pathway collection (no GO Slim)
pw_list <- build_pathway_collection(
  species        = "Homo sapiens",
  min_size       = 15,
  max_size       = 500,
  include_goslim = FALSE,
  exclude_variants = TRUE
)

# --- 3. Run fgseaMultilevel + collapsePathways per contrast
set.seed(42)

fgsea_results <- list()
independent_pw <- list()

for (ctr in contrasts) {
  message(sprintf("\n--- %s ---", ctr))
  ranks <- stats_list[[ctr]]

  res_dt <- fgseaMultilevel(
    pathways    = pw_list,
    stats       = ranks,
    minSize     = 15,
    maxSize     = 500,
    nPermSimple = 10000,
    eps         = 0
  )

  # Classify database
  res_df <- as.data.frame(res_dt)
  res_df$database <- classify_database(res_df$pathway)

  # collapsePathways on significant terms
  sig_dt <- res_dt[!is.na(res_dt$padj) & res_dt$padj < 0.05, ]
  if (nrow(sig_dt) > 0) {
    collapsed <- collapsePathways(
      fgseaRes = sig_dt,
      pathways = pw_list,
      stats    = ranks
    )
    indep <- collapsed$mainPathways
    message(sprintf("collapsePathways: %d sig -> %d independent",
                    nrow(sig_dt), length(indep)))
    independent_pw[[ctr]] <- indep
  } else {
    message("No significant terms for collapsePathways")
    independent_pw[[ctr]] <- character(0)
  }

  fgsea_results[[ctr]] <- tibble::as_tibble(res_df) %>%
    mutate(contrast = ctr)
}

# --- 4. Take union of independent pathways across contrasts
union_pw <- unique(unlist(independent_pw))
message(sprintf("\nUnion of independent pathways: %d (CvH: %d, TR: %d)",
                length(union_pw),
                length(independent_pw[["Cancer_vs_Healthy"]]),
                length(independent_pw[["Training_CR"]])))

# --- 5. Pivot to wide format (union pathways only)
fgsea_long <- bind_rows(fgsea_results) %>%
  filter(pathway %in% union_pw)

fgsea_wide <- fgsea_long %>%
  dplyr::select(pathway, contrast, NES, padj, size, database) %>%
  pivot_wider(id_cols = c(pathway, database), names_from = contrast,
              values_from = c(NES, padj, size)) %>%
  filter(!is.na(NES_Cancer_vs_Healthy), !is.na(NES_Training_CR)) %>%
  mutate(set_size = coalesce(size_Cancer_vs_Healthy, size_Training_CR))

# Classify significance
fgsea_wide <- fgsea_wide %>%
  mutate(
    sig_CvH = !is.na(padj_Cancer_vs_Healthy) & padj_Cancer_vs_Healthy < 0.05,
    sig_TR  = !is.na(padj_Training_CR)       & padj_Training_CR < 0.05,
    significance = case_when(
      sig_CvH & sig_TR ~ "Sig Both",
      sig_CvH          ~ "Sig Cancer only",
      sig_TR           ~ "Sig Training only",
      TRUE             ~ "NS"
    ) %>% factor(levels = names(SIG_COLORS_F4)),
    pathway_label = clean_pathway_name(pathway)
  )

fgsea_sig <- fgsea_wide %>% filter(significance != "NS")

message(sprintf("  %d total pathways | %d significant | by database: %s",
                nrow(fgsea_wide), nrow(fgsea_sig),
                paste(names(table(fgsea_wide$database)),
                      table(fgsea_wide$database), sep = ":", collapse = ", ")))

# --- 6. Spearman correlation
nes_cor_all <- cor.test(fgsea_wide$NES_Cancer_vs_Healthy,
                        fgsea_wide$NES_Training_CR, method = "spearman")
nes_ci_all  <- fisher_z_ci(nes_cor_all$estimate, nrow(fgsea_wide))
nes_cor_sig <- if (nrow(fgsea_sig) >= 3) {
  cor.test(fgsea_sig$NES_Cancer_vs_Healthy,
           fgsea_sig$NES_Training_CR, method = "spearman")
} else NULL

nes_lim <- max(abs(c(fgsea_wide$NES_Cancer_vs_Healthy,
                      fgsea_wide$NES_Training_CR))) * 1.15

# Quadrant counts on sig terms
n_conc_tr  <- sum(fgsea_sig$NES_Cancer_vs_Healthy > 0 & fgsea_sig$NES_Training_CR > 0)
n_conc_bl  <- sum(fgsea_sig$NES_Cancer_vs_Healthy < 0 & fgsea_sig$NES_Training_CR < 0)
n_disc_q2  <- sum(fgsea_sig$NES_Cancer_vs_Healthy < 0 & fgsea_sig$NES_Training_CR > 0)
n_disc_q4  <- sum(fgsea_sig$NES_Cancer_vs_Healthy > 0 & fgsea_sig$NES_Training_CR < 0)

n_conc_pw     <- n_conc_tr + n_conc_bl
n_total_sig   <- nrow(fgsea_sig)
pw_conc_frac  <- if (n_total_sig > 0) n_conc_pw / n_total_sig else 0
pw_conc_binom <- if (n_total_sig > 0) binom.test(n_conc_pw, n_total_sig) else NULL
pw_conc_ci    <- if (!is.null(pw_conc_binom)) pw_conc_binom$conf.int * 100 else c(NA, NA)

message(sprintf("  NES Spearman (all): rho = %.3f [%.3f, %.3f], p = %.2g",
                nes_cor_all$estimate, nes_ci_all[1], nes_ci_all[2],
                nes_cor_all$p.value))
if (!is.null(nes_cor_sig)) {
  nes_ci_sig <- fisher_z_ci(nes_cor_sig$estimate, nrow(fgsea_sig))
  message(sprintf("  NES Spearman (sig): rho = %.3f [%.3f, %.3f], p = %.2g",
                  nes_cor_sig$estimate, nes_ci_sig[1], nes_ci_sig[2],
                  nes_cor_sig$p.value))
}
message(sprintf("  Pathway concordance: %d/%d sig (%.1f%%) [%.1f, %.1f]",
                n_conc_pw, n_total_sig, pw_conc_frac * 100,
                pw_conc_ci[1], pw_conc_ci[2]))

# --- 7. Label top 20 by combined significance
txt_gene <- scale_text(BASE_GENE, PG_W)
txt_quad <- scale_text(BASE_QUADRANT, PG_W)

label_pw <- fgsea_sig %>%
  mutate(max_padj = pmax(
    replace_na(padj_Cancer_vs_Healthy, 1),
    replace_na(padj_Training_CR, 1)
  )) %>%
  arrange(max_padj) %>%
  slice_head(n = 20) %>%
  mutate(
    label_fill     = SIG_LABEL_FILL_F4[as.character(significance)],
    label_text_col = SIG_LABEL_TEXT_F4[as.character(significance)],
    nudge_y = case_when(
      significance == "Sig Both"          ~  0.15,
      significance == "Sig Cancer only"   ~ -0.15,
      significance == "Sig Training only" ~  0.10,
      TRUE ~ 0
    )
  ) %>%
  arrange(significance)

# --- 8. Build plot layers
# NS terms (background)
ns_df <- fgsea_wide %>% filter(significance == "NS")

# Sig terms (foreground) with database-coded border
sig_df <- fgsea_wide %>%
  filter(significance != "NS") %>%
  mutate(
    border_col = case_when(
      database == "Hallmark" ~ "black",
      database == "KEGG"     ~ DB_COLORS[["KEGG"]],
      database == "Reactome" ~ DB_COLORS[["Reactome"]],
      database == "GO:BP"    ~ DB_COLORS[["GO:BP"]],
      TRUE                   ~ "grey75"
    ),
    bubble_alpha = case_when(
      significance == "Sig Both"          ~ 0.75,
      significance == "Sig Cancer only"   ~ 0.85,
      significance == "Sig Training only" ~ 0.85,
      TRUE ~ 0.60
    ),
    draw_order = factor(significance,
      levels = c("Sig Training only", "Sig Cancer only", "Sig Both"))
  ) %>%
  arrange(draw_order)

# Subtitle
rho_sig_str <- if (!is.null(nes_cor_sig)) {
  sprintf(", rho(sig) = %.2f", nes_cor_sig$estimate)
} else ""
subtitle_str <- sprintf(
  "H + KEGG + Reactome + GO:BP (collapsePathways) | %d pathways (%d sig.) | fGSEA on limma t-statistics\n\u03c1(all) = %.2f [%.2f, %.2f], p %s%s | %.0f%% concordant",
  nrow(fgsea_wide), n_total_sig,
  nes_cor_all$estimate, nes_ci_all[1], nes_ci_all[2],
  ifelse(nes_cor_all$p.value < 0.001, "< 0.001",
         sprintf("= %.3f", nes_cor_all$p.value)),
  rho_sig_str, pw_conc_frac * 100
)

# Database border color legend (manual)
db_legend_label <- sprintf(
  "Border: %s",
  paste(c("Hallmark = black", "KEGG = orange", "Reactome = blue", "GO:BP = teal"),
        collapse = ", ")
)

pG <- ggplot(mapping = aes(x = NES_Cancer_vs_Healthy, y = NES_Training_CR)) +
  # Quadrant backgrounds
  annotate("rect", xmin = 0, xmax = Inf,  ymin = 0, ymax = Inf,
           fill = "#FFE0E0", alpha = 0.55) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#FFE0E0", alpha = 0.55) +
  annotate("rect", xmin = 0, xmax = Inf,  ymin = -Inf, ymax = 0,
           fill = "#DCEEFF", alpha = 0.55) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#DCEEFF", alpha = 0.55) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  # NS points
  geom_point(data = ns_df, size = 1.5, fill = "grey70",
             shape = 21, color = "grey55", alpha = 0.40, stroke = 0.4) +
  # Sig points with database border color
  geom_point(data = sig_df, aes(fill = significance, size = set_size),
             shape = 21, color = sig_df$border_col,
             alpha = sig_df$bubble_alpha, stroke = 0.8) +
  scale_fill_manual(values = SIG_COLORS_F4, name = "Significance") +
  scale_size_continuous(range = c(2, 8), name = "Set size",
                        breaks = c(20, 50, 100, 200)) +
  # Labels (top 20 by combined significance)
  geom_label_repel(data = label_pw, aes(label = pathway_label),
                   fill = label_pw$label_fill, color = label_pw$label_text_col,
                   nudge_y = label_pw$nudge_y,
                   size = txt_gene, fontface = "bold",
                   max.overlaps = 40,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.5, force = 3, force_pull = 0.5,
                   label.padding = unit(1.5, "pt"),
                   label.r = unit(1, "pt"),
                   label.size = 0.15, seed = 42) +
  # Quadrant count annotations
  annotate("label", x = Inf, y = Inf,
           label = sprintf("Concordant Up  n = %d", n_conc_tr),
           hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = -Inf,
           label = sprintf("Concordant Down  n = %d", n_conc_bl),
           hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = Inf,
           label = sprintf("Discordant  n = %d", n_disc_q2),
           hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = Inf, y = -Inf,
           label = sprintf("Discordant  n = %d", n_disc_q4),
           hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  scale_x_continuous(expand = expansion(0, 0)) +
  scale_y_continuous(expand = expansion(0, 0)) +
  coord_cartesian(xlim = c(-nes_lim, nes_lim), ylim = c(-nes_lim, nes_lim)) +
  labs(
    title    = "Pathway-Level Concordance (collapsePathways)",
    subtitle = subtitle_str,
    x = "NES (Cancer vs Healthy)",
    y = "NES (Training CR)",
    caption = db_legend_label
  ) +
  FIG_THEME +
  theme(legend.position = "none",
        plot.caption = element_text(size = 7, color = "grey40", hjust = 0))

# --- 9. Save outputs
ggsave(file.path(RPT, "supp_nes_scatter_collapsed_SUPP.pdf"), pG,
       width = PG_W, height = PG_W, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "supp_nes_scatter_collapsed_SUPP.png"), pG,
       width = PG_W, height = PG_W, units = "mm", dpi = 300)

# Export all terms (not just significant)
fgsea_wide %>%
  transmute(
    pathway, pathway_label, database,
    NES_Cancer_vs_Healthy = round(NES_Cancer_vs_Healthy, 3),
    NES_Training_CR       = round(NES_Training_CR, 3),
    padj_Cancer_vs_Healthy = signif(padj_Cancer_vs_Healthy, 4),
    padj_Training_CR       = signif(padj_Training_CR, 4),
    significance           = as.character(significance),
    set_size,
    independent_in_CvH = pathway %in% independent_pw[["Cancer_vs_Healthy"]],
    independent_in_TR  = pathway %in% independent_pw[["Training_CR"]]
  ) %>%
  arrange(significance, desc(abs(NES_Cancer_vs_Healthy) + abs(NES_Training_CR))) %>%
  write_csv(file.path(DAT, "nes_scatter_collapsed.csv"))

cat("F04 CRvH Supp NES scatter (collapsePathways) done\n")
