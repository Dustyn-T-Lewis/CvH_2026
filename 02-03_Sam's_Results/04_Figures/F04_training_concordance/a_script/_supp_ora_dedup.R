#!/usr/bin/env Rscript
# SUPP Panel A: ORA Dedup Sensitivity
# Defends main Panel A — shows ORA results are robust to Jaccard dedup cutoff.
# Sourced by 02_supp_panels.R — expects style.R + pathway_utils.R already loaded.
# Exports: pS_ora_dedup (ggplot)

BASE <- "02-03_Sam's_Results/04_Figures/F04_training_concordance"
DAT  <- file.path(BASE, "c_data", "panel_supp")
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

# Data
dep_df <- read_csv("02-03_Sam's_Results/03_DEP/c_data/F04_combined_CREvPLA.csv",
                   show_col_types = FALSE)

scatter_df <- dep_df |>
  transmute(gene,
            logFC_TY = logFC_Training_CRE,
            logFC_TO = logFC_Training_PLA) |>
  filter(!is.na(logFC_TY), !is.na(logFC_TO)) |>
  mutate(quadrant = case_when(
    logFC_TY > 0 & logFC_TO > 0 ~ "Concordant Up",
    logFC_TY < 0 & logFC_TO < 0 ~ "Concordant Down",
    logFC_TY > 0 & logFC_TO < 0 ~ "Discordant (CRE Up / PLA Down)",
    TRUE                         ~ "Discordant (CRE Down / PLA Up)"))

universe <- scatter_df$gene

# Pathway collection (same as _panel_A_ORA.R)
pw_collection <- build_pathway_collection(min_size = 15, max_size = 500,
                                           include_goslim = FALSE,
                                           exclude_variants = TRUE)

# Run ORA at multiple Jaccard cutoffs
cutoffs <- c(0.3, 0.5, 0.7, 1.0)
target_quads <- c("Concordant Up", "Concordant Down")

results <- list()
for (quad in target_quads) {
  genes <- scatter_df$gene[scatter_df$quadrant == quad]
  if (length(genes) < 5) next

  for (jc in cutoffs) {
    message(sprintf("  ORA: %s | Jaccard = %.1f", quad, jc))
    ora_res <- tryCatch({
      run_ora_deduplicated(
        genes = genes, universe = universe, pathways = pw_collection,
        jaccard_cutoff = jc, min_size = 15, max_size = 500, padj_cutoff = 0.05
      )
    }, error = function(e) { message("    error: ", e$message); tibble() })

    n_sig <- if (nrow(ora_res) > 0) sum(ora_res$padj < 0.05) else 0L
    results[[length(results) + 1]] <- tibble(
      quadrant = quad, jaccard_cutoff = jc,
      n_enriched = n_sig, n_genes = length(genes))
  }
}

sens_df <- bind_rows(results)
write_csv(sens_df, file.path(DAT, "SUPP_ora_dedup_sensitivity.csv"))
message("Sensitivity data:\n", paste(capture.output(print(sens_df)), collapse = "\n"))

# Plot
sens_df <- sens_df |>
  mutate(cutoff_label = factor(sprintf("J = %.1f", jaccard_cutoff),
                               levels = sprintf("J = %.1f", cutoffs)))

pS_ora_dedup <- ggplot(sens_df, aes(x = quadrant, y = n_enriched,
                                     fill = cutoff_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = n_enriched),
            position = position_dodge(width = 0.7), vjust = -0.3,
            size = BASE_COUNT, fontface = "bold") +
  scale_fill_manual(values = c("J = 0.3" = "grey80", "J = 0.5" = "grey60",
                               "J = 0.7" = "grey40", "J = 1.0" = "grey20"),
                    name = "Jaccard cutoff") +
  labs(title = "ORA Dedup Sensitivity",
       subtitle = "# enriched pathways (FDR < 0.05) by Jaccard dedup cutoff",
       x = NULL, y = "Enriched pathways (FDR < 0.05)") +
  FIG_THEME +
  theme(legend.key.size = unit(1.8, "mm"))

# Standalone save
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)

PW <- 89; PH <- 70
ggsave(file.path(RPT_PNG, "SUPP_ora_dedup.png"), pS_ora_dedup,
       width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_ora_dedup.pdf"), pS_ora_dedup,
       width = PW, height = PH, units = "mm", device = pdf_device)

pS_ora_dedup <- strip_for_composite(pS_ora_dedup)

message("SUPP Panel A (ORA dedup) done")
