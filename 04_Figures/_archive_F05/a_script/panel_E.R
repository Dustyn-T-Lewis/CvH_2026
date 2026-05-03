# Figure 5 — Panel E: NES Scatter (Training CRE × Training PLA)
# De novo fGSEA on Hallmark + GO Slim (a priori collection).
# X = NES_Training_CRE, Y = NES_Training_PLA
# Identity line (slope=1, concordant response), significance coloring,
# Spearman correlation, pathway labels, quadrant counts.
# Outputs: panel_E_nes_supp.pdf/png, panel_E/nes_scatter.csv

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(stringr)
})

PE_W <- 200; PE_H <- 200

RPT <- "04_Figures/F05/b_reports"
DAT <- "04_Figures/F05/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_E"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# --- Build Hallmark + GO Slim pathway collection ---
pw_collection <- build_hallmark_goslim_collection(min_size = 10, max_size = 500)

# --- Load DEP results (t-statistics for ranking) ---
dep_cr <- read_csv("03_DEP/c_data/03_combined_results_CR.csv",
                    show_col_types = FALSE)

# --- Run fGSEA per contrast on Hallmark + GO Slim ---
set.seed(42)
fgsea_list <- list()
contrasts_supp <- c("Training_CRE", "Training_PLA")

for (ctr in contrasts_supp) {
  tcol  <- paste0("t_", ctr)
  stats <- setNames(dep_cr[[tcol]], dep_cr$gene)
  stats <- stats[!is.na(stats) & is.finite(stats)]

  if (anyDuplicated(names(stats))) {
    dup_df <- tibble(gene = names(stats), t = unname(stats)) |>
      group_by(gene) |>
      slice_max(abs(t), n = 1, with_ties = FALSE) |>
      ungroup()
    stats <- setNames(dup_df$t, dup_df$gene)
  }
  stats <- sort(stats, decreasing = TRUE)

  res <- run_fgsea_perdb(
    ranks    = stats,
    pathways = pw_collection,
    nperm    = 10000,
    min_size = 10,
    max_size = 500
  )
  res$contrast <- ctr
  fgsea_list[[ctr]] <- res
}

fgsea_hg <- bind_rows(fgsea_list)

# --- Pivot to wide ---
fgsea_wide <- fgsea_hg |>
  select(pathway, contrast, NES, padj, size, database) |>
  pivot_wider(id_cols = c(pathway, database),
              names_from = contrast,
              values_from = c(NES, padj, size)) |>
  filter(!is.na(NES_Training_CRE), !is.na(NES_Training_PLA)) |>
  mutate(set_size = coalesce(size_Training_CRE, size_Training_PLA))

# --- Significance categories (custom supplement palette) ---
SIG_COLORS_SUPP <- c(
  "Sig Both"       = "#2E7D32",
  "Sig CRE only"   = "#2166AC",
  "Sig PLA only"   = "#D6604D",
  "NS"             = "grey70"
)
SIG_LABEL_FILL_SUPP <- c(
  "Sig Both"     = scales::alpha("#2E7D32", 0.75),
  "Sig CRE only" = scales::alpha("#2166AC", 0.75),
  "Sig PLA only" = scales::alpha("#D6604D", 0.75),
  "NS"           = scales::alpha("grey70",  0.75)
)
SIG_LABEL_TEXT_SUPP <- c(
  "Sig Both"     = "white",
  "Sig CRE only" = "white",
  "Sig PLA only" = "white",
  "NS"           = "white"
)

fgsea_wide <- fgsea_wide |>
  mutate(
    sig_CRE = !is.na(padj_Training_CRE) & padj_Training_CRE < 0.05,
    sig_PLA = !is.na(padj_Training_PLA) & padj_Training_PLA < 0.05,
    significance = case_when(
      sig_CRE & sig_PLA ~ "Sig Both",
      sig_CRE           ~ "Sig CRE only",
      sig_PLA           ~ "Sig PLA only",
      TRUE              ~ "NS"
    ) |> factor(levels = names(SIG_COLORS_SUPP)),
    pathway_label = clean_pathway_name(pathway)
  )

fgsea_sig <- fgsea_wide |> filter(significance != "NS")

message(sprintf("  %d total pathways (Hallmark: %d, GO Slim: %d) | %d significant",
                nrow(fgsea_wide),
                sum(fgsea_wide$database == "Hallmark"),
                sum(fgsea_wide$database == "GOSlim"),
                nrow(fgsea_sig)))

# --- Spearman correlation (all terms + sig-only) ---
nes_cor_all <- cor.test(fgsea_wide$NES_Training_CRE,
                        fgsea_wide$NES_Training_PLA, method = "spearman")
nes_ci_all  <- fisher_z_ci(nes_cor_all$estimate, nrow(fgsea_wide))
nes_cor_sig <- if (nrow(fgsea_sig) >= 3) {
  cor.test(fgsea_sig$NES_Training_CRE,
           fgsea_sig$NES_Training_PLA, method = "spearman")
} else NULL

nes_lim <- max(abs(c(fgsea_wide$NES_Training_CRE,
                      fgsea_wide$NES_Training_PLA))) * 1.15

# --- Quadrant counts (sig terms only) ---
n_q1 <- sum(fgsea_sig$NES_Training_CRE > 0 & fgsea_sig$NES_Training_PLA > 0)
n_q2 <- sum(fgsea_sig$NES_Training_CRE < 0 & fgsea_sig$NES_Training_PLA > 0)
n_q3 <- sum(fgsea_sig$NES_Training_CRE < 0 & fgsea_sig$NES_Training_PLA < 0)
n_q4 <- sum(fgsea_sig$NES_Training_CRE > 0 & fgsea_sig$NES_Training_PLA < 0)

n_concordant <- n_q1 + n_q3
n_discordant <- n_q2 + n_q4

message(sprintf("  NES Spearman (all): rho = %.3f [%.3f, %.3f], p = %.2g",
                nes_cor_all$estimate, nes_ci_all[1], nes_ci_all[2],
                nes_cor_all$p.value))
if (!is.null(nes_cor_sig)) {
  nes_ci_sig <- fisher_z_ci(nes_cor_sig$estimate, nrow(fgsea_sig))
  message(sprintf("  NES Spearman (sig): rho = %.3f [%.3f, %.3f], p = %.2g",
                  nes_cor_sig$estimate, nes_ci_sig[1], nes_ci_sig[2],
                  nes_cor_sig$p.value))
}
message(sprintf("  Concordant: %d/%d sig, Discordant: %d/%d sig",
                n_concordant, nrow(fgsea_sig), n_discordant, nrow(fgsea_sig)))

# --- Scaled text ---
txt_gene <- scale_text(BASE_GENE, PE_W)
txt_quad <- scale_text(BASE_QUADRANT, PE_W)

# --- Labels: all sig terms ---
label_pw <- fgsea_sig |>
  mutate(
    label_fill     = SIG_LABEL_FILL_SUPP[as.character(significance)],
    label_text_col = SIG_LABEL_TEXT_SUPP[as.character(significance)]
  ) |>
  mutate(pathway_label = pathway_label |>
    str_replace("Amino Acid Metabolic.*", "Amino Acid Metabolism") |>
    str_replace("Muscle System.*", "Muscle System") |>
    str_replace("Ketone Metabolic.*", "Ketone Metabolism")
  ) |>
  mutate(nudge_y = case_when(
    significance == "Sig Both"     ~  0.15,
    significance == "Sig CRE only" ~ -0.15,
    significance == "Sig PLA only" ~  0.10,
    TRUE ~ 0
  )) |>
  arrange(significance)

# --- Split data for layered plotting ---
ns_df  <- fgsea_wide |> filter(significance == "NS")
sig_df <- fgsea_wide |> filter(significance != "NS") |>
  mutate(
    border_col   = ifelse(database == "Hallmark", "black", "grey75"),
    bubble_alpha = case_when(
      significance == "Sig Both"     ~ 0.75,
      significance == "Sig CRE only" ~ 0.85,
      significance == "Sig PLA only" ~ 0.85,
      TRUE ~ 0.60
    ),
    draw_order = factor(significance,
      levels = c("Sig PLA only", "Sig CRE only", "Sig Both"))
  ) |>
  arrange(draw_order)

# --- Subtitle ---
rho_sig_str <- if (!is.null(nes_cor_sig)) sprintf(", rho(sig) = %.2f", nes_cor_sig$estimate) else ""
subtitle_str <- sprintf(
  "GO Slim + Hallmark (a priori) | %d pathways (%d sig.) | fGSEA on limma t-statistics\n\u03c1(all) = %.2f [%.2f, %.2f], p %s%s | %d concordant, %d discordant",
  nrow(fgsea_wide), nrow(fgsea_sig),
  nes_cor_all$estimate, nes_ci_all[1], nes_ci_all[2],
  ifelse(nes_cor_all$p.value < 0.001, "< 0.001", sprintf("= %.3f", nes_cor_all$p.value)),
  rho_sig_str, n_concordant, n_discordant
)

# --- Plot ---
pE <- ggplot(mapping = aes(x = NES_Training_CRE, y = NES_Training_PLA)) +
  # Reference lines
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  # Identity line (concordant response)
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  # NS points (background)
  geom_point(data = ns_df, size = 1.5, fill = "grey70",
             shape = 21, color = "grey55", alpha = 0.40, stroke = 0.4) +
  # Sig points (foreground)
  geom_point(data = sig_df, aes(fill = significance, size = set_size),
             shape = 21, color = sig_df$border_col,
             alpha = sig_df$bubble_alpha, stroke = 0.8) +
  scale_fill_manual(values = SIG_COLORS_SUPP, name = "Significance") +
  scale_size_continuous(range = c(2, 8), name = "Set size",
                        breaks = c(20, 50, 100, 200)) +
  # Labels
  geom_label_repel(data = label_pw, aes(label = pathway_label),
                   fill = label_pw$label_fill, color = label_pw$label_text_col,
                   nudge_y = label_pw$nudge_y,
                   size = txt_gene, fontface = "bold",
                   max.overlaps = 50,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.35, force = 5, force_pull = 0.3,
                   label.padding = unit(1.2, "pt"),
                   label.r = unit(1, "pt"),
                   label.size = 0.15, seed = 42) +
  # Quadrant count annotations
  annotate("label", x = Inf, y = Inf,
           label = sprintf("Concordant Up  n = %d", n_q1),
           hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
           color = "#2E7D32", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = -Inf,
           label = sprintf("Concordant Down  n = %d", n_q3),
           hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
           color = "#2E7D32", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = Inf,
           label = sprintf("Discordant  n = %d", n_q2),
           hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
           color = "#B71C1C", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = Inf, y = -Inf,
           label = sprintf("Discordant  n = %d", n_q4),
           hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
           color = "#B71C1C", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  scale_x_continuous(expand = expansion(0, 0)) +
  scale_y_continuous(expand = expansion(0, 0)) +
  coord_cartesian(xlim = c(-nes_lim, nes_lim), ylim = c(-nes_lim, nes_lim)) +
  labs(
    title    = "Pathway NES: Creatine vs Placebo Training",
    subtitle = subtitle_str,
    x = paste0("NES \u2014 ", CTR_SHORT["Training_CRE"]),
    y = paste0("NES \u2014 ", CTR_SHORT["Training_PLA"]),
    tag = "E"
  ) +
  FIG_THEME +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 8, face = "bold"),
    legend.text     = element_text(size = 7),
    legend.key.size = unit(3, "mm"),
    legend.margin   = margin(0, 0, 0, 0)
  ) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 3, alpha = 0.8)),
         size = guide_legend(nrow = 1))

# --- Save ---
ggsave(file.path(RPT, "panel_E_nes_supp.pdf"), pE,
       width = PE_W, height = PE_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_E_nes_supp.png"), pE,
       width = PE_W, height = PE_H, units = "mm", dpi = 300)

# --- Export ---
fgsea_wide |>
  transmute(
    pathway, pathway_label, database,
    NES_Training_CRE  = round(NES_Training_CRE, 3),
    NES_Training_PLA  = round(NES_Training_PLA, 3),
    padj_Training_CRE = signif(padj_Training_CRE, 4),
    padj_Training_PLA = signif(padj_Training_PLA, 4),
    significance      = as.character(significance),
    set_size
  ) |>
  arrange(significance, desc(abs(NES_Training_CRE) + abs(NES_Training_PLA))) |>
  write_csv(file.path(DAT, "panel_E", "nes_scatter.csv"))

cat("Panel E (NES scatter CRE vs PLA, Hallmark + GO Slim) done.\n")
