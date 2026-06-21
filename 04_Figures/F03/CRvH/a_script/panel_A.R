# F03/CRvH — Panel A: DEPs per Contrast (Pseudo-log Stacked Bar)
# Contrasts: Cancer_vs_Healthy, Training_CR
# Outputs: pA (ggplot object), audit_panel_A_dep_fraction_ci.csv
#
# CONTRACT — downstream scripts depend on these objects remaining in the
# environment: sig_sets, dir_map, all_genes, SET_LABELS, SET_DISPLAY_COLORS,
# pi_total, fdr_total.  Do NOT wrap in a function or local scope.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F03/a_script/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(tidyr)
})

DEP_FILE <- "03_DEP/a_non_imputed/c_data/combined_results_pi.csv"
RPT      <- "04_Figures/F03/CRvH/b_reports"
RPT_PDF       <- file.path(RPT, "main", "pdf")
RPT_PNG       <- file.path(RPT, "main", "png")
RPT_SUPP_PDF  <- file.path(RPT, "supp", "pdf")
RPT_SUPP_PNG  <- file.path(RPT, "supp", "png")
DAT      <- "04_Figures/F03/CRvH/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

CONTRASTS <- c("Cancer_vs_Healthy", "Training_CR")
dep_df    <- read_csv(DEP_FILE, show_col_types = FALSE) |>
  dplyr::mutate(contrast = dplyr::recode(contrast,
                                         CRvH_Baseline = "Cancer_vs_Healthy",
                                         CR_Training   = "Training_CR")) |>
  tidyr::pivot_wider(id_cols = c(uniprot_id, gene, protein, description),
                     names_from = contrast,
                     values_from = c(logFC, t, P.Value, adj.P.Val, pi_score, sig_pi),
                     names_glue = "{.value}_{contrast}")
pdf_device <- get_pdf_device()
PA_W <- 170
PA_H <- 60

SET_LABELS <- CTR_SHORT[CONTRASTS]
all_genes  <- unique(dep_df$gene[!is.na(dep_df$gene)])

sig_sets <- list()
dir_map  <- list()

for (ctr in CONTRASTS) {
  pi_vals  <- dep_df[[paste0("pi_score_", ctr)]]
  lfc_vals <- dep_df[[paste0("logFC_", ctr)]]
  is_sig   <- !is.na(pi_vals) & pi_vals < 0.05
  sig_sets[[ctr]] <- dep_df$gene[is_sig]
  dir_map[[ctr]]  <- setNames(ifelse(lfc_vals[is_sig] > 0, "Up", "Down"),
                               dep_df$gene[is_sig])
}

n_total <- length(all_genes)
pi_ci <- data.frame(
  contrast = CONTRASTS,
  n_sig    = sapply(sig_sets, length),
  n_total  = n_total,
  pct      = 100 * sapply(sig_sets, length) / n_total,
  ci_lo    = sapply(sig_sets, function(s)
    100 * binom.test(length(s), n_total)$conf.int[1]),
  ci_hi    = sapply(sig_sets, function(s)
    100 * binom.test(length(s), n_total)$conf.int[2])
)

pi_total  <- sum(sapply(sig_sets, length))
fdr_total <- sum(sapply(CONTRASTS, function(ctr) {
  fdr_col <- paste0("adj.P.Val_", ctr)
  if (fdr_col %in% names(dep_df))
    sum(dep_df[[fdr_col]] < 0.10, na.rm = TRUE) else 0
}))

frac_list <- lapply(CONTRASTS, function(ctr) {
  fdr_col <- paste0("adj.P.Val_", ctr)
  p_col   <- paste0("P.Value_", ctr)
  tibble(
    contrast  = SET_LABELS[ctr],
    threshold = c("p < 0.05", "q < 0.10", "\u03A0 < 0.05"),
    n = c(sum(!is.na(dep_df[[p_col]])   & dep_df[[p_col]]   < 0.05),
          sum(!is.na(dep_df[[fdr_col]]) & dep_df[[fdr_col]] < 0.10),
          length(sig_sets[[ctr]]))
  )
})

display_levels <- rev(unname(SET_LABELS[CONTRASTS]))

frac_df <- bind_rows(frac_list) |>
  mutate(
    contrast  = factor(contrast, levels = display_levels),
    threshold = factor(threshold, levels = c("p < 0.05", "q < 0.10",
                                              "\u03A0 < 0.05")),
    pct       = 100 * n / length(all_genes),
    fill_key  = paste(contrast, threshold, sep = "___")
  ) |>
  filter(n > 1)

SET_DISPLAY_COLORS <- setNames(
  unname(CONTRAST_COLORS[CONTRASTS]),
  unname(SET_LABELS[CONTRASTS])
)

FRAC_FILL <- c()
for (cname in names(SET_DISPLAY_COLORS)) {
  col <- unname(SET_DISPLAY_COLORS[cname])
  FRAC_FILL[paste(cname, "p < 0.05",      sep = "___")] <- adjustcolor(col, alpha.f = 0.25)
  FRAC_FILL[paste(cname, "q < 0.10",      sep = "___")] <- adjustcolor(col, alpha.f = 0.55)
  FRAC_FILL[paste(cname, "\u03A0 < 0.05", sep = "___")] <- col
}

THRESH_LABEL <- c("p < 0.05" = "p \u2264 0.05", "q < 0.10" = "FDR \u2264 0.10",
                  "\u03A0 < 0.05" = "\u03A0 \u2264 0.05")

label_df <- frac_df |>
  group_by(contrast) |> arrange(contrast, threshold) |>
  mutate(label     = THRESH_LABEL[as.character(threshold)],
         next_pct  = lead(pct, default = 0),
         seg_width = pct - next_pct,
         label_y   = (next_pct + pct) / 2,
         text_col  = if_else(threshold == "p < 0.05", "grey20", "white")) |>
  filter(seg_width > 0.3) |>
  ungroup()

# Count labels (n = X) at bar end
count_labels <- frac_df |>
  filter(threshold == "p < 0.05") |>
  mutate(label = sprintf("n = %s", format(n, big.mark = ",")))

# Background annotation rectangles (one per contrast, bottom to top)
bg_rects <- lapply(seq_along(CONTRASTS), function(i) {
  annotate("rect",
           xmin = length(CONTRASTS) - i + 0.5,
           xmax = length(CONTRASTS) - i + 1.5,
           ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS[CONTRASTS[i]], alpha = 0.20,
           color = "grey70", linewidth = 0.2)
})

pA <- ggplot(frac_df, aes(x = contrast, y = pct, fill = fill_key)) +
  bg_rects +
  geom_col(position = "identity", width = 0.75, color = "black", linewidth = 0.3) +
  geom_text(data = label_df,
            aes(x = contrast, y = label_y, label = label, color = I(text_col)),
            inherit.aes = FALSE, hjust = 0.5,
            size = scale_text(BASE_COUNT - 1.5, PA_W), fontface = "bold") +
  geom_text(data = count_labels,
            aes(x = contrast, y = pct, label = label),
            inherit.aes = FALSE, hjust = -0.1, vjust = 0.5,
            size = 2.5, fontface = "bold", color = "grey30") +
  scale_fill_manual(values = FRAC_FILL) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 5, base = 10),
                     expand = expansion(mult = c(0, 0.08)),
                     breaks = c(0, 2, 5, 10, 20, 35)) +
  coord_flip() +
  labs(title = "DEPs per Contrast (CRvH Model)",
       subtitle = sprintf("Primary threshold: \u03A0 < 0.05 | %s proteins",
                          format(length(all_genes), big.mark = ",")),
       x = NULL, y = "% of proteome",
       tag = "A") +
  FIG_THEME + theme(legend.position = "none",
                    axis.text.y = element_text(face = "bold",
                                               size = FIG_AXIS_TEXT - 0.5))

write.csv(pi_ci, file.path(DAT, "audit_panel_A_dep_fraction_ci.csv"), row.names = FALSE)

ggsave(file.path(RPT_PDF, "panel_A_dep_counts_MAIN.pdf"), pA,
       width = PA_W, height = PA_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_PNG, "panel_A_dep_counts_MAIN.png"), pA,
       width = PA_W, height = PA_H, units = "mm", dpi = 300)

cat("F03/CRvH Panel A done.\n")
