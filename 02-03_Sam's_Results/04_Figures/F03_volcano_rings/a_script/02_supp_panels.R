#!/usr/bin/env Rscript
# F03 Supp — 5-panel diagnostic composite (2×2 + full-row)
# A: p-value histograms  B: Pi-score distributions  C: FDR distributions
# D: MA plots  E: outlier sensitivity (full bottom row)

setwd(rprojroot::find_rstudio_root_file())

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(patchwork)
library(cowplot)

source("02-03_Sam's_Results/04_Figures/shared/style.R")

# cairo_pdf requires X11 (unavailable on this system); fall back to base pdf()
get_pdf_device <- function() grDevices::pdf

# Helpers present in YvO style.R but absent from CvH style.R
if (!exists("strip_for_composite")) {
  strip_for_composite <- function(p) {
    p + labs(title = NULL, subtitle = NULL, tag = NULL) +
      theme(legend.position = "none")
  }
}
if (!exists("composite_text_sizes")) {
  composite_text_sizes <- function(comp_h_mm) {
    list(title    = pmax(6, pmin(8, round(5 + comp_h_mm / 80))),
         subtitle = pmax(4, pmin(6, round(3 + comp_h_mm / 100))),
         tag      = 8)
  }
}

BASE    <- "02-03_Sam's_Results/04_Figures/F03_volcano_rings"
RPT_PNG <- file.path(BASE, "b_reports", "supp", "png", "panels")
RPT_PDF <- file.path(BASE, "b_reports", "supp", "pdf", "panels")
DAT     <- file.path(BASE, "c_data", "supp")
for (d in c(RPT_PNG, RPT_PDF, DAT)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

pdf_dev  <- get_pdf_device()
DEP_DIR  <- "02-03_Sam's_Results/03_DEP/c_data/04_per_contrast_results"
CTRS     <- c("Cancer_vs_Healthy", "Training_CR", "Training_CRE", "Training_PLA")

per_contrast <- lapply(CTRS, function(ctr) {
  read_csv(file.path(DEP_DIR, paste0(ctr, ".csv")), show_col_types = FALSE) |>
    as.data.frame()
})
names(per_contrast) <- CTRS

# ── Parametric histogram builder (collapses A–C) ────────────────────────────

make_dist_panel <- function(col, fill, xlab, vline = NULL, stat_fmt, title, tag) {
  hist_df <- bind_rows(lapply(CTRS, function(ctr) {
    tibble(contrast = ctr, value = per_contrast[[ctr]][[col]])
  })) |>
    filter(!is.na(value)) |>
    mutate(contrast = factor(contrast, levels = CTRS))

  n_sig <- hist_df |>
    summarise(n_sig = sum(value < 0.05), .by = contrast)

  n_bins <- 20

  p <- ggplot(hist_df, aes(value)) +
    geom_histogram(breaks = seq(0, 1, length.out = n_bins + 1),
                   fill = fill, color = "white", linewidth = 0.3)

  if (!is.null(vline))
    p <- p + geom_vline(xintercept = vline, linetype = "dashed",
                        color = "grey30", linewidth = 0.4)

  if (col == "P.Value") {
    uniform_ref <- hist_df |> summarise(h = n() / n_bins, .by = contrast)
    p <- p + geom_hline(data = uniform_ref, aes(yintercept = h),
                        linetype = "dashed", color = "grey40", linewidth = 0.4)
  }

  p <- p +
    geom_text(data = n_sig,
              aes(x = 0.5, y = Inf, label = CTR_SHORT[as.character(contrast)]),
              inherit.aes = FALSE, hjust = 0.5, vjust = 1.2,
              size = 2.0, fontface = "bold", color = "grey20") +
    geom_text(data = n_sig,
              aes(x = 0.5, y = Inf, label = sprintf(stat_fmt, n_sig)),
              inherit.aes = FALSE, hjust = 0.5, vjust = 3.8,
              size = 2.8, fontface = "bold", color = "grey40") +
    facet_wrap(~contrast, ncol = 2, scales = "free_y",
               labeller = labeller(contrast = CTR_SHORT)) +
    labs(title = title,
         subtitle = sprintf("%s proteins | 20 bins",
                            format(round(nrow(hist_df) / length(CTRS)), big.mark = ",")),
         x = xlab, y = "Proteins", tag = tag) +
    FIG_THEME + theme(strip.background = element_blank(), strip.text = element_blank())

  write_csv(hist_df, file.path(DAT, sprintf("panel_%s.csv", gsub("[. ]", "_", col))))
  ggsave(file.path(RPT_PNG, sprintf("SUPP_%s.png", gsub("[. ]", "_", title))), p,
         width = 89, height = 75, units = "mm", dpi = 300)
  p
}

# ── Panels A–C: histogram variants ──────────────────────────────────────────

pA_supp <- make_dist_panel("P.Value", "#5DA5DA", "p-value", NULL,
                            "p < 0.05: %d", "Raw p-value distribution", "a")
pA_title <- "Raw p-value distribution"

pB_supp <- make_dist_panel("pi_score", "#E05A4E", "Pi-score", 0.05,
                             "Pi < 0.05: %d", "Pi-score distribution", "b")
pB_title <- "Pi-score distribution"

pC_supp <- make_dist_panel("adj.P.Val", "#9B7FBF", "FDR (BH)", 0.05,
                             "FDR < 0.10: %d", "FDR distribution", "c")
pC_title <- "FDR distribution"

pA_supp <- strip_for_composite(pA_supp)
pB_supp <- strip_for_composite(pB_supp)
pC_supp <- strip_for_composite(pC_supp)

# ── Panel D: MA plots ────────────────────────────────────────────────────────

ma_df <- bind_rows(lapply(CTRS, function(ctr) {
  per_contrast[[ctr]] |>
    transmute(contrast = ctr, gene, average_intensity, logFC,
              sig_pi, direction = case_when(
                sig_pi ==  1 ~ "Up",
                sig_pi == -1 ~ "Down",
                TRUE         ~ "NS"))
})) |>
  mutate(contrast  = factor(contrast, levels = CTRS),
         direction = factor(direction, levels = c("Up", "Down", "NS")))

n_dep_ma <- ma_df |> filter(direction != "NS") |> count(contrast, name = "n_dep")

pD_supp <- ggplot(ma_df, aes(average_intensity, logFC, color = direction)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_point(data = \(d) filter(d, direction == "NS"), alpha = 0.25, size = 0.6) +
  geom_point(data = \(d) filter(d, direction != "NS"), alpha = 0.85, size = 0.9) +
  geom_text(data = n_dep_ma,
            aes(x = Inf, y = Inf, label = sprintf("Pi: %d", n_dep)),
            inherit.aes = FALSE, hjust = 1.05, vjust = 1.5, size = 2.2,
            fontface = "bold", color = "grey20") +
  scale_color_manual(values = DIR_COLORS, name = "Pi < 0.05") +
  facet_wrap(~contrast, ncol = 2, labeller = labeller(contrast = CTR_SHORT)) +
  labs(title = "MA plots", x = "Mean log2 intensity", y = "logFC", tag = "d") +
  FIG_THEME + theme(strip.background = element_blank(), strip.text = element_blank(),
                    legend.position = "top", legend.key.size = unit(3, "mm"))
write_csv(ma_df, file.path(DAT, "panel_D_ma.csv"))
ggsave(file.path(RPT_PNG, "SUPP_panel_D_ma.png"), pD_supp,
       width = 89, height = 75, units = "mm", dpi = 300)
pD_title  <- "MA plots"
pD_supp   <- strip_for_composite(pD_supp)

# ── Panel E: Outlier sensitivity ─────────────────────────────────────────────
# Uses main CvH pipeline's 11_outlier_sensitivity.csv (covers all 6 contrasts)

OUT_SENS_PATH <- "03_DEP/c_data/11_outlier_sensitivity.csv"

out_sens <- tryCatch(
  read_csv(OUT_SENS_PATH, show_col_types = FALSE) |>
    filter(Contrast %in% CTRS),
  error = function(e) NULL)

if (!is.null(out_sens) && nrow(out_sens) > 0) {
  long_df <- out_sens |>
    mutate(Contrast = factor(Contrast, levels = CTRS)) |>
    pivot_longer(c(FDR_full, FDR_reduced, Pi_full, Pi_reduced),
                 names_to = "metric_cohort", values_to = "n") |>
    mutate(
      metric = factor(sub("_.*", "", metric_cohort),
                      levels = c("FDR", "Pi"),
                      labels = c("FDR < 0.05", "Pi < 0.05")),
      cohort = factor(sub(".*_", "", metric_cohort), levels = c("full", "reduced"))
    )

  pE_supp <- ggplot(long_df, aes(cohort, n, fill = cohort)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = n), vjust = -0.3, size = 2.8, fontface = "bold") +
    facet_grid(metric ~ Contrast, scales = "free_y", switch = "y",
               labeller = labeller(Contrast = CTR_SHORT)) +
    scale_fill_manual(values = c(full = "#2166AC", reduced = "#B2182B"),
                      labels = c("Full", "Outlier-removed"), name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(title = "DEP retention (outlier removal)", x = NULL, y = "DEPs", tag = "e") +
    FIG_THEME + theme(strip.text.x = element_blank(),
                      strip.text.y = element_text(face = "bold", size = FIG_STRIP_SIZE - 1),
                      legend.position = "top", legend.key.size = unit(3, "mm"))
  write_csv(long_df, file.path(DAT, "panel_E_outlier_sensitivity.csv"))
  ggsave(file.path(RPT_PNG, "SUPP_panel_E_outlier.png"), pE_supp,
         width = 178, height = 60, units = "mm", dpi = 300)
} else {
  message("WARNING: outlier sensitivity data not found at ", OUT_SENS_PATH,
          " — panel E will be a placeholder")
  pE_supp <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = "Outlier sensitivity data unavailable\n(run 03_DEP pipeline first)",
             size = 3, color = "grey40", hjust = 0.5) +
    theme_void() +
    labs(tag = "e")
}
pE_title <- "DEP retention (outlier removal)"
pE_supp  <- strip_for_composite(pE_supp)

# ── Supp composite (3×2 grid) ────────────────────────────────────────────────

SUPP_PDF <- file.path(BASE, "b_reports", "supp", "pdf")
SUPP_PNG <- file.path(BASE, "b_reports", "supp", "png")
dir.create(SUPP_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_PNG, recursive = TRUE, showWarnings = FALSE)

COMP_W <- 178; COMP_H <- 225
txt  <- composite_text_sizes(COMP_H)
grid <- (pA_supp | pB_supp) / (pC_supp | pD_supp) / pE_supp &
  theme(plot.margin = margin(9, 4, 4, 4))

X_L <- 0.012; X_R <- 0.512; X_TTL <- 0.029; SUB_OFF <- 0.014
Y_R1 <- 0.974; Y_R2 <- 0.651; Y_R3 <- 0.321

all_titles <- list(pA_title, pB_title, pC_title, pD_title, pE_title)
all_tags   <- LETTERS[1:5]
all_xs     <- c(X_L, X_R, X_L, X_R, X_L)
all_ys     <- c(Y_R1, Y_R1, Y_R2, Y_R2, Y_R3)

composite_supp <- ggdraw(grid)
for (i in seq_along(all_tags)) {
  composite_supp <- composite_supp +
    draw_label(all_tags[i], x = all_xs[i], y = all_ys[i],
               fontface = "bold", size = txt$tag, hjust = 0, vjust = 1) +
    draw_label(all_titles[[i]], x = all_xs[i] + X_TTL, y = all_ys[i],
               fontface = "bold", size = txt$title, hjust = 0, vjust = 1)
}

ggsave(file.path(SUPP_PDF, "SUPP_F03_composite.pdf"), composite_supp,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_dev)
ggsave(file.path(SUPP_PNG, "SUPP_F03_composite.png"), composite_supp,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message("F03 supp composite done")
