# F03/CRvH — Panel B: DEP Rank Location (Barcode Plot) — Side-by-side variant
# Left: Cancer_vs_Healthy | Right: Training_CR
# Outputs: panel_B_barcode_2x2_SUPP.pdf/.png

setwd(here::here())
source("04_Figures/F03/a_script/style.R")

library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(patchwork)

DEP_FILE <- "03_DEP/a_non_imputed/c_data/combined_results_pi.csv"
RPT      <- "04_Figures/F03/CRvH/b_reports"
RPT_PDF       <- file.path(RPT, "main", "pdf")
RPT_PNG       <- file.path(RPT, "main", "png")
RPT_SUPP_PDF  <- file.path(RPT, "supp", "pdf")
RPT_SUPP_PNG  <- file.path(RPT, "supp", "png")
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

dir.create(RPT_PDF,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PNG,      recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_SUPP_PNG, recursive = TRUE, showWarnings = FALSE)
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

PD_W <- 210
PD_H <- 80

# --- Build long-form data ---
rank_list <- lapply(CONTRASTS, function(ctr) {
  t_col   <- paste0("t_", ctr)
  pi_col  <- paste0("pi_score_", ctr)
  lfc_col <- paste0("logFC_", ctr)

  dep_df |>
    filter(!is.na(.data[[t_col]])) |>
    arrange(.data[[t_col]]) |>
    mutate(
      rank_frac = seq_len(n()) / n(),
      is_dep    = !is.na(.data[[pi_col]]) & .data[[pi_col]] < 0.05,
      direction = case_when(
        !is_dep              ~ NA_character_,
        .data[[lfc_col]] > 0 ~ "Up",
        TRUE                 ~ "Down"
      ),
      contrast = ctr
    ) |>
    select(gene, contrast, rank_frac, is_dep, direction)
})
rank_df <- bind_rows(rank_list)

dep_only_all <- rank_df |> filter(is_dep)
dep_only_all$direction <- factor(dep_only_all$direction, levels = c("Up", "Down"))

# --- Density ---
DENS_PAD <- 0.06
dens_all <- lapply(split(dep_only_all, dep_only_all$contrast), function(ctr_df) {
  lapply(split(ctr_df, ctr_df$direction, drop = TRUE), function(dir_df) {
    if (nrow(dir_df) < 2) return(NULL)
    d <- density(dir_df$rank_frac, adjust = 1.8,
                 from = -DENS_PAD, to = 1 + DENS_PAD, n = 512)
    tibble(x = d$x, y = d$y, direction = dir_df$direction[1],
           contrast = dir_df$contrast[1])
  }) |> bind_rows()
}) |> bind_rows()

dens_all <- dens_all |>
  group_by(contrast) |>
  mutate(y_norm = y / max(y)) |>
  ungroup()
dens_all$direction <- factor(dens_all$direction, levels = c("Up", "Down"))

TICK_DEPTH <- -0.40
ANNOT_SZ <- scale_text(BASE_STAT - 1.2, PD_W / 2)

n_down <- dep_only_all |> filter(direction == "Down") |> count(contrast) |> tibble::deframe()
n_up   <- dep_only_all |> filter(direction == "Up")   |> count(contrast) |> tibble::deframe()
DESC_DOWN <- c(Cancer_vs_Healthy = "proteins lower in cancer",
               Training_CR       = "proteins dec. with training")
DESC_UP   <- c(Cancer_vs_Healthy = "proteins higher in cancer",
               Training_CR       = "proteins inc. with training")
LABELS_DOWN <- setNames(paste(n_down[CONTRASTS], DESC_DOWN[CONTRASTS]), CONTRASTS)
LABELS_UP   <- setNames(paste(n_up[CONTRASTS],   DESC_UP[CONTRASTS]),   CONTRASTS)

# --- Helper: build one panel ---
make_panel <- function(ctr) {
  dep_ctr  <- dep_only_all |> filter(contrast == ctr)
  dens_ctr <- dens_all |> filter(contrast == ctr)
  bg_fill  <- unname(CONTRAST_COLORS[ctr])

  peaks <- dens_ctr |>
    group_by(direction) |>
    slice_max(y_norm, n = 1, with_ties = FALSE) |>
    ungroup()

  pk_down <- peaks |> filter(direction == "Down")
  pk_up   <- peaks |> filter(direction == "Up")

  p <- ggplot() +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
             fill = bg_fill, alpha = 0.20) +
    geom_ribbon(data = dens_ctr,
                aes(x = x, ymin = 0, ymax = y_norm, fill = direction),
                alpha = 0.30, outline.type = "upper") +
    geom_line(data = dens_ctr,
              aes(x = x, y = y_norm, color = direction),
              linewidth = 0.5) +
    geom_segment(data = dep_ctr,
                 aes(x = rank_frac, xend = rank_frac,
                     y = 0, yend = TICK_DEPTH, color = direction),
                 linewidth = 0.35, alpha = 0.8) +
    geom_hline(yintercept = 0, linewidth = 0.25, color = "grey50")

  if (nrow(pk_down) > 0) {
    p <- p +
      geom_label(aes(x = pk_down$x + 0.04, y = pk_down$y_norm,
                      label = LABELS_DOWN[ctr]),
                 hjust = 0, vjust = 0.5, size = ANNOT_SZ,
                 color = "white", fontface = "bold",
                 fill = scales::alpha(unname(DIR_COLORS["Down"]), 0.85),
                 linewidth = 0, label.padding = unit(1, "mm"),
                 label.r = unit(0.6, "mm"))
  }

  if (nrow(pk_up) > 0) {
    p <- p +
      geom_label(aes(x = pk_up$x - 0.04, y = pk_up$y_norm,
                      label = LABELS_UP[ctr]),
                 hjust = 1, vjust = 0.5, size = ANNOT_SZ,
                 color = "white", fontface = "bold",
                 fill = scales::alpha(unname(DIR_COLORS["Up"]), 0.85),
                 linewidth = 0, label.padding = unit(1, "mm"),
                 label.r = unit(0.6, "mm"))
  }

  p +
    scale_fill_manual(values = c(Up = unname(DIR_COLORS["Up"]),
                                  Down = unname(DIR_COLORS["Down"]))) +
    scale_color_manual(values = c(Up = unname(DIR_COLORS["Up"]),
                                   Down = unname(DIR_COLORS["Down"]))) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    coord_cartesian(xlim = c(-DENS_PAD, 1 + DENS_PAD),
                    ylim = c(TICK_DEPTH, 1.15), clip = "off") +
    labs(title = CTR_FACET[ctr], x = NULL, y = NULL) +
    FIG_THEME +
    theme(
      legend.position    = "none",
      axis.text.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(face = "bold", size = FIG_STRIP_SIZE,
                                         hjust = 0.5),
      plot.margin        = margin(2, 4, 2, 4)
    )
}

# --- Build 1x2 side-by-side ---
p_cvh <- make_panel("Cancer_vs_Healthy")
p_tr  <- make_panel("Training_CR")

pB_2x2 <- (p_cvh | p_tr) +
  plot_annotation(
    title    = "DEP Rank Location",
    subtitle = "\u03A0 \u2264 0.05 DEPs in t-stat-ranked proteome",
    tag_levels = list(c("B", "")),
    theme = theme(
      plot.title    = element_text(face = "bold", size = FIG_TITLE_SIZE),
      plot.subtitle = element_text(face = "bold.italic", size = FIG_SUBTITLE_SIZE,
                                    color = "grey30"),
      plot.tag      = element_text(face = "bold", size = 15)
    )
  ) &
  theme(plot.margin = margin(2, 4, 2, 4))

pB_2x2 <- pB_2x2 + plot_annotation(
  caption = "Rank position (by t-statistic)",
  theme = theme(plot.caption = element_text(face = "plain", size = 9,
                                             hjust = 0.5))
)

ggsave(file.path(RPT_SUPP_PDF, "panel_B_barcode_2x2_SUPP.pdf"), pB_2x2,
       width = PD_W, height = PD_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT_SUPP_PNG, "panel_B_barcode_2x2_SUPP.png"), pB_2x2,
       width = PD_W, height = PD_H, units = "mm", dpi = 300)

message("F03/CRvH Panel B (barcode side-by-side SUPP) done")
