# F04 CRvH: All Volcano Rings (1x2: Cancer_vs_Healthy, Training_CR)
# Self-contained — generates both volcano rings + 1x2 composite
# Outputs: panel_A/B individual PDFs/PNGs + F04_CRvH_volcanoes_MAIN composite
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F04/a_script/style.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

VW <- 190; VH <- 180
RPT <- "04_Figures/F04/CRvH/b_reports"
DAT <- "04_Figures/F04/CRvH/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# DEP results: new proteoDA long output -> wide per-contrast columns the panel expects.
# Contrast labels recoded to the figure's display names (Cancer_vs_Healthy, Training_CR).
dep_df <- read_csv("03_DEP/a_non_imputed/c_data/combined_results_pi.csv",
                   show_col_types = FALSE) |>
  mutate(contrast = recode(contrast,
                           CRvH_Baseline = "Cancer_vs_Healthy",
                           CR_Training   = "Training_CR")) |>
  pivot_wider(id_cols = c(uniprot_id, gene, protein, description),
              names_from = contrast,
              values_from = c(logFC, t, P.Value, adj.P.Val, pi_score, sig_pi),
              names_glue = "{.value}_{contrast}")

# --- fGSEA cache (shared, regenerated from current DEP t-stats) ---
fgsea_all <- read_csv("04_Figures/shared/fgsea_CRvH.csv", show_col_types = FALSE)

# --- Contrast definitions ---
volcano_specs <- list(
  list(contrast = "Cancer_vs_Healthy",
       title    = "Cancer vs Healthy",
       subtitle = "CR_T1 - H_T1",
       tag      = "A"),
  list(contrast = "Training_CR",
       title    = "Training Response (CR)",
       subtitle = "CR_T2 - CR_T1",
       tag      = "B")
)

# --- Generate each volcano ring ---
volcano_plots <- lapply(volcano_specs, function(spec) {
  message(sprintf("  Building volcano: %s", spec$contrast))

  # Enrich subtitle with DEP + pathway counts
  pi_col <- paste0("pi_score_", spec$contrast)
  n_dep  <- if (pi_col %in% names(dep_df)) sum(dep_df[[pi_col]] < 0.05, na.rm = TRUE) else 0
  n_path <- sum(!is.na(fgsea_all$padj) & fgsea_all$padj < 0.05 &
                fgsea_all$contrast == spec$contrast)
  enriched_sub <- sprintf("%s | %d DEPs, %d pathways", spec$subtitle, n_dep, n_path)

  top_terms <- select_ring_terms(fgsea_all, spec$contrast)
  ring_data <- build_ring_with_gaps(top_terms, spec$contrast, fgsea_all)

  p <- make_volcano_ring(
    de_df              = dep_df,
    go_df              = fgsea_all,
    contrast           = spec$contrast,
    title              = NULL,
    contrast_title     = spec$title,
    contrast_subtitle  = enriched_sub,
    ring_data_override = ring_data,
    label_size         = scale_text(BASE_PATHWAY, VW),
    title_size         = scale_text(BASE_TAG, VW),
    point_size         = 1.2,
    point_alpha        = 0.55,
    count_label_size   = scale_text(BASE_COUNT, VW)
  )

  # Save individual panel
  fname <- tolower(gsub("[^a-zA-Z0-9]", "_", spec$title))
  ggsave(file.path(RPT, sprintf("panel_%s_%s_MAIN.pdf", spec$tag, fname)),
         p, width = VW, height = VH, units = "mm", device = pdf_device)
  ggsave(file.path(RPT, sprintf("panel_%s_%s_MAIN.png", spec$tag, fname)),
         p, width = VW, height = VH, units = "mm", dpi = 300)

  # Save ring data
  ring_out <- attr(p, "ring_data")
  if (!is.null(ring_out) && nrow(ring_out) > 0) {
    dir.create(file.path(DAT, paste0("panel_", spec$tag)), showWarnings = FALSE, recursive = TRUE)
    write_csv(ring_out %>% dplyr::select(-gene_list),
              file.path(DAT, paste0("panel_", spec$tag), "ring_terms.csv"))
  }

  p
})

# --- 1x2 composite ---
composite <- volcano_plots[[1]] | volcano_plots[[2]]
composite <- composite +
  plot_annotation(
    tag_levels = list(c("A", "B")),
    theme = theme(plot.margin = margin(2, 2, 2, 2, "mm"))
  )

COMP_W <- VW * 2 + 15  # ~395mm
COMP_H <- VH + 10      # ~190mm

ggsave(file.path(RPT, "F04_CRvH_volcanoes_MAIN.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "F04_CRvH_volcanoes_MAIN.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message("F04 CRvH volcanoes done (1x2)")
