# F04 CR: All Volcano Rings (2x2: Training_CRE, Training_PLA, Supplement_Interaction, Baseline_Supplement)
# Self-contained — generates all 4 volcano rings + 2x2 composite
# Outputs: panel_A-D individual PDFs/PNGs + F04_CR_volcanoes_MAIN composite
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F04/a_script/style.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(patchwork)
})

VW <- 190; VH <- 180
RPT <- "04_Figures/F04/CR/b_reports"
DAT <- "04_Figures/F04/CR/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

dep_df <- read_csv("03_DEP/c_data/03_combined_results_CR.csv", show_col_types = FALSE)

# --- fGSEA cache (from F03/CR) ---
fgsea_path <- "04_Figures/F03/CR/c_data/01_panel_C_fgsea_results.csv"
if (!file.exists(fgsea_path)) stop("fGSEA cache not found - run F03/CR first")
fgsea_all <- read_csv(fgsea_path, show_col_types = FALSE)

# Database selection for rings
avail_dbs <- unique(fgsea_all$database)
ring_dbs  <- intersect(c("Hallmark", "GO Slim", "GO:BP"), avail_dbs)
if (length(ring_dbs) == 0) ring_dbs <- avail_dbs[1:2]

# --- Contrast definitions (PLA = reference like Young, CRE = comparison like Old) ---
volcano_specs <- list(
  list(contrast = "Training_CRE",
       title    = "Training Response (Creatine)",
       subtitle = "CRE_T2 - CRE_T1",
       tag      = "A"),
  list(contrast = "Training_PLA",
       title    = "Training Response (Placebo)",
       subtitle = "PLA_T2 - PLA_T1",
       tag      = "B"),
  list(contrast = "Supplement_Interaction",
       title    = "Supplement x Training Interaction",
       subtitle = "Training_CRE - Training_PLA",
       tag      = "C"),
  list(contrast = "Baseline_Supplement",
       title    = "Baseline Supplement Effect",
       subtitle = "CRE_T1 - PLA_T1",
       tag      = "D")
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

  top_terms <- select_ring_terms(fgsea_all, spec$contrast, databases = ring_dbs)
  ring_data <- build_ring_with_gaps(top_terms, spec$contrast, fgsea_all,
                                     databases = ring_dbs)

  p <- make_volcano_ring(
    de_df              = dep_df,
    go_df              = fgsea_all,
    contrast           = spec$contrast,
    title              = NULL,
    contrast_title     = spec$title,
    contrast_subtitle  = enriched_sub,
    ring_data_override = ring_data,
    databases          = ring_dbs,
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

# --- 2x2 composite ---
composite <- (volcano_plots[[1]] | volcano_plots[[2]]) /
             (volcano_plots[[3]] | volcano_plots[[4]]) +
  plot_annotation(
    tag_levels = list(c("A", "B", "C", "D")),
    theme = theme(plot.margin = margin(2, 2, 2, 2, "mm"))
  )

COMP_W <- VW * 2 + 20  # ~400mm
COMP_H <- VH * 2 + 20  # ~380mm

ggsave(file.path(RPT, "F04_CR_volcanoes_MAIN.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "F04_CR_volcanoes_MAIN.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)

message("F04 CR volcanoes done (2x2)")
