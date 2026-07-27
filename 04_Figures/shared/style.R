# CvH Shared Figure Style
# Single source of truth: palettes, themes, sizing constants, helpers.

pacman::p_load(ggplot2, scales)

# ── Palettes ──
DIR_COLORS <- c(Up = "#D6604D", Down = "#4393C3", NS = "grey70")

GROUP_FILL <- c(
  CRE_T1 = scales::alpha("#2166AC", 0.7),
  CRE_T2 = scales::alpha("#67A9CF", 0.7),
  PLA_T1 = scales::alpha("#D6604D", 0.7),
  PLA_T2 = scales::alpha("#F4A582", 0.7),
  H_T1   = scales::alpha("#4DAF4A", 0.7)
)

# Phenotype (F01) bar fills. The T1/T2 alpha ramp encodes pre vs post, so this is
# a deliberate ramp, not a duplicate of GROUP_FILL; SUPP_COLORS is the base hue.
SUPP_FILL <- c(
  CRE_T1 = scales::alpha("#2166AC", 0.5), CRE_T2 = scales::alpha("#2166AC", 0.9),
  PLA_T1 = scales::alpha("#D6604D", 0.5), PLA_T2 = scales::alpha("#D6604D", 0.9),
  H_T1   = scales::alpha("#4DAF4A", 0.7)
)
SUPP_COLORS <- c(CRE = "#2166AC", PLA = "#D6604D", H = "#4DAF4A")

# Database palettes, keyed to the fgsea `database` values. DB_COLORS is the dark
# canonical set (axis/label text, white in-bar fit-text); ORA_DB_COLORS is the
# light tint for bar fills that sit under dark labels.
DB_COLORS <- c(
  Hallmark = "#E41A1C", Reactome = "#4DAF4A", KEGG = "#377EB8",
  `GO:BP` = "#FF7F00", `GO Slim` = "#1B9E77",
  WikiPathways = "#984EA3", Other = "grey60"
)
DB_ORDER <- c("Hallmark", "Reactome", "KEGG", "GO:BP", "GO Slim")
ORA_DB_COLORS <- c(
  Hallmark = "#F4A7A6", Reactome = "#B4DDB2", KEGG = "#A9C4E0",
  `GO:BP` = "#FFD199", `GO Slim` = "#A6D9C6",
  WikiPathways = "#D3BCE0", Other = "grey80"
)

# ── Sizing ──
PANEL_MD <- 180
BASE_GENE <- 3.2
BASE_STAT <- 3.5
BASE_COUNT <- 3.5

scale_text <- function(base_size, panel_width_mm, ref_width = PANEL_MD) {
  base_size * sqrt(panel_width_mm / ref_width)
}

# ── Font sizing constants ──
FIG_TITLE_SIZE <- 12
FIG_SUBTITLE_SIZE <- 9
FIG_STRIP_SIZE <- 10
FIG_AXIS_TEXT <- 8.5
FIG_LEGEND_TITLE <- 9.5
FIG_LEGEND_TEXT <- 8.5

# ── Theme ──
# Sizes scale by panel width (scale_text): at the reference width (PANEL_MD) the
# function reproduces the base hierarchy, so a narrower panel gets smaller type.
# FIG_THEME is the default-width object every panel appends with + FIG_THEME.
theme_cvh <- function(base_size = 10, panel_width_mm = PANEL_MD) {
  s <- function(pt) scale_text(pt, panel_width_mm)
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = s(FIG_TITLE_SIZE)
      ),
      plot.subtitle = ggplot2::element_text(
        face = "bold.italic", size = s(FIG_SUBTITLE_SIZE),
        colour = "grey30"
      ),
      plot.tag = ggplot2::element_text(face = "bold", size = s(15)),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold", size = s(FIG_STRIP_SIZE)
      ),
      axis.title = ggplot2::element_text(face = "bold", size = s(10)),
      axis.text = ggplot2::element_text(size = s(FIG_AXIS_TEXT)),
      legend.title = ggplot2::element_text(
        face = "bold", size = s(FIG_LEGEND_TITLE)
      ),
      legend.text = ggplot2::element_text(size = s(FIG_LEGEND_TEXT)),
      legend.key.size = grid::unit(3, "mm"),
      panel.grid.minor = ggplot2::element_blank()
    )
}

FIG_THEME <- theme_cvh()

# Panel letter baked into the title at a constant gap, so every panel reads
# uniformly regardless of its y-axis width (steadier than a floating plot.tag).
add_tag <- function(p, tag) {
  cur <- p$labels$title
  p + ggplot2::labs(
    tag = NULL,
    title = paste0(tag, "  ", if (is.null(cur)) "" else cur)
  )
}

# ── Utility functions ──
get_pdf_device <- function() {
  cairo_ok <- capabilities("cairo") && tryCatch(
    {
      f <- tempfile(fileext = ".pdf")
      grDevices::cairo_pdf(f)
      grDevices::dev.off()
      unlink(f)
      TRUE
    },
    error = function(e) FALSE
  )
  if (cairo_ok) grDevices::cairo_pdf else grDevices::pdf
}

fmt_p <- function(p) {
  ifelse(p < 0.001, "< 0.001",
    ifelse(p < 0.01, sprintf("= %.3f", p),
      sprintf("= %.2f", p)
    )
  )
}





clean_pathway_name <- function(name, max_chars = NULL) {
  name |>
    stringr::str_remove("^HALLMARK_") |>
    stringr::str_remove("^GOSLIM_") |>
    stringr::str_remove("^GOBP_") |>
    stringr::str_remove("^GOCC_") |>
    stringr::str_remove("^GOMF_") |>
    stringr::str_remove("^REACTOME_") |>
    stringr::str_remove("^KEGG_MEDICUS_") |>
    stringr::str_remove("^KEGG_") |>
    stringr::str_replace_all("_", " ") |>
    stringr::str_to_title() |>
    stringr::str_remove("^Reference ") |>
    stringr::str_replace("Mtorc1", "mTORC1") |>
    stringr::str_replace("Myc ", "MYC ") |>
    stringr::str_replace("E2f ", "E2F ") |>
    stringr::str_replace("Dna ", "DNA ") |>
    stringr::str_replace("Rna ", "RNA ") |>
    stringr::str_replace("Tnfa ", "TNFa ") |>
    stringr::str_replace("Uv ", "UV ") |>
    stringr::str_replace("G2m ", "G2M ") |>
    stringr::str_replace("Il6 ", "IL6 ") |>
    stringr::str_replace("Il2 ", "IL2 ") |>
    stringr::str_replace("Kras ", "KRAS ") |>
    stringr::str_replace("P53 ", "p53 ") |>
    stringr::str_replace("Tgf ", "TGF ") |>
    stringr::str_replace("Nf Kb", "NF-kB") |>
    stringr::str_replace("Atp ", "ATP ") |>
    stringr::str_replace("Nadh ", "NADH ") |>
    stringr::str_replace("Oxidative Phosphorylation", "OXPHOS") |>
    stringr::str_replace("\\bTrna\\b", "tRNA") |>
    stringr::str_replace("\\bMhc\\b", "MHC") |>
    stringr::str_replace("\\bIqgaps?\\b", "IQGAPs") |>
    stringr::str_replace("\\bIv\\b", "IV") |>
    stringr::str_replace("\\bIii\\b", "III") |>
    stringr::str_replace("\\bIi\\b", "II")
}

make_sigmoid_ribbon <- function(x0, x1, y0_top, y0_bot, y1_top, y1_bot,
                                n_pts = 50, ribbon_id) {
  t <- seq(0, 1, length.out = n_pts)
  blend <- (1 - cos(pi * t)) / 2
  tibble::tibble(
    x = c(x0 + (x1 - x0) * t, rev(x0 + (x1 - x0) * t)),
    y = c(
      y0_top + (y1_top - y0_top) * blend,
      rev(y0_bot + (y1_bot - y0_bot) * blend)
    ),
    ribbon_id = ribbon_id
  )
}

# ── F04 Reversal (Cancer Recovery: CRvH_Baseline vs CR_Training) ──
SIG_COLORS_F4 <- c(
  "Sig Both"          = "#2E7D32",
  "Sig Cancer only"   = "#4CAF50",
  "Sig Training only" = "#9C27B0",
  "NS"                = "grey70"
)

SIG_LABEL_FILL_F4 <- c(
  "Sig Both"          = scales::alpha("#2E7D32", 0.75),
  "Sig Cancer only"   = scales::alpha("#4CAF50", 0.75),
  "Sig Training only" = scales::alpha("#9C27B0", 0.75),
  "NS"                = scales::alpha("grey70", 0.75)
)
SIG_LABEL_TEXT_F4 <- c(
  "Sig Both"          = "white",
  "Sig Cancer only"   = "white",
  "Sig Training only" = "white",
  "NS"                = "white"
)


# ── Contrast palette (all 6 contrasts, both models) ──
CONTRAST_COLORS <- c(
  CRvH_Baseline          = "#4CAF50",
  CR_Training            = "#9C27B0",
  Resid                  = "#795548",
  Baseline_Supplement    = "#00897B",
  Training_CRE           = "#2166AC",
  Training_PLA           = "#D6604D",
  Supplement_Interaction = "#FF8F00"
)

# ── Contrast labels ──


# ── PCA palette (5 groups, F01) ──
PCA_COLORS <- c(
  CRE_T1 = "#2166AC", CRE_T2 = "#67A9CF",
  PLA_T1 = "#D6604D", PLA_T2 = "#F4A582",
  H_T1   = "#4DAF4A"
)

PCA_SHAPES <- c(
  CRE_T1 = 16, CRE_T2 = 17,
  PLA_T1 = 16, PLA_T2 = 17,
  H_T1   = 15
)

# ── Supplement group labels (F01) ──

# ── F04 Reversal pattern colors (Cancer Recovery) ──


# ── F05 model-specific contrast groupings ──

# ── F04 CRvH concordance quadrant colors ──

# ── F05 cancer-direction colors ──
CANCER_DIR_COLORS <- c(
  "Cancer Up"   = "#E57373",
  "Cancer Down" = "#64B5F6"
)
