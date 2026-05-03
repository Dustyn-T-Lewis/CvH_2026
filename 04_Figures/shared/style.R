# CvH Shared Figure Style
# Single source of truth: palettes, themes, sizing constants, helpers.

suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
})

# ── Palettes ──
GROUP_COLORS <- c(
  CRE_T1 = "#2166AC", CRE_T2 = "#67A9CF",
  PLA_T1 = "#D6604D", PLA_T2 = "#F4A582",
  H_T1   = "#4DAF4A")

DIR_COLORS <- c(Up = "#D6604D", Down = "#4393C3", NS = "grey70")

GROUP_FILL <- c(
  CRE_T1 = scales::alpha("#2166AC", 0.7),
  CRE_T2 = scales::alpha("#67A9CF", 0.7),
  PLA_T1 = scales::alpha("#D6604D", 0.7),
  PLA_T2 = scales::alpha("#F4A582", 0.7),
  H_T1   = scales::alpha("#4DAF4A", 0.7))

DB_COLORS <- c(
  Hallmark    = "#E41A1C", KEGG      = "#377EB8",
  Reactome    = "#4DAF4A", WikiPathways = "#984EA3",
  `GO:BP`     = "#FF7F00", BioCarta  = "#A65628",
  PID         = "#F781BF", GOSlim    = "#1B9E77")

# ── Sizing ──
PANEL_MD      <- 180
BASE_PATHWAY  <- 4.0
BASE_GENE     <- 3.2
BASE_STAT     <- 3.5
BASE_QUADRANT <- 4.0
BASE_COUNT    <- 3.5
BASE_TAG      <- 18

scale_text <- function(base_size, panel_width_mm, ref_width = PANEL_MD)
  base_size * sqrt(panel_width_mm / ref_width)

# ── Font sizing constants ──
FIG_TITLE_SIZE    <- 12
FIG_SUBTITLE_SIZE <- 9
FIG_STRIP_SIZE    <- 10
FIG_AXIS_TEXT     <- 8.5
FIG_LEGEND_TITLE  <- 9.5
FIG_LEGEND_TEXT   <- 8.5

# ── Theme ──
FIG_THEME <- theme_bw(base_size = 10) +
  theme(
    plot.title       = element_text(face = "bold", size = FIG_TITLE_SIZE),
    plot.subtitle    = element_text(face = "bold.italic", size = FIG_SUBTITLE_SIZE,
                                    colour = "grey30"),
    plot.tag         = element_text(face = "bold", size = 15),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = FIG_STRIP_SIZE),
    axis.title       = element_text(face = "bold", size = 10),
    axis.text        = element_text(size = FIG_AXIS_TEXT),
    legend.title     = element_text(face = "bold", size = FIG_LEGEND_TITLE),
    legend.text      = element_text(size = FIG_LEGEND_TEXT),
    legend.key.size  = unit(3, "mm"),
    panel.grid.minor = element_blank())

# ── Utility functions ──
get_pdf_device <- function() {
  if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
}

fmt_p <- function(p) {
  ifelse(p < 0.001, "< 0.001",
    ifelse(p < 0.01, sprintf("= %.3f", p),
      sprintf("= %.2f", p)))
}

sig_stars <- function(padj) {
  ifelse(padj < 0.001, "***",
    ifelse(padj < 0.01, "**",
      ifelse(padj < 0.05, "*", "ns")))
}

reorder_within <- function(x, by, within, fun = mean, sep = "___") {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

scale_y_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_y_discrete(labels = function(x) gsub(reg, "", x), ...)
}

darken_color <- function(col, factor = 0.7) {
  r <- grDevices::col2rgb(col)
  grDevices::rgb(r[1]*factor, r[2]*factor, r[3]*factor, maxColorValue = 255)
}

fisher_z_ci <- function(r, n, level = 0.95) {
  z  <- atanh(r)
  se <- 1 / sqrt(n - 3)
  q  <- qnorm((1 + level) / 2)
  lo <- tanh(z - q * se)
  hi <- tanh(z + q * se)
  c(lo = lo, hi = hi)
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
    stringr::str_replace("Oxidative Phosphorylation", "OXPHOS")
}

make_sigmoid_ribbon <- function(x0, x1, y0_top, y0_bot, y1_top, y1_bot,
                                n_pts = 50, ribbon_id) {
  t <- seq(0, 1, length.out = n_pts)
  blend <- (1 - cos(pi * t)) / 2
  tibble::tibble(
    x = c(x0 + (x1 - x0) * t, rev(x0 + (x1 - x0) * t)),
    y = c(y0_top + (y1_top - y0_top) * blend,
          rev(y0_bot + (y1_bot - y0_bot) * blend)),
    ribbon_id = ribbon_id
  )
}

# ── F03 (Supplement Concordance: CRE vs PLA) ──
classify_proteins_f3 <- function(pi_CRE, pi_PLA, pi_int, threshold = 0.05) {
  dplyr::case_when(
    pi_int < threshold                      ~ "Interaction",
    pi_CRE < threshold & pi_PLA < threshold ~ "Sig Both",
    pi_CRE < threshold                      ~ "Sig CRE only",
    pi_PLA < threshold                      ~ "Sig PLA only",
    TRUE                                    ~ "NS"
  ) |>
    factor(levels = c("Interaction", "Sig Both",
                       "Sig CRE only", "Sig PLA only", "NS"))
}

SIG_COLORS_F3 <- c(
  "Interaction"  = "#FF8F00",
  "Sig Both"     = "#2E7D32",
  "Sig CRE only" = "#2166AC",
  "Sig PLA only" = "#D6604D",
  "NS"           = "grey70"
)

SIG_LABEL_FILL_F3 <- c(
  "Interaction"  = scales::alpha("#FF8F00", 0.75),
  "Sig Both"     = scales::alpha("#2E7D32", 0.75),
  "Sig CRE only" = scales::alpha("#2166AC", 0.75),
  "Sig PLA only" = scales::alpha("#D6604D", 0.75),
  "NS"           = scales::alpha("grey70",  0.75)
)
SIG_LABEL_TEXT_F3 <- c(
  "Interaction"  = "white",
  "Sig Both"     = "white",
  "Sig CRE only" = "white",
  "Sig PLA only" = "white",
  "NS"           = "white"
)

ORA_QUAD_COLORS_F3 <- c(
  "Concordant Up"   = "#E57373",
  "Concordant Down" = "#64B5F6",
  "Discordant (CRE Up / PLA Down)" = "#FFB74D",
  "Discordant (CRE Down / PLA Up)" = "#81C784"
)

# ── F04 (Cancer Recovery: Cancer_vs_Healthy vs Training_CR) ──
classify_proteins_f4 <- function(pi_CvH, pi_TR, threshold = 0.05) {
  dplyr::case_when(
    pi_CvH < threshold & pi_TR < threshold ~ "Sig Both",
    pi_CvH < threshold                     ~ "Sig Cancer only",
    pi_TR < threshold                      ~ "Sig Training only",
    TRUE                                   ~ "NS"
  ) |>
    factor(levels = c("Sig Both",
                       "Sig Cancer only", "Sig Training only", "NS"))
}

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
  "NS"                = scales::alpha("grey70",  0.75)
)
SIG_LABEL_TEXT_F4 <- c(
  "Sig Both"          = "white",
  "Sig Cancer only"   = "white",
  "Sig Training only" = "white",
  "NS"                = "white"
)

ORA_QUAD_COLORS_F4 <- c(
  "Reversed (Cancer Up)"   = "#64B5F6",
  "Reversed (Cancer Down)" = "#E57373",
  "Exacerbated Up"         = "#FFB74D",
  "Exacerbated Down"       = "#81C784"
)

# ── Contrast palette (all 6 contrasts, both models) ──
CONTRAST_COLORS <- c(
  Cancer_vs_Healthy      = "#4CAF50",
  Training_CR            = "#9C27B0",
  Baseline_Supplement    = "#00897B",
  Training_CRE           = "#2166AC",
  Training_PLA           = "#D6604D",
  Supplement_Interaction = "#FF8F00")

# ── Contrast labels ──
CTR_SHORT <- c(
  Cancer_vs_Healthy      = "CR vs H",
  Training_CR            = "Tr.(CR)",
  Baseline_Supplement    = "BL(CRE\u2013PLA)",
  Training_CRE           = "Tr.(CRE)",
  Training_PLA           = "Tr.(PLA)",
  Supplement_Interaction = "CRE\u2013PLA")

CTR_FACET <- CTR_SHORT
CTR_AXIS  <- CTR_SHORT

# ── PCA palette (5 groups, F01) ──
PCA_COLORS <- c(
  CRE_T1 = "#2166AC", CRE_T2 = "#67A9CF",
  PLA_T1 = "#D6604D", PLA_T2 = "#F4A582",
  H_T1   = "#4DAF4A")

PCA_SHAPES <- c(
  CRE_T1 = 16, CRE_T2 = 17,
  PLA_T1 = 16, PLA_T2 = 17,
  H_T1   = 15)

# ── Supplement group labels (F01) ──
SUPP_LABELS <- c(CRE = "Creatine", PLA = "Placebo", H = "Healthy",
                 CR = "Cancer Recovery")

# ── F03 pattern colors (Supplement Concordance) ──
PATTERN_COLS_F3 <- c(
  "Shared"       = "#457B9D",
  "CRE-specific" = "#2166AC",
  "PLA-specific" = "#D6604D",
  "Interaction"  = "#FF8F00"
)
PATTERN_ORDER_F3 <- c("Shared", "CRE-specific", "PLA-specific", "Interaction")

PATTERN_SUBTITLES_F3 <- c(
  "Shared"       = "Sig CRE & PLA",
  "CRE-specific" = "Sig CRE only",
  "PLA-specific" = "Sig PLA only",
  "Interaction"  = "Sig Interaction"
)

# ── F04 pattern colors (Cancer Recovery) ──
PATTERN_COLS_F4 <- c(
  "Reversed"           = "#2563EB",
  "Partially Reversed" = "#64B5F6",
  "Persistent"         = "#FFB74D",
  "Exacerbated"        = "#DC2626")

PATTERN_ORDER_F4 <- c("Reversed", "Partially Reversed",
                       "Persistent", "Exacerbated")

# ── F05 model-specific contrast groupings ──
ALL_CONTRASTS_CRVH <- c("Cancer_vs_Healthy", "Training_CR")
ALL_CONTRASTS_CR   <- c("Baseline_Supplement", "Training_CRE",
                         "Training_PLA", "Supplement_Interaction")
ALL_CONTRASTS      <- c(ALL_CONTRASTS_CRVH, ALL_CONTRASTS_CR)

# ── F04 CRvH concordance quadrant colors ──
ORA_QUAD_COLORS_F4_CONC <- c(
  "Concordant Up"   = "#E57373",
  "Concordant Down" = "#64B5F6",
  "Discordant (Cancer Up / Training Down)" = "#FFB74D",
  "Discordant (Cancer Down / Training Up)" = "#81C784"
)

# ── F05 cancer-direction colors ──
CANCER_DIR_COLORS <- c(
  "Cancer Up"   = "#E57373",
  "Cancer Down" = "#64B5F6")
