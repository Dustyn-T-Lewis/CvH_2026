# F07 CR Training Concordance — Figure-specific style
# Training_CRE (x) vs Training_PLA (y)
# Do CRE and PLA groups show concordant training responses?
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

# ── F07 contrast colors ──
CONTRAST_COLORS_F07 <- c(
  Training_CRE           = "#2166AC",
  Training_PLA           = "#D6604D",
  Supplement_Interaction = "#FF8F00",
  Baseline_Supplement    = "#00897B"
)

# ── Significance classification ──
classify_proteins_f07 <- function(pi_CRE, pi_PLA, pi_int, threshold = 0.05) {
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

SIG_COLORS_F07 <- c(
  "Interaction"  = "#FF8F00",
  "Sig Both"     = "#2E7D32",
  "Sig CRE only" = "#2166AC",
  "Sig PLA only" = "#D6604D",
  "NS"           = "grey70"
)

SIG_LABEL_FILL_F07 <- c(
  "Interaction"  = scales::alpha("#FF8F00", 0.75),
  "Sig Both"     = scales::alpha("#2E7D32", 0.75),
  "Sig CRE only" = scales::alpha("#2166AC", 0.75),
  "Sig PLA only" = scales::alpha("#D6604D", 0.75),
  "NS"           = scales::alpha("grey70",  0.75)
)
SIG_LABEL_TEXT_F07 <- c(
  "Interaction"  = "white",
  "Sig Both"     = "white",
  "Sig CRE only" = "white",
  "Sig PLA only" = "white",
  "NS"           = "white"
)

# ── Concordance quadrant colors ──
ORA_QUAD_COLORS_F07 <- c(
  "Concordant Up"                    = "#E57373",
  "Concordant Down"                  = "#64B5F6",
  "Discordant (CRE Up / PLA Down)"  = "#FFB74D",
  "Discordant (CRE Down / PLA Up)"  = "#81C784"
)

# ── Display labels for ORA bars ──
DISPLAY_LABELS_F07 <- c()
