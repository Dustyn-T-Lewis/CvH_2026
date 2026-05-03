# F06 Supplement Reversal — Figure-specific style (CR model only)
# Baseline_Supplement (x) vs Supplement_Interaction (y)
# Does the differential training response reverse baseline supplement differences?
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

# Significance colors
SIG_COLORS_F06 <- c(
  "Sig Both"             = "#2E7D32",
  "Sig Baseline only"    = "#E6AB02",
  "Sig Interaction only" = "#9B7FBF",
  "NS"                   = "grey70"
)
SIG_LABEL_FILL_F06 <- c(
  "Sig Both"             = scales::alpha("#2E7D32", 0.80),
  "Sig Baseline only"    = scales::alpha("#E6AB02", 0.80),
  "Sig Interaction only" = scales::alpha("#9B7FBF", 0.80),
  "NS"                   = scales::alpha("grey70",  0.75)
)
SIG_LABEL_TEXT_F06 <- c(
  "Sig Both"             = "white",
  "Sig Baseline only"    = "black",
  "Sig Interaction only" = "white",
  "NS"                   = "white"
)

# Reversal quadrant colors
ORA_QUAD_COLORS_F06 <- c(
  "Reversed Up"     = "#4393C3",
  "Reversed Down"   = "#4393C3",
  "Exacerbated Up"  = "#D6604D",
  "Exacerbated Down"= "#D6604D"
)

# Classification function
classify_supplement_reversal <- function(pi_base, pi_int, threshold = 0.05) {
  dplyr::case_when(
    pi_base < threshold & pi_int < threshold ~ "Sig Both",
    pi_base < threshold                      ~ "Sig Baseline only",
    pi_int  < threshold                      ~ "Sig Interaction only",
    TRUE                                     ~ "NS"
  )
}
