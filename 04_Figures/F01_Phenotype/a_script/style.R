# F01 Phenotype — Figure-specific style
source("04_Figures/shared/style.R")

# Supplement colors for phenotype bars
SUPP_FILL <- c(
  CRE_T1 = scales::alpha("#2166AC", 0.5), CRE_T2 = scales::alpha("#2166AC", 0.9),
  PLA_T1 = scales::alpha("#D6604D", 0.5), PLA_T2 = scales::alpha("#D6604D", 0.9),
  H_T1   = scales::alpha("#4DAF4A", 0.7))

SUPP_COLORS <- c(CRE = "#2166AC", PLA = "#D6604D", H = "#4DAF4A")
