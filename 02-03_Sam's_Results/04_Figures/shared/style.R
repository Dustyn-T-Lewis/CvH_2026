# Sam's-Results figure style.
# Sources the canonical CvH pipeline style so palettes, themes, and sizing
# constants stay aligned with the main 04_Figures/ output.
#
# To diverge intentionally (e.g., add a Sam-specific palette), override
# the relevant constant *after* the source() call below.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
