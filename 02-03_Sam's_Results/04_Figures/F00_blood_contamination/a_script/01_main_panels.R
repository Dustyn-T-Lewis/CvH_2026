# F00 main-panels driver. Sources each _panel_*.R in order; each script
# leaves its ggplot object in .GlobalEnv as panel_A/B/C/D.

setwd(rprojroot::find_rstudio_root_file())
base_dir <- "02-03_Sam's_Results/04_Figures/F00_blood_contamination/a_script"

source(file.path(base_dir, "_panel_A_blood_stacked.R"))
source(file.path(base_dir, "_panel_B_bm_ratio_dotplot.R"))
source(file.path(base_dir, "_panel_C_bm_vs_mahalanobis.R"))
source(file.path(base_dir, "_panel_D_filter_venns.R"))

stopifnot(exists("panel_A"), exists("panel_B"),
          exists("panel_C"), exists("panel_D"))
message("F00 main panels: all 4 ggplot objects built.")
