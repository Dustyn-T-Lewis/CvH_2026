# F02 supp-panels driver. Sources each _supp_panel_*.R; objects:
# supp_A / supp_B / supp_C / supp_D.

setwd(rprojroot::find_rstudio_root_file())
base_dir <- "02-03_Sam's_Results/04_Figures/F02_blood_contamination/a_script"

source(file.path(base_dir, "_supp_panel_A_marker_intensity.R"))
source(file.path(base_dir, "_supp_panel_B_cutoff_sensitivity.R"))
source(file.path(base_dir, "_supp_panel_C_corr_heatmap.R"))
source(file.path(base_dir, "_supp_panel_D_blood_markers_in_DEP.R"))

stopifnot(exists("supp_A"), exists("supp_B"),
          exists("supp_C"), exists("supp_D"))
message("F02 supp panels: all 4 ggplot objects built.")
