# Reversal — Master Orchestrator
# Sources supp panels first (enrichment_heatmap not used), then supp stitcher,
# then main stitcher. Final cleanup.
setwd(rprojroot::find_rstudio_root_file())

message("===== Reversal: Master Orchestrator =====")

# Main composite FIRST (before supp loads AnnotationDbi which masks dplyr::select)
message("\n--- Step 1: Main composite ---")
source("04_Figures/Reversal/a_script/main/90_stitch_main.R")

# Supplementary composite (AnnotationDbi masking OK — supp scripts use dplyr:: prefix)
message("\n--- Step 2: Supplementary panels ---")
source("04_Figures/Reversal/a_script/supp/90_stitch_supp.R")

message("\n===== Reversal figure pipeline complete =====")
