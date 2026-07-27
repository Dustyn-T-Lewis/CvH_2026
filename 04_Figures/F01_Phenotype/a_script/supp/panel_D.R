# F01 Supplementary Panel D: LBM (Pre/Post x Supplement)
source("04_Figures/F01_Phenotype/a_script/supp/pre_post_panel.R")

pD <- pre_post_panel(
  pre_col = "pre_LBM_kg", post_col = "post_LBM_kg",
  title = "Lean Body Mass", y_lab = "LBM (kg)",
  tag = "D", stem = "panel_D_lbm"
)
