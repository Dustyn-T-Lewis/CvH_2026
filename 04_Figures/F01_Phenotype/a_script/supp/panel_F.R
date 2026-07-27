# F01 Supplementary Panel F: Leg Extension (Pre/Post x Supplement)
source("04_Figures/F01_Phenotype/a_script/supp/pre_post_panel.R")

pF <- pre_post_panel(
  pre_col = "pre_leg_ext_lbs", post_col = "post_leg_ext_lbs",
  title = "Leg Extension", y_lab = "Leg Extension (lbs)",
  tag = "F", stem = "panel_F_leg_ext"
)
