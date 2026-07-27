# F01 Supplementary Panel G: Grip Strength (Pre/Post x Supplement)
source("04_Figures/F01_Phenotype/a_script/supp/pre_post_panel.R")

pG <- pre_post_panel(
  pre_col = "pre_grip_lbs", post_col = "post_grip_lbs",
  title = "Grip Strength", y_lab = "Grip Strength (lbs)",
  tag = "G", stem = "panel_G_grip"
)
