# F01 Supplementary Panel E: Chest Press (Pre/Post x Supplement)
source("04_Figures/F01_Phenotype/a_script/supp/pre_post_panel.R")

pE <- pre_post_panel(
  pre_col = "pre_chest_press_lbs", post_col = "post_chest_press_lbs",
  title = "Chest Press", y_lab = "Chest Press (lbs)",
  tag = "E", stem = "panel_E_chest_press"
)
