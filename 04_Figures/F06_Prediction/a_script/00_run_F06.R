# Rebuild F06 end to end: feature spaces, prediction engines, circularity ladder, figures.
setwd(here::here())

scripts <- c(
  "00_features.R",
  "01_classify.R",
  "supp/circularity_ladder.R",
  "02_render.R"
)

for (s in scripts) {
  message("=== ", s, " ===")
  source(file.path("04_Figures/F06_Prediction/a_script", s))
}
