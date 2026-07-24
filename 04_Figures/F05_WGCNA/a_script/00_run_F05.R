# Rebuild F05 end to end from the imputed DAList: network + module-trait audit,
# per-module ORA, the two main panels, and the composite.
setwd(here::here())

scripts <- c(
  "00_build_wgcna.R",
  "01_module_stats.R",
  "02_clustering.R"
)

for (s in scripts) {
  message("=== ", s, " ===")
  source(file.path("04_Figures/F05_WGCNA/a_script", s))
}
