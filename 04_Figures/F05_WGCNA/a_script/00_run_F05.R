# Rebuild F05 end to end from the imputed DAList: network + module-trait audit,
# per-module ORA, the two main panels, and the composite.
setwd(here::here())

scripts <- c(
  "CvH_WGCNA_run.R",
  "module_ora_prep.R",
  "02_member_response.R",
  "90_stitch_F05.R"
)

for (s in scripts) {
  message("=== ", s, " ===")
  source(file.path("04_Figures/F05_WGCNA/a_script", s))
}
