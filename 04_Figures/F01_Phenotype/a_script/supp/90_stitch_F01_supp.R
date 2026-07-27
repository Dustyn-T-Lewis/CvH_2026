# F01 phenotype supplement: sources each pre/post outcome panel (which also self-saves
# under b_reports/supp) and assembles them into one supplement composite.

setwd(here::here())
source("04_Figures/shared/style.R")
pacman::p_load(patchwork)

panel_dir <- "04_Figures/F01_Phenotype/a_script/supp"
for (p in c("panel_D.R", "panel_E.R", "panel_F.R", "panel_G.R")) {
  source(file.path(panel_dir, p))
}

# Each panel is already a two-plot patchwork carrying its own plot_layout(widths=);
# wrap_elements keeps each one intact instead of flattening it into this grid.
composite <- (wrap_elements(pD) | wrap_elements(pE)) /
  (wrap_elements(pF) | wrap_elements(pG))

RPT <- "04_Figures/F01_Phenotype/b_reports/supp"
for (sub in c("pdf", "png")) {
  dir.create(file.path(RPT, sub), recursive = TRUE, showWarnings = FALSE)
}
ggsave(file.path(RPT, "png/F01_phenotype_supplement_SUPP.png"), composite,
  width = 340, height = 240, units = "mm", dpi = 300, bg = "white", limitsize = FALSE
)
ggsave(file.path(RPT, "pdf/F01_phenotype_supplement_SUPP.pdf"), composite,
  width = 340, height = 240, units = "mm", device = get_pdf_device(), limitsize = FALSE
)

message("F01 phenotype supplement saved (4 outcome panels)")
