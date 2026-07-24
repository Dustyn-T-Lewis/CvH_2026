# F02 QC supplement: sources each QC panel (which also self-saves under
# b_reports/supp/panels) and assembles them into one supplement composite.

setwd(here::here())
source("04_Figures/shared/style.R")
pacman::p_load(cowplot)

# panel_A_ICC pulls in psych, whose alpha() (Cronbach's) masks scales::alpha;
# bind the scales version so bare alpha() in later panels resolves correctly.
alpha <- scales::alpha

panel_dir <- "04_Figures/F02_Proteome_Overview/a_script/supp"
for (p in c(
  "panel_cv_violin.R", "panel_cv_scatter.R", "panel_icc.R",
  "panel_dbrda.R", "panel_baseline_cv.R"
)) {
  source(file.path(panel_dir, p))
}

composite <- cowplot::plot_grid(
  pA, pB, pA_ICC, pE_dbrda, pF,
  ncol = 2, labels = "AUTO", label_fontface = "bold"
)

RPT <- "04_Figures/F02_Proteome_Overview/b_reports/supp"
for (sub in c("pdf", "png")) {
  dir.create(file.path(RPT, sub), recursive = TRUE, showWarnings = FALSE)
}
ggsave(file.path(RPT, "png/F02_QC_supplement_SUPP.png"), composite,
  width = 300, height = 340, units = "mm", dpi = 300, limitsize = FALSE
)
ggsave(file.path(RPT, "pdf/F02_QC_supplement_SUPP.pdf"), composite,
  width = 300, height = 340, units = "mm", device = get_pdf_device(), limitsize = FALSE
)

message("F02 QC supplement composite saved")
