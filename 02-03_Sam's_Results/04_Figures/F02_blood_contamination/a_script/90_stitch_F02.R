# F02 composite stitcher.
# Composes MAIN (2x2: A,B / C,D) and SUPP (2x2: A,B / C,D) into PDFs + PNGs.

suppressPackageStartupMessages({
  library(cowplot); library(ggplot2)
})

setwd(rprojroot::find_rstudio_root_file())
base_dir <- "02-03_Sam's_Results/04_Figures/F02_blood_contamination/a_script"

source(file.path(base_dir, "01_main_panels.R"))
source(file.path(base_dir, "02_supp_panels.R"))

main_composite <- plot_grid(
  panel_A, panel_B,
  panel_C, panel_D,
  ncol = 2, labels = c("A", "B", "C", "D"),
  label_size = 14, align = "hv"
)

supp_composite <- plot_grid(
  supp_A, supp_B,
  supp_C, supp_D,
  ncol = 2, labels = c("A", "B", "C", "D"),
  label_size = 14, align = "hv"
)

out_main_pdf <- "02-03_Sam's_Results/04_Figures/F02_blood_contamination/b_reports/main/pdf/MAIN_F02_blood_contamination.pdf"
out_main_png <- "02-03_Sam's_Results/04_Figures/F02_blood_contamination/b_reports/main/png/MAIN_F02_blood_contamination.png"
out_supp_pdf <- "02-03_Sam's_Results/04_Figures/F02_blood_contamination/b_reports/supp/pdf/SUPP_F02_blood_contamination.pdf"
out_supp_png <- "02-03_Sam's_Results/04_Figures/F02_blood_contamination/b_reports/supp/png/SUPP_F02_blood_contamination.png"

ggsave(out_main_pdf, main_composite, width = 12, height = 9, units = "in")
ggsave(out_main_png, main_composite, width = 12, height = 9, units = "in", dpi = 300, bg = "white")
ggsave(out_supp_pdf, supp_composite, width = 12, height = 9, units = "in")
ggsave(out_supp_png, supp_composite, width = 12, height = 9, units = "in", dpi = 300, bg = "white")

message("F02 composites written:\n  ", out_main_pdf, "\n  ", out_supp_pdf)
