# F05_WGCNA driver (Mito F03 card layout). One row per module across five aligned
# columns: protein count, module-trait heatmap, member-response fry tiles,
# eigengene trajectory with post-hoc brackets, top-5 ORA pathways. The main card
# holds modules with a nominal group effect (omnibus F p < 0.05); the rest go to
# the supplement.

setwd(here::here())
source("04_Figures/F05_WGCNA/a_script/style.R")
source("04_Figures/F05_WGCNA/a_script/card_panels.R")
pacman::p_load(readr, dplyr, patchwork)

DAT <- "04_Figures/F05_WGCNA/c_data"
RPT_PDF <- "04_Figures/F05_WGCNA/b_reports/main/pdf"
RPT_PNG <- "04_Figures/F05_WGCNA/b_reports/main/png"
SUPP_PDF <- "04_Figures/F05_WGCNA/b_reports/supp/pdf"
SUPP_PNG <- "04_Figures/F05_WGCNA/b_reports/supp/png"
for (d in c(RPT_PDF, RPT_PNG, SUPP_PDF, SUPP_PNG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

mod_bio <- read_csv(file.path(DAT, "mod_bio_labels.csv"), show_col_types = FALSE)
lmm <- read_csv(file.path(DAT, "wgcna_lmm_contrast_audit.csv"), show_col_types = FALSE)
settests <- read_csv(file.path(DAT, "module_set_tests.csv"), show_col_types = FALSE)
omnibus <- read_csv(file.path(DAT, "module_omnibus_F.csv"), show_col_types = FALSE)
traj <- read_csv(file.path(DAT, "trajectory_eigengenes.csv"), show_col_types = FALSE) |>
  filter(module_color != "grey")
ora <- read_csv(file.path(DAT, "module_ora.csv"), show_col_types = FALSE)

mod_order <- mod_bio |>
  arrange(desc(n_proteins)) |>
  pull(module_color)
responsive <- omnibus |>
  filter(p < 0.05) |>
  pull(module_color)
main_mods <- mod_order[mod_order %in% responsive]
supp_mods <- setdiff(mod_order, main_mods)

save_card <- function(modules, title, subtitle, stem, dir_png, dir_pdf) {
  card <- assemble_card(mod_bio, lmm, settests, traj, ora, modules) +
    plot_annotation(
      title = title, subtitle = subtitle,
      theme = theme(
        plot.title = element_text(face = "bold", size = 15, colour = "grey10"),
        plot.subtitle = element_text(face = "italic", size = 8.5, colour = "grey40")
      )
    )
  h <- 44 + 34 * length(modules)
  ggsave(file.path(dir_png, paste0(stem, ".png")), card,
    width = 300, height = h, units = "mm", dpi = 300, bg = "white", limitsize = FALSE
  )
  ggsave(file.path(dir_pdf, paste0(stem, ".pdf")), card,
    width = 300, height = h, units = "mm", device = pdf_device, limitsize = FALSE
  )
}

save_card(
  main_mods,
  "WGCNA co-expression modules and their eigengene response",
  "Each row is one module; main = modules with a nominal group effect (omnibus F p < 0.05), the rest supplementary.",
  "F05_wgcna_MAIN", RPT_PNG, RPT_PDF
)
message(sprintf("F05 main card saved (%d modules)", length(main_mods)))

if (length(supp_mods)) {
  save_card(
    supp_mods,
    "WGCNA modules without a nominal group effect (omnibus F p ≥ 0.05)",
    "Same layout as the main figure. These modules carry no detectable group effect on their eigengene.",
    "SUPP_F05_wgcna_other", SUPP_PNG, SUPP_PDF
  )
  message(sprintf("F05 supplement card saved (%d modules)", length(supp_mods)))
}
