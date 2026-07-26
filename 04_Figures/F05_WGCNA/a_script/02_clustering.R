# F05 stage 02: render driver. Caches per-module ORA, assembles the module card
# (main = omnibus-F p < 0.05, rest supplementary), renders the construction and
# phenotype-method supplements, and writes one supplementary workbook.

setwd(here::here())
source("04_Figures/F05_WGCNA/a_script/style.R")
source("04_Figures/F05_WGCNA/a_script/panels/module_card.R")
source("04_Figures/F05_WGCNA/a_script/panels/construction.R")
source("04_Figures/F05_WGCNA/a_script/supp/construction.R")
source("04_Figures/F05_WGCNA/a_script/supp/supplement_arm.R")
source("04_Figures/shared/pathway_utils.R")
source("04_Figures/shared/figure_supplement_helpers.R")
pacman::p_load(readr, dplyr, tidyr, purrr, patchwork)

DAT <- "04_Figures/F05_WGCNA/c_data"
RPT_PNG <- "04_Figures/F05_WGCNA/b_reports/main/png"
RPT_PDF <- "04_Figures/F05_WGCNA/b_reports/main/pdf"
SUPP_PNG <- "04_Figures/F05_WGCNA/b_reports/supp/png"
SUPP_PDF <- "04_Figures/F05_WGCNA/b_reports/supp/pdf"
for (d in c(RPT_PNG, RPT_PDF, SUPP_PNG, SUPP_PDF)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
pdf_device <- get_pdf_device()

# --- per-module ORA cache ---
module_df <- read_csv(file.path(DAT, "wgcna_module_assignments.csv"), show_col_types = FALSE)
mod_bio <- read_csv(file.path(DAT, "mod_bio_labels.csv"), show_col_types = FALSE)
universe <- unique(na.omit(module_df$gene))
pw <- build_pathway_collection(min_size = 15, max_size = 500, include_goslim = FALSE)
ora <- map_dfr(mod_bio$module_color, function(mod) {
  genes <- unique(na.omit(module_df$gene[module_df$module_color == mod]))
  res <- run_ora_deduplicated(
    genes = genes, universe = universe, pathways = pw,
    jaccard_cutoff = 0.5, min_size = 10, max_size = 500, padj_cutoff = 1
  )
  if (is.null(res) || nrow(res) == 0) {
    return(tibble())
  }
  res |>
    mutate(module_color = mod) |>
    select(module_color, pathway, database, pval, padj, overlap, size, odds_ratio)
})
write_csv(ora, file.path(DAT, "module_ora.csv"))

# --- assemble the card ---
nes <- read_csv(file.path(DAT, "module_fgsea_nes.csv"), show_col_types = FALSE)
settests <- read_csv(file.path(DAT, "module_set_tests.csv"), show_col_types = FALSE)
omnibus <- read_csv(file.path(DAT, "module_omnibus_F.csv"), show_col_types = FALSE)
lmm <- read_csv(file.path(DAT, "wgcna_lmm_contrast_audit.csv"), show_col_types = FALSE)
cor_pheno <- read_csv(file.path(DAT, "module_trait_cor.csv"), show_col_types = FALSE)
traj <- read_csv(file.path(DAT, "trajectory_eigengenes.csv"), show_col_types = FALSE) |>
  filter(module_color != "grey")

mod_order <- mod_bio |>
  arrange(desc(n_proteins)) |>
  pull(module_color)
main_mods <- mod_order
supp_mods <- setdiff(mod_order, main_mods)

save_card <- function(modules, title, subtitle, stem, dir_png, dir_pdf) {
  card <- assemble_card(mod_bio, nes, settests, traj, lmm, ora, modules) +
    plot_annotation(
      title = title, subtitle = subtitle,
      theme = theme(
        plot.title = element_text(face = "bold", size = 15, colour = "grey10"),
        plot.subtitle = element_text(face = "italic", size = 8.5, colour = "grey40")
      )
    )
  h <- 48 + 35 * length(modules)
  ggsave(file.path(dir_png, paste0(stem, ".png")), card,
    width = 300, height = h, units = "mm", dpi = 300, bg = "white", limitsize = FALSE
  )
  ggsave(file.path(dir_pdf, paste0(stem, ".pdf")), card,
    width = 300, height = h, units = "mm", device = pdf_device, limitsize = FALSE
  )
}

if (length(main_mods)) {
  save_card(
    main_mods,
    "WGCNA co-expression modules and their eigengene response across the design",
    "Rows: all co-expression modules, ordered by size. Trajectory brackets = post-hoc contrasts; ✱ = FDR < 0.05.",
    "F05_wgcna_MAIN", RPT_PNG, RPT_PDF
  )
  message(sprintf("F05 main card saved (%d modules)", length(main_mods)))
}
if (length(supp_mods)) {
  save_card(
    supp_mods,
    "WGCNA modules without a nominal group effect (omnibus F p ≥ 0.05)",
    "Same layout as the main figure.",
    "SUPP_F05_wgcna_other", SUPP_PNG, SUPP_PDF
  )
  message(sprintf("F05 supplement card saved (%d modules)", length(supp_mods)))
}

# --- supplements ---
cor_matched <- read_csv(file.path(DAT, "module_trait_cor_matched.csv"), show_col_types = FALSE)
lmm_pheno <- read_csv(file.path(DAT, "module_trait_lmm.csv"), show_col_types = FALSE)
render_construction_supp(DAT, SUPP_PNG, SUPP_PDF, pdf_device)
render_phenotype_supp(cor_pheno, cor_matched, lmm_pheno, SUPP_PNG, SUPP_PDF, pdf_device)
render_supplement_arm(DAT, SUPP_PNG, SUPP_PDF, pdf_device)

# --- one supplementary workbook ---
overview <- tibble::tribble(
  ~Sheet, ~Contents,
  "module_summary", "Module id, colour, protein count.",
  "eigengene_lmm", "Eigengene ME ~ group_time contrasts (r_equiv, Kenward-Roger).",
  "fry_camera", "fry (gate) and camera (reported) per module per contrast.",
  "module_omnibus_F", "Eigengene omnibus moderated-F; main-figure gate.",
  "module_fgsea_nes", "Per-module fGSEA NES, modules as gene sets ranked by DE t.",
  "phenotype_baseline", "Primary: T1 eigengene vs pre_ outcomes, corPvalueStudent.",
  "phenotype_matched", "Sensitivity: all samples, matched-timepoint Pearson.",
  "phenotype_lmm", "Sensitivity: eigengene ~ trait + (1|subject), std beta.",
  "module_ora", "Per-module over-representation, Jaccard-deduplicated."
)
build_workbook(
  file.path(DAT, "F05_supplementary.xlsx"),
  "F05 WGCNA modules — supplementary tables",
  "Signed bicor network on the missForest-imputed matrix. Card: fGSEA NES over fry, baseline phenotype r, eigengene trajectory, top-5 ORA.",
  overview,
  list(
    list(name = "module_summary", df = mod_bio),
    list(name = "eigengene_lmm", df = lmm),
    list(name = "fry_camera", df = settests),
    list(name = "module_omnibus_F", df = omnibus),
    list(name = "module_fgsea_nes", df = nes),
    list(name = "phenotype_baseline", df = cor_pheno),
    list(name = "phenotype_matched", df = cor_matched),
    list(name = "phenotype_lmm", df = lmm_pheno),
    list(name = "module_ora", df = ora)
  )
)

message("F05 render complete")
