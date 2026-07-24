# Global proteome overview (canonical F02): PCA, DEP counts, effect size, DEP
# overlap, direction, pathway enrichment. Main pools CR vs Healthy; the
# supplement variant splits CRE / PLA / Healthy across every panel.

setwd(here::here())
source("04_Figures/shared/style.R")
source("04_Figures/F02_Proteome_Overview/a_script/overview_panels.R")

RPT <- "04_Figures/F02_Proteome_Overview/b_reports"
for (sub in c("main/pdf", "main/png", "supp/pdf", "supp/png")) {
  dir.create(file.path(RPT, sub), recursive = TRUE, showWarnings = FALSE)
}
pdf_device <- get_pdf_device()

add_tag <- function(p, tag) {
  p + labs(tag = tag) +
    theme(plot.tag = element_text(face = "bold", size = 10))
}

# Compact styling for the dense 6-panel grid.
compact <- theme(
  plot.title = element_text(size = 8.5, margin = margin(b = 1)),
  plot.subtitle = element_text(size = 6, margin = margin(b = 2)),
  legend.key.size = unit(2, "mm"),
  legend.text = element_text(size = 6),
  plot.margin = margin(3, 3, 2, 2)
)

build_overview <- function(mode, ctr_map, pal) {
  dep <- load_dep(ctr_map)
  a <- add_tag(panel_pca(mode), "A")
  b <- add_tag(panel_dep_counts(dep, pal), "B")
  c <- add_tag(panel_effect(dep, pal), "C")
  d <- add_tag(panel_overlap(dep, pal), "D")
  e <- add_tag(panel_direction(dep, pal), "E")
  f <- add_tag(panel_pathway(ctr_map), "F")
  a + b + c + d + e + f +
    plot_layout(design = "AABBCC\nDDEEFF", widths = c(1.1, 1, 0.95), heights = c(1, 1)) &
    compact
}

main <- build_overview("main", MAIN_CTR, MAIN_PAL)
supp <- build_overview("supp", SUPP_CTR, SUPP_PAL)

W <- 210
H <- 165
ggsave(file.path(RPT, "main/png/F02_proteome_overview_MAIN.png"), main,
  width = W, height = H, units = "mm", dpi = 300
)
ggsave(file.path(RPT, "main/pdf/F02_proteome_overview_MAIN.pdf"), main,
  width = W, height = H, units = "mm", device = pdf_device
)
ggsave(file.path(RPT, "supp/png/F02_proteome_overview_SUPP.png"), supp,
  width = W, height = H, units = "mm", dpi = 300
)
ggsave(file.path(RPT, "supp/pdf/F02_proteome_overview_SUPP.pdf"), supp,
  width = W, height = H, units = "mm", device = pdf_device
)

message("F02 proteome overview (main + supp) saved")
