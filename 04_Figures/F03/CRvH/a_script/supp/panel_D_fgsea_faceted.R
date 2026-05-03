# F03/CRvH — Supp: fGSEA Faceted Bar Chart (Pathway Enrichment)
# Per-database BH faceted bars (supplementary figure)
# DISPLAY ONLY — reads fGSEA cache from panel_C.R (run that first)
# CRvH model: 2 contrasts
# Outputs: supp/panel_D_fgsea_faceted_SUPP.pdf/.png

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F03/a_script/style.R")

library(dplyr)
library(readr)
library(tidyr)

DAT <- "04_Figures/F03/CRvH/c_data"
RPT <- "04_Figures/F03/CRvH/b_reports/supp"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

fgsea_raw <- read_csv(file.path(DAT, "01_panel_C_fgsea_results.csv"),
                      show_col_types = FALSE)
pdf_device <- get_pdf_device()

DISPLAY_DBS       <- c("Hallmark", "GO Slim", "GO:BP", "KEGG", "Reactome")
DISPLAY_CONTRASTS <- c("Cancer_vs_Healthy", "Training_CR")
PC_W <- 110

db_ctr_totals <- fgsea_raw |>
  filter(database %in% DISPLAY_DBS, contrast %in% DISPLAY_CONTRASTS) |>
  group_by(contrast, database) |>
  summarise(n_total = n(), .groups = "drop")

count_df <- fgsea_raw |>
  filter(!is.na(padj), padj < 0.05,
         database %in% DISPLAY_DBS, contrast %in% DISPLAY_CONTRASTS) |>
  group_by(contrast, database) |>
  summarise(
    Up   = sum(NES > 0),
    Down = sum(NES < 0),
    .groups = "drop"
  ) |>
  pivot_longer(cols = c(Up, Down), names_to = "direction",
               values_to = "count") |>
  left_join(db_ctr_totals, by = c("contrast", "database")) |>
  mutate(fraction = count / n_total)

nonempty_dbs <- count_df |>
  group_by(database) |> filter(sum(count) > 0) |> pull(database) |> unique()
count_df <- count_df |> filter(database %in% nonempty_dbs)

db_label_n <- db_ctr_totals |>
  filter(database %in% nonempty_dbs) |>
  group_by(database) |>
  summarise(n_label = as.integer(median(n_total)), .groups = "drop")
db_labels <- setNames(
  sprintf("%s (n=%d testable)", db_label_n$database, db_label_n$n_label),
  db_label_n$database
)

count_df$contrast  <- factor(count_df$contrast, levels = DISPLAY_CONTRASTS)
count_df$database  <- factor(count_df$database, levels = intersect(DISPLAY_DBS, nonempty_dbs))
count_df$direction <- factor(count_df$direction, levels = c("Up", "Down"))

n_facets <- length(levels(count_df$database))
PC_H <- max(80, n_facets * 55)

lbl_sz <- scale_text(BASE_COUNT, PC_W)

pD <- ggplot(count_df, aes(x = contrast, y = fraction * 100, fill = direction)) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Cancer_vs_Healthy"], alpha = 0.20,
           color = "grey70", linewidth = 0.2) +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = CONTRAST_COLORS["Training_CR"], alpha = 0.20,
           color = "grey70", linewidth = 0.2) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.3) +
  geom_text(aes(y = fraction * 100 / 2,
                label = ifelse(count > 0, count, "")),
            position = position_dodge(width = 0.7),
            vjust = 0.5, hjust = 0.5, size = lbl_sz,
            color = "white", fontface = "bold", show.legend = FALSE) +
  facet_wrap(~ database, ncol = 1, scales = "free_y", strip.position = "top",
             labeller = as_labeller(db_labels)) +
  scale_x_discrete(labels = CTR_SHORT) +
  scale_fill_manual(values = DIR_COLORS) +
  labs(title = "Fraction of Pathways Enriched Up/Down (CRvH)",
       subtitle = "fGSEA per-database BH corrected fraction",
       x = NULL, y = "% of database significant",
       tag = "D") +
  FIG_THEME +
  theme(axis.text.x     = element_text(angle = 35, hjust = 1,
                                       size = FIG_AXIS_TEXT - 0.5),
        legend.position  = "none",
        strip.text       = element_text(size = FIG_STRIP_SIZE - 1, face = "bold",
                                        margin = margin(1, 0, 1, 0)),
        strip.placement  = "outside")

ggsave(file.path(RPT, "panel_D_fgsea_faceted_SUPP.pdf"), pD,
       width = PC_W, height = PC_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_D_fgsea_faceted_SUPP.png"), pD,
       width = PC_W, height = PC_H, units = "mm", dpi = 300)

cat("F03/CRvH Supp: fGSEA faceted done.\n")
