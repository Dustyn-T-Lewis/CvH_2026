# F02 CRvH — Panel G: Proteome Compartment Composition (Intensity-Weighted)
# Stacked bars for CR_T1, CR_T2, H_T1 showing HPA subcellular fractions.
# Companion dot plot: Plasma protein intensity fraction with bootstrap CIs.
# Outputs: panel_G_compartment_composition_SUPP.png, audit CSVs

setwd(here::here())
source("04_Figures/F02/a_script/style.R")

pacman::p_load(dplyr, tidyr, readr, ggplot2, cowplot)

PG_W <- 200; PG_H <- 140

RPT <- "04_Figures/F02/CRvH/b_reports"
DAT <- "04_Figures/F02/CRvH/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

# --- Load data ---
.dal <- readRDS("02_Normalization/c_data/DAList_normalized.rds")
norm_df <- tibble::as_tibble(cbind(
  .dal$annotation[, c("uniprot_id", "protein", "gene", "description")],
  as.data.frame(.dal$data)))
meta    <- read_csv("00_input/CvH_meta.csv", show_col_types = FALSE)
hpa     <- read.delim("00_input/HPA_annotations.tsv",
                       check.names = FALSE, stringsAsFactors = FALSE)

ann_cols   <- c("uniprot_id", "protein", "gene", "description")
samp_names <- setdiff(names(norm_df), ann_cols)
meta       <- meta |> filter(Col_ID %in% samp_names)

# --- Define 3 groups ---
grp_def <- list(
  "CR_T1" = meta$Col_ID[meta$Group %in% c("CR_CRE", "CR_PLA") & meta$Timepoint == "T1"],
  "CR_T2" = meta$Col_ID[meta$Group %in% c("CR_CRE", "CR_PLA") & meta$Timepoint == "T2"],
  "H_T1"  = meta$Col_ID[meta$Group == "PPS"]
)

# --- HPA annotations ---
hpa_sub <- hpa[, c("Gene", "Protein class", "Subcellular main location")]
names(hpa_sub) <- c("gene", "protein_class", "subcellular_main")
hpa_sub <- hpa_sub[!duplicated(hpa_sub$gene), ]

# Primary compartment: first value before comma
hpa_sub$primary_loc <- trimws(sapply(strsplit(hpa_sub$subcellular_main, ","), `[`, 1))
hpa_sub$primary_loc[is.na(hpa_sub$primary_loc) | hpa_sub$primary_loc == ""] <- "Unknown"

# Consolidate to 8 categories
hpa_sub$is_plasma   <- grepl("Plasma proteins", hpa_sub$protein_class, fixed = TRUE)
hpa_sub$is_secreted <- grepl("Predicted secreted", hpa_sub$protein_class, fixed = TRUE)

consolidate_compartment <- function(loc, is_plasma, is_secreted) {
  case_when(
    loc == "Cytosol"                  ~ "Cytosol",
    loc == "Mitochondria"             ~ "Mitochondria",
    loc %in% c("Nucleoplasm", "Nuclear membrane", "Nuclear speckles",
               "Nuclear bodies", "Nucleoli", "Nucleoli fibrillar center",
               "Nucleoli rim")        ~ "Nucleus",
    loc == "Plasma membrane"          ~ "Plasma membrane",
    loc %in% c("Vesicles", "Endosomes", "Lysosomes") ~ "Vesicles",
    loc == "Endoplasmic reticulum"    ~ "ER",
    loc %in% c("Microtubules", "Actin filaments", "Intermediate filaments",
               "Focal adhesion sites", "Centrosome", "Centriolar satellite",
               "Midbody")             ~ "Cytoskeleton",
    is_plasma | is_secreted           ~ "Secreted/Plasma",
    TRUE                              ~ "Unknown"
  )
}

hpa_sub$compartment <- consolidate_compartment(
  hpa_sub$primary_loc, hpa_sub$is_plasma, hpa_sub$is_secreted
)

COMPARTMENT_ORDER <- c("Secreted/Plasma", "Cytosol", "Mitochondria", "Nucleus",
                        "Plasma membrane", "Vesicles", "ER", "Cytoskeleton", "Unknown")

COMPARTMENT_COLORS <- c(
  "Secreted/Plasma" = "#C62828",
  "Cytosol"         = "#66C2A5",
  "Mitochondria"    = "#FC8D62",
  "Nucleus"         = "#8DA0CB",
  "Plasma membrane" = "#E78AC3",
  "Vesicles"        = "#A6D854",
  "ER"              = "#FFD92F",
  "Cytoskeleton"    = "#B3B3B3",
  "Unknown"         = "#D9D9D9"
)

# --- Join annotation to proteins ---
prot_df <- norm_df |>
  dplyr::select(gene) |>
  left_join(hpa_sub |> dplyr::select(gene, compartment, is_plasma) |> distinct(gene, .keep_all = TRUE),
            by = "gene") |>
  mutate(
    compartment = replace_na(compartment, "Other/Unknown"),
    is_plasma   = replace_na(is_plasma, FALSE),
    compartment = factor(compartment, levels = COMPARTMENT_ORDER)
  )

# --- Intensity-weighted compartment fractions ---
lin_mat <- 2^as.matrix(norm_df[, samp_names])

comp_fractions <- lapply(names(grp_def), function(grp_name) {
  ids <- grp_def[[grp_name]]
  mean_int <- rowMeans(lin_mat[, ids, drop = FALSE], na.rm = TRUE)
  total <- sum(mean_int, na.rm = TRUE)

  tibble(compartment = prot_df$compartment, mean_intensity = mean_int) |>
    group_by(compartment) |>
    summarise(
      intensity_sum = sum(mean_intensity, na.rm = TRUE),
      n_proteins    = n(),
      .groups = "drop"
    ) |>
    mutate(fraction = intensity_sum / total, group = grp_name)
}) |> bind_rows()

comp_fractions$group <- factor(comp_fractions$group, levels = c("CR_T1", "CR_T2", "H_T1"))

# --- Plasma protein fraction with bootstrap CI ---
set.seed(42)
N_BOOT <- 1000

plasma_stats <- lapply(names(grp_def), function(grp_name) {
  ids <- grp_def[[grp_name]]
  sub_mat <- lin_mat[, ids, drop = FALSE]

  mean_int <- rowMeans(sub_mat, na.rm = TRUE)
  total    <- sum(mean_int, na.rm = TRUE)
  plasma_sum <- sum(mean_int[prot_df$is_plasma], na.rm = TRUE)
  obs_frac <- plasma_sum / total

  boot_fracs <- replicate(N_BOOT, {
    boot_ids <- sample(ids, replace = TRUE)
    boot_mean <- rowMeans(lin_mat[, boot_ids, drop = FALSE], na.rm = TRUE)
    sum(boot_mean[prot_df$is_plasma], na.rm = TRUE) / sum(boot_mean, na.rm = TRUE)
  })

  tibble(
    group        = grp_name,
    plasma_frac  = obs_frac,
    ci_lo        = quantile(boot_fracs, 0.025),
    ci_hi        = quantile(boot_fracs, 0.975),
    n_plasma     = sum(prot_df$is_plasma),
    n_total      = nrow(prot_df),
    n_samples    = length(ids)
  )
}) |> bind_rows()

plasma_stats$group <- factor(plasma_stats$group, levels = c("CR_T1", "CR_T2", "H_T1"))

# --- Stacked bar plot ---
pct_labels <- comp_fractions |>
  mutate(pct = fraction * 100) |>
  filter(pct >= 3)

GROUP_LABELS <- c(
  CR_T1 = "CR Baseline\n(microneedle)",
  CR_T2 = "CR Post-Trn\n(microneedle)",
  H_T1  = "Healthy\n(Bergstrom)"
)

p_stack <- ggplot(comp_fractions, aes(x = group, y = fraction, fill = compartment)) +
  geom_col(position = "fill", width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(data = pct_labels,
            aes(label = sprintf("%.0f%%", pct)),
            position = position_fill(vjust = 0.5),
            size = scale_text(BASE_STAT, PG_W) * 0.75,
            color = "grey20", fontface = "bold") +
  scale_fill_manual(values = COMPARTMENT_COLORS, name = "Compartment",
                    breaks = COMPARTMENT_ORDER) +
  scale_x_discrete(labels = GROUP_LABELS) +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(
    title = "Proteome Compartment Composition",
    subtitle = sprintf("%s proteins | Intensity-weighted fractions (mean linear)",
                       format(nrow(norm_df), big.mark = ",")),
    y = "Fraction of total intensity",
    x = NULL, tag = "G"
  ) +
  FIG_THEME +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "right",
    legend.key.size = unit(4, "mm")
  )

# --- Plasma fraction dot plot ---
BAR_COLORS <- c(CR_T1 = "#2166AC", CR_T2 = "#67A9CF", H_T1 = "#4DAF4A")

p_plasma <- ggplot(plasma_stats, aes(x = plasma_frac, y = group, color = group)) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.25, linewidth = 0.8) +
  geom_point(size = 3) +
  scale_color_manual(values = BAR_COLORS, guide = "none") +
  scale_y_discrete(labels = GROUP_LABELS, limits = rev(levels(plasma_stats$group))) +
  scale_x_continuous(labels = scales::percent, limits = c(0, NA), expand = expansion(mult = c(0, 0.1))) +
  labs(
    title = "Plasma Protein\nIntensity Fraction",
    subtitle = sprintf("n=%d flagged | Bootstrap 95%% CI", sum(prot_df$is_plasma)),
    x = "Fraction of total intensity",
    y = NULL
  ) +
  FIG_THEME +
  theme(
    panel.grid.major.y = element_blank(),
    plot.title = element_text(size = FIG_TITLE_SIZE * 0.85)
  )

# --- Combine ---
pG <- cowplot::plot_grid(p_stack, p_plasma, nrow = 1,
                          rel_widths = c(3, 1.8), align = "h", axis = "tb")

ggsave(file.path(RPT, "panel_G_compartment_composition_SUPP.png"), pG,
       width = PG_W, height = PG_H, units = "mm", dpi = 300)

# --- Audit exports ---
write_csv(comp_fractions |>
  mutate(fraction = round(fraction, 6), intensity_sum = round(intensity_sum, 2)),
  file.path(DAT, "audit_panel_G_compartment_fractions.csv"))

write_csv(plasma_stats |>
  mutate(across(c(plasma_frac, ci_lo, ci_hi), ~ round(.x, 6))),
  file.path(DAT, "audit_panel_G_plasma_fraction.csv"))

message("F02/CRvH Panel G done")
