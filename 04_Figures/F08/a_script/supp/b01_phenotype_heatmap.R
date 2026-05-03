# Supplementary: Module-Phenotype Correlation Heatmap (CvH)
# Correlates module eigengenes with available sample-level variables
# CvH design: Cancer (CR) vs Healthy, with CRE/PLA supplement arms
# Since CvH lacks continuous phenotypes (BMI, VL, LBM), this panel
# shows eigengene correlations with:
#   - cancer_binary (1 = CR, 0 = H)
#   - time_binary (1 = T2, 0 = T1; CR only)
#   - supplement_binary (1 = CRE, 0 = PLA; CR only)
# Plus group-mean eigengene bar/box plots per module
#
# Generates: b01_phenotype_heatmap_SUPP.pdf/.png, c_data/b01_phenotype_*.csv

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F08/a_script/style.R")

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(patchwork)
})

RPT <- "04_Figures/F08/b_reports/supp"
DAT <- "04_Figures/F08/c_data"
SRC <- "04_Figures/F06/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

message("Supplementary: module-phenotype heatmap (CvH)...")

# --- Load data ---
MEs  <- readRDS(file.path(SRC, "MEs.rds"))
meta <- read_csv(file.path(SRC, "meta.csv"), show_col_types = FALSE)
mod_bio <- read_csv(file.path(SRC, "mod_bio_labels.csv"), show_col_types = FALSE)

if (!"display_label" %in% colnames(mod_bio)) {
  mod_bio <- mod_bio %>%
    mutate(display_label = paste0(module_id, ": ", module_color, " (n=", n_proteins, ")"))
}
mod_display_vec <- setNames(mod_bio$display_label, paste0("ME", mod_bio$module_color))

# Module ordering: by size
mod_order <- mod_bio %>%
  arrange(desc(n_proteins)) %>%
  pull(module_color) %>%
  paste0("ME", .)
mod_order <- intersect(mod_order, colnames(MEs))

# --- Build phenotype matrix ---
# Cancer binary (all samples)
meta$cancer_binary <- ifelse(meta$group %in% c("CR_CRE", "CR_PLA"), 1, 0)

# Time binary (CR samples only: T1 = 0, T2 = 1)
meta$time_binary <- ifelse(meta$time == "T2", 1,
                    ifelse(meta$time == "T1", 0, NA_integer_))

# Supplement binary (CR samples only: CRE = 1, PLA = 0)
meta$supplement_binary <- ifelse(meta$supplement == "CRE", 1,
                          ifelse(meta$supplement == "PLA", 0, NA_integer_))

trait_cols <- c("cancer_binary", "time_binary", "supplement_binary")
trait_labels <- c(cancer_binary = "Cancer Status",
                  time_binary = "Timepoint (T2 vs T1)",
                  supplement_binary = "Supplement (CRE vs PLA)")

# --- Compute correlations: eigengene vs trait ---
cor_list <- list()
for (mod in mod_order) {
  me_vec <- MEs[meta$sample_id, mod]
  for (trait in trait_cols) {
    tv <- meta[[trait]]
    ok <- complete.cases(me_vec, tv)
    n_pair <- sum(ok)
    if (n_pair < 4) {
      cor_list <- c(cor_list, list(tibble(
        module = mod, trait = trait, r = NA_real_,
        p_raw = NA_real_, n = n_pair,
        ci_lo = NA_real_, ci_hi = NA_real_
      )))
      next
    }
    ct <- cor.test(me_vec[ok], tv[ok], method = "pearson")
    cor_list <- c(cor_list, list(tibble(
      module = mod, trait = trait,
      r = round(ct$estimate, 4),
      p_raw = ct$p.value, n = n_pair,
      ci_lo = round(ct$conf.int[1], 4),
      ci_hi = round(ct$conf.int[2], 4)
    )))
  }
}

cor_df <- bind_rows(cor_list)

# BH correction per trait
cor_df <- cor_df %>%
  group_by(trait) %>%
  mutate(p_bh = p.adjust(p_raw, method = "BH")) %>%
  ungroup()

# --- Build heatmap ---
cor_df <- cor_df %>%
  mutate(
    stars = sig_stars(p_bh),
    label = ifelse(is.na(r), "", sprintf("%.2f%s", r,
                   ifelse(stars == "ns", "", paste0("\n", stars)))),
    trait_label = trait_labels[trait],
    module = factor(module, levels = rev(mod_order))
  )

PH_W <- 180
PH_H <- 200

p_heat <- ggplot(cor_df, aes(x = trait_label, y = module, fill = r)) +
  geom_tile(color = "black", linewidth = 0.3) +
  # FDR border
  geom_tile(data = cor_df %>% filter(!is.na(p_bh) & p_bh < 0.05),
            color = "black", linewidth = 1.0, fill = NA) +
  # Nominal dot
  geom_point(data = cor_df %>% filter(!is.na(p_raw) & p_raw < 0.05 &
                                        (is.na(p_bh) | p_bh >= 0.05)),
             shape = 16, size = 0.8, color = "grey30") +
  # Text
  geom_text(data = cor_df %>% filter(is.na(p_bh) | p_bh >= 0.05),
            aes(label = label), size = 2.5, color = "black", lineheight = 0.85) +
  geom_text(data = cor_df %>% filter(!is.na(p_bh) & p_bh < 0.05),
            aes(label = label), size = 2.8, fontface = "bold",
            color = "white", lineheight = 0.85) +
  scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                       midpoint = 0, limits = c(-0.8, 0.8),
                       oob = scales::squish,
                       name = "Pearson r",
                       na.value = "grey90") +
  scale_y_discrete(labels = function(x) {
    lbl <- mod_display_vec[x]
    ifelse(is.na(lbl), x, lbl)
  }) +
  labs(x = NULL, y = NULL,
       title = "Module-Phenotype Correlations (CvH)",
       subtitle = "Pearson r | BH-corrected per trait | Solid = FDR < 0.05 | Dot = nominal p < 0.05",
       caption = paste0("Cancer Status: all ", nrow(meta), " samples | ",
                        "Time/Supplement: CR samples only (n=",
                        sum(!is.na(meta$time_binary)), ")")) +
  FIG_THEME +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, size = 9, face = "bold"),
    axis.text.y = element_text(size = 7.5),
    legend.position = "bottom",
    plot.caption = element_text(size = 7, color = "grey40", hjust = 0)
  )

# --- Eigengene boxplots per group_time (top 4 key modules) ---
km_file <- file.path(SRC, "key_modules.txt")
KEY_MODULES <- if (file.exists(km_file)) {
  readLines(km_file)
} else {
  head(mod_bio$module_color, 4)
}
KEY_MODULES <- KEY_MODULES[nzchar(trimws(KEY_MODULES))]

me_long <- tibble()
for (mod in paste0("ME", KEY_MODULES)) {
  me_long <- bind_rows(me_long, tibble(
    sample_id = meta$sample_id,
    group_time = meta$group_time,
    module = mod,
    eigengene = MEs[meta$sample_id, mod]
  ))
}

me_long$group_time <- factor(me_long$group_time,
                              levels = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "H_T1"))
me_long$module <- factor(me_long$module, levels = paste0("ME", KEY_MODULES))

p_box <- ggplot(me_long, aes(x = group_time, y = eigengene, fill = group_time)) +
  geom_boxplot(width = 0.6, outlier.size = 0.8, linewidth = 0.3) +
  facet_wrap(~ module, scales = "free_y", ncol = 2,
             labeller = labeller(module = function(x) {
               lbl <- mod_display_vec[x]
               ifelse(is.na(lbl), x, lbl)
             })) +
  scale_fill_manual(values = GROUP_COLORS, guide = "none") +
  labs(x = NULL, y = "Module Eigengene",
       title = "Eigengene Distribution by Group") +
  FIG_THEME +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    strip.text  = element_text(size = 8, face = "bold")
  )

# --- Composite ---
composite <- p_heat / p_box + plot_layout(heights = c(1.2, 1))

ggsave(file.path(RPT, "b01_phenotype_heatmap_SUPP.pdf"), composite,
       width = PH_W, height = PH_H * 1.6, units = "mm",
       device = pdf_device, limitsize = FALSE)
ggsave(file.path(RPT, "b01_phenotype_heatmap_SUPP.png"), composite,
       width = PH_W, height = PH_H * 1.6, units = "mm",
       dpi = 300, limitsize = FALSE)

# --- Save audit ---
write_csv(cor_df %>% dplyr::select(module, trait, r, p_raw, p_bh, n, ci_lo, ci_hi),
          file.path(DAT, "b01_phenotype_correlation_audit.csv"))

message("  Supplementary phenotype heatmap saved")
