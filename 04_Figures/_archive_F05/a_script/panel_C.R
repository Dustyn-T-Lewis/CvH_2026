# Figure 5 — Panel C: GSEA Concordance Matrix (6x6 NES Correlation)
# De novo fGSEA on Hallmark + GO Slim for all 6 contrasts.
# Spearman correlation of NES vectors for each pair of contrasts.
# Lower triangle: colored tiles, upper triangle: correlation values.
# Outputs: panel_C_concordance_matrix.pdf/png, concordance_matrix.csv

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})

PC_W <- 200; PC_H <- 200

RPT <- "04_Figures/F05/b_reports"
DAT <- "04_Figures/F05/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# --- Build Hallmark + GO Slim pathway collection ---
pw_collection <- build_hallmark_goslim_collection(min_size = 10, max_size = 500)

# --- Load DEP results (t-statistics for ranking) ---
dep_crvh <- read_csv("03_DEP/c_data/03_combined_results_CRvH.csv",
                      show_col_types = FALSE)
dep_cr   <- read_csv("03_DEP/c_data/03_combined_results_CR.csv",
                      show_col_types = FALSE)

# --- Run fGSEA per contrast on Hallmark + GO Slim ---
set.seed(42)
fgsea_all_list <- list()

for (ctr in ALL_CONTRASTS) {
  tcol    <- paste0("t_", ctr)
  dep_src <- if (ctr %in% ALL_CONTRASTS_CRVH) dep_crvh else dep_cr
  stats   <- setNames(dep_src[[tcol]], dep_src$gene)
  stats   <- stats[!is.na(stats) & is.finite(stats)]

  if (anyDuplicated(names(stats))) {
    dup_df <- tibble(gene = names(stats), t = unname(stats)) |>
      group_by(gene) |>
      slice_max(abs(t), n = 1, with_ties = FALSE) |>
      ungroup()
    stats <- setNames(dup_df$t, dup_df$gene)
  }
  stats <- sort(stats, decreasing = TRUE)

  res <- run_fgsea_perdb(
    ranks    = stats,
    pathways = pw_collection,
    nperm    = 10000,
    min_size = 10,
    max_size = 500
  )
  res$contrast <- ctr
  fgsea_all_list[[ctr]] <- res
}

fgsea_all <- bind_rows(fgsea_all_list) |>
  filter(contrast %in% ALL_CONTRASTS) |>
  select(pathway, contrast, NES)

cat(sprintf("Hallmark + GO Slim concordance: %d rows, %d unique pathways\n",
            nrow(fgsea_all), n_distinct(fgsea_all$pathway)))

# --- Pivot to wide: rows = pathways, columns = contrasts ---
nes_wide <- fgsea_all |>
  pivot_wider(id_cols = pathway, names_from = contrast, values_from = NES)

cat(sprintf("Total pathways in union: %d\n", nrow(nes_wide)))

# --- Compute pairwise Spearman correlations ---
n_ctr <- length(ALL_CONTRASTS)
cor_mat <- matrix(NA_real_, nrow = n_ctr, ncol = n_ctr,
                  dimnames = list(ALL_CONTRASTS, ALL_CONTRASTS))
n_mat   <- matrix(NA_integer_, nrow = n_ctr, ncol = n_ctr,
                  dimnames = list(ALL_CONTRASTS, ALL_CONTRASTS))

for (i in seq_len(n_ctr)) {
  for (j in seq_len(n_ctr)) {
    xi <- nes_wide[[ALL_CONTRASTS[i]]]
    xj <- nes_wide[[ALL_CONTRASTS[j]]]
    complete <- !is.na(xi) & !is.na(xj)
    n_mat[i, j] <- sum(complete)
    if (sum(complete) >= 5) {
      cor_mat[i, j] <- cor(xi[complete], xj[complete], method = "spearman")
    }
  }
}

# --- Build plot data frame ---
plot_rows <- list()

for (i in seq_len(n_ctr)) {
  for (j in seq_len(n_ctr)) {
    row_data <- data.frame(
      x = j, y = n_ctr - i + 1,
      contrast_x = ALL_CONTRASTS[j],
      contrast_y = ALL_CONTRASTS[i],
      rho = cor_mat[i, j],
      n_shared = n_mat[i, j],
      stringsAsFactors = FALSE
    )

    if (i == j) {
      row_data$tile_fill <- NA_real_
      row_data$text_label <- CTR_SHORT[ALL_CONTRASTS[i]]
      row_data$text_size  <- 3.2
      row_data$text_face  <- "bold"
      row_data$region     <- "diag"
    } else if (i > j) {
      row_data$tile_fill  <- cor_mat[i, j]
      row_data$text_label <- NA_character_
      row_data$text_size  <- NA_real_
      row_data$text_face  <- NA_character_
      row_data$region     <- "lower"
    } else {
      row_data$tile_fill  <- NA_real_
      row_data$text_label <- sprintf("%.2f\n(%d)", cor_mat[i, j], n_mat[i, j])
      row_data$text_size  <- 2.8
      row_data$text_face  <- "plain"
      row_data$region     <- "upper"
    }

    plot_rows[[length(plot_rows) + 1]] <- row_data
  }
}

plot_df <- do.call(rbind, plot_rows)

lower_df <- plot_df[plot_df$region == "lower", ]
upper_df <- plot_df[plot_df$region == "upper", ]
diag_df  <- plot_df[plot_df$region == "diag", ]

# --- Plot ---
pC <- ggplot() +
  geom_tile(data = lower_df,
            aes(x = x, y = y, fill = tile_fill),
            color = "white", linewidth = 0.8) +
  geom_text(data = upper_df,
            aes(x = x, y = y, label = text_label),
            size = 2.8, color = "grey20", lineheight = 0.85) +
  geom_tile(data = diag_df,
            aes(x = x, y = y),
            fill = "grey95", color = "white", linewidth = 0.8) +
  geom_text(data = diag_df,
            aes(x = x, y = y, label = text_label),
            size = 3.2, fontface = "bold", color = "grey20") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D6604D",
                       midpoint = 0, limits = c(-1, 1),
                       name = "Spearman\nrho", na.value = "grey95") +
  coord_fixed() +
  labs(title = "GSEA Concordance Across Contrasts",
       subtitle = "Hallmark + GO Slim | Pairwise Spearman correlation of NES vectors",
       tag = "C") +
  FIG_THEME +
  theme(
    axis.text  = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    legend.position = "right"
  )

# --- Save ---
ggsave(file.path(RPT, "panel_C_concordance_matrix.pdf"), pC,
       width = PC_W, height = PC_H, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_C_concordance_matrix.png"), pC,
       width = PC_W, height = PC_H, units = "mm", dpi = 300)

# --- Export ---
cor_export <- as.data.frame(cor_mat)
cor_export$contrast <- rownames(cor_export)
cor_export <- cor_export[, c("contrast", ALL_CONTRASTS)]

write.csv(cor_export,
          file.path(DAT, "concordance_matrix.csv"), row.names = FALSE)

cat(sprintf("Panel C (concordance matrix, Hallmark + GO Slim) done: %d x %d contrasts.\n",
            n_ctr, n_ctr))
