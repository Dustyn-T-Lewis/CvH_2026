# F04 CRvH Supplementary: RRHO with Marginal Signal Profiles (Concordance)
# Non-stratified RRHO (Plaisier 2010) for smooth continuous visualization.
# Signed -log10(p): positive = over-enrichment (concordance), negative = under-enrichment.
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F04/a_script/style.R")

library(tidyverse)
library(patchwork)

RPT <- "04_Figures/F04/CRvH/b_reports/supp"
DAT <- "04_Figures/F04/CRvH/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

# --- Data prep ----------------------------------------------------------------
dep_df <- read_csv("03_DEP/c_data/03_combined_results_CRvH.csv", show_col_types = FALSE)

rr_df <- dep_df %>%
  transmute(gene, t_cvh = t_Cancer_vs_Healthy, t_tr = t_Training_CR) %>%
  filter(!is.na(t_cvh), !is.na(t_tr)) %>%
  distinct(gene, .keep_all = TRUE)

N <- nrow(rr_df)
scores1 <- rr_df$t_cvh
scores2 <- rr_df$t_tr

# --- Non-stratified RRHO (Plaisier 2010, pure-R phyper) -----------------------
stepsize <- 20
rank1 <- rank(-scores1, ties.method = "first")   # 1 = most upregulated
rank2 <- rank(-scores2, ties.method = "first")

# 2D prefix sum: count_mat[i,j] = |{genes with rank1 <= i AND rank2 <= j}|
count_mat <- matrix(0L, N, N)
for (g in seq_len(N)) count_mat[rank1[g], rank2[g]] <- 1L
count_mat <- apply(count_mat, 2, cumsum)
count_mat <- t(apply(count_mat, 1, cumsum))

steps <- seq(stepsize, N, by = stepsize)
n_steps <- length(steps)

hmat <- matrix(0, n_steps, n_steps)
for (i in seq_along(steps)) {
  s <- steps[i]
  for (j in seq_along(steps)) {
    M <- steps[j]
    k <- count_mat[s, M]
    expected <- s * M / N
    if (k > expected) {
      hmat[i, j] <- -log10(pmax(phyper(k - 1, M, N - M, s, lower.tail = FALSE), 1e-320))
    } else if (k < expected) {
      hmat[i, j] <- log10(pmax(phyper(k, M, N - M, s, lower.tail = TRUE), 1e-320))
    }
  }
}
nr <- nrow(hmat); nc <- ncol(hmat)
message(sprintf("  RRHO matrix: %d x %d, range [%.1f, %.1f]", nr, nc, min(hmat), max(hmat)))

# Zero-crossing positions (where t-statistic flips sign)
zero_idx1 <- which.min(abs(steps - sum(scores1 > 0)))
zero_idx2 <- which.min(abs(steps - sum(scores2 > 0)))

# Reference gene counts from RRHO2 panel_E
rrho2_file <- file.path(DAT, "panel_E/rrho2_summary.csv")
if (file.exists(rrho2_file)) {
  rrho2_sum <- read_csv(rrho2_file, show_col_types = FALSE)
  n_UU <- rrho2_sum$n_hotspot_genes[rrho2_sum$quadrant == "Concordant_Up"]
  n_DD <- rrho2_sum$n_hotspot_genes[rrho2_sum$quadrant == "Concordant_Down"]
} else {
  n_UU <- NA
  n_DD <- NA
  message("RRHO2 summary not found -- run panel_E.R first for quadrant counts")
}

# Asymmetry p-value
asym_file <- file.path(DAT, "supp/h_directional_asymmetry_tests.csv")
if (file.exists(asym_file)) {
  asym_df <- read_csv(asym_file, show_col_types = FALSE)
  perm_p <- asym_df$p_value[asym_df$test == "Permutation (10K)"]
} else {
  perm_p <- NA
  message("Asymmetry tests not found -- run directional_asymmetry.R first")
}

# --- Marginal profiles (max over-enrichment along each axis) ------------------
top_max   <- apply(hmat, 1, max)
right_max <- apply(hmat, 2, max)

COL_LEFT  <- "#E57373"   # Concordant Up
COL_RIGHT <- "#64B5F6"   # Concordant Down

# --- Heatmap (jet colormap, matching RRHO2 panels) ----------------------------
JET_COLORS <- c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F",
                "yellow", "#FF7F00", "red", "#7F0000")

hmat_df <- expand.grid(row = 1:nr, col = 1:nc) %>%
  mutate(value = pmax(as.vector(hmat), 0))   # clamp negligible negatives to 0

max_val <- max(hmat_df$value)

uu_label <- if (!is.na(n_UU)) sprintf("Conc. Up\nn = %d", n_UU) else "Conc. Up"
dd_label <- if (!is.na(n_DD)) sprintf("Conc. Down\nn = %d", n_DD) else "Conc. Down"

p_heat <- ggplot(hmat_df, aes(x = row, y = col, fill = value)) +
  geom_raster() +
  scale_fill_gradientn(
    colors = JET_COLORS, limits = c(0, max_val), na.value = "#00007F",
    name = expression(-log[10](P)),
    guide = guide_colorbar(barwidth = unit(30, "mm"), barheight = unit(2.5, "mm"),
                           title.position = "top", title.hjust = 0.5,
                           title.theme = element_text(size = 6, face = "bold"))
  ) +
  geom_vline(xintercept = zero_idx1, linetype = "dashed", color = "grey70",
             linewidth = 0.3) +
  geom_hline(yintercept = zero_idx2, linetype = "dashed", color = "grey70",
             linewidth = 0.3) +
  annotate("text", x = zero_idx1 * 0.4, y = zero_idx2 * 0.15,
           label = uu_label,
           color = "white", fontface = "bold", size = 2.8, lineheight = 0.85) +
  annotate("text", x = zero_idx1 + (nr - zero_idx1) * 0.6,
           y = zero_idx2 + (nc - zero_idx2) * 0.88,
           label = dd_label,
           color = "white", fontface = "bold", size = 2.8, lineheight = 0.85) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = expression("Cancer vs Healthy rank" ~ (Up %->% Down)),
       y = expression("Training CR rank" ~ (Up %->% Down))) +
  FIG_THEME +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.border = element_blank(), panel.grid = element_blank(),
        legend.position = "bottom", legend.margin = margin(0, 0, 0, 0),
        plot.margin = margin(1, 1, 2, 2, "mm"))

# --- Top marginal -------------------------------------------------------------
top_df <- tibble(idx = seq_along(top_max), value = top_max)

asym_str <- if (!is.na(perm_p)) {
  sprintf(" | Asymmetry: n = %s/%s, perm %s",
          ifelse(is.na(n_UU), "?", n_UU),
          ifelse(is.na(n_DD), "?", n_DD),
          fmt_p(perm_p))
} else ""

p_top <- ggplot() +
  geom_ribbon(data = filter(top_df, idx <= zero_idx1),
              aes(x = idx, ymin = 0, ymax = value), fill = COL_LEFT, alpha = 0.7) +
  geom_ribbon(data = filter(top_df, idx >= zero_idx1),
              aes(x = idx, ymin = 0, ymax = value), fill = COL_RIGHT, alpha = 0.7) +
  geom_line(data = top_df, aes(x = idx, y = value), linewidth = 0.3, color = "grey30") +
  geom_vline(xintercept = zero_idx1, linetype = "dashed", color = "grey30",
             linewidth = 0.3, alpha = 0.7) +
  scale_x_continuous(limits = c(1, nr), expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(title = "RRHO with Marginal Signal Profiles (Plaisier 2010)",
       subtitle = sprintf("Concordance | %d genes%s", N, asym_str),
       y = expression(max ~ -log[10](P))) +
  FIG_THEME +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(), panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(2, 1, 0, 2, "mm"))

# --- Right marginal -----------------------------------------------------------
right_df <- tibble(idx = seq_along(right_max), value = right_max)

p_right <- ggplot() +
  geom_ribbon(data = filter(right_df, idx <= zero_idx2),
              aes(y = idx, xmin = 0, xmax = value), fill = COL_LEFT, alpha = 0.7) +
  geom_ribbon(data = filter(right_df, idx >= zero_idx2),
              aes(y = idx, xmin = 0, xmax = value), fill = COL_RIGHT, alpha = 0.7) +
  geom_line(data = right_df, aes(y = idx, x = value), linewidth = 0.3, color = "grey30") +
  geom_hline(yintercept = zero_idx2, linetype = "dashed", color = "grey30",
             linewidth = 0.3, alpha = 0.7) +
  scale_y_continuous(limits = c(1, nc), expand = c(0, 0)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = expression(max ~ -log[10](P))) +
  FIG_THEME +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.title.y = element_blank(), panel.grid = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(1, 2, 2, 0, "mm"))

# --- Composite layout ---------------------------------------------------------
composite <- p_top + plot_spacer() + p_heat + p_right +
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))

ggsave(file.path(RPT, "supp_rrho_marginals_SUPP.pdf"), composite,
       width = 240, height = 260, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "supp_rrho_marginals_SUPP.png"), composite,
       width = 240, height = 260, units = "mm", dpi = 300)

message("F04 CRvH RRHO-M supplementary saved (Plaisier 2010)")
message(sprintf("  Matrix: %dx%d | Range: [%.1f, %.1f]", nr, nc, min(hmat), max(hmat)))
