# F06 CR Panel B: RRHO2 — Supplement Reversal
# Baseline_Supplement (x) vs Supplement_Interaction (y)
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F06/a_script/style.R")
suppressPackageStartupMessages({ library(tidyverse) })

PG_W <- 180
RPT <- "04_Figures/F06/CR/b_reports"
DAT <- "04_Figures/F06/CR/c_data"
dir.create(file.path(DAT, "panel_B"), recursive = TRUE, showWarnings = FALSE)
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

CTR_X <- "Baseline_Supplement"; CTR_Y <- "Supplement_Interaction"
dep_df <- read_csv("03_DEP/c_data/03_combined_results_CR.csv", show_col_types = FALSE)

rank_df <- dep_df %>%
  filter(!is.na(.data[[paste0("t_", CTR_X)]]), !is.na(.data[[paste0("t_", CTR_Y)]])) %>%
  mutate(rank_X = rank(-.data[[paste0("t_", CTR_X)]], ties.method = "average"),
         rank_Y = rank(-.data[[paste0("t_", CTR_Y)]], ties.method = "average"))

# Pure-R RRHO2 (phyper-based)
n <- nrow(rank_df)
step <- max(1L, floor(n / 200))
breaks <- seq(step, n, by = step)
nb <- length(breaks)
mat <- matrix(0, nrow = nb, ncol = nb)

gene_x <- rank_df$gene[order(rank_df$rank_X)]
gene_y <- rank_df$gene[order(rank_df$rank_Y)]

for (i in seq_len(nb)) {
  set_x <- gene_x[seq_len(breaks[i])]
  for (j in seq_len(nb)) {
    set_y <- gene_y[seq_len(breaks[j])]
    overlap <- length(intersect(set_x, set_y))
    pval <- phyper(overlap - 1, breaks[i], n - breaks[i], breaks[j], lower.tail = FALSE)
    mat[i, j] <- -log10(max(pval, .Machine$double.xmin))
  }
}

# Sign by concordance
mid <- ceiling(nb / 2)
sign_mat <- matrix(1, nrow = nb, ncol = nb)
sign_mat[seq_len(mid), (mid+1):nb] <- -1
sign_mat[(mid+1):nb, seq_len(mid)] <- -1
mat <- mat * sign_mat

grid_df <- expand.grid(i = seq_len(nb), j = seq_len(nb)) %>%
  mutate(x = breaks[i], y = breaks[j], signed_score = mat[cbind(i, j)])
max_abs <- max(abs(grid_df$signed_score))

r_val <- cor(rank_df[[paste0("t_", CTR_X)]], rank_df[[paste0("t_", CTR_Y)]], method = "pearson")
rho_val <- cor(rank_df[[paste0("t_", CTR_X)]], rank_df[[paste0("t_", CTR_Y)]], method = "spearman")
r_ci <- fisher_z_ci(r_val, n)

# For reversal: concordant quadrants = exacerbated, discordant = reversed
pB <- ggplot(grid_df, aes(x = x, y = y, fill = signed_score)) +
  geom_tile() +
  scale_fill_gradient2(low = "#4393C3", mid = "white", high = "#D6604D",
                        midpoint = 0, limits = c(-max_abs, max_abs),
                        name = expression(-log[10](p) %*% sign)) +
  geom_hline(yintercept = n/2, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  geom_vline(xintercept = n/2, linetype = "dashed", color = "grey40", linewidth = 0.3) +
  annotate("text", x = nb*step*0.25, y = nb*step*0.75, label = "Exacerbated\n(same sign)",
           fontface = "bold", size = 3.5, color = "#D6604D") +
  annotate("text", x = nb*step*0.75, y = nb*step*0.25, label = "Exacerbated\n(same sign)",
           fontface = "bold", size = 3.5, color = "#D6604D") +
  annotate("text", x = nb*step*0.75, y = nb*step*0.75, label = "Reversed",
           fontface = "bold", size = 3.5, color = "#4393C3") +
  annotate("text", x = nb*step*0.25, y = nb*step*0.25, label = "Reversed",
           fontface = "bold", size = 3.5, color = "#4393C3") +
  coord_fixed() +
  labs(title = "RRHO2: Supplement Reversal",
       subtitle = sprintf("N = %d | r = %.3f [%.3f, %.3f] | rho = %.3f",
                            n, r_val, r_ci[1], r_ci[2], rho_val),
       x = "Rank by t(Baseline Supplement) -> most up",
       y = "Rank by t(Supplement Interaction) -> most up") +
  FIG_THEME +
  theme(legend.position = "right", legend.key.width = unit(3, "mm"),
        legend.key.height = unit(15, "mm"), axis.text = element_text(size = 8))

ggsave(file.path(RPT, "panel_B_RRHO2_MAIN.pdf"), pB,
       width = PG_W, height = PG_W, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_B_RRHO2_MAIN.png"), pB,
       width = PG_W, height = PG_W, units = "mm", dpi = 300)
message("F06 CR Panel B done")
