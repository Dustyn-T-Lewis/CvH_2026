# Supp C — 6x6 Spearman correlation heatmap of blood markers + B_M_ratio.

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(ggplot2)
})

setwd(rprojroot::find_rstudio_root_file())
source("02-03_Sam's_Results/04_Figures/build_data_index.R")
source("02-03_Sam's_Results/04_Figures/shared/style.R")

sam <- readRDS(sam_idx$sam$dalist_rds)
meta <- as.data.frame(sam$metadata)

vars <- c("HBB_pct", "HBA1_pct", "MB_pct", "ALB_pct", "CKM_pct", "B_M_ratio")
m <- meta[, vars]

cor_mat <- cor(m, method = "spearman", use = "pairwise.complete.obs")
p_mat <- outer(seq_along(vars), seq_along(vars), Vectorize(function(i, j) {
  if (i == j) return(NA_real_)
  suppressWarnings(cor.test(m[[i]], m[[j]], method = "spearman", exact = FALSE)$p.value)
}))
rownames(p_mat) <- colnames(p_mat) <- vars

cor_long <- as.data.frame(as.table(cor_mat)) |>
  rename(var1 = Var1, var2 = Var2, rho = Freq) |>
  mutate(p = as.vector(p_mat))

write_csv(cor_long,
          "02-03_Sam's_Results/04_Figures/F00_blood_contamination/c_data/supp_C_corr_matrix.csv")

cor_long$var1 <- factor(cor_long$var1, levels = vars)
cor_long$var2 <- factor(cor_long$var2, levels = rev(vars))

p_supp_C <- ggplot(cor_long, aes(var1, var2, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman ρ") +
  labs(x = NULL, y = NULL,
       title = "Inter-marker + B_M_ratio correlation (Spearman, N=35)") +
  coord_equal() +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

ggsave("02-03_Sam's_Results/04_Figures/F00_blood_contamination/b_reports/supp/png/panels/supp_C.png",
       p_supp_C, width = 5, height = 4.5, dpi = 300, bg = "white")

supp_C <- p_supp_C
