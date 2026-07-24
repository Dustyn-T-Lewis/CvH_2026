# SUPP: Bootstrap CI for Reversal Pearson r
# 1000 bootstrap replicates of r(logFC_CvH, logFC_TR)
setwd(here::here())
source("04_Figures/shared/style.R")
pacman::p_load(tidyverse)

RPT_PNG <- "04_Figures/F04_Reversal/b_reports/supp/png/panels"
RPT_PDF <- "04_Figures/F04_Reversal/b_reports/supp/pdf/panels"
DAT     <- "04_Figures/F04_Reversal/c_data/panel_supp"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT,     recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

# ── Data ─────────────────────────────────────────────────────────────────────
source("04_Figures/F04_Reversal/a_script/f04_data.R")
dep_df <- dep_df %>%
  transmute(gene,
            logFC_CvH = logFC_CRvH_Baseline,
            logFC_TR  = logFC_CR_Training) %>%
  filter(!is.na(logFC_CvH), !is.na(logFC_TR))

observed_r <- cor(dep_df$logFC_CvH, dep_df$logFC_TR, method = "pearson")
n <- nrow(dep_df)

# ── Bootstrap ────────────────────────────────────────────────────────────────
set.seed(42)
B <- 1000
boot_r <- replicate(B, {
  idx <- sample.int(n, replace = TRUE)
  cor(dep_df$logFC_CvH[idx], dep_df$logFC_TR[idx], method = "pearson")
})

ci <- quantile(boot_r, probs = c(0.025, 0.975))
boot_out <- tibble(replicate = seq_len(B), r_boot = boot_r)
write_csv(boot_out, file.path(DAT, "SUPP_r_bootstrap.csv"))

# ── Plot ─────────────────────────────────────────────────────────────────────
pS_r_boot <- ggplot(boot_out, aes(x = r_boot)) +
  geom_histogram(bins = 40, fill = "grey70", colour = "grey40", linewidth = 0.3) +
  annotate("rect", xmin = ci[1], xmax = ci[2], ymin = -Inf, ymax = Inf,
           fill = DIR_COLORS["Down"], alpha = 0.15) +
  geom_vline(xintercept = observed_r, colour = DIR_COLORS["Up"],
             linewidth = 0.8, linetype = "solid") +
  annotate("text",
           x = observed_r, y = Inf, vjust = -0.5, hjust = -0.1,
           label = sprintf("r = %.3f", observed_r),
           colour = DIR_COLORS["Up"], fontface = "bold", size = 3.5) +
  annotate("text",
           x = mean(ci), y = Inf, vjust = -2,
           label = sprintf("95%% CI [%.3f, %.3f]", ci[1], ci[2]),
           colour = "grey30", size = 3) +
  labs(title    = "Bootstrap Distribution of Reversal r",
       subtitle = sprintf("1000 replicates, n = %d proteins", n),
       x = "Pearson r (logFC Cancer vs logFC Training)",
       y = "Count") +
  FIG_THEME

ggsave(file.path(RPT_PNG, "SUPP_r_bootstrap.png"), pS_r_boot,
       width = 140, height = 100, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_r_bootstrap.pdf"), pS_r_boot,
       width = 140, height = 100, units = "mm", device = pdf_device)

message("Done: SUPP_r_bootstrap")
