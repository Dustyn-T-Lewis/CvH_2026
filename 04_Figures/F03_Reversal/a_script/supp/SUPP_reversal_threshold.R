# SUPP Reversal Threshold Sensitivity
# Line plot: % Reversed / Exacerbated / Negligible across |logFC| thresholds
setwd(here::here())
source("04_Figures/shared/style.R")

pacman::p_load(tidyverse)

RPT_PNG <- "04_Figures/F03_Reversal/b_reports/supp/png/panels"
RPT_PDF <- "04_Figures/F03_Reversal/b_reports/supp/pdf/panels"
DAT     <- "04_Figures/F03_Reversal/c_data/panel_supp"

dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT,     recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# -- Load data -----------------------------------------------------------------
source("04_Figures/F03_Reversal/a_script/f03_data.R")
dep <- dep_df %>%
  dplyr::select(gene, logFC_CvH = logFC_CRvH_Baseline,
         logFC_TR = logFC_CR_Training)

# -- Sweep thresholds ---------------------------------------------------------
thresholds <- seq(0.05, 0.30, by = 0.01)
n_total    <- nrow(dep)

sweep_df <- map_dfr(thresholds, function(thr) {
  dep %>%
    mutate(
      opposite  = sign(logFC_CvH) != sign(logFC_TR),
      both_pass = abs(logFC_CvH) >= thr & abs(logFC_TR) >= thr,
      class = case_when(
        opposite & both_pass  ~ "Reversed",
        !opposite & both_pass ~ "Exacerbated",
        TRUE                  ~ "Negligible"
      )
    ) %>%
    count(class) %>%
    mutate(pct = 100 * n / n_total, threshold = thr)
})

sweep_df$class <- factor(sweep_df$class,
                         levels = c("Reversed", "Exacerbated", "Negligible"))

# -- Export CSV ----------------------------------------------------------------
write.csv(sweep_df, file.path(DAT, "SUPP_reversal_threshold.csv"),
          row.names = FALSE)

# -- Plot ----------------------------------------------------------------------
line_cols <- c(Reversed = "#2563EB", Exacerbated = "#DC2626",
               Negligible = "grey50")

pS_threshold <- ggplot(sweep_df, aes(threshold, pct, colour = class)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.2) +
  scale_colour_manual(values = line_cols) +
  scale_x_continuous(breaks = seq(0.05, 0.30, 0.05),
                     labels = sprintf("%.2f", seq(0.05, 0.30, 0.05))) +
  labs(x = "|logFC| threshold", y = "Proteins (%)",
       title = "Reversal classification sensitivity",
       subtitle = "Across dual |logFC| thresholds",
       colour = "Class") +
  FIG_THEME +
  theme(legend.position = "bottom")

ggsave(file.path(RPT_PNG, "SUPP_reversal_threshold.png"), pS_threshold,
       width = 140, height = 100, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(RPT_PDF, "SUPP_reversal_threshold.pdf"), pS_threshold,
       width = 140, height = 100, units = "mm", device = pdf_device)

message("Done: SUPP_reversal_threshold  [", nrow(sweep_df), " rows exported]")
