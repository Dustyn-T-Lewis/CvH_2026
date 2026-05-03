# F01 CRvH Supplementary Panel E: Cancer Clinical Characteristics
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F01/a_script/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

PW <- 170; PH <- 80
RPT <- "04_Figures/F01/CRvH/b_reports/supp"
DAT <- "04_Figures/F01/CRvH/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

meta <- read.csv("00_input/CRm_meta.csv", stringsAsFactors = FALSE)

# CR subjects only, one row per subject (T1)
cr <- meta %>%
  filter(timepoint == "T1", cancer == "SURV") %>%
  select(pid, stage, yrs_since_dx, yrs_since_chemo)

# --- Left: Cancer stage distribution
stage_df <- cr %>%
  count(stage) %>%
  mutate(stage = factor(stage))

pE_left <- ggplot(stage_df, aes(x = stage, y = n)) +
  geom_bar(stat = "identity", width = 0.6, fill = SUPP_COLORS["CRE"],
           color = "grey30", linewidth = 0.3) +
  geom_text(aes(label = n), vjust = -0.5, size = 3, fontface = "bold") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20)),
                     breaks = scales::pretty_breaks(n = 5)) +
  labs(title = "Cancer Stage Distribution",
       subtitle = sprintf("CR subjects (n = %d)", nrow(cr)),
       y = "Count", x = "Cancer Stage", tag = "E") +
  FIG_THEME +
  theme(plot.subtitle = element_text(size = 8, color = "grey30",
                                     face = "bold.italic"))

# --- Middle: Years since diagnosis
pE_mid <- ggplot(cr, aes(x = "CR", y = yrs_since_dx)) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           fill = SUPP_COLORS["CRE"], alpha = 0.08,
           color = "grey85", linewidth = 0.15) +
  geom_bar(stat = "summary", fun = mean, width = 0.55,
           fill = SUPP_COLORS["CRE"], color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se,
                width = 0.15, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.5, alpha = 0.35,
              shape = 16, color = "grey30") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Years Since Dx", y = "Years", x = NULL) +
  FIG_THEME

# --- Right: Years since chemo
pE_right <- ggplot(cr, aes(x = "CR", y = yrs_since_chemo)) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           fill = SUPP_COLORS["CRE"], alpha = 0.08,
           color = "grey85", linewidth = 0.15) +
  geom_bar(stat = "summary", fun = mean, width = 0.55,
           fill = SUPP_COLORS["CRE"], color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se,
                width = 0.15, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.5, alpha = 0.35,
              shape = 16, color = "grey30") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Years Since Chemo", y = "Years", x = NULL) +
  FIG_THEME

pE <- pE_left | pE_mid | pE_right

# --- Audit CSV
audit_E <- cr %>%
  summarise(
    n = n(),
    stage_1 = sum(stage == 1, na.rm = TRUE),
    stage_2 = sum(stage == 2, na.rm = TRUE),
    stage_3 = sum(stage == 3, na.rm = TRUE),
    stage_4 = sum(stage == 4, na.rm = TRUE),
    yrs_dx_mean = mean(yrs_since_dx, na.rm = TRUE),
    yrs_dx_sd   = sd(yrs_since_dx, na.rm = TRUE),
    yrs_dx_median = median(yrs_since_dx, na.rm = TRUE),
    yrs_chemo_mean = mean(yrs_since_chemo, na.rm = TRUE),
    yrs_chemo_sd   = sd(yrs_since_chemo, na.rm = TRUE),
    yrs_chemo_median = median(yrs_since_chemo, na.rm = TRUE)
  )
write.csv(audit_E, file.path(DAT, "panel_E_clinical.csv"), row.names = FALSE)

ggsave(file.path(RPT, "panel_E_clinical.pdf"), pE,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT, "panel_E_clinical.png"), pE,
       width = PW, height = PH, units = "mm", dpi = 300)
cat("F01 CRvH Supp Panel E done\n")
