# F01 CRvH Panel A: Age Distribution (CR vs Healthy)
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F01/a_script/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggsignif)
})

PW <- 90; PH <- 100
RPT <- "04_Figures/F01/CRvH/b_reports"
DAT <- "04_Figures/F01/CRvH/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

meta <- read.csv("00_input/CRm_meta.csv", stringsAsFactors = FALSE)

# One row per subject: T1 only
subj <- meta %>%
  filter(timepoint == "T1") %>%
  mutate(group = ifelse(cancer == "SURV", "CR", "H"),
         group = factor(group, levels = c("CR", "H")))

# --- Statistics
stats_A <- t.test(age ~ group, data = subj)

sw_cr <- shapiro.test(subj$age[subj$group == "CR"])
sw_h  <- shapiro.test(subj$age[subj$group == "H"])
norm_sub <- sprintf("Shapiro-Wilk: CR %s, H %s",
                    fmt_p(sw_cr$p.value), fmt_p(sw_h$p.value))

# --- Audit CSV
audit_A <- subj %>%
  group_by(group) %>%
  summarise(n = n(), mean = mean(age), sd = sd(age),
            sem = sd(age) / sqrt(n()), .groups = "drop") %>%
  mutate(shapiro_p = c(sw_cr$p.value, sw_h$p.value),
         t_test_p  = stats_A$p.value)
write.csv(audit_A, file.path(DAT, "panel_A_age.csv"), row.names = FALSE)

bar_colors <- c(CR = unname(SUPP_COLORS["CRE"]),
                H  = unname(SUPP_COLORS["H"]))

pA <- ggplot(subj, aes(x = group, y = age, fill = group)) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           fill = SUPP_COLORS["CRE"], alpha = 0.08,
           color = "grey85", linewidth = 0.15) +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = SUPP_COLORS["H"], alpha = 0.08,
           color = "grey85", linewidth = 0.15) +
  geom_bar(stat = "summary", fun = mean, width = 0.6,
           color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se,
                width = 0.2, linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.35,
              shape = 16, color = "grey30") +
  geom_signif(
    comparisons = list(c("CR", "H")),
    annotations = fmt_p(stats_A$p.value),
    textsize = 3, tip_length = 0.02,
    y_position = max(subj$age) * 1.15
  ) +
  scale_fill_manual(values = bar_colors) +
  scale_x_discrete(labels = c(CR = "Cancer Recovery", H = "Healthy")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Age", subtitle = norm_sub,
       y = "Age (years)", x = NULL, tag = "A") +
  FIG_THEME + theme(legend.position = "none",
                    plot.subtitle = element_text(size = 8, color = "grey30",
                                                face = "bold.italic"))

ggsave(file.path(RPT, "panel_A_age.pdf"), pA,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT, "panel_A_age.png"), pA,
       width = PW, height = PH, units = "mm", dpi = 300)
cat("F01 CRvH Panel A done\n")
