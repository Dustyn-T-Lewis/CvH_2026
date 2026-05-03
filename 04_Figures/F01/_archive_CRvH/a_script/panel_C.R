# F01 CRvH Panel C: Baseline Weight + CR Training Change
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F01/a_script/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(patchwork)
  library(ggsignif)
})

PW <- 170; PH <- 80
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

# --- Left panel: Baseline weight (CR vs H)
baseline <- subj %>%
  filter(!is.na(pre_weight_kg)) %>%
  select(pid, group, pre_weight_kg)

stats_baseline <- t.test(pre_weight_kg ~ group, data = baseline)

sw_cr_bl <- shapiro.test(baseline$pre_weight_kg[baseline$group == "CR"])
sw_h_bl  <- shapiro.test(baseline$pre_weight_kg[baseline$group == "H"])

# --- Right panel: Delta weight (CR only, complete cases)
cr_delta <- subj %>%
  filter(group == "CR", !is.na(pre_weight_kg), !is.na(post_weight_kg)) %>%
  mutate(delta_wt = post_weight_kg - pre_weight_kg)

stats_delta <- t.test(cr_delta$delta_wt, mu = 0)

sw_delta <- shapiro.test(cr_delta$delta_wt)

norm_sub <- sprintf("Shapiro-Wilk: CR %s, H %s, delta %s",
                    fmt_p(sw_cr_bl$p.value), fmt_p(sw_h_bl$p.value),
                    fmt_p(sw_delta$p.value))
full_sub <- sprintf("Baseline: %s | Delta (CR): %s\n%s",
                    fmt_p(stats_baseline$p.value),
                    fmt_p(stats_delta$p.value), norm_sub)

# --- Audit CSV
audit_C <- data.frame(
  test = c("unpaired_t_baseline", "one_sample_t_delta"),
  comparison = c("CR vs H", "CR delta vs 0"),
  statistic = c(stats_baseline$statistic, stats_delta$statistic),
  p_value = c(stats_baseline$p.value, stats_delta$p.value),
  df = c(stats_baseline$parameter, stats_delta$parameter),
  mean_diff = c(diff(rev(stats_baseline$estimate)), stats_delta$estimate),
  ci_lo = c(stats_baseline$conf.int[1], stats_delta$conf.int[1]),
  ci_hi = c(stats_baseline$conf.int[2], stats_delta$conf.int[2]),
  shapiro_p = c(NA, sw_delta$p.value),
  n_CR = c(sum(baseline$group == "CR"), nrow(cr_delta)),
  n_H  = c(sum(baseline$group == "H"), NA)
)
write.csv(audit_C, file.path(DAT, "panel_C_weight.csv"), row.names = FALSE)

# --- Left plot: Baseline bars
bar_colors_bl <- c(CR = unname(SUPP_COLORS["CRE"]),
                   H  = unname(SUPP_COLORS["H"]))

y_max_left <- max(baseline$pre_weight_kg, na.rm = TRUE)

pC_left <- ggplot(baseline, aes(x = group, y = pre_weight_kg, fill = group)) +
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
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.35,
              shape = 21, color = "black", stroke = 0.3) +
  geom_signif(
    comparisons = list(c("CR", "H")),
    annotations = fmt_p(stats_baseline$p.value),
    y_position = y_max_left * 1.10, textsize = 2.5, tip_length = 0.01
  ) +
  scale_fill_manual(values = bar_colors_bl) +
  scale_x_discrete(labels = c(CR = "Cancer\nRecovery", H = "Healthy")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(title = "Body Weight", subtitle = full_sub,
       y = "Weight (kg)", x = NULL, tag = "C") +
  FIG_THEME +
  theme(plot.subtitle = element_text(size = 7, color = "grey30",
                                     face = "bold.italic"),
        legend.position = "none")

# --- Right plot: Delta weight (CR only)
pC_right <- ggplot(cr_delta, aes(x = "CR", y = delta_wt)) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           fill = SUPP_COLORS["CRE"], alpha = 0.08,
           color = "grey85", linewidth = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
             linewidth = 0.3) +
  geom_bar(stat = "summary", fun = mean, width = 0.55,
           fill = SUPP_COLORS["CRE"], color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se,
                width = 0.15, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.35,
              shape = 21, color = "black", stroke = 0.3) +
  annotate("text", x = 1, y = max(abs(cr_delta$delta_wt)) * 1.25,
           label = paste("p", fmt_p(stats_delta$p.value)),
           size = 2.5, color = "grey30") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(y = expression(Delta * " Weight (kg)"), x = NULL) +
  FIG_THEME + theme(legend.position = "none")

pC <- (pC_left | pC_right) + plot_layout(widths = c(0.6, 0.4))

ggsave(file.path(RPT, "panel_C_weight.pdf"), pC,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT, "panel_C_weight.png"), pC,
       width = PW, height = PH, units = "mm", dpi = 300)
cat("F01 CRvH Panel C done\n")
