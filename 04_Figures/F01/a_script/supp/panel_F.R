# F01 Supplementary Panel F: Leg Extension (Pre/Post x Supplement)
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F01/a_script/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ggsignif)
  library(rstatix)
})

PW <- 170; PH <- 80
RPT <- "04_Figures/F01/b_reports/supp"
DAT <- "04_Figures/F01/c_data/supp"

meta <- read.csv("00_input/CvH_meta.csv", stringsAsFactors = FALSE) |>
  dplyr::rename(pid = Subject_ID, timepoint = Timepoint, supp = Supplement)

subj <- meta %>%
  filter(timepoint == "T1", cancer == "SURV",
         !is.na(pre_leg_ext_lbs), !is.na(post_leg_ext_lbs)) %>%
  mutate(supp = factor(supp, levels = c("CRE", "PLA")),
         delta_le = post_leg_ext_lbs - pre_leg_ext_lbs)

subj_long <- subj %>%
  select(pid, supp, pre_leg_ext_lbs, post_leg_ext_lbs) %>%
  pivot_longer(cols = c(pre_leg_ext_lbs, post_leg_ext_lbs),
               names_to = "time", values_to = "le") %>%
  mutate(time = ifelse(grepl("^pre_", time), "Pre", "Post"),
         time = factor(time, levels = c("Pre", "Post")),
         supp_time = factor(paste0(supp, "_T",
                                   ifelse(time == "Pre", "1", "2")),
                            levels = c("CRE_T1", "CRE_T2",
                                       "PLA_T1", "PLA_T2")))

stats_anova <- rstatix::anova_test(data = subj_long, dv = le,
                                    wid = pid,
                                    between = supp, within = time)

cre_subj <- subj %>% filter(supp == "CRE")
pla_subj <- subj %>% filter(supp == "PLA")
stats_paired_cre <- t.test(cre_subj$post_leg_ext_lbs,
                           cre_subj$pre_leg_ext_lbs, paired = TRUE)
stats_paired_pla <- t.test(pla_subj$post_leg_ext_lbs,
                           pla_subj$pre_leg_ext_lbs, paired = TRUE)
stats_delta      <- t.test(delta_le ~ supp, data = subj)

anova_tbl <- as.data.frame(stats_anova)
anova_sub <- sprintf("Supp %s   Time %s   Interaction %s",
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "supp"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "time"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "supp:time"]))

sw_cre <- shapiro.test(cre_subj$delta_le)
sw_pla <- shapiro.test(pla_subj$delta_le)
n_cre <- nrow(cre_subj); n_pla <- nrow(pla_subj)
full_sub <- paste0(anova_sub, "\n",
                   sprintf("Shapiro-Wilk (delta): CRE %s, PLA %s | CRE n=%d, PLA n=%d",
                           fmt_p(sw_cre$p.value), fmt_p(sw_pla$p.value), n_cre, n_pla))

audit_F <- data.frame(
  test = c("paired_t_CRE", "paired_t_PLA", "unpaired_t_delta"),
  group = c("CRE", "PLA", "CRE vs PLA"),
  statistic = c(stats_paired_cre$statistic, stats_paired_pla$statistic, stats_delta$statistic),
  p_value = c(stats_paired_cre$p.value, stats_paired_pla$p.value, stats_delta$p.value),
  df = c(stats_paired_cre$parameter, stats_paired_pla$parameter, stats_delta$parameter),
  mean_diff = c(stats_paired_cre$estimate, stats_paired_pla$estimate, diff(stats_delta$estimate)),
  ci_lo = c(stats_paired_cre$conf.int[1], stats_paired_pla$conf.int[1], stats_delta$conf.int[1]),
  ci_hi = c(stats_paired_cre$conf.int[2], stats_paired_pla$conf.int[2], stats_delta$conf.int[2]),
  shapiro_p = c(sw_cre$p.value, sw_pla$p.value, NA)
)
write.csv(audit_F, file.path(DAT, "panel_F_leg_ext.csv"), row.names = FALSE)

y_max_left <- max(subj_long$le, na.rm = TRUE)

pF_left <- ggplot(subj_long, aes(x = supp_time, y = le, fill = supp_time)) +
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = SUPP_COLORS["CRE"], alpha = 0.08, color = "grey85", linewidth = 0.15) +
  annotate("rect", xmin = 2.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
           fill = SUPP_COLORS["PLA"], alpha = 0.08, color = "grey85", linewidth = 0.15) +
  geom_bar(stat = "summary", fun = mean, width = 0.65, color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.5, shape = 21, color = "black", stroke = 0.3) +
  geom_signif(comparisons = list(c("CRE_T1", "CRE_T2")),
              annotations = fmt_p(stats_paired_cre$p.value),
              y_position = y_max_left * 1.05, textsize = 2.5, tip_length = 0.01) +
  geom_signif(comparisons = list(c("PLA_T1", "PLA_T2")),
              annotations = fmt_p(stats_paired_pla$p.value),
              y_position = y_max_left * 1.05, textsize = 2.5, tip_length = 0.01) +
  annotate("text", x = 1.5, y = -Inf, label = "Creatine",
           vjust = 4.2, fontface = "bold", size = 3.2, color = "grey25") +
  annotate("text", x = 3.5, y = -Inf, label = "Placebo",
           vjust = 4.2, fontface = "bold", size = 3.2, color = "grey25") +
  scale_fill_manual(values = SUPP_FILL) +
  scale_x_discrete(labels = c(CRE_T1 = "Pre", CRE_T2 = "Post",
                               PLA_T1 = "Pre", PLA_T2 = "Post")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  coord_cartesian(clip = "off") +
  labs(title = "Leg Extension", subtitle = full_sub,
       y = "Leg Extension (lbs)", x = NULL, tag = "F") +
  FIG_THEME +
  theme(plot.subtitle = element_text(size = 7, color = "grey40", face = "italic"),
        plot.margin = margin(5, 5, 20, 5), legend.position = "none")

delta_bar_colors <- c(CRE = unname(SUPP_COLORS["CRE"]), PLA = unname(SUPP_COLORS["PLA"]))
y_max_right <- max(subj$delta_le, na.rm = TRUE)

pF_right <- ggplot(subj, aes(x = supp, y = delta_le, fill = supp)) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           fill = SUPP_COLORS["CRE"], alpha = 0.08, color = "grey85", linewidth = 0.15) +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = SUPP_COLORS["PLA"], alpha = 0.08, color = "grey85", linewidth = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_bar(stat = "summary", fun = mean, width = 0.55, color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.15, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.5, shape = 21, color = "black", stroke = 0.3) +
  geom_signif(comparisons = list(c("CRE", "PLA")),
              annotations = fmt_p(stats_delta$p.value),
              textsize = 2.5, tip_length = 0.02, y_position = y_max_right * 1.10) +
  scale_fill_manual(values = delta_bar_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(y = expression(Delta * " Leg Extension (lbs)"), x = NULL) +
  FIG_THEME + theme(legend.position = "none")

pF <- (pF_left | pF_right) + plot_layout(widths = c(0.65, 0.35))

ggsave(file.path(RPT, "panel_F_leg_ext_SUPP.pdf"), pF,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT, "panel_F_leg_ext_SUPP.png"), pF,
       width = PW, height = PH, units = "mm", dpi = 300)
cat("F01 Supp Panel F done\n")
