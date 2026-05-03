# F01 CR Panel B: LBM (pre/post + change by supplement)
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
RPT <- "04_Figures/F01/CR/b_reports"
DAT <- "04_Figures/F01/CR/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

meta <- read.csv("00_input/CRm_meta.csv", stringsAsFactors = FALSE)

# CR subjects only, T1 row, complete pre+post for LBM
subj <- meta %>%
  filter(timepoint == "T1", cancer == "SURV",
         !is.na(pre_LBM_kg), !is.na(post_LBM_kg)) %>%
  mutate(supp = factor(supp, levels = c("CRE", "PLA")),
         delta_lbm = post_LBM_kg - pre_LBM_kg)

# Long form for ANOVA + bar plots
subj_long <- subj %>%
  select(pid, supp, pre_LBM_kg, post_LBM_kg) %>%
  pivot_longer(cols = c(pre_LBM_kg, post_LBM_kg),
               names_to = "time", values_to = "lbm") %>%
  mutate(time = ifelse(grepl("^pre_", time), "Pre", "Post"),
         time = factor(time, levels = c("Pre", "Post")),
         supp_time = factor(paste0(supp, "_T",
                                   ifelse(time == "Pre", "1", "2")),
                            levels = c("CRE_T1", "CRE_T2",
                                       "PLA_T1", "PLA_T2")))

# --- Mixed ANOVA: Supplement (between) x Time (within)
stats_anova <- rstatix::anova_test(data = subj_long, dv = lbm,
                                    wid = pid,
                                    between = supp, within = time)

cre_subj <- subj %>% filter(supp == "CRE")
pla_subj <- subj %>% filter(supp == "PLA")
stats_paired_cre <- t.test(cre_subj$post_LBM_kg,
                           cre_subj$pre_LBM_kg, paired = TRUE)
stats_paired_pla <- t.test(pla_subj$post_LBM_kg,
                           pla_subj$pre_LBM_kg, paired = TRUE)
stats_delta      <- t.test(delta_lbm ~ supp, data = subj)

anova_tbl <- as.data.frame(stats_anova)
anova_sub <- sprintf("Supp %s   Time %s   Interaction %s",
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "supp"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "time"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "supp:time"]))

sw_cre <- shapiro.test(cre_subj$delta_lbm)
sw_pla <- shapiro.test(pla_subj$delta_lbm)
norm_sub <- sprintf("Shapiro-Wilk (delta): CRE %s, PLA %s",
                    fmt_p(sw_cre$p.value), fmt_p(sw_pla$p.value))

n_cre <- nrow(cre_subj); n_pla <- nrow(pla_subj)
full_sub <- paste0(anova_sub, "\n", norm_sub,
                   sprintf(" | CRE n=%d, PLA n=%d", n_cre, n_pla))

# --- Audit CSV
audit_B <- data.frame(
  test = c("paired_t_CRE", "paired_t_PLA", "unpaired_t_delta"),
  group = c("CRE", "PLA", "CRE vs PLA"),
  statistic = c(stats_paired_cre$statistic, stats_paired_pla$statistic,
                stats_delta$statistic),
  p_value = c(stats_paired_cre$p.value, stats_paired_pla$p.value,
              stats_delta$p.value),
  df = c(stats_paired_cre$parameter, stats_paired_pla$parameter,
         stats_delta$parameter),
  mean_diff = c(stats_paired_cre$estimate, stats_paired_pla$estimate,
                diff(stats_delta$estimate)),
  ci_lo = c(stats_paired_cre$conf.int[1], stats_paired_pla$conf.int[1],
            stats_delta$conf.int[1]),
  ci_hi = c(stats_paired_cre$conf.int[2], stats_paired_pla$conf.int[2],
            stats_delta$conf.int[2]),
  shapiro_p = c(sw_cre$p.value, sw_pla$p.value, NA)
)
write.csv(audit_B, file.path(DAT, "panel_B_lbm.csv"), row.names = FALSE)

# --- Left plot: Absolute pre/post by supplement
y_max_left <- max(subj_long$lbm, na.rm = TRUE)

pB_left <- ggplot(subj_long, aes(x = supp_time, y = lbm, fill = supp_time)) +
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = SUPP_COLORS["CRE"], alpha = 0.08,
           color = "grey85", linewidth = 0.15) +
  annotate("rect", xmin = 2.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
           fill = SUPP_COLORS["PLA"], alpha = 0.08,
           color = "grey85", linewidth = 0.15) +
  geom_bar(stat = "summary", fun = mean, width = 0.65,
           color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se,
                width = 0.2, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.35,
              shape = 21, color = "black", stroke = 0.3) +
  geom_signif(comparisons = list(c("CRE_T1", "CRE_T2")),
              annotations = fmt_p(stats_paired_cre$p.value),
              y_position = y_max_left * 1.05, textsize = 2.5,
              tip_length = 0.01) +
  geom_signif(comparisons = list(c("PLA_T1", "PLA_T2")),
              annotations = fmt_p(stats_paired_pla$p.value),
              y_position = y_max_left * 1.05, textsize = 2.5,
              tip_length = 0.01) +
  annotate("text", x = 1.5, y = -Inf, label = "Creatine",
           vjust = 4.2, fontface = "bold", size = 3.2, color = "grey25") +
  annotate("text", x = 3.5, y = -Inf, label = "Placebo",
           vjust = 4.2, fontface = "bold", size = 3.2, color = "grey25") +
  scale_fill_manual(values = SUPP_FILL) +
  scale_x_discrete(labels = c(CRE_T1 = "Pre", CRE_T2 = "Post",
                               PLA_T1 = "Pre", PLA_T2 = "Post")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  coord_cartesian(clip = "off") +
  labs(title = "Lean Body Mass", subtitle = full_sub,
       y = "LBM (kg)", x = NULL, tag = "B") +
  FIG_THEME +
  theme(plot.subtitle = element_text(size = 7, color = "grey30",
                                     face = "bold.italic"),
        plot.margin = margin(5, 5, 20, 5), legend.position = "none")

# --- Right plot: Delta by supplement
delta_bar_colors <- c(CRE = unname(SUPP_COLORS["CRE"]),
                      PLA = unname(SUPP_COLORS["PLA"]))

y_max_right <- max(subj$delta_lbm, na.rm = TRUE)

pB_right <- ggplot(subj, aes(x = supp, y = delta_lbm, fill = supp)) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           fill = SUPP_COLORS["CRE"], alpha = 0.08,
           color = "grey85", linewidth = 0.15) +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = SUPP_COLORS["PLA"], alpha = 0.08,
           color = "grey85", linewidth = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50",
             linewidth = 0.3) +
  geom_bar(stat = "summary", fun = mean, width = 0.55,
           color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se,
                width = 0.15, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.35,
              shape = 21, color = "black", stroke = 0.3) +
  geom_signif(comparisons = list(c("CRE", "PLA")),
              annotations = fmt_p(stats_delta$p.value),
              textsize = 2.5, tip_length = 0.02,
              y_position = y_max_right * 1.10) +
  scale_fill_manual(values = delta_bar_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(y = expression(Delta * " LBM (kg)"), x = NULL) +
  FIG_THEME + theme(legend.position = "none")

pB <- (pB_left | pB_right) + plot_layout(widths = c(0.65, 0.35))

ggsave(file.path(RPT, "panel_B_lbm.pdf"), pB,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT, "panel_B_lbm.png"), pB,
       width = PW, height = PH, units = "mm", dpi = 300)
cat("F01 CR Panel B done\n")
