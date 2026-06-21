# F01 Panel B: ALM (Pre/Post x Supplement + Healthy reference)
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
RPT <- "04_Figures/F01/b_reports"
DAT <- "04_Figures/F01/c_data"

meta <- read.csv("00_input/CvH_meta.csv", stringsAsFactors = FALSE) |>
  dplyr::rename(pid = Subject_ID, timepoint = Timepoint, supp = Supplement)

# --- CR subjects: complete pre+post ALM (one row per subject = T1 row)
cr_subj <- meta %>%
  filter(timepoint == "T1", cancer == "SURV",
         !is.na(pre_ALM_kg), !is.na(post_ALM_kg)) %>%
  mutate(supp = factor(supp, levels = c("CRE", "PLA")),
         delta_alm = post_ALM_kg - pre_ALM_kg)

# --- Healthy subjects: baseline ALM only
h_subj <- meta %>%
  filter(timepoint == "T1", cancer == "CTL", !is.na(pre_ALM_kg))

# --- Long form for CR (ANOVA + bars)
cr_long <- cr_subj %>%
  select(pid, supp, pre_ALM_kg, post_ALM_kg) %>%
  pivot_longer(cols = c(pre_ALM_kg, post_ALM_kg),
               names_to = "time", values_to = "alm") %>%
  mutate(time = ifelse(grepl("^pre_", time), "Pre", "Post"),
         time = factor(time, levels = c("Pre", "Post")),
         supp_time = factor(paste0(supp, "_T",
                                   ifelse(time == "Pre", "1", "2")),
                            levels = c("CRE_T1", "CRE_T2",
                                       "PLA_T1", "PLA_T2")))

# --- Combined long form for 5-bar plot
all_long <- bind_rows(
  h_subj %>% transmute(pid, supp_time = "H_T1", alm = pre_ALM_kg),
  cr_long %>% select(pid, supp_time, alm)
) %>%
  mutate(supp_time = factor(supp_time,
    levels = c("H_T1", "CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2")))

# --- Statistics ---

# 1. Mixed ANOVA: Supplement x Time (CR only)
stats_anova <- rstatix::anova_test(data = cr_long, dv = alm,
                                    wid = pid,
                                    between = supp, within = time)

# 2. Baseline: pooled CR_T1 vs H_T1
baseline_df <- bind_rows(
  cr_subj %>% transmute(pid, group = "CR", pre_ALM_kg),
  h_subj  %>% transmute(pid, group = "H",  pre_ALM_kg)
)
stats_baseline <- t.test(pre_ALM_kg ~ group, data = baseline_df)

# 3. Paired t-tests within supplement groups
cre_subj <- cr_subj %>% filter(supp == "CRE")
pla_subj <- cr_subj %>% filter(supp == "PLA")
stats_paired_cre <- t.test(cre_subj$post_ALM_kg,
                           cre_subj$pre_ALM_kg, paired = TRUE)
stats_paired_pla <- t.test(pla_subj$post_ALM_kg,
                           pla_subj$pre_ALM_kg, paired = TRUE)

# 4. Delta comparison: CRE vs PLA
stats_delta <- t.test(delta_alm ~ supp, data = cr_subj)

# 5. Shapiro-Wilk on deltas
sw_cre <- shapiro.test(cre_subj$delta_alm)
sw_pla <- shapiro.test(pla_subj$delta_alm)

# --- Subtitle
anova_tbl <- as.data.frame(stats_anova)
anova_sub <- sprintf("Supp %s   Time %s   Interaction %s",
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "supp"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "time"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "supp:time"]))

n_cre <- nrow(cre_subj); n_pla <- nrow(pla_subj); n_h <- nrow(h_subj)
norm_sub <- sprintf("H vs CR(BL) %s | Shapiro-Wilk (delta): CRE %s, PLA %s | CRE n=%d, PLA n=%d, H n=%d",
                    fmt_p(stats_baseline$p.value),
                    fmt_p(sw_cre$p.value), fmt_p(sw_pla$p.value),
                    n_cre, n_pla, n_h)
full_sub <- paste0(anova_sub, "\n", norm_sub)

# --- Audit CSV
audit_B <- data.frame(
  test = c("mixed_anova_supp", "mixed_anova_time", "mixed_anova_int",
           "baseline_CR_vs_H", "paired_t_CRE", "paired_t_PLA",
           "unpaired_t_delta"),
  comparison = c("Supplement", "Time", "Supp:Time",
                 "CR(BL) vs H", "CRE Pre vs Post", "PLA Pre vs Post",
                 "CRE delta vs PLA delta"),
  statistic = c(anova_tbl$F[anova_tbl$Effect == "supp"],
                anova_tbl$F[anova_tbl$Effect == "time"],
                anova_tbl$F[anova_tbl$Effect == "supp:time"],
                stats_baseline$statistic,
                stats_paired_cre$statistic, stats_paired_pla$statistic,
                stats_delta$statistic),
  p_value = c(anova_tbl$p[anova_tbl$Effect == "supp"],
              anova_tbl$p[anova_tbl$Effect == "time"],
              anova_tbl$p[anova_tbl$Effect == "supp:time"],
              stats_baseline$p.value,
              stats_paired_cre$p.value, stats_paired_pla$p.value,
              stats_delta$p.value),
  shapiro_p = c(NA, NA, NA, NA, sw_cre$p.value, sw_pla$p.value, NA),
  n = c(n_cre + n_pla, n_cre + n_pla, n_cre + n_pla,
        n_cre + n_pla + n_h, n_cre, n_pla, n_cre + n_pla)
)
write.csv(audit_B, file.path(DAT, "panel_B_alm.csv"), row.names = FALSE)

# --- Left plot: 5-bar (H_T1 | CRE_T1, CRE_T2 | PLA_T1, PLA_T2)
y_max_left <- max(all_long$alm, na.rm = TRUE)

pB_left <- ggplot(all_long, aes(x = supp_time, y = alm, fill = supp_time)) +
  # Lane shadings: Healthy | Creatine | Placebo
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           fill = SUPP_COLORS["H"], alpha = 0.08,
           color = "grey85", linewidth = 0.15) +
  annotate("rect", xmin = 1.5, xmax = 3.5, ymin = -Inf, ymax = Inf,
           fill = SUPP_COLORS["CRE"], alpha = 0.08,
           color = "grey85", linewidth = 0.15) +
  annotate("rect", xmin = 3.5, xmax = 5.5, ymin = -Inf, ymax = Inf,
           fill = SUPP_COLORS["PLA"], alpha = 0.08,
           color = "grey85", linewidth = 0.15) +
  geom_bar(stat = "summary", fun = mean, width = 0.65,
           color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se,
                width = 0.2, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.35,
              shape = 21, color = "black", stroke = 0.3) +
  # Paired t-test brackets within CRE and PLA
  geom_signif(comparisons = list(c("CRE_T1", "CRE_T2")),
              annotations = fmt_p(stats_paired_cre$p.value),
              y_position = y_max_left * 1.05, textsize = 2.5,
              tip_length = 0.01) +
  geom_signif(comparisons = list(c("PLA_T1", "PLA_T2")),
              annotations = fmt_p(stats_paired_pla$p.value),
              y_position = y_max_left * 1.05, textsize = 2.5,
              tip_length = 0.01) +
  # Group labels below x-axis
  annotate("text", x = 1.0, y = -Inf, label = "Healthy",
           vjust = 4.2, fontface = "bold", size = 3.0, color = "grey25") +
  annotate("text", x = 2.5, y = -Inf, label = "Creatine",
           vjust = 4.2, fontface = "bold", size = 3.0, color = "grey25") +
  annotate("text", x = 4.5, y = -Inf, label = "Placebo",
           vjust = 4.2, fontface = "bold", size = 3.0, color = "grey25") +
  scale_fill_manual(values = SUPP_FILL) +
  scale_x_discrete(labels = c(H_T1 = "BL", CRE_T1 = "Pre", CRE_T2 = "Post",
                               PLA_T1 = "Pre", PLA_T2 = "Post")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  coord_cartesian(clip = "off") +
  labs(title = "Appendicular Lean Mass", subtitle = full_sub,
       y = "ALM (kg)", x = NULL, tag = "B") +
  FIG_THEME +
  theme(plot.subtitle = element_text(size = 7, color = "grey30",
                                     face = "bold.italic"),
        plot.margin = margin(5, 5, 20, 5), legend.position = "none")

# --- Right plot: Delta by supplement (CRE vs PLA)
delta_bar_colors <- c(CRE = unname(SUPP_COLORS["CRE"]),
                      PLA = unname(SUPP_COLORS["PLA"]))

y_max_right <- max(abs(cr_subj$delta_alm), na.rm = TRUE)

pB_right <- ggplot(cr_subj, aes(x = supp, y = delta_alm, fill = supp)) +
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
  labs(y = expression(Delta * " ALM (kg)"), x = NULL) +
  FIG_THEME + theme(legend.position = "none")

pB <- (pB_left | pB_right) + plot_layout(widths = c(0.70, 0.30))

ggsave(file.path(RPT, "panel_B_alm_MAIN.pdf"), pB,
       width = PW, height = PH, units = "mm", device = get_pdf_device())
ggsave(file.path(RPT, "panel_B_alm_MAIN.png"), pB,
       width = PW, height = PH, units = "mm", dpi = 300)
cat("F01 Panel B done\n")
