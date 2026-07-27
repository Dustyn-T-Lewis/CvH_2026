# One pre/post x supplement outcome panel: absolute means on the left, the
# within-subject delta on the right. Panels D-G differ only in which outcome
# column they read and how they are labelled.

setwd(here::here())
source("04_Figures/shared/style.R")

pacman::p_load(dplyr, tidyr, patchwork, ggsignif, rstatix)

PW <- 170
PH <- 80
RPT <- "04_Figures/F01_Phenotype/b_reports/supp"
for (sub in c("pdf/panels", "png/panels")) {
  dir.create(file.path(RPT, sub), recursive = TRUE, showWarnings = FALSE)
}
DAT <- "04_Figures/F01_Phenotype/c_data/supp"

pre_post_panel <- function(pre_col, post_col, title, y_lab, tag, stem) {
  meta <- read.csv("00_input/CvH_meta.csv", stringsAsFactors = FALSE) |>
    dplyr::rename(pid = Subject_ID, timepoint = Timepoint, supp = Supplement)

  subj <- meta |>
    filter(
      timepoint == "T1", cancer == "SURV",
      !is.na(.data[[pre_col]]), !is.na(.data[[post_col]])
    ) |>
    mutate(
      supp = factor(supp, levels = c("CRE", "PLA")),
      delta = .data[[post_col]] - .data[[pre_col]]
    )

  subj_long <- subj |>
    select(pid, supp, all_of(c(pre_col, post_col))) |>
    pivot_longer(all_of(c(pre_col, post_col)),
      names_to = "time", values_to = "value"
    ) |>
    mutate(
      time = factor(ifelse(grepl("^pre_", time), "Pre", "Post"),
        levels = c("Pre", "Post")
      ),
      supp_time = factor(
        paste0(supp, "_T", ifelse(time == "Pre", "1", "2")),
        levels = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2")
      )
    )

  anova_tbl <- as.data.frame(rstatix::anova_test(
    data = subj_long, dv = value, wid = pid, between = supp, within = time
  ))

  cre_subj <- subj |> filter(supp == "CRE")
  pla_subj <- subj |> filter(supp == "PLA")
  t_cre <- t.test(cre_subj[[post_col]], cre_subj[[pre_col]], paired = TRUE)
  t_pla <- t.test(pla_subj[[post_col]], pla_subj[[pre_col]], paired = TRUE)
  t_delta <- t.test(delta ~ supp, data = subj)
  sw_cre <- shapiro.test(cre_subj$delta)
  sw_pla <- shapiro.test(pla_subj$delta)

  full_sub <- paste0(
    sprintf(
      "Supp %s   Time %s   Interaction %s",
      fmt_p(anova_tbl$p[anova_tbl$Effect == "supp"]),
      fmt_p(anova_tbl$p[anova_tbl$Effect == "time"]),
      fmt_p(anova_tbl$p[anova_tbl$Effect == "supp:time"])
    ), "\n",
    sprintf(
      "Shapiro-Wilk (delta): CRE %s, PLA %s | CRE n=%d, PLA n=%d",
      fmt_p(sw_cre$p.value), fmt_p(sw_pla$p.value),
      nrow(cre_subj), nrow(pla_subj)
    )
  )

  audit <- data.frame(
    test = c("paired_t_CRE", "paired_t_PLA", "unpaired_t_delta"),
    group = c("CRE", "PLA", "CRE vs PLA"),
    statistic = c(t_cre$statistic, t_pla$statistic, t_delta$statistic),
    p_value = c(t_cre$p.value, t_pla$p.value, t_delta$p.value),
    df = c(t_cre$parameter, t_pla$parameter, t_delta$parameter),
    mean_diff = c(t_cre$estimate, t_pla$estimate, diff(t_delta$estimate)),
    ci_lo = c(t_cre$conf.int[1], t_pla$conf.int[1], t_delta$conf.int[1]),
    ci_hi = c(t_cre$conf.int[2], t_pla$conf.int[2], t_delta$conf.int[2]),
    shapiro_p = c(sw_cre$p.value, sw_pla$p.value, NA)
  )
  write.csv(audit, file.path(DAT, paste0(stem, ".csv")), row.names = FALSE)

  lane <- function(xmin, xmax, key) {
    annotate("rect",
      xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf,
      fill = SUPP_COLORS[key], alpha = 0.08, color = "grey85", linewidth = 0.15
    )
  }
  jitter_pts <- geom_point(
    position = position_jitter(width = 0.12, seed = 42),
    size = 1.2, alpha = 0.35, shape = 21, color = "black", stroke = 0.3
  )

  y_max_left <- max(subj_long$value, na.rm = TRUE)
  left <- ggplot(subj_long, aes(x = supp_time, y = value, fill = supp_time)) +
    lane(0.5, 2.5, "CRE") +
    lane(2.5, 4.5, "PLA") +
    geom_bar(stat = "summary", fun = mean, width = 0.65, color = "grey30", linewidth = 0.3) +
    geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.2, linewidth = 0.4) +
    jitter_pts +
    geom_signif(
      comparisons = list(c("CRE_T1", "CRE_T2")),
      annotations = fmt_p(t_cre$p.value),
      y_position = y_max_left * 1.05, textsize = 2.5, tip_length = 0.01
    ) +
    geom_signif(
      comparisons = list(c("PLA_T1", "PLA_T2")),
      annotations = fmt_p(t_pla$p.value),
      y_position = y_max_left * 1.05, textsize = 2.5, tip_length = 0.01
    ) +
    annotate("text",
      x = 1.5, y = -Inf, label = "Creatine",
      vjust = 4.2, fontface = "bold", size = 3.2, color = "grey25"
    ) +
    annotate("text",
      x = 3.5, y = -Inf, label = "Placebo",
      vjust = 4.2, fontface = "bold", size = 3.2, color = "grey25"
    ) +
    scale_fill_manual(values = SUPP_FILL) +
    scale_x_discrete(labels = c(
      CRE_T1 = "Pre", CRE_T2 = "Post",
      PLA_T1 = "Pre", PLA_T2 = "Post"
    )) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    coord_cartesian(clip = "off") +
    labs(title = title, subtitle = full_sub, y = y_lab, x = NULL, tag = tag) +
    FIG_THEME +
    theme(
      plot.subtitle = element_text(size = 7, color = "grey40", face = "italic"),
      plot.margin = margin(5, 5, 20, 5), legend.position = "none"
    )

  # abs() so the bracket sits above the bars even when every delta is negative
  y_max_right <- max(abs(subj$delta), na.rm = TRUE)
  right <- ggplot(subj, aes(x = supp, y = delta, fill = supp)) +
    lane(0.5, 1.5, "CRE") +
    lane(1.5, 2.5, "PLA") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
    geom_bar(stat = "summary", fun = mean, width = 0.55, color = "grey30", linewidth = 0.3) +
    geom_errorbar(stat = "summary", fun.data = mean_se, width = 0.15, linewidth = 0.4) +
    jitter_pts +
    geom_signif(
      comparisons = list(c("CRE", "PLA")),
      annotations = fmt_p(t_delta$p.value),
      textsize = 2.5, tip_length = 0.02, y_position = y_max_right * 1.10
    ) +
    scale_fill_manual(values = c(
      CRE = unname(SUPP_COLORS["CRE"]), PLA = unname(SUPP_COLORS["PLA"])
    )) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    labs(y = bquote(Delta * .(paste0(" ", y_lab))), x = NULL) +
    FIG_THEME +
    theme(legend.position = "none")

  panel <- (left | right) + plot_layout(widths = c(0.65, 0.35))

  ggsave(file.path(RPT, "pdf/panels", paste0(stem, "_SUPP.pdf")), panel,
    width = PW, height = PH, units = "mm", device = get_pdf_device()
  )
  ggsave(file.path(RPT, "png/panels", paste0(stem, "_SUPP.png")), panel,
    width = PW, height = PH, units = "mm", dpi = 300
  )
  panel
}
