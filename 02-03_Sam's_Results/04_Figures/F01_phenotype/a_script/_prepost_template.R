# F01 pre/post panel template — sourced by scripts that define a `cfg` list.
# Sam CvH adaptation of YvO F01 _prepost_template.R.
#
# Required cfg: pre_col, post_col, y_label, delta_label, title, tag,
#   output_prefix, file_tag, audit_file, rpt_png, rpt_pdf, dat, file_prefix
# Optional cfg: y_breaks, y_labels, left_margin, right_margin, left_x_expand
#
# Data layout:
#   Sam's metadata is WIDE — pre_ and post_ values sit in the same row.
#   We deduplicate to one row per subject (prefer T1; fall back to T2 for
#   T2-only subjects), then compute delta = post - pre.
#   Statistical model: paired t-test within each arm (pre vs post) +
#   unpaired t-test on deltas (CRE vs PLA).

stopifnot(exists("cfg"), is.list(cfg))

# ── helpers ─────────────────────────────────────────────────────────────────
# These mirror the YvO style.R helpers; defined here so Sam's style.R
# (which sources the CvH shared style) doesn't need to re-export them.
.fmt_p_plot <- function(p, threshold = 0.05) {
  label <- if (p < 0.001) "p < 0.001" else if (p < 0.01) sprintf("p = %.3f", p) else sprintf("p = %.2f", p)
  if (p < threshold) paste0('bold("', label, '")') else paste0('"', label, '"')
}

.strip_for_composite <- function(p) {
  p + labs(title = NULL, subtitle = NULL, tag = NULL) +
    theme(legend.position = "none")
}

.composite_text_sizes <- function(comp_h_mm) {
  list(
    title    = pmax(6, pmin(8, round(5 + comp_h_mm / 80))),
    subtitle = pmax(4, pmin(6, round(3 + comp_h_mm / 100))),
    tag      = 8
  )
}

# ── libraries ────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ggsignif)
})

source("02-03_Sam's_Results/04_Figures/shared/style.R")

# cairo_pdf requires X11 (unavailable on this system); fall back to base pdf().
get_pdf_device <- function() grDevices::pdf

# ── defaults ─────────────────────────────────────────────────────────────────
if (is.null(cfg$left_margin))   cfg$left_margin   <- margin(2, 2, 2, 2)
if (is.null(cfg$right_margin))  cfg$right_margin  <- margin(2, 2, 2, 2)
if (is.null(cfg$left_x_expand)) cfg$left_x_expand <- expansion(add = 0.3)

PW <- 170; PH <- 80
for (d in c(cfg$rpt_png, cfg$rpt_pdf, cfg$dat))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ── load data ────────────────────────────────────────────────────────────────
sam <- readRDS("02-03_Sam's_Results/00_input/01_normalized_DAList_SURV_stringent_muscle.RDS")
meta_full <- as.data.frame(sam$metadata)

# SURV subjects only, one row per subject (prefer T1, fall back to T2)
surv_raw <- meta_full[meta_full$supp != "" & !is.na(meta_full$supp), ]
surv_raw  <- surv_raw[order(surv_raw$pid, surv_raw$timepoint), ]
meta      <- surv_raw[!duplicated(surv_raw$pid), ]

meta <- meta |>
  mutate(
    Supp_Group = factor(supp, levels = c("CRE", "PLA")),
    Supp_Time_Pre  = paste0(supp, "_T1"),
    Supp_Time_Post = paste0(supp, "_T2")
  )

PRE  <- cfg$pre_col
POST <- cfg$post_col

# ── reshape to paired-delta frame ────────────────────────────────────────────
pheno <- meta |>
  select(pid, Supp_Group, DV_Pre = all_of(PRE), DV_Post = all_of(POST)) |>
  filter(!is.na(DV_Pre), !is.na(DV_Post)) |>
  mutate(delta_DV = DV_Post - DV_Pre)

grp_cre <- pheno |> filter(Supp_Group == "CRE")
grp_pla <- pheno |> filter(Supp_Group == "PLA")

# ── statistics ───────────────────────────────────────────────────────────────
stats_paired_cre <- t.test(grp_cre$DV_Post, grp_cre$DV_Pre, paired = TRUE)
stats_paired_pla <- t.test(grp_pla$DV_Post, grp_pla$DV_Pre, paired = TRUE)
stats_delta      <- if (length(unique(pheno$Supp_Group)) == 2) {
  t.test(delta_DV ~ Supp_Group, data = pheno)
} else {
  list(p.value = NA, statistic = NA, parameter = NA,
       estimate = c(0, 0), conf.int = c(NA, NA))
}

sw_dc <- shapiro.test(grp_cre$delta_DV)
sw_dp <- if (nrow(grp_pla) >= 3) shapiro.test(grp_pla$delta_DV) else list(p.value = NA)

n_cre <- nrow(grp_cre); n_pla <- nrow(grp_pla)
anova_sub <- sprintf(
  "CRE: %s   PLA: %s   CRE vs PLA delta: %s",
  if (n_cre >= 2) {
    p <- stats_paired_cre$p.value
    if (p < 0.001) "p < 0.001" else if (p < 0.01) sprintf("p = %.3f", p) else sprintf("p = %.2f", p)
  } else "n < 2",
  if (n_pla >= 2) {
    p <- stats_paired_pla$p.value
    if (p < 0.001) "p < 0.001" else if (p < 0.01) sprintf("p = %.3f", p) else sprintf("p = %.2f", p)
  } else "n < 2",
  if (!is.na(stats_delta$p.value)) {
    p <- stats_delta$p.value
    if (p < 0.001) "p < 0.001" else if (p < 0.01) sprintf("p = %.3f", p) else sprintf("p = %.2f", p)
  } else "NA"
)

# ── audit CSV ────────────────────────────────────────────────────────────────
audit_df <- data.frame(
  test      = c("paired_t_CRE", "paired_t_PLA", "unpaired_t_delta"),
  Group     = c("CRE", "PLA", "CRE vs PLA"),
  n         = c(n_cre, n_pla, nrow(pheno)),
  statistic = c(stats_paired_cre$statistic,
                stats_paired_pla$statistic,
                if (is.list(stats_delta)) stats_delta$statistic else NA),
  p_value   = c(stats_paired_cre$p.value,
                stats_paired_pla$p.value,
                if (is.list(stats_delta)) stats_delta$p.value else NA),
  df        = c(stats_paired_cre$parameter,
                stats_paired_pla$parameter,
                if (is.list(stats_delta)) stats_delta$parameter else NA),
  mean_diff = c(stats_paired_cre$estimate,
                stats_paired_pla$estimate,
                if (is.list(stats_delta)) diff(rev(stats_delta$estimate)) else NA),
  ci_lo     = c(stats_paired_cre$conf.int[1],
                stats_paired_pla$conf.int[1],
                if (is.list(stats_delta)) stats_delta$conf.int[1] else NA),
  ci_hi     = c(stats_paired_cre$conf.int[2],
                stats_paired_pla$conf.int[2],
                if (is.list(stats_delta)) stats_delta$conf.int[2] else NA),
  shapiro_p = c(sw_dc$p.value, sw_dp$p.value, NA))
write.csv(audit_df, file.path(cfg$dat, cfg$audit_file), row.names = FALSE)

# ── long frame for left panel (pre/post bars) ────────────────────────────────
pheno_long <- pheno |>
  pivot_longer(cols = c(DV_Pre, DV_Post),
               names_to = "Timepoint", values_to = "value") |>
  mutate(
    Timepoint  = factor(sub("DV_", "", Timepoint), levels = c("Pre", "Post")),
    Supp_Time  = factor(paste0(as.character(Supp_Group), "_",
                               ifelse(Timepoint == "Pre", "T1", "T2")),
                        levels = c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2"))
  )

# Supp background colors
SUPP_BG <- c(CRE = scales::alpha(GROUP_COLORS["CRE_T1"], 0.12),
             PLA = scales::alpha(GROUP_COLORS["PLA_T1"], 0.12))

y_max_left <- max(pheno_long$value, na.rm = TRUE)

p_left <- ggplot(pheno_long, aes(Supp_Time, value, fill = Supp_Time)) +
  annotate("rect", xmin = 0.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = SUPP_BG["CRE"], color = "grey85", linewidth = 0.15) +
  annotate("rect", xmin = 2.5, xmax = 4.5, ymin = -Inf, ymax = Inf,
           fill = SUPP_BG["PLA"], color = "grey85", linewidth = 0.15) +
  geom_bar(stat = "summary", fun = mean, width = 0.65,
           color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se,
                width = 0.2, linewidth = 0.4) +
  geom_point(size = 0.9, alpha = 0.45, shape = 21,
             color = "black", stroke = 0.2,
             position = position_jitter(width = 0.12, seed = 42)) +
  {
    if (n_cre >= 2)
      geom_signif(comparisons = list(c("CRE_T1", "CRE_T2")),
                  annotations = .fmt_p_plot(stats_paired_cre$p.value),
                  parse = TRUE, y_position = y_max_left * 1.05,
                  textsize = 1.5, size = 0.3, tip_length = 0.01)
  } +
  {
    if (n_pla >= 2)
      geom_signif(comparisons = list(c("PLA_T1", "PLA_T2")),
                  annotations = .fmt_p_plot(stats_paired_pla$p.value),
                  parse = TRUE, y_position = y_max_left * 1.05,
                  textsize = 1.5, size = 0.3, tip_length = 0.01)
  } +
  scale_fill_manual(values = GROUP_FILL[c("CRE_T1","CRE_T2","PLA_T1","PLA_T2")]) +
  scale_x_discrete(
    labels = c(CRE_T1 = "Pre", CRE_T2 = "Post", PLA_T1 = "Pre", PLA_T2 = "Post"),
    expand = cfg$left_x_expand) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = cfg$title, subtitle = anova_sub,
       y = cfg$y_label, x = NULL, tag = cfg$tag) +
  FIG_THEME + theme(plot.margin = cfg$left_margin, legend.position = "none")

# ── right panel: delta bars ───────────────────────────────────────────────────
delta_bar_colors <- c(CRE = unname(GROUP_FILL["CRE_T2"]),
                      PLA = unname(GROUP_FILL["PLA_T2"]))
y_max_right <- max(abs(pheno$delta_DV), na.rm = TRUE)
y_lim_right <- max(pheno$delta_DV, na.rm = TRUE)

p_right <- ggplot(pheno, aes(Supp_Group, delta_DV, fill = Supp_Group)) +
  annotate("rect", xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf,
           fill = SUPP_BG["CRE"], color = "grey85", linewidth = 0.15) +
  annotate("rect", xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf,
           fill = SUPP_BG["PLA"], color = "grey85", linewidth = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_bar(stat = "summary", fun = mean, width = 0.55,
           color = "grey30", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se,
                width = 0.15, linewidth = 0.4) +
  geom_point(size = 0.9, alpha = 0.45, shape = 21,
             color = "black", stroke = 0.2,
             position = position_jitter(width = 0.12, seed = 42)) +
  {
    if (!is.na(stats_delta$p.value) && n_cre >= 2 && n_pla >= 2)
      geom_signif(comparisons = list(c("CRE", "PLA")),
                  annotations = .fmt_p_plot(stats_delta$p.value),
                  parse = TRUE, textsize = 1.5, size = 0.3, tip_length = 0.02,
                  y_position = y_lim_right * 1.10)
  } +
  scale_fill_manual(values = delta_bar_colors) +
  scale_x_discrete(expand = expansion(add = 0.3)) +
  labs(y = cfg$delta_label, x = NULL) +
  FIG_THEME + theme(legend.position = "none",
                    axis.title.y = element_text(margin = margin(r = 1)),
                    plot.margin = cfg$right_margin)

if (!is.null(cfg$y_breaks)) {
  p_right <- p_right +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.22)),
                       breaks = cfg$y_breaks, labels = cfg$y_labels)
} else {
  p_right <- p_right + scale_y_continuous(expand = expansion(mult = c(0.08, 0.22)))
}

# ── per-panel PNG + PDF ───────────────────────────────────────────────────────
p_combo <- (p_left | p_right) + plot_layout(widths = c(0.65, 0.35))
ggsave(file.path(cfg$rpt_png, paste0(cfg$file_prefix, "_", cfg$file_tag, ".png")),
       p_combo, width = PW, height = PH, units = "mm", dpi = 300)
ggsave(file.path(cfg$rpt_pdf, paste0(cfg$file_prefix, "_", cfg$file_tag, ".pdf")),
       p_combo, width = PW, height = PH, units = "mm", device = get_pdf_device())
message(sprintf("F01 %s done", cfg$title))

# ── export stripped versions for composite ────────────────────────────────────
pfx <- cfg$output_prefix
assign(paste0(pfx, "_title"),    cfg$title,    envir = .GlobalEnv)
assign(paste0(pfx, "_subtitle"), anova_sub,    envir = .GlobalEnv)
assign(paste0(pfx, "_left"),  .strip_for_composite(p_left),  envir = .GlobalEnv)
assign(paste0(pfx, "_right"), .strip_for_composite(p_right), envir = .GlobalEnv)
