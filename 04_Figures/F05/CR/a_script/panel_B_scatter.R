# F05/CR Panel B: Concordance Scatter -- logFC Training_CRE vs Training_PLA
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/F05/a_script/style.R")
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(boot)
})

PD_W <- 200

RPT <- "04_Figures/F05/CR/b_reports"
DAT <- "04_Figures/F05/CR/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(DAT, "panel_B"), recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

dep_df <- read_csv("03_DEP/c_data/03_combined_results_CR.csv", show_col_types = FALSE)
imputation_df <- read_csv("02_Imputation/c_data/02_mar_mnar_classification.csv",
                           show_col_types = FALSE) %>%
  transmute(gene, imputed = classification != "Complete")

scatter_df <- dep_df %>%
  transmute(
    gene,
    logFC_Training_CRE, logFC_Training_PLA,
    pi_CRE = pi_score_Training_CRE,
    pi_PLA = pi_score_Training_PLA,
    pi_Int = pi_score_Supplement_Interaction
  ) %>%
  filter(!is.na(logFC_Training_CRE), !is.na(logFC_Training_PLA)) %>%
  left_join(imputation_df, by = "gene") %>%
  mutate(
    imputed = replace_na(imputed, FALSE),
    significance = classify_proteins_f3(pi_CRE, pi_PLA, pi_Int),
    quadrant = case_when(
      logFC_Training_CRE > 0 & logFC_Training_PLA > 0 ~ "Concordant Up",
      logFC_Training_CRE < 0 & logFC_Training_PLA < 0 ~ "Concordant Down",
      logFC_Training_CRE > 0 & logFC_Training_PLA < 0 ~ "Discordant (CRE Up / PLA Down)",
      TRUE ~ "Discordant (CRE Down / PLA Up)"
    ),
    border_col = ifelse(imputed, "grey50", "grey75"),
    point_size = ifelse(significance == "NS", 1.8, 2.3),
    point_stroke = ifelse(significance == "NS", 0.6, 0.9),
    bubble_alpha = case_when(
      significance == "NS"          ~ 0.30,
      significance == "Interaction" ~ 0.55,
      significance == "Sig Both"    ~ 0.75,
      TRUE ~ 0.85
    )
  )

# Correlation with 95% CIs
cor_r   <- cor.test(scatter_df$logFC_Training_CRE, scatter_df$logFC_Training_PLA,
                    method = "pearson",  conf.level = 0.95)
cor_rho <- cor.test(scatter_df$logFC_Training_CRE, scatter_df$logFC_Training_PLA,
                    method = "spearman", conf.level = 0.95)
n_obs <- nrow(scatter_df)
rho_z <- atanh(cor_rho$estimate)
rho_se <- 1 / sqrt(n_obs - 3)
rho_ci <- tanh(rho_z + c(-1, 1) * qnorm(0.975) * rho_se)

# Sign concordance on ALL proteins
sign_concordance <- mean(sign(scatter_df$logFC_Training_CRE) ==
                         sign(scatter_df$logFC_Training_PLA)) * 100

# Bootstrap 95% CI for sign concordance (BCa, 10000 replicates)
set.seed(42)
boot_sign_conc <- boot::boot(
  data = scatter_df,
  statistic = function(d, i) {
    mean(sign(d$logFC_Training_CRE[i]) == sign(d$logFC_Training_PLA[i])) * 100
  },
  R = 10000
)
boot_ci <- tryCatch(
  boot::boot.ci(boot_sign_conc, type = "bca", conf = 0.95)$bca[4:5],
  error = function(e) quantile(boot_sign_conc$t, c(0.025, 0.975))
)

sig_mask <- scatter_df$significance != "NS"
conc_sig <- mean(sign(scatter_df$logFC_Training_CRE[sig_mask]) ==
                 sign(scatter_df$logFC_Training_PLA[sig_mask])) * 100
n_sig <- sum(sig_mask)

# Correlations on significant proteins only
cor_r_sig   <- cor.test(scatter_df$logFC_Training_CRE[sig_mask],
                        scatter_df$logFC_Training_PLA[sig_mask],
                        method = "pearson", conf.level = 0.95)
cor_rho_sig <- cor.test(scatter_df$logFC_Training_CRE[sig_mask],
                        scatter_df$logFC_Training_PLA[sig_mask],
                        method = "spearman", conf.level = 0.95)
# Bootstrap CI for Spearman rho (sig-only)
set.seed(43)
boot_rho_sig <- boot::boot(
  data = scatter_df[sig_mask, ],
  statistic = function(d, i)
    cor(d$logFC_Training_CRE[i], d$logFC_Training_PLA[i], method = "spearman"),
  R = 10000
)
rho_sig_ci <- tryCatch(
  boot::boot.ci(boot_rho_sig, type = "bca", conf = 0.95)$bca[4:5],
  error = function(e) quantile(boot_rho_sig$t, c(0.025, 0.975))
)

concordance_stats <- tibble(
  metric = c("Pearson_r", "Spearman_rho", "Sign_concordance_pct",
             "Pearson_r_sig", "Spearman_rho_sig", "Sign_concordance_sig_pct"),
  estimate = c(cor_r$estimate, cor_rho$estimate, sign_concordance,
               cor_r_sig$estimate, cor_rho_sig$estimate, conc_sig),
  ci_lower = c(cor_r$conf.int[1], rho_ci[1], boot_ci[1],
               cor_r_sig$conf.int[1], rho_sig_ci[1], NA_real_),
  ci_upper = c(cor_r$conf.int[2], rho_ci[2], boot_ci[2],
               cor_r_sig$conf.int[2], rho_sig_ci[2], NA_real_),
  p_value  = c(cor_r$p.value, cor_rho$p.value, NA_real_,
               cor_r_sig$p.value, cor_rho_sig$p.value, NA_real_),
  n        = c(n_obs, n_obs, n_obs, n_sig, n_sig, n_sig),
  note     = c("95% CI from cor.test()",
               "95% CI via Fisher z-transformation",
               "95% BCa bootstrap CI (10000 replicates, all proteins)",
               "Sig proteins only -- 95% CI from cor.test()",
               "Sig proteins only -- 95% BCa bootstrap CI (10000 replicates)",
               "Sign concordance among significant proteins only")
)
write_csv(concordance_stats, file.path(DAT, "panel_B", "concordance_stats.csv"))

x_range <- range(scatter_df$logFC_Training_CRE, na.rm = TRUE)
y_range <- range(scatter_df$logFC_Training_PLA, na.rm = TRUE)
x_pad   <- diff(x_range) * 0.10
y_pad   <- diff(y_range) * 0.10
xlim_range <- c(x_range[1] - x_pad, x_range[2] + x_pad)
ylim_range <- c(y_range[1] - y_pad, y_range[2] + y_pad)

q_df <- scatter_df %>%
  mutate(q = case_when(
    logFC_Training_CRE > 0 & logFC_Training_PLA > 0 ~ "Q1",
    logFC_Training_CRE < 0 & logFC_Training_PLA < 0 ~ "Q3",
    logFC_Training_CRE > 0 & logFC_Training_PLA < 0 ~ "Q4",
    TRUE ~ "Q2"
  ))
q_counts <- q_df %>% count(q) %>% deframe()
q_sig    <- q_df %>% filter(significance != "NS") %>% count(q) %>% deframe()
for (qq in c("Q1","Q2","Q3","Q4")) { if (is.na(q_sig[qq])) q_sig[qq] <- 0 }

label_df <- scatter_df %>%
  filter(significance != "NS") %>%
  group_by(significance) %>%
  arrange(desc(abs(logFC_Training_CRE) + abs(logFC_Training_PLA))) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  mutate(
    label_fill     = SIG_LABEL_FILL_F3[as.character(significance)],
    label_text_col = SIG_LABEL_TEXT_F3[as.character(significance)],
    nudge_y = case_when(
      significance == "Interaction"    ~ -0.03,
      significance == "Sig CRE only"   ~  0.03,
      significance == "Sig Both"       ~  0.04,
      significance == "Sig PLA only"   ~ -0.04,
      TRUE ~ 0
    )
  )

ns_df  <- scatter_df %>% filter(significance == "NS")
sig_df <- scatter_df %>% filter(significance != "NS")

txt_gene <- scale_text(BASE_GENE, PD_W)
txt_quad <- scale_text(BASE_QUADRANT, PD_W)

pD <- ggplot(mapping = aes(x = logFC_Training_CRE, y = logFC_Training_PLA)) +
  annotate("rect", xmin = 0, xmax = Inf,  ymin = 0, ymax = Inf,
           fill = "#FFE0E0", alpha = 0.55) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
           fill = "#FFE0E0", alpha = 0.55) +
  annotate("rect", xmin = 0, xmax = Inf,  ymin = -Inf, ymax = 0,
           fill = "#DCEEFF", alpha = 0.55) +
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
           fill = "#DCEEFF", alpha = 0.55) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "black", linewidth = 0.3) +
  geom_point(data = ns_df,
             aes(x = logFC_Training_CRE, y = logFC_Training_PLA),
             color = "grey80", fill = "grey85", shape = 21,
             size = 1.0, alpha = 0.3, stroke = 0.2) +
  geom_point(data = sig_df,
             aes(fill = significance),
             shape = 21,
             size = sig_df$point_size,
             color = sig_df$border_col,
             alpha = sig_df$bubble_alpha,
             stroke = sig_df$point_stroke) +
  scale_fill_manual(values = SIG_COLORS_F3, name = "Significance") +
  geom_label_repel(data = label_df, aes(label = gene),
                   fill = label_df$label_fill, color = label_df$label_text_col,
                   nudge_y = label_df$nudge_y,
                   size = txt_gene, fontface = "italic",
                   max.overlaps = 40,
                   segment.size = 0.2, segment.color = "grey50",
                   min.segment.length = 0, show.legend = FALSE,
                   box.padding = 0.6, force = 3, force_pull = 0.5,
                   label.padding = unit(1.5, "pt"),
                   label.r = unit(1, "pt"),
                   label.size = 0.15, seed = 42,
                   xlim = xlim_range * 0.9, ylim = ylim_range * 0.9) +
  annotate("label", x = Inf, y = Inf,
           label = sprintf("Concordant Up\u2002n = %s/%s", q_sig["Q1"], q_counts["Q1"]),
           hjust = 1, vjust = 1, size = txt_quad, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = -Inf,
           label = sprintf("Concordant Down\u2002n = %s/%s", q_sig["Q3"], q_counts["Q3"]),
           hjust = 0, vjust = 0, size = txt_quad, fontface = "bold",
           color = "#DC2626", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = Inf,
           label = sprintf("Discordant\u2002n = %s/%s", q_sig["Q2"], q_counts["Q2"]),
           hjust = 0, vjust = 1, size = txt_quad, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  annotate("label", x = Inf, y = -Inf,
           label = sprintf("Discordant\u2002n = %s/%s", q_sig["Q4"], q_counts["Q4"]),
           hjust = 1, vjust = 0, size = txt_quad, fontface = "bold",
           color = "#2563EB", fill = alpha("white", 0.9),
           label.padding = unit(2.5, "pt")) +
  coord_fixed(ratio = 1, xlim = xlim_range, ylim = ylim_range, expand = FALSE) +
  labs(
    title = "Protein-Level Concordance  |  limma + dupCor, missForest-imputed",
    subtitle = sprintf("All (n = %s): r = %.2f [%.2f, %.2f], \u03c1 = %.2f [%.2f, %.2f] | concordance = %.0f%%\nSig. (n = %d): r = %.2f [%.2f, %.2f], \u03c1 = %.2f [%.2f, %.2f] | concordance = %.0f%%",
                       format(n_obs, big.mark = ","),
                       cor_r$estimate, cor_r$conf.int[1], cor_r$conf.int[2],
                       cor_rho$estimate, rho_ci[1], rho_ci[2],
                       sign_concordance,
                       n_sig,
                       cor_r_sig$estimate, cor_r_sig$conf.int[1], cor_r_sig$conf.int[2],
                       cor_rho_sig$estimate, rho_sig_ci[1], rho_sig_ci[2],
                       conc_sig),
    x = expression(log[2]*FC ~ "(Training CRE)"),
    y = expression(log[2]*FC ~ "(Training PLA)")
  ) +
  FIG_THEME +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 8, face = "bold"),
    legend.text     = element_text(size = 7),
    legend.key.size = unit(3, "mm"),
    legend.margin   = margin(0, 0, 0, 0)
  ) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 3, alpha = 0.8)))

ggsave(file.path(RPT, "panel_B_concordance_MAIN.pdf"), pD,
       width = PD_W, height = PD_W, units = "mm", device = pdf_device)
ggsave(file.path(RPT, "panel_B_concordance_MAIN.png"), pD,
       width = PD_W, height = PD_W, units = "mm", dpi = 300)

scatter_df %>%
  transmute(
    gene,
    logFC_Training_CRE = round(logFC_Training_CRE, 4),
    logFC_Training_PLA = round(logFC_Training_PLA, 4),
    pi_score_CRE       = round(pi_CRE, 6),
    pi_score_PLA       = round(pi_PLA, 6),
    pi_score_Interaction = round(pi_Int, 6),
    significance         = as.character(significance),
    quadrant, imputed
  ) %>%
  arrange(significance, desc(abs(logFC_Training_CRE) + abs(logFC_Training_PLA))) %>%
  write_csv(file.path(DAT, "panel_B", "concordance.csv"))

message("F05/CR Panel B done")
