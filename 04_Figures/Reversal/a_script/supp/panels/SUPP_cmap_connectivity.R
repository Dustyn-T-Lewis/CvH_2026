# SUPP: CMap Connectivity Score for Reversal
# Lamb et al. 2006, Science (PMID 17008526): Connectivity Map / KS enrichment
# Subramanian et al. 2017, Cell (PMID 29195078): L1000 / next-gen CMap
# Samart et al. 2021, Brief Bioinform (PMID 34013329): score reconciliation
#
# Weighted KS enrichment score:
#   Query = cancer-up / cancer-down gene sets (Pi < 0.05)
#   Reference = all proteins ranked by Training_CR t-statistic
#   Connectivity = (ES_up - ES_down) / 2
#   Negative connectivity = reversal; positive = exacerbation
#   Null: 10,000 random gene-set permutations
setwd(here::here())
source("04_Figures/shared/style.R")
pacman::p_load(tidyverse, patchwork)

RPT_PNG <- "04_Figures/Reversal/b_reports/supp/png/panels"
RPT_PDF <- "04_Figures/Reversal/b_reports/supp/pdf/panels"
DAT     <- "04_Figures/Reversal/c_data/panel_supp"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT,     recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()

# ── Data ─────────────────────────────────────────────────────────────────────
dep_df <- read_csv("03_DEP/c_data/03_combined_results_CRvH.csv",
                   show_col_types = FALSE) %>%
  filter(!is.na(logFC_Cancer_vs_Healthy), !is.na(logFC_Training_CR),
         !is.na(t_Training_CR), !is.na(t_Cancer_vs_Healthy))

# Rank all proteins by training t-statistic (descending)
all_ranked <- dep_df %>%
  arrange(desc(t_Training_CR)) %>%
  mutate(rank = row_number())
n_all <- nrow(all_ranked)

# Cancer signature sets
cancer_sig <- dep_df %>% filter(pi_score_Cancer_vs_Healthy < 0.05)
up_genes <- cancer_sig$gene[cancer_sig$logFC_Cancer_vs_Healthy > 0]
dn_genes <- cancer_sig$gene[cancer_sig$logFC_Cancer_vs_Healthy < 0]
n_up <- length(up_genes)
n_dn <- length(dn_genes)

up_ranks   <- all_ranked$rank[all_ranked$gene %in% up_genes]
dn_ranks   <- all_ranked$rank[all_ranked$gene %in% dn_genes]
up_weights <- abs(all_ranked$t_Cancer_vs_Healthy[up_ranks])
dn_weights <- abs(all_ranked$t_Cancer_vs_Healthy[dn_ranks])

message(sprintf("  Cancer signature: %d up, %d down", n_up, n_dn))

# ── Weighted KS enrichment score ─────────────────────────────────────────────
compute_es <- function(ranks_in_set, n_total, weights = NULL) {
  n_set <- length(ranks_in_set)
  if (n_set == 0) return(list(es = 0, running = rep(0, n_total)))
  hit_pos <- sort(ranks_in_set)
  if (is.null(weights)) weights <- rep(1, n_set)
  weights <- abs(weights[order(match(ranks_in_set, hit_pos))])

  hit_indicator <- rep(FALSE, n_total)
  hit_indicator[hit_pos] <- TRUE
  hit_weights <- numeric(n_total)
  hit_weights[hit_pos] <- weights
  hit_cum  <- cumsum(hit_weights * hit_indicator) / sum(weights)
  miss_cum <- cumsum(!hit_indicator) / (n_total - n_set)
  running  <- hit_cum - miss_cum
  es <- if (abs(max(running)) > abs(min(running))) max(running) else min(running)
  list(es = es, running = running)
}

compute_connectivity <- function(up_r, dn_r, n_total, up_w = NULL, dn_w = NULL) {
  es_u <- compute_es(up_r, n_total, up_w)$es
  es_d <- compute_es(dn_r, n_total, dn_w)$es
  if (sign(es_u) == sign(es_d)) return(0)
  (es_u - es_d) / 2
}

es_up_obj <- compute_es(up_ranks, n_all, up_weights)
es_dn_obj <- compute_es(dn_ranks, n_all, dn_weights)

if (sign(es_up_obj$es) == sign(es_dn_obj$es)) {
  connectivity <- 0; coherent <- FALSE
} else {
  connectivity <- (es_up_obj$es - es_dn_obj$es) / 2; coherent <- TRUE
}

message(sprintf("  ES(up) = %.3f, ES(down) = %.3f, Connectivity = %.3f (%s)",
                es_up_obj$es, es_dn_obj$es, connectivity,
                ifelse(connectivity < 0, "REVERSAL", "EXACERBATION")))

# ── Permutation test (10,000×) ──────────────────────────────────────────────
set.seed(42)
B <- 10000
null_conn <- numeric(B)
for (i in seq_len(B)) {
  perm_up <- sample(n_all, n_up)
  perm_dn <- sample(n_all, n_dn)
  null_conn[i] <- compute_connectivity(perm_up, perm_dn, n_all)
}

# One-sided: is observed connectivity more negative than null?
perm_p <- (sum(null_conn <= connectivity) + 1) / (B + 1)
message(sprintf("  Permutation p = %.5f (n_perm = %d)", perm_p, B))

# ── Export CSV ───────────────────────────────────────────────────────────────
cmap_summary <- tibble(
  ES_up = es_up_obj$es, ES_down = es_dn_obj$es,
  connectivity = connectivity, coherent = coherent,
  n_up = n_up, n_dn = n_dn, n_total = n_all,
  perm_p = perm_p, n_perm = B
)
write_csv(cmap_summary, file.path(DAT, "SUPP_cmap_connectivity.csv"))
write_csv(tibble(replicate = seq_len(B), null_connectivity = null_conn),
          file.path(DAT, "SUPP_cmap_null_dist.csv"))

# ── Visualization ────────────────────────────────────────────────────────────
# Left: running enrichment curves
run_df <- tibble(
  rank = rep(1:n_all, 2),
  ES   = c(es_up_obj$running, es_dn_obj$running),
  Set  = rep(c("Cancer Up", "Cancer Down"), each = n_all)
)

marks_up <- tibble(rank = up_ranks)
marks_dn <- tibble(rank = dn_ranks)

p_curves <- ggplot(run_df, aes(x = rank, y = ES, color = Set)) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60",
             linewidth = 0.3) +
  geom_rug(data = marks_up, aes(x = rank), sides = "b",
           color = "#E57373", alpha = 0.2, length = unit(1.5, "mm"),
           inherit.aes = FALSE) +
  geom_rug(data = marks_dn, aes(x = rank), sides = "t",
           color = "#64B5F6", alpha = 0.2, length = unit(1.5, "mm"),
           inherit.aes = FALSE) +
  scale_color_manual(values = c("Cancer Up" = "#E57373", "Cancer Down" = "#64B5F6"),
                     name = NULL) +
  annotate("text", x = n_all * 0.02, y = min(run_df$ES) * 0.85,
           label = sprintf("ES(up) = %.3f\nES(dn) = %.3f\nCS = %.3f",
                            es_up_obj$es, es_dn_obj$es, connectivity),
           hjust = 0, size = 3, fontface = "bold", color = "grey25") +
  labs(x = sprintf("Rank by Training CR t-stat (n = %d)", n_all),
       y = "Enrichment Score") +
  FIG_THEME +
  theme(legend.position = c(0.82, 0.15),
        legend.background = element_rect(fill = alpha("white", 0.85),
                                          linewidth = 0),
        legend.key.size = unit(3, "mm"),
        legend.text = element_text(size = 7))

# Right: null distribution
p_null <- ggplot(tibble(x = null_conn), aes(x = x)) +
  geom_histogram(bins = 60, fill = "grey70", color = "grey40", linewidth = 0.2) +
  geom_vline(xintercept = connectivity, color = "#2563EB", linewidth = 0.8) +
  annotate("text", x = connectivity, y = Inf,
           label = sprintf("Observed = %.3f\nperm p %s", connectivity, fmt_p(perm_p)),
           hjust = -0.1, vjust = 1.5, size = 3, fontface = "bold", color = "#2563EB") +
  labs(x = "Connectivity Score (null)", y = "Count") +
  FIG_THEME

# Combine
pS_cmap <- p_curves + p_null +
  plot_layout(widths = c(1.3, 1)) +
  plot_annotation(
    title = "CMap Connectivity Score",
    subtitle = sprintf("Cancer sig. (%d up, %d dn) on Training CR rank | CS = %.3f (%s) | perm p %s",
                        n_up, n_dn, connectivity,
                        ifelse(connectivity < 0, "reversal", "exacerbation"),
                        fmt_p(perm_p)),
    theme = theme(plot.title = element_text(size = FIG_TITLE_SIZE, face = "bold"),
                  plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, color = "grey30")))

ggsave(file.path(RPT_PNG, "SUPP_cmap_connectivity.png"), pS_cmap,
       width = 260, height = 100, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "SUPP_cmap_connectivity.pdf"), pS_cmap,
       width = 260, height = 100, units = "mm", device = pdf_device)

message("Done: SUPP_cmap_connectivity")
