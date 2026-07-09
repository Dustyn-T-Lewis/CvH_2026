# Reversal Panel F: Trajectory clustering of the disease signature
# Soft (fuzzy c-means) clustering of each disease-signature protein's standardized
# abundance profile across the ordered axis H_pre -> CR_pre -> CR_post, then a
# reversal-class label per cluster + per-cluster ORA. Pattern-first complement to
# the per-contrast scatter (A) and NES scatter (B).
#   Engine = e1071::cmeans (the fuzzy c-means Mfuzz wraps; Futschik & Carlisle 2005).
#   Descriptive only -- inference stays with fry / RRHO2 / permutation null.
setwd(here::here())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")
pacman::p_load(tidyverse, e1071, ggrepel, patchwork)

RPT_PNG <- "04_Figures/Reversal/b_reports/main/png/panels"
RPT_PDF <- "04_Figures/Reversal/b_reports/main/pdf/panels"
DAT     <- "04_Figures/Reversal/c_data/panel_F"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()
set.seed(42)

ORDER   <- c("H_pre", "CR_pre", "CR_post")          # healthy -> diseased -> trained
AXISLAB <- c(H_pre = "Healthy", CR_pre = "CR pre", CR_post = "CR post")
K       <- 6                                         # soft clusters
M_FUZZ  <- 1.25                                      # fuzzifier (Mfuzz-style)
PHI_BAND <- 0.25                                     # engine band (REVERSAL_PHI_BAND)
CLASS_COLORS <- c(Normalized = "#1B7837", Persistent = "#878787",
                  Exacerbated = "#D6604D")

# ── Data: group-mean trajectories of the disease signature ───────────────────
source("04_Figures/Reversal/a_script/reversal_inputs.R")   # dep_df, dal
stopifnot(all(ORDER %in% dal$metadata$group_time))

gmean <- sapply(ORDER, function(g)
  rowMeans(dal$data[, dal$metadata$Col_ID[dal$metadata$group_time == g], drop = FALSE]))
rownames(gmean) <- rownames(dal$data)                      # uniprot_id

# disease signature = Pi < 0.05 on Cancer_vs_Healthy. Per-protein reversal class
# from the engine band on the actual D/T logFC (phi = -T/D), so the colour is the
# rigorous reversal call (consistent with fry / permutation null), not a cluster guess.
sig <- dep_df |>
  filter(pi_score_Cancer_vs_Healthy < 0.05) |>
  mutate(phi = -logFC_Training_CR / logFC_Cancer_vs_Healthy,
         reversal_class = case_when(phi >=  PHI_BAND ~ "Normalized",
                                    phi <= -PHI_BAND ~ "Exacerbated",
                                    TRUE             ~ "Persistent")) |>
  distinct(uniprot_id, gene, phi, reversal_class)
gmean <- gmean[rownames(gmean) %in% sig$uniprot_id, , drop = FALSE]

# z-score each protein's 3-point profile (drop constant rows)
z <- t(apply(gmean, 1, function(r) (r - mean(r)) / sd(r)))
z <- z[is.finite(rowSums(z)), , drop = FALSE]
message(sprintf("  Trajectory clustering: %d disease-signature proteins", nrow(z)))

# ── Fuzzy c-means soft clustering ────────────────────────────────────────────
fcm <- cmeans(z, centers = K, m = M_FUZZ, iter.max = 200)
memb_max <- apply(fcm$membership, 1, max)
assign <- tibble(uniprot_id = rownames(z),
                 cluster = fcm$cluster, membership = memb_max) |>
  left_join(sig, by = "uniprot_id")             # per-protein phi + reversal_class
write_csv(assign, file.path(DAT, "trajectory_clusters.csv"))

# centroid shape per cluster -> dominant-pattern label for the facet strip only
cent <- as_tibble(fcm$centers, .name_repair = "minimal")
names(cent) <- ORDER
cent <- cent |>
  mutate(cluster = row_number(),
         d_pre = CR_pre - H_pre,
         phi_c = -(CR_post - CR_pre) / d_pre,
         pattern = case_when(abs(d_pre) < 0.4   ~ "training-emergent",
                             phi_c >=  PHI_BAND ~ "normalizing",
                             phi_c <= -PHI_BAND ~ "exacerbating",
                             TRUE               ~ "persistent"))

clab <- cent |>
  mutate(n = as.integer(table(factor(assign$cluster, levels = seq_len(K)))[cluster]),
         lab = sprintf("Cluster %d · %s (n=%d)", cluster, pattern, n))
class_n <- assign |> count(reversal_class)
rev_dir <- mean(assign$phi > 0, na.rm = TRUE)   # directional reverse fraction (phi > 0)

# ── Per-cluster ORA (top pathway per cluster) ────────────────────────────────
universe <- unique(dep_df$gene)
pw <- build_pathway_collection(min_size = 15, max_size = 500,
                               include_goslim = FALSE, exclude_variants = TRUE)
top_path <- map_dfr(seq_len(K), function(k) {
  g <- assign$gene[assign$cluster == k]
  if (length(g) < 5) return(tibble(cluster = k, pathway_label = NA_character_, padj = NA_real_))
  r <- tryCatch(run_ora_deduplicated(genes = g, universe = universe, pathways = pw,
                                     jaccard_cutoff = 0.5, min_size = 15, max_size = 500,
                                     padj_cutoff = 1), error = function(e) tibble())
  if (nrow(r) == 0) return(tibble(cluster = k, pathway_label = NA_character_, padj = NA_real_))
  r |> arrange(padj) |> slice_head(n = 1) |>
    transmute(cluster = k, pathway_label = clean_pathway_name(pathway), padj)
})
if (nrow(top_path)) write_csv(top_path, file.path(DAT, "cluster_top_pathway.csv"))
clab <- clab |> left_join(top_path, by = "cluster") |>
  mutate(path_lab = ifelse(is.na(pathway_label), "",
                           sprintf("%s (FDR %s)", pathway_label, fmt_p(padj))))

# ── Long frame for trajectory plotting ───────────────────────────────────────
long <- as_tibble(z, rownames = "uniprot_id") |>
  pivot_longer(all_of(ORDER), names_to = "cond", values_to = "zval") |>
  left_join(select(assign, uniprot_id, cluster, membership, reversal_class), by = "uniprot_id") |>
  mutate(cond = factor(cond, levels = ORDER))
cent_long <- cent |>
  select(cluster, all_of(ORDER)) |>
  pivot_longer(all_of(ORDER), names_to = "cond", values_to = "zval") |>
  mutate(cond = factor(cond, levels = ORDER))
facet_lab <- setNames(clab$lab, clab$cluster)

p_traj <- ggplot(long, aes(cond, zval, group = uniprot_id)) +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
  geom_line(aes(color = reversal_class, alpha = membership), linewidth = 0.25) +
  geom_line(data = cent_long, aes(cond, zval, group = cluster),
            inherit.aes = FALSE, color = "black", linewidth = 1.1) +
  geom_text(data = clab, aes(x = 2, y = Inf, label = path_lab), inherit.aes = FALSE,
            vjust = 1.4, size = 2.4, fontface = "italic", color = "grey25", lineheight = 0.9) +
  facet_wrap(~ cluster, ncol = 3, labeller = labeller(cluster = facet_lab)) +
  scale_color_manual(values = CLASS_COLORS, name = "Reversal class") +
  scale_alpha_continuous(range = c(0.05, 0.5), guide = "none") +
  scale_x_discrete(labels = AXISLAB) +
  labs(x = NULL, y = "standardized abundance (z)") +
  FIG_THEME +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 7.5, face = "bold"))

composite <- p_traj +
  plot_annotation(
    title = "Cancer Recovery Reversal: trajectory clustering of the disease signature",
    subtitle = sprintf("Fuzzy c-means k=%d | %d disease DEPs | %.1f%% reverse | Normalized %d · Persistent %d · Exacerbated %d",
                       K, nrow(z), 100 * rev_dir,
                       sum(class_n$n[class_n$reversal_class == "Normalized"]),
                       sum(class_n$n[class_n$reversal_class == "Persistent"]),
                       sum(class_n$n[class_n$reversal_class == "Exacerbated"])),
    theme = theme(plot.title = element_text(size = FIG_TITLE_SIZE, face = "bold"),
                  plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, color = "grey30")))

COMP_W <- 200; COMP_H <- 150
ggsave(file.path(RPT_PNG, "MAIN_panel_F_trajectory_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "MAIN_panel_F_trajectory_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
message("Reversal Panel F trajectory composite done")
