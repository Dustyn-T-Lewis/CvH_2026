# Reversal Panel G: enriched volcano of the residual contrast (Resid = CR_post - H)
# "What's still off after training" -- the per-protein residual gap from healthy,
# with fgsea pathway enrichment on the residual t-statistic alongside. Pairs with
# the trajectory clusters (F, "what reversed").
setwd(here::here())
source("04_Figures/shared/style.R")
source("04_Figures/shared/pathway_utils.R")
pacman::p_load(tidyverse, fgsea, ggrepel, patchwork)

RPT_PNG <- "04_Figures/F04_Reversal/b_reports/main/png/panels"
RPT_PDF <- "04_Figures/F04_Reversal/b_reports/main/pdf/panels"
DAT     <- "04_Figures/F04_Reversal/c_data/panel_G"
dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)
pdf_device <- get_pdf_device()
set.seed(42)

source("04_Figures/F04_Reversal/a_script/f04_data.R")   # dep_df

# ── Volcano frame: residual contrast ─────────────────────────────────────────
v <- dep_df |>
  transmute(gene,
            logFC = logFC_Resid,
            P     = pmax(P.Value_Resid, 1e-300),
            pi    = pi_score_Resid) |>
  filter(!is.na(logFC), !is.na(P)) |>
  mutate(neglogP = -log10(P),
         sig = pi < 0.05,
         dir = case_when(!sig ~ "NS", logFC > 0 ~ "Up", TRUE ~ "Down"))
n_sig <- sum(v$sig)
write_csv(v, file.path(DAT, "resid_volcano.csv"))

lab <- v |> filter(sig) |> arrange(pi) |> slice_head(n = 18)

p_volc <- ggplot(v, aes(logFC, neglogP)) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_point(data = filter(v, !sig), color = "grey80", size = 0.5, alpha = 0.4) +
  geom_point(data = filter(v, sig), aes(color = dir), size = 0.9, alpha = 0.8) +
  geom_text_repel(data = lab, aes(label = gene), size = 2.4, fontface = "italic",
                  max.overlaps = 20, segment.size = 0.2, segment.color = "grey50",
                  min.segment.length = 0, box.padding = 0.3, seed = 42) +
  scale_color_manual(values = c(Up = unname(DIR_COLORS["Up"]), Down = unname(DIR_COLORS["Down"])),
                     name = NULL, labels = c(Down = "lower than healthy", Up = "higher than healthy")) +
  labs(x = expression(log[2]*FC~"(CR post − Healthy)"), y = expression(-log[10]~italic(P))) +
  FIG_THEME +
  theme(legend.position = "top", panel.grid.minor = element_blank())

# ── Pathway enrichment on the residual t-statistic (fgsea) ───────────────────
rr <- dep_df |> select(gene, t_Resid) |> filter(is.finite(t_Resid)) |> distinct(gene, .keep_all = TRUE)
ranks <- setNames(rr$t_Resid, rr$gene)
pw <- build_pathway_collection(min_size = 15, max_size = 500,
                               include_goslim = FALSE, exclude_variants = TRUE)
gs <- fgsea(pathways = pw, stats = sort(ranks), eps = 0)

top <- gs |> as_tibble() |> filter(padj < 0.10) |>
  mutate(pathway_label = clean_pathway_name(pathway),
         side = ifelse(NES > 0, "Higher than healthy", "Lower than healthy")) |>
  group_by(side) |> arrange(padj) |> slice_head(n = 8) |> ungroup() |>
  arrange(NES) |> mutate(pathway_label = factor(pathway_label, levels = pathway_label))
write_csv(select(top, pathway, pathway_label, NES, pval, padj, side), file.path(DAT, "resid_fgsea_top.csv"))

p_path <- ggplot(top, aes(NES, pathway_label, color = NES > 0)) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = NES, yend = pathway_label), linewidth = 0.5) +
  geom_point(aes(size = -log10(padj))) +
  scale_color_manual(values = c(`TRUE` = unname(DIR_COLORS["Up"]), `FALSE` = unname(DIR_COLORS["Down"])),
                     guide = "none") +
  scale_size_continuous(name = expression(-log[10]~FDR), range = c(1.5, 4)) +
  labs(x = "NES (residual)", y = NULL,
       title = "Pathways still off after training") +
  FIG_THEME +
  theme(plot.title = element_text(size = 9, face = "bold"),
        axis.text.y = element_text(size = 7), legend.position = "right",
        panel.grid.minor = element_blank())

composite <- (p_volc | p_path) +
  plot_layout(widths = c(1, 1.15)) +
  plot_annotation(
    title = "Cancer Recovery Reversal: residual contrast (what training did NOT fix)",
    subtitle = sprintf("Resid = CR_post − Healthy | %d proteins still differ (Π < 0.05) | fgsea FDR < 0.10",
                       n_sig),
    theme = theme(plot.title = element_text(size = FIG_TITLE_SIZE, face = "bold"),
                  plot.subtitle = element_text(size = FIG_SUBTITLE_SIZE, color = "grey30")))

COMP_W <- 220; COMP_H <- 120
ggsave(file.path(RPT_PNG, "MAIN_panel_G_resid_volcano_composite.png"), composite,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)
ggsave(file.path(RPT_PDF, "MAIN_panel_G_resid_volcano_composite.pdf"), composite,
       width = COMP_W, height = COMP_H, units = "mm", device = pdf_device)
message(sprintf("Reversal Panel G residual volcano done | %d residual DEPs", n_sig))
