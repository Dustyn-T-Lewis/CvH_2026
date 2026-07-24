# SUPP GO Slim Bars
# Stacked horizontal bars: GO Slim consolidated category distribution by reversal quadrant
setwd(here::here())
source("04_Figures/shared/style.R")
source("04_Figures/shared/go_slim_categories.R")

pacman::p_load(tidyverse)

RPT_PNG <- "04_Figures/F04_Reversal/b_reports/supp/png/panels"
RPT_PDF <- "04_Figures/F04_Reversal/b_reports/supp/pdf/panels"
DAT     <- "04_Figures/F04_Reversal/c_data/panel_supp"

dir.create(RPT_PNG, recursive = TRUE, showWarnings = FALSE)
dir.create(RPT_PDF, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT,     recursive = TRUE, showWarnings = FALSE)

pdf_device <- get_pdf_device()

# -- Load data -----------------------------------------------------------------
source("04_Figures/F04_Reversal/a_script/f04_data.R")
dep <- dep_df

# -- Filter Pi-significant cancer proteins & classify --------------------------
sig <- dep %>%
  filter(pi_score_CRvH_Baseline < 0.05) %>%
  mutate(quadrant = case_when(
    logFC_CRvH_Baseline > 0 & logFC_CR_Training < 0 ~ "Reversed Up",
    logFC_CRvH_Baseline < 0 & logFC_CR_Training > 0 ~ "Reversed Down",
    TRUE ~ "Non-reversed"
  ))

# -- Assign GO Slim categories -------------------------------------------------
all_genes <- unique(dep$gene)
fg_genes  <- unique(sig$gene)

slim_map <- assign_go_slim_consolidated(fg_genes, all_genes)

sig_slim <- sig %>%
  left_join(slim_map, by = "gene") %>%
  filter(!is.na(consolidated))

# -- Count per quadrant x category ---------------------------------------------
count_df <- sig_slim %>%
  count(quadrant, consolidated, name = "n") %>%
  group_by(quadrant) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup() %>%
  complete(quadrant, consolidated, fill = list(n = 0L, pct = 0))

count_df$quadrant <- factor(count_df$quadrant,
                            levels = c("Reversed Up", "Reversed Down",
                                       "Non-reversed"))

write.csv(count_df, file.path(DAT, "SUPP_goslim_distribution.csv"),
          row.names = FALSE)

# -- Plot ----------------------------------------------------------------------
pS_goslim <- ggplot(count_df, aes(x = n, y = consolidated, fill = consolidated)) +
  geom_col(show.legend = FALSE) +
  facet_wrap(~ quadrant, ncol = 1, scales = "free_x") +
  scale_fill_manual(values = CONSOLIDATED_COLORS, drop = FALSE) +
  scale_y_discrete(limits = rev(CONSOLIDATED_PATHWAY_ORDER)) +
  labs(x = "Number of proteins", y = NULL,
       title = "GO Slim category distribution",
       subtitle = "Pi-significant cancer proteins by reversal quadrant") +
  FIG_THEME +
  theme(strip.text = element_text(face = "bold", size = FIG_STRIP_SIZE),
        axis.text.y = element_text(size = 7))

ggsave(file.path(RPT_PNG, "SUPP_goslim_bars.png"), pS_goslim,
       width = 170, height = 140, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(RPT_PDF, "SUPP_goslim_bars.pdf"), pS_goslim,
       width = 170, height = 140, units = "mm", device = pdf_device)

message("Done: SUPP_goslim_bars  [", n_distinct(sig_slim$gene), " genes, ",
        n_distinct(count_df$consolidated), " categories]")
