# F03 DEP Overview Keys — Contrast colors, direction, threshold legend
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/keys/extract_keys.R")

# --- Contrast colors (all 6) ---
ctr_df <- tibble::tibble(
  contrast = names(CONTRAST_COLORS),
  y        = seq_along(CONTRAST_COLORS)
)
p_ctr <- ggplot(ctr_df, aes(x = 1, y = y, fill = contrast)) +
  geom_tile() +
  scale_fill_manual(values = CONTRAST_COLORS, name = "Contrast") +
  FIG_THEME + theme(legend.position = "bottom")
k_ctr <- wrap_key(extract_key(p_ctr))

# --- Direction legend ---
dir_df <- tibble::tibble(
  direction = names(DIR_COLORS),
  y         = seq_along(DIR_COLORS)
)
p_dir <- ggplot(dir_df, aes(x = 1, y = y, fill = direction)) +
  geom_tile() +
  scale_fill_manual(values = DIR_COLORS, name = "Direction") +
  FIG_THEME + theme(legend.position = "bottom")
k_dir <- wrap_key(extract_key(p_dir))

# --- Threshold legend (text) ---
thresh_df <- tibble::tibble(
  threshold = c("p < 0.05", "FDR < 0.10", "Pi < 0.05"),
  y = 1:3
)
p_thresh <- ggplot(thresh_df, aes(x = 1, y = y, label = threshold)) +
  geom_text(size = 3) +
  labs(title = "Significance\nThresholds") +
  theme_void() + theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
k_thresh <- wrap_elements(p_thresh)

composite <- k_ctr | k_dir | k_thresh
save_key(composite, "F3_keys", width = 300, height = 60)
