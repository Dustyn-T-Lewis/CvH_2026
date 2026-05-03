# F02 QC Keys — Group fill, PCA shapes, contrast colors
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/keys/extract_keys.R")

# --- Group fill legend (CV violins) ---
group_df <- tibble::tibble(
  group = names(GROUP_FILL),
  fill  = unname(GROUP_FILL),
  y     = seq_along(GROUP_FILL)
)
p_group <- ggplot(group_df, aes(x = 1, y = y, fill = group)) +
  geom_tile() +
  scale_fill_manual(values = GROUP_FILL, name = "Group") +
  FIG_THEME + theme(legend.position = "bottom")
k_group <- wrap_key(extract_key(p_group))

# --- PCA color/shape legend ---
pca_df <- tibble::tibble(
  group = names(PCA_COLORS),
  color = unname(PCA_COLORS),
  shape = unname(PCA_SHAPES),
  y     = seq_along(PCA_COLORS)
)
p_pca <- ggplot(pca_df, aes(x = 1, y = y, color = group, shape = group)) +
  geom_point(size = 3) +
  scale_color_manual(values = PCA_COLORS, name = "Group") +
  scale_shape_manual(values = PCA_SHAPES, name = "Group") +
  FIG_THEME + theme(legend.position = "bottom")
k_pca <- wrap_key(extract_key(p_pca))

# --- Contrast colors legend ---
ctr_df <- tibble::tibble(
  contrast = names(CONTRAST_COLORS),
  color    = unname(CONTRAST_COLORS),
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
  fill      = unname(DIR_COLORS),
  y         = seq_along(DIR_COLORS)
)
p_dir <- ggplot(dir_df, aes(x = 1, y = y, fill = direction)) +
  geom_tile() +
  scale_fill_manual(values = DIR_COLORS, name = "Direction") +
  FIG_THEME + theme(legend.position = "bottom")
k_dir <- wrap_key(extract_key(p_dir))

# --- Compose and save ---
composite <- k_group | k_pca | k_ctr | k_dir
save_key(composite, "F2_keys", width = 350, height = 60)
