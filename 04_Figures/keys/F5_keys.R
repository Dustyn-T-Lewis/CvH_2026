# F05 Reversal Keys — Pattern, quadrant, cancer direction, significance legends
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/keys/extract_keys.R")

# --- Reversal pattern colors ---
pat_df <- tibble::tibble(
  pattern = names(PATTERN_COLS_F4),
  y       = seq_along(PATTERN_COLS_F4)
)
p_pat <- ggplot(pat_df, aes(x = 1, y = y, fill = pattern)) +
  geom_tile() +
  scale_fill_manual(values = PATTERN_COLS_F4, name = "Reversal Pattern") +
  FIG_THEME + theme(legend.position = "bottom")
k_pat <- wrap_key(extract_key(p_pat))

# --- Reversal quadrant colors ---
quad_df <- tibble::tibble(
  quadrant = names(ORA_QUAD_COLORS_F4),
  y        = seq_along(ORA_QUAD_COLORS_F4)
)
p_quad <- ggplot(quad_df, aes(x = 1, y = y, fill = quadrant)) +
  geom_tile() +
  scale_fill_manual(values = ORA_QUAD_COLORS_F4, name = "Quadrant") +
  FIG_THEME + theme(legend.position = "bottom")
k_quad <- wrap_key(extract_key(p_quad))

# --- Cancer direction colors ---
cdir_df <- tibble::tibble(
  direction = names(CANCER_DIR_COLORS),
  y         = seq_along(CANCER_DIR_COLORS)
)
p_cdir <- ggplot(cdir_df, aes(x = 1, y = y, fill = direction)) +
  geom_tile() +
  scale_fill_manual(values = CANCER_DIR_COLORS, name = "Cancer Direction") +
  FIG_THEME + theme(legend.position = "bottom")
k_cdir <- wrap_key(extract_key(p_cdir))

# --- CRvH significance ---
sig_df <- tibble::tibble(
  sig_class = names(SIG_COLORS_F4),
  y         = seq_along(SIG_COLORS_F4)
)
p_sig <- ggplot(sig_df, aes(x = 1, y = y, fill = sig_class)) +
  geom_tile() +
  scale_fill_manual(values = SIG_COLORS_F4, name = "Significance") +
  FIG_THEME + theme(legend.position = "bottom")
k_sig <- wrap_key(extract_key(p_sig))

composite <- k_pat | k_quad | k_cdir | k_sig
save_key(composite, "F5_keys", width = 350, height = 70)
