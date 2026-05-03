# F04 Concordance Keys — Significance, direction, database, quadrant legends
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/keys/extract_keys.R")

# --- CRvH significance colors ---
sig4_df <- tibble::tibble(
  sig_class = names(SIG_COLORS_F4),
  y         = seq_along(SIG_COLORS_F4)
)
p_sig4 <- ggplot(sig4_df, aes(x = 1, y = y, fill = sig_class)) +
  geom_tile() +
  scale_fill_manual(values = SIG_COLORS_F4, name = "CRvH Significance") +
  FIG_THEME + theme(legend.position = "bottom")
k_sig4 <- wrap_key(extract_key(p_sig4))

# --- CR significance colors ---
sig3_df <- tibble::tibble(
  sig_class = names(SIG_COLORS_F3),
  y         = seq_along(SIG_COLORS_F3)
)
p_sig3 <- ggplot(sig3_df, aes(x = 1, y = y, fill = sig_class)) +
  geom_tile() +
  scale_fill_manual(values = SIG_COLORS_F3, name = "CR Significance") +
  FIG_THEME + theme(legend.position = "bottom")
k_sig3 <- wrap_key(extract_key(p_sig3))

# --- Database colors ---
db_df <- tibble::tibble(
  database = names(DB_COLORS),
  y        = seq_along(DB_COLORS)
)
p_db <- ggplot(db_df, aes(x = 1, y = y, fill = database)) +
  geom_tile() +
  scale_fill_manual(values = DB_COLORS, name = "Database") +
  FIG_THEME + theme(legend.position = "bottom")
k_db <- wrap_key(extract_key(p_db))

composite <- k_sig4 | k_sig3 | k_db
save_key(composite, "F4_keys", width = 350, height = 70)
