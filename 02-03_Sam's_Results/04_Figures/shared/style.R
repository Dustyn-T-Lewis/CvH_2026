# Sam's-Results figure style.
# Sources the canonical CvH pipeline style so palettes, themes, and sizing
# constants stay aligned with the main 04_Figures/ output.
#
# Overrides below bring Sam figures into alignment with YvO style constants
# (J Physiol double-column spec: base_size = 6, Helvetica, tighter annotation sizes).
# The main CvH 04_Figures/ pipeline keeps base_size = 10 for its own figures;
# only Sam's sub-pipeline adopts the YvO sizing.

setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

# ── YvO-aligned size overrides ───────────────────────────────────────────────
# Sizing constants (match YvO 04_Figures/shared/style.R exactly)
PANEL_MD      <- 178    # J Physiol double-column width (was 180)
BASE_PATHWAY  <- 2.8    # ~8pt pathway labels            (was 4.0)
BASE_GENE     <- 2.5    # ~7pt gene labels               (was 3.2)
BASE_STAT     <- 2.5    # ~7pt stat annotations          (was 3.5)
BASE_QUADRANT <- 2.8    # ~8pt quadrant labels           (was 4.0)
BASE_COUNT    <- 2.5    # ~7pt bar counts                (was 3.5)
BASE_TAG      <- 8      # panel tag pt                   (was 18)

# Font size hierarchy (J Physiol spec)
FIG_TITLE_SIZE    <- 7   # (was 12)
FIG_SUBTITLE_SIZE <- 4   # (was 9)
FIG_STRIP_SIZE    <- 5   # (was 10)
FIG_AXIS_TEXT     <- 5   # (was 8.5)
FIG_LEGEND_TITLE  <- 5   # (was 9.5)
FIG_LEGEND_TEXT   <- 4   # (was 8.5)

# Rebuild FIG_THEME at base_size = 6 with Helvetica (YvO spec)
FIG_THEME <- theme_bw(base_size = 6, base_family = "Helvetica") +
  theme(
    plot.title         = element_text(face = "bold", size = FIG_TITLE_SIZE,
                                      margin = margin(b = 1)),
    plot.subtitle      = element_text(face = "bold.italic", size = FIG_SUBTITLE_SIZE,
                                      colour = "grey30", margin = margin(t = 0, b = 2)),
    plot.tag           = element_text(face = "bold", size = BASE_TAG),
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = FIG_STRIP_SIZE),
    axis.title.x       = element_text(face = "bold", size = 5,
                                      margin = margin(t = 0)),
    axis.title.y       = element_text(face = "bold", size = 5,
                                      margin = margin(r = -1)),
    axis.text          = element_text(size = FIG_AXIS_TEXT, color = "grey15"),
    legend.title       = element_text(face = "bold", size = FIG_LEGEND_TITLE,
                                      color = "grey20"),
    legend.text        = element_text(size = FIG_LEGEND_TEXT, color = "grey15"),
    legend.key.size    = unit(2.5, "mm"),
    panel.grid.minor   = element_blank()
  )

# composite_text_sizes() helper (same formula as YvO)
composite_text_sizes <- function(comp_h_mm) {
  list(
    title    = pmax(6, pmin(8, round(5 + comp_h_mm / 80))),
    subtitle = pmax(4, pmin(6, round(3 + comp_h_mm / 100))),
    tag      = 8
  )
}

# strip_for_composite() helper (same as YvO)
strip_for_composite <- function(p) {
  p + labs(title = NULL, subtitle = NULL, tag = NULL) +
    theme(legend.position = "none")
}
