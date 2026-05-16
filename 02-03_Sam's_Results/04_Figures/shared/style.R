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

# ── Protein classification helpers (CvH / Sam versions) ─────────────────────
# F04 (concordance): classifies by pi_CRE + pi_PLA + optional pi_interaction.
# Returns factor with levels matching key_lvls used in F04 _panel_A_ORA.R.
classify_proteins_f2 <- function(pi_Y, pi_O, pi_int, threshold = 0.05) {
  # pi_int may be NULL or all NA — treat as absent interaction contrast.
  # Construct a safe logical vector for the interaction condition.
  n <- length(pi_Y)
  int_sig <- if (!is.null(pi_int) && length(pi_int) == n) {
    !is.na(pi_int) & pi_int < threshold
  } else {
    rep(FALSE, n)
  }
  dplyr::case_when(
    int_sig                                      ~ "Interaction",
    pi_Y < threshold & pi_O < threshold          ~ "Sig Both",
    pi_Y < threshold                             ~ "Sig CRE only",
    pi_O < threshold                             ~ "Sig PLA only",
    TRUE                                         ~ "NS"
  ) |>
    factor(levels = c("Interaction", "Sig Both", "Sig CRE only", "Sig PLA only", "NS"))
}

# F05 (reversal): classifies by pi_CvH + pi_Training_CR.
# Returns factor with levels matching key_lvls used in F05 _panel_A_ORA.R.
classify_proteins_f3 <- function(pi_aging, pi_training_old, threshold = 0.05) {
  dplyr::case_when(
    pi_aging < threshold & pi_training_old < threshold ~ "Sig Both",
    pi_aging < threshold                               ~ "Sig CvH only",
    pi_training_old < threshold                        ~ "Sig CR only",
    TRUE                                               ~ "NS"
  ) |>
    factor(levels = c("Sig Both", "Sig CvH only", "Sig CR only", "NS"))
}

# Fisher Z CI for Pearson r (Bonett & Wright 2000)
fisher_z_ci <- function(r, n, k = 0, level = 0.95) {
  n_eff <- n - k
  if (n_eff < 4 || is.na(r)) return(c(lo = NA_real_, hi = NA_real_))
  z    <- atanh(r)
  se   <- 1 / sqrt(n_eff - 3)
  crit <- qnorm(1 - (1 - level) / 2)
  c(lo = tanh(z - crit * se), hi = tanh(z + crit * se))
}

sig_stars <- function(padj) {
  dplyr::case_when(
    padj < 0.001 ~ "***",
    padj < 0.01  ~ "**",
    padj < 0.05  ~ "*",
    TRUE         ~ ""
  )
}

# ── Color constants for Sam's F04/F05 (CvH naming) ──────────────────────────
# F04: Training concordance (CRE vs PLA)
SIG_COLORS_F2 <- c(
  "Interaction"  = "#7B5EA7",
  "Sig Both"     = "#2E7D32",
  "Sig CRE only" = "#E05A4E",
  "Sig PLA only" = "#5DA5DA",
  "NS"           = "grey70"
)
SIG_LABEL_FILL_F2 <- c(
  "Interaction"  = scales::alpha("#7B5EA7", 0.75),
  "Sig Both"     = scales::alpha("#2E7D32", 0.75),
  "Sig CRE only" = scales::alpha("#E05A4E", 0.75),
  "Sig PLA only" = scales::alpha("#5DA5DA", 0.75),
  "NS"           = scales::alpha("grey70",  0.75)
)
SIG_LABEL_TEXT_F2 <- setNames(rep("white", 5), names(SIG_LABEL_FILL_F2))

ORA_QUAD_COLORS_F2 <- c(
  "Concordant Up"                 = "#E57373",
  "Concordant Down"               = "#64B5F6",
  "Discordant (CRE Up / PLA Down)" = "#FFB74D",
  "Discordant (CRE Down / PLA Up)" = "#81C784"
)

# F05: Recovery reversal (Cancer_vs_Healthy vs Training_CR)
SIG_COLORS_F3 <- c(
  "Sig Both"    = "#2E7D32",
  "Sig CvH only" = "#E05A4E",
  "Sig CR only" = "#5DA5DA",
  "NS"          = "grey70"
)
SIG_LABEL_FILL_F3 <- c(
  "Sig Both"    = scales::alpha("#2E7D32", 0.75),
  "Sig CvH only" = scales::alpha("#E05A4E", 0.75),
  "Sig CR only" = scales::alpha("#5DA5DA", 0.75),
  "NS"          = scales::alpha("grey70",  0.75)
)
SIG_LABEL_TEXT_F3 <- setNames(rep("white", 4), names(SIG_LABEL_FILL_F3))

ORA_QUAD_COLORS_F3 <- c(
  "Reversed (CvH Up / CR Down)"  = "#E57373",
  "Reversed (CvH Down / CR Up)"  = "#64B5F6",
  "Exacerbated Up"               = "#FFB74D",
  "Exacerbated Down"             = "#81C784"
)

# AGE_COLORS: used by panel_D_nes_scatter.R for quadrant background shading.
# "Old" = warm/red = concordant direction; "Young" = cool/blue = discordant.
AGE_COLORS <- c(Young = "#4393C3", Old = "#D6604D")
