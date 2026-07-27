# Statistics the figure panels share. Style lives in style.R; this file holds
# the calculations, so they can be tested without loading ggplot.

# Fisher z-transform CI for a Pearson r. Undefined below n = 4, where the
# z-standard-error 1/sqrt(n - 3) blows up or goes imaginary.
fisher_z_ci <- function(r, n, level = 0.95) {
  z <- atanh(r)
  se <- 1 / sqrt(n - 3)
  q <- qnorm((1 + level) / 2)
  c(lo = tanh(z - q * se), hi = tanh(z + q * se))
}

# F04 protein classes: which of the two Pi-score axes clears threshold.
classify_proteins_f4 <- function(pi_CvH, pi_TR, threshold = 0.05) {
  dplyr::case_when(
    pi_CvH < threshold & pi_TR < threshold ~ "Sig Both",
    pi_CvH < threshold ~ "Sig Cancer only",
    pi_TR < threshold ~ "Sig Training only",
    TRUE ~ "NS"
  ) |>
    factor(levels = c("Sig Both", "Sig Cancer only", "Sig Training only", "NS"))
}
