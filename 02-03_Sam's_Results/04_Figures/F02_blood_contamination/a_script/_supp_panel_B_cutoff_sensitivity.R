# Supp B — B_M_ratio cutoff sensitivity for Cancer_vs_Healthy DEP count.
# Sweeps top-X% B_M_ratio sample drops and refits limma via proteoDA.
#
# Design: ~0 + cancer_time + sex (no random effect for sweep; repeated measures
# would require ≥2 observations per subject which breaks at high cutoffs).
# Contrast: Cancer_vs_Healthy = (SURV_T1 + SURV_T2)/2 - CTL_T1
# FDR threshold: 0.10 (BH correction via p.adjust).
#
# Outputs:
#   c_data/supp_B_cutoff_sensitivity.csv
#   b_reports/supp/png/panels/supp_B.png

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(limma)
  library(proteoDA)
})

setwd(rprojroot::find_rstudio_root_file())
source("02-03_Sam's_Results/04_Figures/build_data_index.R")
source("02-03_Sam's_Results/04_Figures/shared/style.R")

# ---------------------------------------------------------------------------
# Load Sam's DAList
# ---------------------------------------------------------------------------

sam   <- readRDS(sam_idx$sam$dalist_rds)
meta  <- as.data.frame(sam$metadata)

fdr_thresh <- 0.10

# ---------------------------------------------------------------------------
# Inspect cancer_time levels and build contrast string
# ---------------------------------------------------------------------------

ct_levels   <- sort(unique(as.character(meta$cancer_time)))
surv_levels <- grep("^SURV", ct_levels, value = TRUE)
ctl_levels  <- grep("^CTL",  ct_levels, value = TRUE)

stopifnot("Need at least one SURV level" = length(surv_levels) >= 1,
          "Need at least one CTL level"  = length(ctl_levels)  >= 1)

surv_expr <- if (length(surv_levels) == 1) {
  surv_levels
} else {
  sprintf("(%s)/%d", paste(surv_levels, collapse = " + "), length(surv_levels))
}
ctl_expr <- if (length(ctl_levels) == 1) {
  ctl_levels
} else {
  sprintf("(%s)/%d", paste(ctl_levels, collapse = " + "), length(ctl_levels))
}

contrast_str <- sprintf("Cancer_vs_Healthy = %s - %s", surv_expr, ctl_expr)
message("Supp B contrast: ", contrast_str)

# ---------------------------------------------------------------------------
# Helper: refit limma on a DAList subset and return DEP count
# ---------------------------------------------------------------------------

fit_for_subset <- function(dal_sub) {
  # Coerce cancer_time and sex to factor with no NA levels
  dal_sub$metadata$cancer_time <- droplevels(factor(dal_sub$metadata$cancer_time))
  dal_sub$metadata$sex         <- droplevels(factor(dal_sub$metadata$sex))

  dal_sub <- add_design(dal_sub, design_formula = ~0 + cancer_time + sex)
  dal_sub <- add_contrasts(dal_sub, contrasts_vector = contrast_str)
  dal_sub <- fit_limma_model(dal_sub)

  # Access results via limma::topTable (proteoDA::results is populated only
  # after extract_DA_results; we bypass that to keep the sweep lean)
  fit <- dal_sub$eBayes_fit
  tt  <- limma::topTable(fit, coef = "Cancer_vs_Healthy", number = Inf,
                         sort.by = "none")
  tt$FDR <- p.adjust(tt$P.Value, method = "BH")
  sum(tt$FDR < fdr_thresh, na.rm = TRUE)
}

# ---------------------------------------------------------------------------
# Sweep cutoffs: 0%, 5%, 10%, 15%, 20% of top-B_M_ratio samples dropped
# ---------------------------------------------------------------------------

cutoffs <- c(0, 5, 10, 15, 20)

sweep_results <- lapply(cutoffs, function(p) {
  n_drop <- ceiling(nrow(meta) * p / 100)

  drop_ids <- if (n_drop == 0) {
    character(0)
  } else {
    meta |>
      arrange(desc(B_M_ratio)) |>
      slice_head(n = n_drop) |>
      pull(sample_id)
  }

  keep_ids <- setdiff(meta$sample_id, drop_ids)

  # Build subset DAList: strip fitted objects, keep data + metadata + annotation.
  # Use `[<-`(list(NULL)) rather than `$<- NULL` to preserve slot names —
  # validate_DAList requires all 7 named slots to exist (even if NULL).
  dal_sub          <- sam
  dal_sub$data     <- dal_sub$data[, keep_ids, drop = FALSE]
  dal_sub$metadata <- dal_sub$metadata[keep_ids, , drop = FALSE]
  dal_sub["design"]     <- list(NULL)
  dal_sub["eBayes_fit"] <- list(NULL)
  dal_sub["results"]    <- list(NULL)

  message(sprintf("  cutoff=%2d%%  drop=%d  keep=%d samples",
                  p, n_drop, length(keep_ids)))

  n_DEP <- tryCatch(
    fit_for_subset(dal_sub),
    error = function(e) {
      warning(sprintf("cutoff=%d%% failed: %s", p, conditionMessage(e)))
      NA_integer_
    }
  )

  data.frame(
    cutoff_pct         = p,
    n_samples_dropped  = n_drop,
    n_samples_kept     = length(keep_ids),
    n_proteins_used    = nrow(dal_sub$data),
    n_DEP_fdr10        = n_DEP,
    dropped_sample_ids = paste(drop_ids, collapse = ";"),
    stringsAsFactors   = FALSE
  )
})

supp_B_data <- do.call(rbind, sweep_results)

message("\nSweep complete:")
print(supp_B_data[, c("cutoff_pct", "n_samples_dropped", "n_samples_kept",
                       "n_DEP_fdr10")])

# ---------------------------------------------------------------------------
# Write CSV
# ---------------------------------------------------------------------------

write_csv(
  supp_B_data,
  "02-03_Sam's_Results/04_Figures/F02_blood_contamination/c_data/supp_B_cutoff_sensitivity.csv"
)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

baseline <- supp_B_data$n_DEP_fdr10[supp_B_data$cutoff_pct == 0]

p_supp_B <- ggplot(supp_B_data, aes(cutoff_pct, n_DEP_fdr10)) +
  geom_hline(yintercept = baseline, linetype = "dashed", color = "grey50") +
  geom_line(color = "steelblue") +
  geom_point(color = "steelblue", size = 2.5) +
  geom_text(aes(label = n_DEP_fdr10), vjust = -1, size = 3) +
  labs(
    x       = "% of top-B_M_ratio samples dropped",
    y       = "Cancer_vs_Healthy DEPs at FDR < 0.10",
    title   = "DEP count sensitivity to blood-contamination cutoff",
    caption = paste0(
      "Dashed line = baseline (0% drop, n = ", baseline, " DEPs). ",
      "Refitted via proteoDA::fit_limma_model (robust eBayes, ~0 + cancer_time + sex). ",
      "n = ", nrow(meta), " total samples; top-B_M_ratio samples dropped sequentially."
    )
  ) +
  scale_x_continuous(breaks = c(0, 5, 10, 15, 20)) +
  theme_minimal(base_size = 10) +
  theme(plot.caption = element_text(size = 7, color = "grey40", hjust = 0))

ggsave(
  "02-03_Sam's_Results/04_Figures/F02_blood_contamination/b_reports/supp/png/panels/supp_B.png",
  p_supp_B,
  width  = 5.5,
  height = 4,
  dpi    = 300,
  bg     = "white"
)

message("Supp B done.")

# Export for driver script
supp_B <- p_supp_B
