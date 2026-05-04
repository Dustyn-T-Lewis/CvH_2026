# Three-way DEP comparison:
#   (1) Sam's limma  - his methodology on his filtered/normalized data
#   (2) Ours on his  - our methodology on Sam's filtered/normalized data
#   (3) Ours on ours - our methodology on our filtered/normalized data
#
# Outputs:
#   comparison_3way.xlsx   one sheet per matched contrast + a summary sheet
#   comparison_3way_summary.csv   per-contrast Spearman rho, N significant, overlap

library(dplyr)
library(readr)
library(readxl)
library(openxlsx)
library(tibble)

setwd(rprojroot::find_rstudio_root_file())

cfg <- list(
  sam_xlsx     = "02-03_Sam's_Results/limma_results_all_muscle_str.xlsx",
  ours_on_his  = "02-03_Sam's_Results/our_rerun/c_data/04_per_contrast_results",
  ours_on_ours = "03_DEP/c_data/04_per_contrast_results",
  out_xlsx     = "02-03_Sam's_Results/our_rerun/c_data/comparison_3way.xlsx",
  out_csv      = "02-03_Sam's_Results/our_rerun/c_data/comparison_3way_summary.csv",
  fdr_thresh   = 0.10
)

# Sam contrast -> our contrast (exact-math only; close matches noted below)
contrast_map <- list(
  Baseline_SURVvCTL = "Cancer_vs_Healthy",        # CR vs H at T1
  Baseline_CREvPLA  = "Baseline_Supplement",      # CRE_T1 vs PLA_T1
  Training_CRE      = "Training_CRE",             # T2-T1 in CRE arm
  Training_PLA      = "Training_PLA",             # T2-T1 in PLA arm
  Interaction_supp  = "Supplement_Interaction"    # (CRE_T2-T1) - (PLA_T2-T1)
)

# Helper -- read one of our per-contrast CSVs and standardize
read_ours <- function(dir, cname) {
  f <- file.path(dir, paste0(cname, ".csv"))
  if (!file.exists(f)) {
    warning("Missing: ", f)
    return(NULL)
  }
  out <- read_csv(f, show_col_types = FALSE) |>
    distinct(uniprot_id, .keep_all = TRUE)
  if (!"gene" %in% names(out)) out$gene <- NA_character_
  out |> select(uniprot_id, gene, logFC, t, P.Value, adj.P.Val)
}

# Helper -- read Sam's limma sheet
read_sam <- function(sheet) {
  out <- read_xlsx(cfg$sam_xlsx, sheet = sheet) |>
    distinct(uniprot_id, .keep_all = TRUE)
  if (!"gene" %in% names(out)) out$gene <- NA_character_
  out |> select(uniprot_id, gene, logFC, t, P.Value, adj.P.Val)
}

# Sheet-name uniqueness guard (substr 31 is Excel's hard limit)
stopifnot(!any(duplicated(substr(names(contrast_map), 1, 31))))

cat("Sheets in Sam's xlsx:\n")
print(excel_sheets(cfg$sam_xlsx))

summary_rows <- list()
wb <- createWorkbook()

for (sam_c in names(contrast_map)) {
  ours_c <- contrast_map[[sam_c]]
  cat(sprintf("\n--- %s  (Sam)  <->  %s  (Ours)\n", sam_c, ours_c))

  sam <- read_sam(sam_c)
  oh  <- read_ours(cfg$ours_on_his,  ours_c)
  oo  <- read_ours(cfg$ours_on_ours, ours_c)

  if (is.null(sam) || is.null(oh) || is.null(oo)) next

  # Three-way join on uniprot_id; gene from Sam (most complete)
  three <- sam |>
    rename(sam_logFC = logFC, sam_t = t, sam_P = P.Value, sam_adj = adj.P.Val) |>
    full_join(oh |> rename(oh_logFC = logFC, oh_t = t, oh_P = P.Value, oh_adj = adj.P.Val) |>
                select(-any_of("gene")),
              by = "uniprot_id") |>
    full_join(oo |> rename(oo_logFC = logFC, oo_t = t, oo_P = P.Value, oo_adj = adj.P.Val) |>
                select(-any_of("gene")),
              by = "uniprot_id")

  n_sam_total <- sum(!is.na(three$sam_logFC))
  n_oh_total  <- sum(!is.na(three$oh_logFC))
  n_oo_total  <- sum(!is.na(three$oo_logFC))
  matched <- three |> filter(!is.na(sam_logFC) & !is.na(oh_logFC) & !is.na(oo_logFC))

  # Significance flags at FDR < 0.10
  sam_sig <- which(matched$sam_adj < cfg$fdr_thresh & !is.na(matched$sam_adj))
  oh_sig  <- which(matched$oh_adj  < cfg$fdr_thresh & !is.na(matched$oh_adj))
  oo_sig  <- which(matched$oo_adj  < cfg$fdr_thresh & !is.na(matched$oo_adj))

  rho_sam_oh <- cor(matched$sam_logFC, matched$oh_logFC, method = "spearman", use = "complete.obs")
  rho_sam_oo <- cor(matched$sam_logFC, matched$oo_logFC, method = "spearman", use = "complete.obs")
  rho_oh_oo  <- cor(matched$oh_logFC,  matched$oo_logFC, method = "spearman", use = "complete.obs")

  summary_rows[[sam_c]] <- tibble(
    sam_contrast = sam_c,
    our_contrast = ours_c,
    n_sam_total  = n_sam_total,
    n_oh_total   = n_oh_total,
    n_oo_total   = n_oo_total,
    n_matched    = nrow(matched),
    n_sig_sam    = length(sam_sig),
    n_sig_oh     = length(oh_sig),
    n_sig_oo     = length(oo_sig),
    overlap_sam_oh = length(intersect(sam_sig, oh_sig)),
    overlap_sam_oo = length(intersect(sam_sig, oo_sig)),
    overlap_oh_oo  = length(intersect(oh_sig,  oo_sig)),
    rho_sam_oh = round(rho_sam_oh, 3),
    rho_sam_oo = round(rho_sam_oo, 3),
    rho_oh_oo  = round(rho_oh_oo,  3)
  )

  cat(sprintf("  Totals: Sam=%d, OurOnHis=%d, OurOnOurs=%d | matched=%d (dropped Sam=%d, OnHis=%d, OnOurs=%d)\n",
              n_sam_total, n_oh_total, n_oo_total, nrow(matched),
              n_sam_total - nrow(matched),
              n_oh_total  - nrow(matched),
              n_oo_total  - nrow(matched)))
  cat(sprintf("  Sig (FDR<%.2f): Sam=%d, OurOnHis=%d, OurOnOurs=%d\n",
              cfg$fdr_thresh, length(sam_sig), length(oh_sig), length(oo_sig)))
  cat(sprintf("  rho: Sam-OurOnHis = %.3f, Sam-OurOnOurs = %.3f, OurOnHis-OurOnOurs = %.3f\n",
              rho_sam_oh, rho_sam_oo, rho_oh_oo))

  # Sheet: per-protein 3-way table, sorted by min adj.P
  sheet_df <- matched |>
    mutate(min_adj = suppressWarnings(pmin(sam_adj, oh_adj, oo_adj, na.rm = TRUE)),
           min_adj = ifelse(is.infinite(min_adj), NA_real_, min_adj)) |>
    arrange(min_adj) |>
    select(uniprot_id, gene,
           sam_logFC, sam_adj,
           oh_logFC,  oh_adj,
           oo_logFC,  oo_adj,
           min_adj) |>
    mutate(across(c(sam_logFC, oh_logFC, oo_logFC), \(x) round(x, 3)),
           across(c(sam_adj, oh_adj, oo_adj, min_adj), \(x) signif(x, 3)))

  sname <- substr(sam_c, 1, 31)
  addWorksheet(wb, sname)
  writeData(wb, sname,
            sprintf("%s (Sam)  <->  %s (Ours).  Three-way comparison.  FDR threshold = %.2f",
                    sam_c, ours_c, cfg$fdr_thresh),
            startRow = 1)
  addStyle(wb, sname,
           createStyle(textDecoration = "bold", fontSize = 11),
           rows = 1, cols = 1)
  writeData(wb, sname, sheet_df, startRow = 3,
            headerStyle = createStyle(textDecoration = "bold", fgFill = "#DCE6F1"))
  freezePane(wb, sname, firstActiveRow = 4, firstActiveCol = 3)
  setColWidths(wb, sname, cols = seq_len(ncol(sheet_df)), widths = "auto")
}

summary_df <- bind_rows(summary_rows)
addWorksheet(wb, "Summary")
writeData(wb, "Summary",
          "Three-way comparison summary.  rho = Spearman of logFC.  Sig = adj.P.Val < 0.10.",
          startRow = 1)
addStyle(wb, "Summary",
         createStyle(textDecoration = "bold", fontSize = 12),
         rows = 1, cols = 1)
writeData(wb, "Summary", summary_df, startRow = 3,
          headerStyle = createStyle(textDecoration = "bold", fgFill = "#DCE6F1"))
freezePane(wb, "Summary", firstActiveRow = 4, firstActiveCol = 3)
setColWidths(wb, "Summary", cols = seq_len(ncol(summary_df)), widths = "auto")

saveWorkbook(wb, cfg$out_xlsx, overwrite = TRUE)
write_csv(summary_df, cfg$out_csv)

cat("\nWrote: ", cfg$out_xlsx, "\n")
cat("Wrote: ", cfg$out_csv, "\n")
print(as.data.frame(summary_df))
