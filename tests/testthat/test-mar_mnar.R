# Real-data tests for the MAR/MNAR consensus classifier.
library(testthat)
source(here::here("R", "mar_mnar.R"))

norm_csv <- here::here("02_Normalization", "c_data", "02_normalized.csv")

test_that("classify_mar_mnar returns the expected structure on the normalized matrix", {
  skip_if_not(file.exists(norm_csv), "Stage 02 output not present")
  df  <- readr::read_csv(norm_csv, show_col_types = FALSE)
  ann <- dplyr::select(df, uniprot_id, gene, protein, description)
  mat <- as.matrix(df[, setdiff(names(df), names(ann))]); rownames(mat) <- ann$uniprot_id

  set.seed(42)
  cl <- classify_mar_mnar(mat)

  expect_equal(nrow(cl), nrow(mat))
  expect_setequal(unique(cl$classification), c("Complete", "MAR", "MNAR"))
  # Complete proteins have zero missing; nothing else does
  expect_true(all(cl$n_miss[cl$classification == "Complete"] == 0))
  expect_true(all(cl$n_miss[cl$classification != "Complete"] > 0))
  # reliability flag matches its definition
  expect_equal(cl$imputation_reliable,
               cl$classification == "Complete" | cl$pct_miss < 50)
})
